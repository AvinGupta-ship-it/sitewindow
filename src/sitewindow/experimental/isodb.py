"""Loader for the NIST/ARPA-E ISODB GitHub mirror and carbon H2 isotherm selection.

The mirror stores one JSON file per isotherm under ``Library/<doi-derived>/``. This module
parses those files, selects single-component hydrogen isotherms on carbon-like adsorbents,
normalizes units, and writes a harmonized Parquet table with per-record provenance.
"""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from sitewindow.units import UnitError, loading_to_mol_per_kg, pressure_to_bar

ISODB_GIT_URL = "https://github.com/NIST-ISODB/isodb-library"

_CARBON_PATTERNS = [
    r"\bcarbon\b",
    r"graph",
    r"\bcnt\b",
    r"nanotube",
    r"\bswnt",
    r"\bmwnt",
    r"swcnt",
    r"mwcnt",
    r"\bac\b",
    r"\bac-",
    r"activated",
    r"\bmsc",
    r"ax-?21",
    r"maxsorb",
    r"\bcmk",
    r"\bztc",
    r"templated",
    r"aerogel",
    r"fiber",
    r"fibre",
    r"\bcdc\b",
    r"carbide[- ]derived",
    r"charcoal",
    r"\bchar\b",
    r"norit",
    r"nanofib",
    r"carbon black",
    r"anthracite",
    r"coke",
    r"biochar",
    r"hydrochar",
    r"xerogel",
    r"ordered mesoporous carbon",
]
_NON_CARBON_PATTERNS = [
    r"mof",
    r"zif",
    r"uio",
    r"hkust",
    r"irmof",
    r"cu-btc",
    r"zeolite",
    r"silica",
    r"\bsba",
    r"\bmcm",
    r"alumina",
    r"mg\b",
    r"lani",
    r"tio2",
    r"\bpd\b",
    r"\bpt\b",
    r"\bni\b",
    r"platinum",
    r"palladium",
    r"nickel",
    r"cof-",
    r"pim-",
    r"polymer",
    r"hcp",
    r"boron",
    r"\bbn\b",
    r"nitride",
    r"\bcu\b",
    r"\bag\b",
    r"\bli\b",
    r"lithium",
    r"\bna\b",
    r"\bk\b",
    r"potassium",
    r"sodium",
    r"titan",
    r"zirc",
    r"\bzn\b",
    r"\bco\b",
    r"cobalt",
    r"\bfe\b",
    r"iron",
    r"ruthen",
    r"spillover",
]


@dataclass
class SelectionResult:
    table: pd.DataFrame
    excluded: pd.DataFrame
    notes: dict = field(default_factory=dict)


def _head_commit(repo: Path) -> str:
    out = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    )
    return out.stdout.strip()


def clone_mirror(dest: Path, depth: int = 1, commit: str | None = None) -> str:
    """Clone the ISODB mirror if absent and return the checked-out commit hash.

    When ``commit`` is given (the hash recorded in ``data/sources.lock``), that commit is
    fetched and checked out so that the corpus is pinned independently of the mirror's
    current HEAD (manual Section 6.2). Without ``commit`` the mirror's HEAD is used.
    """
    dest = Path(dest)
    if not (dest / ".git").exists():
        subprocess.run(
            ["git", "clone", "--depth", str(depth), ISODB_GIT_URL, str(dest)], check=True
        )
    if commit is not None and _head_commit(dest) != commit:
        subprocess.run(
            ["git", "-C", str(dest), "fetch", "--depth", str(depth), "origin", commit], check=True
        )
        subprocess.run(["git", "-C", str(dest), "checkout", "--quiet", commit], check=True)
    return _head_commit(dest)


def _is_h2(d: dict) -> bool:
    ads = d.get("adsorbates", [])
    if len(ads) != 1:
        return False
    name = str(ads[0].get("name", "")).strip().lower()
    return name in {"hydrogen", "h2", "dihydrogen"}


def classify_material(name: str) -> tuple[str, bool]:
    """Return (material_class, is_carbon_only) from an adsorbent name."""
    n = name.lower()
    carbon = any(re.search(p, n) for p in _CARBON_PATTERNS)
    non = any(re.search(p, n) for p in _NON_CARBON_PATTERNS)
    if not carbon:
        return "non_carbon", False
    if non:
        return "carbon_composite_or_doped", False
    if re.search(r"swcnt|swnt|single[- ]wall", n):
        return "SWCNT", True
    if re.search(r"mwcnt|mwnt|multi[- ]wall|nanotube|cnt", n):
        return "MWCNT", True
    if re.search(r"reduced graphene|rgo|r-go", n):
        return "rGO", True
    if re.search(r"graphene oxide|\bgo\b", n):
        return "GO", True
    if re.search(r"graphene|graphite|graphit", n):
        return "graphene", True
    if re.search(r"templat|ztc|cmk|ordered mesoporous", n):
        return "templated_carbon", True
    if re.search(r"activated|\bac\b|\bac-|maxsorb|ax-?21|norit|msc", n):
        return "activated_carbon", True
    return "carbon_other", True


def _parse_isotherm(path: Path) -> dict | None:
    try:
        d = json.loads(path.read_text())
    except Exception:
        return None
    if not _is_h2(d):
        return None
    data = d.get("isotherm_data", [])
    if not data:
        return None
    P = np.array([pt.get("pressure", np.nan) for pt in data], dtype=float)
    n = np.array([pt.get("total_adsorption", np.nan) for pt in data], dtype=float)
    return {
        "isodb_filename": d.get("filename", path.name),
        "doi": d.get("DOI", ""),
        "adsorbent_name": d.get("adsorbent", {}).get("name", ""),
        "adsorbent_hashkey": d.get("adsorbent", {}).get("hashkey", ""),
        "temperature_K": float(d.get("temperature", np.nan)),
        "pressure_units": d.get("pressureUnits", ""),
        "adsorption_units": d.get("adsorptionUnits", ""),
        "category": d.get("category", ""),
        "digitized": bool(d.get("digitizer")) and not bool(d.get("tabular_data")),
        "tabular": bool(d.get("tabular_data")),
        "n_points": int(len(data)),
        "pressure_raw": P,
        "loading_raw": n,
    }


def load_carbon_h2(mirror_root: Path) -> SelectionResult:
    """Parse the mirror, select carbon-only H2 isotherms, normalize units."""
    files = sorted(Path(mirror_root, "Library").glob("*/*.json"))
    rows, excluded = [], []
    for f in files:
        rec = _parse_isotherm(f)
        if rec is None:
            continue
        mclass, carbon_only = classify_material(rec["adsorbent_name"])
        rec["material_class"] = mclass
        if not carbon_only:
            excluded.append(
                {
                    **{k: v for k, v in rec.items() if k not in ("pressure_raw", "loading_raw")},
                    "reason": f"not carbon-only ({mclass})",
                }
            )
            continue
        try:
            P_bar = pressure_to_bar(rec["pressure_raw"], rec["pressure_units"])
            n_mol_kg = loading_to_mol_per_kg(rec["loading_raw"], rec["adsorption_units"])
        except UnitError as e:
            excluded.append(
                {
                    **{k: v for k, v in rec.items() if k not in ("pressure_raw", "loading_raw")},
                    "reason": f"unit: {e}",
                }
            )
            continue
        ok = np.isfinite(P_bar) & np.isfinite(n_mol_kg)
        if ok.sum() < 3:
            excluded.append(
                {
                    **{k: v for k, v in rec.items() if k not in ("pressure_raw", "loading_raw")},
                    "reason": "fewer than 3 finite points",
                }
            )
            continue
        order = np.argsort(P_bar[ok])
        P_bar, n_mol_kg = P_bar[ok][order], n_mol_kg[ok][order]
        flags = []
        if np.any(np.diff(n_mol_kg) < -0.05 * np.nanmax(np.abs(n_mol_kg))):
            flags.append("nonmonotonic")
        if np.any(np.diff(P_bar) == 0):
            flags.append("duplicate_pressure")
        if np.any(n_mol_kg < 0):
            flags.append("negative_loading")
        rec.update(
            {
                "pressure_bar": P_bar,
                "loading_mol_kg": n_mol_kg,
                "quality_flags": ",".join(flags),
                "loading_basis": "excess" if "excess" in rec["category"].lower() else "unknown",
            }
        )
        rec.pop("pressure_raw")
        rec.pop("loading_raw")
        rows.append(rec)
    table = pd.DataFrame(rows)
    exc = pd.DataFrame(excluded)
    notes = {"n_files": len(files), "n_selected": len(table), "n_excluded": len(exc)}
    return SelectionResult(table, exc, notes)


def write_parquet(result: SelectionResult, out_dir: Path, commit: str) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t = result.table.copy()
    t["isodb_commit"] = commit
    t["pressure_bar"] = t["pressure_bar"].apply(lambda a: a.tolist())
    t["loading_mol_kg"] = t["loading_mol_kg"].apply(lambda a: a.tolist())
    p = out_dir / "isodb_carbon_h2.parquet"
    t.to_parquet(p, index=False)
    result.excluded.to_csv(out_dir / "isodb_carbon_h2_excluded.csv", index=False)
    return p


__all__ = [
    "clone_mirror",
    "load_carbon_h2",
    "classify_material",
    "write_parquet",
    "SelectionResult",
]
