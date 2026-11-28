"""Material sets of the harmonized carbon H2 corpus and the admissibility rules of C6.

A material is one adsorbent record of one publication, keyed by the ISODB adsorbent hashkey
together with the DOI: the database assigns generic names such as "Carbon" to the comparison
carbons of many unrelated papers under one hashkey, so the hashkey alone does not identify a
sample, while within one paper the digitized isotherms of one adsorbent name belong to one
sample. Every isotherm that the inference does not use is logged with its reason (acceptance
criterion v); points are dropped only when they cannot enter a logarithmic pressure axis or a
loading model (non-positive pressure, negative loading).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from sitewindow.experimental.eos import T_CRIT_H2_K, h2_saturation_pressure_bar

PRESSURE_MIN_BAR = 1.0e-12
PRESSURE_MAX_BAR = 1000.0
MIN_POINTS = 3
GRAPHITIC_CLASSES = ("graphene", "SWCNT", "MWCNT")


def material_id(hashkey: str, doi: str) -> str:
    """Stable identifier of a (adsorbent, publication) pair usable as a file name."""
    key = hashkey.replace("NIST-MATDB-", "")[:12]
    stub = re.sub(r"[^A-Za-z0-9.]+", "_", doi).strip("_")
    return f"{key}__{stub}"


@dataclass
class IsothermData:
    filename: str
    temperature_K: float
    pressure_bar: np.ndarray
    loading_mol_kg: np.ndarray
    digitized: bool
    tabular: bool
    quality_flags: str = ""
    points_dropped: int = 0

    @property
    def n_points(self) -> int:
        return int(len(self.pressure_bar))


@dataclass
class MaterialSet:
    material_id: str
    doi: str
    adsorbent_hashkey: str
    adsorbent_name: str
    material_class: str
    isotherms: list[IsothermData] = field(default_factory=list)

    @property
    def temperatures(self) -> list[float]:
        return sorted({float(i.temperature_K) for i in self.isotherms})

    @property
    def n_temperatures(self) -> int:
        return len(self.temperatures)

    @property
    def generic_name(self) -> bool:
        return self.adsorbent_name.strip().lower() in {"carbon", "activated carbon", "ac"}

    @property
    def graphitic(self) -> bool:
        return self.material_class in GRAPHITIC_CLASSES


def load_corpus(path: Path | str = "data/processed/isodb_carbon_h2.parquet") -> pd.DataFrame:
    """The harmonized corpus with the pressure and loading columns as arrays."""
    t = pd.read_parquet(path)
    t["pressure_bar"] = t["pressure_bar"].apply(lambda a: np.asarray(a, dtype=float))
    t["loading_mol_kg"] = t["loading_mol_kg"].apply(lambda a: np.asarray(a, dtype=float))
    return t


def _admit(row: pd.Series) -> tuple[IsothermData | None, str | None]:
    if str(row.get("category", "") or "").strip().lower() == "mod":
        return None, "modelled isotherm (ISODB category mod), not a measurement"
    P = np.asarray(row["pressure_bar"], dtype=float)
    n = np.asarray(row["loading_mol_kg"], dtype=float)
    T = float(row["temperature_K"])
    keep = np.isfinite(P) & np.isfinite(n) & (P > 0.0) & (n >= 0.0)
    dropped = int((~keep).sum())
    P, n = P[keep], n[keep]
    if len(P) < MIN_POINTS:
        return None, f"fewer than {MIN_POINTS} admissible points after dropping {dropped}"
    if P.min() < PRESSURE_MIN_BAR or P.max() > PRESSURE_MAX_BAR:
        return None, (
            f"pressure axis outside [{PRESSURE_MIN_BAR:g}, {PRESSURE_MAX_BAR:g}] bar "
            f"({P.min():.3g} to {P.max():.3g})"
        )
    if T < T_CRIT_H2_K and P.max() > h2_saturation_pressure_bar(T):
        return None, (
            f"subcritical temperature {T:g} K with pressures above saturation "
            f"({h2_saturation_pressure_bar(T):.3g} bar)"
        )
    order = np.argsort(P, kind="stable")
    iso = IsothermData(
        filename=str(row["isodb_filename"]),
        temperature_K=T,
        pressure_bar=P[order],
        loading_mol_kg=n[order],
        digitized=bool(row["digitized"]),
        tabular=bool(row["tabular"]),
        quality_flags=str(row.get("quality_flags", "") or ""),
        points_dropped=dropped,
    )
    return iso, None


def material_sets(table: pd.DataFrame) -> tuple[list[MaterialSet], pd.DataFrame]:
    """Group the corpus into material sets and log the isotherms the rules exclude."""
    sets: dict[str, MaterialSet] = {}
    excluded: list[dict] = []
    for _, row in table.iterrows():
        mid = material_id(str(row["adsorbent_hashkey"]), str(row["doi"]))
        iso, reason = _admit(row)
        if iso is None:
            excluded.append(
                {
                    "isodb_filename": row["isodb_filename"],
                    "material_id": mid,
                    "doi": row["doi"],
                    "adsorbent_name": row["adsorbent_name"],
                    "temperature_K": row["temperature_K"],
                    "reason": reason,
                }
            )
            continue
        ms = sets.get(mid)
        if ms is None:
            ms = MaterialSet(
                material_id=mid,
                doi=str(row["doi"]),
                adsorbent_hashkey=str(row["adsorbent_hashkey"]),
                adsorbent_name=str(row["adsorbent_name"]),
                material_class=str(row["material_class"]),
            )
            sets[mid] = ms
        ms.isotherms.append(iso)
    cols = ["isodb_filename", "material_id", "doi", "adsorbent_name", "temperature_K", "reason"]
    return list(sets.values()), pd.DataFrame(excluded, columns=cols)


def multi_temperature(sets: list[MaterialSet]) -> list[MaterialSet]:
    return [s for s in sets if s.n_temperatures >= 2]


__all__ = [
    "GRAPHITIC_CLASSES",
    "IsothermData",
    "MaterialSet",
    "load_corpus",
    "material_id",
    "material_sets",
    "multi_temperature",
]


def apply_metadata(
    sets: list[MaterialSet], excluded: pd.DataFrame, simulated: set[str]
) -> tuple[list[MaterialSet], pd.DataFrame]:
    """Remove the material sets the metadata registry marks as simulated, logging their isotherms."""
    keep, rows = [], []
    for s in sets:
        if s.material_id in simulated:
            for iso in s.isotherms:
                rows.append(
                    {
                        "isodb_filename": iso.filename,
                        "material_id": s.material_id,
                        "doi": s.doi,
                        "adsorbent_name": s.adsorbent_name,
                        "temperature_K": iso.temperature_K,
                        "reason": "modelled isotherm (material metadata registry: simulated), not a measurement",
                    }
                )
        else:
            keep.append(s)
    log = pd.concat([excluded, pd.DataFrame(rows, columns=excluded.columns)], ignore_index=True)
    return keep, log


__all__.append("apply_metadata")
