"""Published-value registry (manual C3; source registry S4).

Loads ``data/registry/published_values.yaml``, validates the site-energy records against the
``PublishedValue`` schema and presents them as the frame the meta-analysis consumes. Records whose
quantity is not a site energy (model parameters of benchmarks B9 and B10, experimental
observations) stay in the raw list and are not part of the frame. Compilation entries that report
a range of earlier values carry a uniform distribution and are excluded from the meta-analysis
unless requested, because their studies are not independent of the single-study entries.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from sitewindow.schema import (
    Distribution,
    Extraction,
    PublishedValue,
    Quantity,
    ReferenceState,
    SiteClassName,
    UncertainScalar,
)
from sitewindow.units import energy_to_ev

REGISTRY = Path("data/registry/published_values.yaml")
DEFAULT_SIGMA_EV = 0.05
"""Default standard error of a published value without a stated uncertainty (manual C3): the
convergence practice of plane-wave chemisorption energies."""
TIER0_LEVELS = frozenset({"DMC", "CCSD(T)", "CCSD(T)/CBS", "QMC"})


def load_registry(path: Path = REGISTRY) -> list[dict[str, Any]]:
    """Every raw entry of the registry file, in order."""
    entries = yaml.safe_load(Path(path).read_text()) or []
    if not isinstance(entries, list):
        raise ValueError(f"{path} must hold a list of records")
    return entries


def _scalar(entry: dict[str, Any]) -> tuple[UncertainScalar, bool] | None:
    """The value of an entry as an uncertain scalar in the entry's unit, with a flag marking a
    range (compilation) entry; ``None`` for structured values."""
    value = entry.get("value")
    unc = entry.get("uncertainty")
    if isinstance(value, bool):
        return None
    if isinstance(value, int | float):
        sigma = float(unc) if isinstance(unc, int | float) and unc is not None else 0.0
        return UncertainScalar(value=float(value), sigma=sigma), False
    if (
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(v, int | float) for v in value)
    ):
        lo, hi = sorted(float(v) for v in value)
        return (
            UncertainScalar(
                value=0.5 * (lo + hi),
                sigma=0.5 * (hi - lo) / math.sqrt(3.0),
                distribution=Distribution.uniform,
            ),
            True,
        )
    return None


def is_site_energy(entry: dict[str, Any]) -> bool:
    """Whether a raw entry is a scalar site energy (a quantity and reference state of the
    schema with a numeric value or range), as opposed to a parameter or observation record."""
    return (
        entry.get("quantity") in Quantity.__members__
        and entry.get("reference_state") in ReferenceState.__members__
        and _scalar(entry) is not None
    )


def published_value(entry: dict[str, Any]) -> PublishedValue:
    """One raw site-energy entry as a validated ``PublishedValue`` in eV."""
    scalar = _scalar(entry)
    if scalar is None:
        raise ValueError(f"{entry.get('id')}: value is not a scalar or a range")
    us, compilation = scalar
    unit = str(entry.get("unit", "eV"))
    value_eV = float(energy_to_ev(us.value, unit))
    sigma_eV = float(energy_to_ev(us.sigma, unit))
    cls = entry.get("site_class")
    notes = entry.get("notes")
    if compilation:
        notes = (notes + "; " if notes else "") + "range of earlier values (compilation)"
    arxiv = entry.get("arxiv")
    return PublishedValue(
        id=str(entry["id"]),
        doi=entry.get("doi"),
        arxiv=None if arxiv is None else str(arxiv),
        system=str(entry.get("system", "")),
        site_class=SiteClassName(cls) if cls in SiteClassName.__members__ else None,
        site_id=entry.get("site_id"),
        quantity=Quantity(str(entry["quantity"])),
        value=UncertainScalar(value=value_eV, sigma=sigma_eV, distribution=us.distribution),
        unit="eV",
        reference_state=ReferenceState(str(entry["reference_state"])),
        method=dict(entry.get("method") or {}),
        extraction=Extraction(entry.get("extraction", "tabular_text")),
        source_id=str(entry.get("source_id", "S4")),
        notes=notes,
    )


def published_values(path: Path = REGISTRY, strict: bool = True) -> list[PublishedValue]:
    """The site-energy records of the registry as validated ``PublishedValue`` models, values in
    eV; a record whose value, quantity or reference state is not a site energy is left out. With
    ``strict`` an invalid site-energy entry raises; otherwise it is skipped."""
    out: list[PublishedValue] = []
    for entry in load_registry(path):
        if not is_site_energy(entry):
            continue
        try:
            out.append(published_value(entry))
        except ValueError:
            if strict:
                raise
    return out


_SUPERCELL = re.compile(r"(\d+)\s*[x×]\s*(\d+)")


def cell_atoms(supercell: Any) -> float:
    """Carbon atoms of a graphene ``n x m`` supercell named in a method record; ``nan`` when the
    cell is a tube, a cluster or unstated, so that the cell-size term of the meta-analysis is
    switched off for it."""
    if supercell is None:
        return float("nan")
    m = _SUPERCELL.search(str(supercell))
    if m is None:
        return float("nan")
    return 2.0 * int(m.group(1)) * int(m.group(2))


def registry_frame(
    values: list[PublishedValue],
    default_sigma_eV: float = DEFAULT_SIGMA_EV,
    include_compilations: bool = False,
) -> pd.DataFrame:
    """The meta-analysis frame: one row per published site energy with its class, method
    descriptors, standard error (stated, else ``default_sigma_eV``) and study identifier."""
    rows: list[dict[str, Any]] = []
    for v in values:
        compilation = v.value.distribution == Distribution.uniform
        if compilation and not include_compilations:
            continue
        method = v.method
        level = str(method.get("level") or "DFT")
        functional = str(method.get("functional") or level)
        n_cell = cell_atoms(method.get("supercell"))
        rows.append(
            {
                "id": v.id,
                "study": v.doi or f"arXiv:{v.arxiv}",
                "site_class": None if v.site_class is None else v.site_class.value,
                "site_id": v.site_id,
                "quantity": v.quantity.value,
                "reference_state": v.reference_state.value,
                "value_eV": v.value.value,
                "sigma_eV": v.value.sigma if v.value.sigma > 0 else default_sigma_eV,
                "stated_sigma": bool(v.value.sigma > 0),
                "level": level,
                "tier0": level in TIER0_LEVELS,
                "functional": functional,
                "dispersion": str(method.get("dispersion") or "none"),
                "code": method.get("code"),
                "supercell": method.get("supercell"),
                "n_cell": n_cell,
                "inv_n_cell": 0.0 if math.isnan(n_cell) else 1.0 / n_cell,
                "spin": method.get("spin"),
                "extraction": v.extraction.value,
                "compilation": compilation,
            }
        )
    return pd.DataFrame(rows)


def method_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Counts of registry values by functional, level and class."""
    if frame.empty:
        return pd.DataFrame(columns=["functional", "level", "site_class", "n"])
    return (
        frame.groupby(["functional", "level", "site_class"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )


def check_registry(path: Path = REGISTRY) -> dict[str, Any]:
    """Validation summary of the registry: every entry needs an id, a DOI or arXiv identifier and
    an extraction field; site-energy entries need a quantity, a reference state and a unit."""
    problems: list[str] = []
    warnings: list[str] = []
    entries = load_registry(path)
    seen: set[str] = set()
    values: list[PublishedValue] = []
    for entry in entries:
        ident = str(entry.get("id", "<missing id>"))
        if ident in seen:
            problems.append(f"{ident}: duplicate id")
        seen.add(ident)
        site_energy = is_site_energy(entry)
        if not (entry.get("doi") or entry.get("arxiv")):
            (problems if site_energy else warnings).append(f"{ident}: no DOI or arXiv identifier")
        if entry.get("extraction") not in Extraction.__members__:
            problems.append(f"{ident}: extraction field missing or unknown")
        if entry.get("quantity") in Quantity.__members__:
            if entry.get("reference_state") not in ReferenceState.__members__:
                problems.append(f"{ident}: site energy without a reference state")
            if not entry.get("unit"):
                problems.append(f"{ident}: site energy without a unit")
        if site_energy:
            try:
                values.append(published_value(entry))
            except ValueError as exc:
                problems.append(f"{ident}: {str(exc).splitlines()[0]}")
    frame = registry_frame(values, include_compilations=True)
    return {
        "n_entries": len(entries),
        "n_site_energies": len(values),
        "n_compilations": int(frame["compilation"].sum()) if not frame.empty else 0,
        "classes": sorted(c for c in frame["site_class"].dropna().unique())
        if not frame.empty
        else [],
        "problems": problems,
        "warnings": warnings,
    }


__all__ = [
    "DEFAULT_SIGMA_EV",
    "REGISTRY",
    "TIER0_LEVELS",
    "cell_atoms",
    "check_registry",
    "load_registry",
    "method_table",
    "published_values",
    "registry_frame",
]
