"""Material metadata registry (S2 role): BET areas, pore volumes and source types of corpus
materials, read from ``data/registry/material_metadata.yaml``.

The ISODB records carry no surface areas or pore volumes, so every value here was read in the
source publication (abstract or full text) and is stored with the place where it was read. A
material is keyed by DOI and adsorbent hashkey like the corpus material sets. Entries may also
mark a record as simulated rather than measured (``simulated: true``), which excludes it from
the experimental inference, or as ambiguous (``sample_identity: ambiguous``) when the paper
measured several carbons and the database name does not say which one the isotherm belongs to;
ambiguous materials keep their isotherms but receive no surface area. A value read only in
another paper's table that cites the source (``value_source: secondary``) is kept for the record
but does not enter benchmark B7, following the rule of the published-value registry that cited
values are not registered as the source's own.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from sitewindow.experimental.corpus import material_id

REGISTRY = Path("data/registry/material_metadata.yaml")


@dataclass(frozen=True)
class MaterialMetadata:
    material_id: str
    doi: str
    adsorbent_hashkey: str
    adsorbent_name: str
    bet_area_m2_g: float | None = None
    bet_area_sd_m2_g: float | None = None
    pore_volume_cm3_g: float | None = None
    sample_identity: str = "unambiguous"
    simulated: bool = False
    value_source: str = "primary"
    where_read: str = ""
    access: str = ""
    notes: str = ""

    @property
    def usable_area(self) -> bool:
        """A surface area enters benchmark B7 only when it was read in the paper itself (or a
        copy of it), belongs to an identified sample and the isotherms are measurements."""
        return (
            self.bet_area_m2_g is not None
            and self.sample_identity == "unambiguous"
            and not self.simulated
            and self.value_source == "primary"
        )


def load_metadata(path: Path | str = REGISTRY) -> dict[str, MaterialMetadata]:
    """Metadata keyed by material id; an absent registry gives an empty mapping."""
    path = Path(path)
    if not path.exists():
        return {}
    doc = yaml.safe_load(path.read_text()) or {}
    out: dict[str, MaterialMetadata] = {}
    for e in doc.get("materials", []) or []:
        mid = material_id(str(e["adsorbent_hashkey"]), str(e["doi"]))
        out[mid] = MaterialMetadata(
            material_id=mid,
            doi=str(e["doi"]),
            adsorbent_hashkey=str(e["adsorbent_hashkey"]),
            adsorbent_name=str(e.get("adsorbent_name", "")),
            bet_area_m2_g=_opt_float(e.get("bet_area_m2_g")),
            bet_area_sd_m2_g=_opt_float(e.get("bet_area_sd_m2_g")),
            pore_volume_cm3_g=_opt_float(e.get("pore_volume_cm3_g")),
            sample_identity=str(e.get("sample_identity", "unambiguous")),
            simulated=bool(e.get("simulated", False)),
            value_source=str(e.get("value_source", "primary")),
            where_read=str(e.get("where_read", "")),
            access=str(e.get("access", "")),
            notes=str(e.get("notes", "")),
        )
    return out


def _opt_float(x) -> float | None:
    return None if x is None else float(x)


def pore_volumes(meta: dict[str, MaterialMetadata]) -> dict[str, float]:
    """Reported pore volumes by material id (primary sources, identified samples), for the
    excess-to-absolute prior."""
    return {
        k: m.pore_volume_cm3_g
        for k, m in meta.items()
        if m.pore_volume_cm3_g is not None
        and m.sample_identity == "unambiguous"
        and m.value_source == "primary"
    }


def simulated_ids(meta: dict[str, MaterialMetadata]) -> set[str]:
    return {k for k, m in meta.items() if m.simulated}


__all__ = ["MaterialMetadata", "REGISTRY", "load_metadata", "pore_volumes", "simulated_ids"]
