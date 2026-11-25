"""Flat tabular form of :class:`SiteEnergyRecord` for Parquet storage.

The table columns are the record fields with the nested provenance and descriptor blocks
flattened; :func:`frame_to_records` inverts :func:`records_to_frame` exactly, so the Parquet
landscape tables validate against the schema on load.
"""

from __future__ import annotations

import math
from datetime import datetime
from typing import Any

import pandas as pd

from sitewindow.schema.models import (
    Distribution,
    Extraction,
    Provenance,
    ProvenanceKind,
    Quantity,
    ReferenceState,
    SiteClassName,
    SiteDescriptors,
    SiteEnergyRecord,
    Tier,
    UncertainScalar,
)

DESCRIPTOR_COLUMNS = (
    "pyramidalization_deg",
    "coordination",
    "n_neighbor_H",
    "distance_to_defect_nm",
    "curvature_inv_nm",
)

RECORD_COLUMNS = (
    "config_hash",
    "site_class",
    "site_id",
    "n_H",
    "tier",
    "engine",
    "quantity",
    "value_eV",
    "sigma_eV",
    "distribution",
    "reference_state",
    "temperature_K",
    "converged",
    "n_atoms",
    "supercell",
    "settings_hash",
    "provenance_kind",
    "source_id",
    "identifier",
    "extraction",
    "package_version",
    "seed",
    "created",
    *DESCRIPTOR_COLUMNS,
)


def _none_if_nan(x: Any) -> Any:
    if x is None:
        return None
    if isinstance(x, float) and math.isnan(x):
        return None
    try:
        if pd.isna(x):
            return None
    except (TypeError, ValueError):
        pass
    return x


def record_to_row(r: SiteEnergyRecord) -> dict[str, Any]:
    d = r.descriptors or SiteDescriptors()
    return {
        "config_hash": r.config_hash,
        "site_class": None if r.site_class is None else r.site_class.value,
        "site_id": r.site_id,
        "n_H": r.n_H,
        "tier": r.tier.value,
        "engine": r.engine,
        "quantity": r.quantity.value,
        "value_eV": r.value.value,
        "sigma_eV": r.value.sigma,
        "distribution": r.value.distribution.value,
        "reference_state": r.reference_state.value,
        "temperature_K": r.temperature_K,
        "converged": r.converged,
        "n_atoms": r.n_atoms,
        "supercell": r.supercell,
        "settings_hash": r.provenance.settings_hash,
        "provenance_kind": r.provenance.kind.value,
        "source_id": r.provenance.source_id,
        "identifier": r.provenance.identifier,
        "extraction": None if r.provenance.extraction is None else r.provenance.extraction.value,
        "package_version": r.provenance.package_version,
        "seed": r.provenance.seed,
        "created": r.provenance.created.isoformat(),
        "pyramidalization_deg": d.pyramidalization_deg,
        "coordination": d.coordination,
        "n_neighbor_H": d.n_neighbor_H,
        "distance_to_defect_nm": d.distance_to_defect_nm,
        "curvature_inv_nm": d.curvature_inv_nm,
    }


def row_to_record(row: dict[str, Any]) -> SiteEnergyRecord:
    g = {k: _none_if_nan(v) for k, v in row.items()}
    descriptors = SiteDescriptors(
        pyramidalization_deg=g.get("pyramidalization_deg"),
        coordination=None if g.get("coordination") is None else int(g["coordination"]),
        n_neighbor_H=None if g.get("n_neighbor_H") is None else int(g["n_neighbor_H"]),
        distance_to_defect_nm=g.get("distance_to_defect_nm"),
        curvature_inv_nm=g.get("curvature_inv_nm"),
    )
    created = g.get("created")
    prov = Provenance(
        kind=ProvenanceKind(g["provenance_kind"]),
        source_id=str(g["source_id"]),
        identifier=str(g["identifier"]),
        extraction=None if g.get("extraction") is None else Extraction(g["extraction"]),
        settings_hash=g.get("settings_hash"),
        package_version=str(g.get("package_version") or ""),
        seed=None if g.get("seed") is None else int(g["seed"]),
        created=datetime.fromisoformat(created) if created else datetime.now().astimezone(),
    )
    return SiteEnergyRecord(
        config_hash=str(g["config_hash"]),
        tier=Tier(g["tier"]),
        engine=str(g["engine"]),
        quantity=Quantity(g["quantity"]),
        value=UncertainScalar(
            value=float(g["value_eV"]),
            sigma=float(g.get("sigma_eV") or 0.0),
            distribution=Distribution(g.get("distribution") or "normal"),
        ),
        reference_state=ReferenceState(g["reference_state"]),
        temperature_K=g.get("temperature_K"),
        converged=bool(g["converged"]),
        provenance=prov,
        site_class=None if g.get("site_class") is None else SiteClassName(g["site_class"]),
        site_id=g.get("site_id"),
        n_H=None if g.get("n_H") is None else int(g["n_H"]),
        descriptors=descriptors,
        n_atoms=None if g.get("n_atoms") is None else int(g["n_atoms"]),
        supercell=g.get("supercell"),
    )


def records_to_frame(records: list[SiteEnergyRecord]) -> pd.DataFrame:
    """Flatten records into a table with the columns of :data:`RECORD_COLUMNS`."""
    if not records:
        return pd.DataFrame(columns=list(RECORD_COLUMNS))
    return pd.DataFrame([record_to_row(r) for r in records], columns=list(RECORD_COLUMNS))


def frame_to_records(frame: pd.DataFrame) -> list[SiteEnergyRecord]:
    """Validate every row of a landscape table back into records."""
    return [row_to_record(dict(row)) for row in frame.to_dict(orient="records")]


__all__ = [
    "RECORD_COLUMNS",
    "DESCRIPTOR_COLUMNS",
    "record_to_row",
    "row_to_record",
    "records_to_frame",
    "frame_to_records",
]
