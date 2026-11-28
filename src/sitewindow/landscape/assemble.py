"""Assembly of the calibrated landscape (manual C3): ``landscape_calibrated.parquet`` and the
calibration report, from whatever Tier-1 anchors exist.

Every Tier-2 (and Tier-3) record of the landscape receives a calibrated mean and standard
deviation on the Tier-1 scale from the per-class transfer and discrepancy process of
:mod:`sitewindow.landscape.calibration`, and on the Tier-0 scale after the class-independent
offset estimated from the DMC and CCSD(T) anchors. Binding energies against the free atom are
calibrated directly; those against half H₂ follow by the Tier-1 atomization constant taken from
the anchor table; sequential energies receive the class transfer as a difference; defect formation
energies receive an offset-only transfer from the anchor of their class. Physisorption records of
the potential tier pass through with the B3 calibration status of their parameter set. The report
lists the calibration status of every class, the transfers, the discrepancy hyperparameters, the
leave-one-anchor-out coverage, the held-out-class errors, the synthetic recovery test, the Tier-0
offset and, when PyMC is available and requested, the meta-analysis of the published values.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from sitewindow import __version__
from sitewindow.io import sha256_file, write_table
from sitewindow.landscape.calibration import (
    DESCRIPTORS,
    T1_ENGINE,
    TIER_SIGMA_EV,
    Calibration,
    anchor_pairs,
    configuration_name,
    held_out_class_errors,
    loo_coverage,
    synthetic_recovery,
    tier0_offset,
)
from sitewindow.landscape.registry import REGISTRY, published_values, registry_frame

log = logging.getLogger(__name__)

CALIBRATED_COLUMNS = (
    "calibrated_mean_T1_eV",
    "calibrated_sd_T1_eV",
    "calibrated_mean_T0_eV",
    "calibrated_sd_T0_eV",
    "calibration_status",
    "calibration_note",
    "n_anchors_class",
    "transfer_a_eV",
    "transfer_b",
    "published_mean_eV",
    "published_sd_eV",
    "n_published",
)
CHEMISORPTION_QUANTITIES = ("E_bind_H_atom", "E_bind_H_half_H2", "E_sequential")


def atomization_constant(anchors: pd.DataFrame, engine: str = T1_ENGINE) -> float | None:
    """E_bind(half H₂) − E_bind(H atom) at Tier 1, the same number for every configuration
    (the Tier-1 energy of H minus half that of H₂), read from the anchor table."""
    a = anchors[anchors["engine"] == engine]
    if a.empty:
        return None
    atom = a[(a["quantity"] == "E_bind_H_atom")].set_index("identifier")["value_eV"]
    half = a[(a["quantity"] == "E_bind_H_half_H2")].set_index("identifier")["value_eV"]
    common = atom.index.intersection(half.index)
    if len(common) == 0:
        return None
    return float((half.loc[common] - atom.loc[common]).mean())


def companion_tables(processed_dir: Path, tier: str, sub: str | None = None) -> list[pd.DataFrame]:
    """Landscape tables of the other tiers of the same calculator (the verify set beside the full
    set and the reverse), read from the sibling directories of ``processed_dir``; they supply the
    tier energy of anchor configurations the primary table lacks, such as the (5,5) tube."""
    out: list[pd.DataFrame] = []
    for other in ("verify", "full"):
        if other == tier:
            continue
        d = processed_dir.parent / other
        if sub is not None:
            d = d / sub
        p = d / f"landscape_raw_{other}.parquet"
        if p.exists():
            out.append(pd.read_parquet(p))
    return out


def _empty_calibration(frame: pd.DataFrame, status: str, note: str) -> pd.DataFrame:
    out = frame.copy()
    out["calibrated_mean_T1_eV"] = np.nan
    out["calibrated_sd_T1_eV"] = np.nan
    out["calibrated_mean_T0_eV"] = np.nan
    out["calibrated_sd_T0_eV"] = np.nan
    out["calibration_status"] = status
    out["calibration_note"] = note
    out["n_anchors_class"] = 0
    out["transfer_a_eV"] = np.nan
    out["transfer_b"] = np.nan
    return out


def calibrate_table(
    table: pd.DataFrame,
    anchors: pd.DataFrame,
    tier: str,
    offset: dict[str, Any],
    kernel: str = "matern52",
    sigma_t: float | None = None,
    companions: list[pd.DataFrame] | None = None,
) -> tuple[pd.DataFrame, Calibration, pd.DataFrame]:
    """Calibrated copy of one landscape table, the fitted calibration and its anchor pairs; the
    ``companions`` (other landscape tables of the same tier, such as the verify set) supply the
    tier energy of an anchor configuration the table itself lacks."""
    sigma_t = TIER_SIGMA_EV.get(tier, 0.05) if sigma_t is None else sigma_t
    pairs = anchor_pairs(anchors, [table, *(companions or [])], "E_bind_H_atom", "H_atom")
    cal = Calibration(tier, kernel, sigma_t=sigma_t).fit(pairs)
    parts: list[pd.DataFrame] = []

    atom = table[(table["quantity"] == "E_bind_H_atom") & (table["reference_state"] == "H_atom")]
    pred = cal.predict(atom)
    pred = pred.rename(
        columns={
            "calibrated_mean_eV": "calibrated_mean_T1_eV",
            "calibrated_sd_eV": "calibrated_sd_T1_eV",
        }
    )
    pred["calibration_note"] = "class transfer with discrepancy process"
    parts.append(pred)

    const = atomization_constant(anchors)
    half = table[
        (table["quantity"] == "E_bind_H_half_H2") & (table["reference_state"] == "half_H2")
    ]
    if const is None:
        parts.append(
            _empty_calibration(
                half, "not_calibrated", "no Tier-1 atomization constant: anchors absent"
            )
        )
    else:
        base = pred.set_index("identifier")
        h = half.copy()
        idx = h["identifier"]
        h["calibrated_mean_T1_eV"] = idx.map(base["calibrated_mean_T1_eV"]).to_numpy() + const
        h["calibrated_sd_T1_eV"] = idx.map(base["calibrated_sd_T1_eV"]).to_numpy()
        h["calibration_status"] = idx.map(base["calibration_status"]).to_numpy()
        h["n_anchors_class"] = idx.map(base["n_anchors_class"]).to_numpy()
        h["transfer_a_eV"] = idx.map(base["transfer_a_eV"]).to_numpy()
        h["transfer_b"] = idx.map(base["transfer_b"]).to_numpy()
        h["calibration_note"] = (
            "atom-referenced calibration shifted by the Tier-1 atomization constant"
        )
        parts.append(h)

    seq_atom = table[(table["quantity"] == "E_sequential") & (table["reference_state"] == "H_atom")]
    seq_half = table[
        (table["quantity"] == "E_sequential") & (table["reference_state"] == "half_H2")
    ]
    if not seq_atom.empty:
        s = cal.predict(seq_atom)
        s = s.rename(
            columns={
                "calibrated_mean_eV": "calibrated_mean_T1_eV",
                "calibrated_sd_eV": "calibrated_sd_T1_eV",
            }
        )
        s["calibrated_sd_T1_eV"] = s["calibrated_sd_T1_eV"] * np.sqrt(2.0)
        s["calibration_note"] = (
            "class transfer applied to the energy difference; uncertainty of two configurations"
        )
        parts.append(s)
        if not seq_half.empty:
            if const is None:
                parts.append(
                    _empty_calibration(
                        seq_half, "not_calibrated", "no Tier-1 atomization constant: anchors absent"
                    )
                )
            else:
                base = s.set_index("identifier")
                h = seq_half.copy()
                idx = h["identifier"]
                h["calibrated_mean_T1_eV"] = (
                    idx.map(base["calibrated_mean_T1_eV"]).to_numpy() + const
                )
                h["calibrated_sd_T1_eV"] = idx.map(base["calibrated_sd_T1_eV"]).to_numpy()
                h["calibration_status"] = idx.map(base["calibration_status"]).to_numpy()
                h["n_anchors_class"] = idx.map(base["n_anchors_class"]).to_numpy()
                h["transfer_a_eV"] = idx.map(base["transfer_a_eV"]).to_numpy()
                h["transfer_b"] = idx.map(base["transfer_b"]).to_numpy()
                h["calibration_note"] = (
                    "atom-referenced sequential calibration shifted by the Tier-1 atomization constant"
                )
                parts.append(h)
    elif not seq_half.empty:
        parts.append(
            _empty_calibration(
                seq_half,
                "not_calibrated",
                "no atom-referenced sequential energies to calibrate from",
            )
        )

    form = table[table["quantity"] == "E_formation_defect"]
    if not form.empty:
        f = form.copy()
        f["calibrated_mean_T1_eV"] = np.nan
        f["calibrated_sd_T1_eV"] = np.nan
        f["calibration_status"] = "awaiting_anchors"
        f["calibration_note"] = "offset from the Tier-1 formation energy of the class"
        f["n_anchors_class"] = 0
        f["transfer_a_eV"] = np.nan
        f["transfer_b"] = np.nan
        a_form = anchors[
            (anchors["quantity"] == "E_formation_defect") & (anchors["engine"] == T1_ENGINE)
        ]
        for cls in f["site_class"].unique():
            sel = a_form[a_form["site_class"] == cls]
            rows = f["site_class"] == cls
            if sel.empty:
                f.loc[rows, "calibrated_mean_T1_eV"] = f.loc[rows, "value_eV"]
                f.loc[rows, "calibrated_sd_T1_eV"] = np.hypot(0.5, sigma_t)
                f.loc[rows, "calibration_note"] = (
                    "no Tier-1 formation energy for this class; tier value carried with a 0.5 eV prior width"
                )
                continue
            e_t1 = float(sel["value_eV"].iloc[0])
            e_t = float(f.loc[rows, "value_eV"].iloc[0])
            f.loc[rows, "calibrated_mean_T1_eV"] = e_t1
            f.loc[rows, "calibrated_sd_T1_eV"] = np.hypot(TIER_SIGMA_EV["T1"], sigma_t)
            f.loc[rows, "calibration_status"] = "offset_calibrated"
            f.loc[rows, "n_anchors_class"] = 1
            f.loc[rows, "transfer_a_eV"] = e_t1 - e_t
            f.loc[rows, "transfer_b"] = 1.0
            f.loc[rows, "calibration_note"] = (
                f"Tier-1 formation energy of the anchor cell ({sel['supercell'].iloc[0]}) replaces the tier value; cell sizes differ"
            )
        parts.append(f)

    covered = {"E_bind_H_atom", "E_bind_H_half_H2", "E_sequential", "E_formation_defect"}
    rest = table[~table["quantity"].isin(covered)]
    if not rest.empty:
        parts.append(_empty_calibration(rest, "not_calibrated", "quantity outside the calibration"))

    out = pd.concat(parts, ignore_index=True)
    chem = out["quantity"].isin(CHEMISORPTION_QUANTITIES).to_numpy()
    out["calibrated_mean_T0_eV"] = out["calibrated_mean_T1_eV"]
    out["calibrated_sd_T0_eV"] = out["calibrated_sd_T1_eV"]
    out.loc[chem, "calibrated_mean_T0_eV"] = (
        out.loc[chem, "calibrated_mean_T1_eV"] - offset["offset_eV"]
    )
    out.loc[chem, "calibrated_sd_T0_eV"] = np.hypot(
        out.loc[chem, "calibrated_sd_T1_eV"], offset["sigma_carried_eV"]
    )
    out["calibration_tier_source"] = tier
    return out, cal, pairs


REGISTRY_SITE_TO_CONFIGURATION: dict[str, str] = {
    "monomer": "monomer",
    "ortho": "ortho",
    "para": "para",
    "meta": "meta",
    "graphane": "graphane",
    "dangling_bond": "seq1",
    "rotated_bond": "monomer_rotated_bond",
    "rotated_pair": "dimer_rotated_bond",
    "edge_CH_bare": "edge_CH",
    "edge_CH_average": "edge_CH",
    "exterior_top": "monomer_exterior",
}
"""Site labels of the published-value registry mapped to the catalogue configuration they
describe, for joining the meta-analysis site means to the landscape rows."""


def attach_published(out: pd.DataFrame, meta_results: list[Any]) -> pd.DataFrame:
    """Posterior site means of the meta-analysis joined to the rows of the same quantity,
    reference state, class and configuration (several registry sites may map to one
    configuration; their values are pooled by count-weighted mean and the larger spread kept)."""
    out = out.copy()
    out["published_mean_eV"] = np.nan
    out["published_sd_eV"] = np.nan
    out["n_published"] = 0
    configuration = out["identifier"].map(configuration_name)
    for res in meta_results:
        sm = res.site_means.copy()
        sm["configuration"] = sm["site_id"].map(REGISTRY_SITE_TO_CONFIGURATION)
        sm = sm[sm["configuration"].notna()]
        if sm.empty:
            continue
        grouped = sm.groupby(["site_class", "configuration"]).apply(
            lambda g: pd.Series(
                {
                    "mean_eV": float(np.average(g["mean_eV"], weights=g["n_values"])),
                    "sd_eV": float(
                        np.sqrt(
                            np.max(
                                g["sd_eV"] ** 2
                                + (g["mean_eV"] - np.average(g["mean_eV"], weights=g["n_values"]))
                                ** 2
                            )
                        )
                    ),
                    "n_values": int(g["n_values"].sum()),
                }
            ),
            include_groups=False,
        )
        mask = (out["quantity"] == res.quantity) & (out["reference_state"] == res.reference_state)
        keys = list(zip(out.loc[mask, "site_class"], configuration[mask], strict=True))
        means = [grouped["mean_eV"].get(k, np.nan) for k in keys]
        sds = [grouped["sd_eV"].get(k, np.nan) for k in keys]
        ns = [int(grouped["n_values"].get(k, 0)) for k in keys]
        out.loc[mask, "published_mean_eV"] = means
        out.loc[mask, "published_sd_eV"] = sds
        out.loc[mask, "n_published"] = ns
    return out


def class_status_table(out: pd.DataFrame, cal: Calibration, registry: pd.DataFrame) -> pd.DataFrame:
    rows = []
    atom = out[(out["quantity"] == "E_bind_H_atom") & (out["reference_state"] == "H_atom")]
    for cls in sorted(atom["site_class"].unique()):
        sub = atom[atom["site_class"] == cls]
        n_pub = int((registry["site_class"] == cls).sum()) if not registry.empty else 0
        rows.append(
            {
                "site_class": cls,
                "status": cal.status(str(cls)),
                "n_anchors": cal.n_anchors(str(cls)),
                "n_configurations": int(len(sub)),
                "mean_shift_T1_minus_tier_eV": float(
                    (sub["calibrated_mean_T1_eV"] - sub["value_eV"]).mean()
                ),
                "mean_sd_T1_eV": float(sub["calibrated_sd_T1_eV"].mean()),
                "n_published": n_pub,
            }
        )
    return pd.DataFrame(rows)


def assemble_calibrated(
    processed_dir: Path,
    tier: str = "full",
    anchors_path: Path = Path("data/processed/tier1/anchors_tier1.parquet"),
    registry_path: Path = REGISTRY,
    out_dir: Path | None = None,
    tier3_dir: Path | None = None,
    meta: bool = True,
    meta_draws: int = 1000,
    meta_chains: int = 4,
    recovery_replicates: int = 100,
    kernel: str = "matern52",
    sigma_t2: float | None = None,
    seed: int = 0,
) -> dict[str, Any]:
    """Build ``landscape_calibrated.parquet`` and ``calibration_report.json`` under ``out_dir``
    (default: ``processed_dir``)."""
    t0 = time.perf_counter()
    processed_dir = Path(processed_dir)
    out_dir = Path(out_dir) if out_dir is not None else processed_dir
    t2_path = processed_dir / f"landscape_raw_{tier}.parquet"
    if not t2_path.exists():
        raise FileNotFoundError(f"no landscape table at {t2_path}")
    t2 = pd.read_parquet(t2_path)
    anchors = (
        pd.read_parquet(anchors_path)
        if Path(anchors_path).exists()
        else pd.DataFrame(columns=t2.columns)
    )
    values = published_values(registry_path)
    registry = registry_frame(values)
    offset = (
        tier0_offset(anchors, registry)
        if not registry.empty
        else tier0_offset(anchors, pd.DataFrame(columns=["tier0"]))
    )

    companions2 = companion_tables(processed_dir, tier)
    out2, cal2, pairs2 = calibrate_table(t2, anchors, "T2", offset, kernel, sigma_t2, companions2)
    tables = [out2]
    report: dict[str, Any] = {
        "tier": tier,
        "package_version": __version__,
        "anchors": {
            "path": str(anchors_path),
            "n_anchor_configurations": int(
                anchors[anchors["quantity"] == "E_bind_H_atom"]["identifier"].nunique()
            )
            if not anchors.empty
            else 0,
        },
        "tier0_offset": offset,
        "T2": {
            "calibration": cal2.report(),
            "loo": loo_coverage(pairs2, "T2", kernel=kernel, sigma_t=cal2.sigma_t)
            if len(pairs2) >= 2
            else {"n": int(len(pairs2)), "note": "fewer than two anchor pairs"},
            "held_out_classes": held_out_class_errors(pairs2, "T2", kernel=kernel)
            if len(pairs2) >= 2
            else {},
            "anchor_pairs": pairs2.drop(columns=list(DESCRIPTORS), errors="ignore").to_dict(
                orient="records"
            ),
        },
    }
    t3_dir = Path(tier3_dir) if tier3_dir is not None else processed_dir / "tier3"
    t3_path = t3_dir / f"landscape_raw_{tier}.parquet"
    if t3_path.exists():
        t3 = pd.read_parquet(t3_path)
        companions3 = companion_tables(processed_dir, tier, "tier3")
        out3, cal3, pairs3 = calibrate_table(t3, anchors, "T3", offset, kernel, None, companions3)
        tables.append(out3)
        report["T3"] = {
            "calibration": cal3.report(),
            "loo": loo_coverage(pairs3, "T3", kernel=kernel, sigma_t=cal3.sigma_t)
            if len(pairs3) >= 2
            else {"n": int(len(pairs3)), "note": "fewer than two anchor pairs"},
            "held_out_classes": held_out_class_errors(pairs3, "T3", kernel=kernel)
            if len(pairs3) >= 2
            else {},
            "anchor_pairs": pairs3.drop(columns=list(DESCRIPTORS), errors="ignore").to_dict(
                orient="records"
            ),
        }
    phys_path = processed_dir / f"physisorption_{tier}.parquet"
    if phys_path.exists():
        phys = pd.read_parquet(phys_path)
        p = _empty_calibration(
            phys,
            "potential_B3_calibrated",
            "pairwise potential; the B3-calibrated parameter set reproduces the CCSD(T) coronene anchor by construction (manual C2)",
        )
        p["calibrated_mean_T1_eV"] = p["value_eV"]
        p["calibrated_sd_T1_eV"] = p["sigma_eV"]
        p["calibrated_mean_T0_eV"] = p["value_eV"]
        p["calibrated_sd_T0_eV"] = p["sigma_eV"]
        p["calibration_tier_source"] = "potential"
        tables.append(p)

    meta_results: list[Any] = []
    meta_report: list[dict[str, Any]] = []
    if meta and not registry.empty:
        try:
            from sitewindow.landscape.meta_analysis import (
                MetaAnalysisUnavailable,
                fit_meta_analysis,
            )

            for (q, r), group in registry.groupby(["quantity", "reference_state"]):
                group = group[group["site_class"].notna()]
                if len(group) < 3:
                    meta_report.append(
                        {
                            "quantity": q,
                            "reference_state": r,
                            "n_values": int(len(group)),
                            "note": "fewer than three values; not fitted",
                        }
                    )
                    continue
                res = fit_meta_analysis(group, draws=meta_draws, chains=meta_chains, seed=seed)
                meta_results.append(res)
                meta_report.append(res.as_dict())
        except MetaAnalysisUnavailable as exc:
            meta_report.append({"note": str(exc)})
    report["meta_analysis"] = meta_report
    report["registry"] = {
        "path": str(registry_path),
        "n_site_energies": int(len(registry)),
        "by_class": registry.groupby("site_class").size().to_dict() if not registry.empty else {},
    }
    if recovery_replicates > 0:
        report["synthetic_recovery_calibration"] = synthetic_recovery(
            recovery_replicates, seed=seed
        )

    out = pd.concat(tables, ignore_index=True)
    out = attach_published(out, meta_results)
    status = class_status_table(out2, cal2, registry)
    report["class_status"] = status.to_dict(orient="records")
    report["wall_s"] = time.perf_counter() - t0

    out_dir.mkdir(parents=True, exist_ok=True)
    p_table = write_table(
        out, out_dir / "landscape_calibrated.parquet", {"tier": tier, "calibration": "C3"}
    )
    p_report = out_dir / "calibration_report.json"
    p_report.write_text(json.dumps(report, indent=1, default=_json_default) + "\n")
    report["files"] = {p_table.name: sha256_file(p_table), p_report.name: sha256_file(p_report)}
    p_report.write_text(json.dumps(report, indent=1, default=_json_default) + "\n")
    log.info("calibrated landscape: %d records, %d classes", len(out), status.shape[0])
    return {"table": out, "report": report, "status": status, "calibration_T2": cal2}


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.bool_):
        return bool(obj)
    return str(obj)


__all__ = [
    "CALIBRATED_COLUMNS",
    "assemble_calibrated",
    "atomization_constant",
    "attach_published",
    "calibrate_table",
    "class_status_table",
    "configuration_name",
]
