"""Prior and kernel sensitivity of the C3 meta-analysis and calibration (manual C3,
Convergence and sensitivity).

Prior sensitivity: every registry group (quantity and reference state) is refitted with the
widths of all priors multiplied by the scales 0.5, 1 and 2; the site means, method effects and
heterogeneities are compared across scales, and the largest shift of any site mean is reported
against its posterior standard deviation at the reference scale. Kernel sensitivity: the
calibration of each tier is refitted with the Matérn-3/2 kernel in place of the Matérn-5/2
kernel; class offsets, predictive standard deviations, leave-one-anchor-out coverage and the
held-out-class errors are compared. Both are diagnostics of the reported model, not model
selection: the reference scale 1 and the Matérn-5/2 kernel remain the reported configuration.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from sitewindow.landscape.assemble import companion_tables
from sitewindow.landscape.calibration import (
    Calibration,
    anchor_pairs,
    held_out_class_errors,
    loo_coverage,
)
from sitewindow.landscape.registry import published_values, registry_frame

log = logging.getLogger(__name__)

PRIOR_SCALES = (0.5, 1.0, 2.0)
KERNELS = ("matern52", "matern32")


def prior_sensitivity(
    registry: pd.DataFrame,
    scales: tuple[float, ...] = PRIOR_SCALES,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    seed: int = 0,
    cores: int | None = None,
) -> dict[str, Any]:
    """Refit every registry group at each prior scale and compare the posteriors."""
    from sitewindow.landscape.meta_analysis import fit_meta_analysis

    groups: list[dict[str, Any]] = []
    for (q, r), group in registry.groupby(["quantity", "reference_state"]):
        group = group[group["site_class"].notna()]
        if len(group) < 3:
            continue
        fits = {}
        for s in scales:
            t0 = time.perf_counter()
            res = fit_meta_analysis(
                group, draws=draws, tune=tune, chains=chains, seed=seed, prior_scale=s, cores=cores
            )
            fits[s] = res
            log.info("prior scale %.1f %s/%s: %d values, rhat %.3f, %.0f s", s, q, r, len(group),
                     res.diagnostics.get("rhat_max", float("nan")), time.perf_counter() - t0)
        ref = fits[1.0]
        site_rows = []
        for _, row in ref.site_means.iterrows():
            entry = {"site": row["site"], "site_class": row["site_class"], "mean_scale_1": row["mean_eV"], "sd_scale_1": row["sd_eV"]}
            shifts = []
            for s in scales:
                sm = fits[s].site_means
                m = sm.loc[sm["site"] == row["site"], "mean_eV"]
                if len(m):
                    entry[f"mean_scale_{s:g}"] = float(m.iloc[0])
                    shifts.append(abs(float(m.iloc[0]) - row["mean_eV"]))
            entry["max_shift_eV"] = max(shifts) if shifts else 0.0
            entry["max_shift_over_sd"] = entry["max_shift_eV"] / row["sd_eV"] if row["sd_eV"] > 0 else None
            site_rows.append(entry)
        method_rows = []

        def _keys(frame: pd.DataFrame) -> pd.Series:
            level = frame["level"].astype(str) if "level" in frame else pd.Series("", index=frame.index)
            return frame["effect"].astype(str) + ":" + level.fillna("")

        ref_keys = _keys(ref.method_effects)
        for i, row in ref.method_effects.iterrows():
            entry = {"effect": ref_keys.loc[i], "mean_scale_1": row["mean_eV"], "sd_scale_1": row["sd_eV"]}
            for s in scales:
                me = fits[s].method_effects
                m = me.loc[_keys(me) == ref_keys.loc[i], "mean_eV"]
                if len(m):
                    entry[f"mean_scale_{s:g}"] = float(m.iloc[0])
            method_rows.append(entry)
        shifts = [r_["max_shift_eV"] for r_ in site_rows]
        groups.append(
            {
                "quantity": q,
                "reference_state": r,
                "n_values": int(len(group)),
                "heterogeneity": {f"{s:g}": fits[s].heterogeneity for s in scales},
                "diagnostics": {f"{s:g}": fits[s].diagnostics for s in scales},
                "site_means": site_rows,
                "method_effects": method_rows,
                "max_site_shift_eV": float(max(shifts)) if shifts else 0.0,
                "max_site_shift_over_sd": float(max(r_["max_shift_over_sd"] or 0.0 for r_ in site_rows)) if site_rows else 0.0,
            }
        )
    return {
        "scales": list(scales),
        "sampler": {"draws": draws, "tune": tune, "chains": chains},
        "groups": groups,
        "max_site_shift_eV": float(max((g["max_site_shift_eV"] for g in groups), default=0.0)),
        "max_site_shift_over_sd": float(max((g["max_site_shift_over_sd"] for g in groups), default=0.0)),
    }


def _tier_pairs(processed_dir: Path, tier: str, anchors: pd.DataFrame, sub: str | None) -> pd.DataFrame | None:
    name = f"landscape_raw_{tier}.parquet"
    path = processed_dir / sub / name if sub else processed_dir / name
    if not path.exists():
        return None
    table = pd.read_parquet(path)
    companions = companion_tables(processed_dir, tier, sub)
    return anchor_pairs(anchors, [table, *companions])


def kernel_sensitivity(
    processed_dir: Path,
    anchors_path: Path,
    tier: str = "full",
    kernels: tuple[str, ...] = KERNELS,
) -> dict[str, Any]:
    """Refit the calibration of each tier with every kernel and compare the results."""
    processed_dir = Path(processed_dir)
    if not Path(anchors_path).exists():
        return {"note": "no anchor table"}
    anchors = pd.read_parquet(anchors_path)
    out: dict[str, Any] = {"kernels": list(kernels), "tiers": {}}
    for label, sub in (("T2", None), ("T3", "tier3")):
        pairs = _tier_pairs(processed_dir, tier, anchors, sub)
        if pairs is None or len(pairs) < 2:
            out["tiers"][label] = {"note": "fewer than two anchor pairs"}
            continue
        per_kernel = {}
        for k in kernels:
            cal = Calibration(label, k).fit(pairs)
            rep = cal.report()
            pred = cal.predict(pairs.rename(columns={"E_t": "value_eV"}))
            per_kernel[k] = {
                "report": rep,
                "anchor_predictive_sd_eV": pred["calibrated_sd_eV"].astype(float).tolist(),
                "loo": loo_coverage(pairs, label, kernel=k),
                "held_out": held_out_class_errors(pairs, label, kernel=k),
            }
        ref, alt = per_kernel[kernels[0]], per_kernel[kernels[-1]]
        offs_ref = {c: v.get("class_offset_eV") for c, v in ref["report"].get("classes", {}).items()}
        offs_alt = {c: v.get("class_offset_eV") for c, v in alt["report"].get("classes", {}).items()}
        offset_shift = {c: (offs_alt.get(c) or 0.0) - (offs_ref.get(c) or 0.0) for c in offs_ref}
        sd_ref = np.asarray(ref["anchor_predictive_sd_eV"])
        sd_alt = np.asarray(alt["anchor_predictive_sd_eV"])
        out["tiers"][label] = {
            "n_anchor_pairs": int(len(pairs)),
            "per_kernel": per_kernel,
            "class_offset_shift_eV": offset_shift,
            "max_class_offset_shift_eV": float(max((abs(v) for v in offset_shift.values()), default=0.0)),
            "predictive_sd_ratio_mean": float(np.mean(sd_alt / sd_ref)) if len(sd_ref) else None,
            "loo_coverage": {k: per_kernel[k]["loo"].get("coverage") for k in kernels},
            "loo_rmse_eV": {k: per_kernel[k]["loo"].get("rmse_eV") for k in kernels},
        }
    return out


def run_sensitivity(
    processed_dir: Path,
    anchors_path: Path,
    registry_path: Path,
    out_path: Path,
    *,
    tier: str = "full",
    prior: bool = True,
    kernel: bool = True,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    seed: int = 0,
    cores: int | None = None,
) -> dict[str, Any]:
    """Run the requested studies and write ``calibration_sensitivity.json``."""
    t0 = time.perf_counter()
    out_path = Path(out_path)
    report: dict[str, Any] = json.loads(out_path.read_text()) if out_path.exists() else {}
    if kernel:
        report["kernel"] = kernel_sensitivity(processed_dir, anchors_path, tier)
    if prior:
        registry = registry_frame(published_values(registry_path))
        report["prior"] = prior_sensitivity(
            registry, draws=draws, tune=tune, chains=chains, seed=seed, cores=cores
        )
    report["wall_s"] = time.perf_counter() - t0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, default=float))
    return report


__all__ = ["KERNELS", "PRIOR_SCALES", "kernel_sensitivity", "prior_sensitivity", "run_sensitivity"]
