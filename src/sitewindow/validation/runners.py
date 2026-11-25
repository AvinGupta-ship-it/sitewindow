"""Benchmark runners. Thresholds are read from ``data/registry/benchmarks.yaml``.

Each runner returns a dict with ``id``, ``passed`` (True/False/None when not executable at the
requested tier), ``predicted``, ``reference``, ``metric`` and ``summary``.
"""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

import numpy as np
import yaml

from sitewindow.constants import EV_TO_KJ_PER_MOL


def _registry() -> dict:
    p = Path("data/registry/benchmarks.yaml")
    if not p.exists():
        p = Path(str(files("sitewindow").joinpath("data/benchmarks.yaml")))
    return yaml.safe_load(p.read_text())


def run_B8(reg: dict) -> dict:
    from sitewindow.theory.entropy import ConstantEntropy
    from sitewindow.theory.window import BHATIA_MYERS, OperatingWindow, optimum_energy, swing

    thr = reg["B8"]["thresholds"]
    eps_star, _ = optimum_energy(BHATIA_MYERS, ConstantEntropy())
    dH = eps_star * EV_TO_KJ_PER_MOL
    Ts = np.linspace(60.0, 300.0, 2401)
    eps = -5.8 / EV_TO_KJ_PER_MOL
    vals = [
        float(swing(eps, OperatingWindow("t", T, 30.0, T, 1.5, fugacity=False), ConstantEntropy()))
        for T in Ts
    ]
    T_opt = float(Ts[int(np.argmax(vals))])
    ok = abs(dH - (-15.1)) < thr["dH_kJ_mol"] and abs(T_opt - 115.0) < thr["T_opt_K"]
    return {
        "id": "B8",
        "passed": bool(ok),
        "predicted": {"dH_opt_kJ_mol": dH, "T_opt_K": T_opt},
        "reference": {"dH_opt_kJ_mol": -15.1, "T_opt_K": 115.0},
        "metric": "absolute difference",
        "summary": f"dH*={dH:.2f} kJ/mol (ref -15.1), T_opt={T_opt:.0f} K (ref 115)",
    }


def run_B9_endpoints(reg: dict) -> dict:
    from scipy.optimize import brentq

    from sitewindow.defects.knockon import (
        debye_mean_square_velocity,
        t_max_eV,
        thermally_smeared_cross_section_barn,
    )

    thr = reg["B9"]["thresholds"]
    v2 = debye_mean_square_velocity(1287.0, 300.0)
    s100 = thermally_smeared_cross_section_barn(100.0, 22.0, v2)
    s80 = thermally_smeared_cross_section_barn(80.0, 22.0, v2)
    onset = float(brentq(lambda E: float(t_max_eV(E)) - 22.0, 50.0, 200.0))
    v2_13 = debye_mean_square_velocity(1287.0 * np.sqrt(12 / 13), 300.0, 13.003355)
    s100_13 = thermally_smeared_cross_section_barn(100.0, 22.0, v2_13, mass_amu=13.003355)
    checks = {
        "onset": abs(onset - 108.0) < thr["onset_keV"],
        "100keV_within_factor": 1 / thr["factor_100keV"] < s100 / 0.2 < thr["factor_100keV"],
        "80keV_within_decade": 0.1 < s80 / 1e-4 < 10.0,
        "13C_below_12C": s100_13 < s100,
    }
    return {
        "id": "B9",
        "passed": bool(all(checks.values())),
        "predicted": {
            "onset_keV": onset,
            "sigma_100keV_barn": s100,
            "sigma_80keV_barn": s80,
            "sigma_100keV_13C_barn": s100_13,
        },
        "reference": {"onset_keV": 108.0, "sigma_100keV_barn": 0.2, "sigma_80keV_barn": 1e-4},
        "metric": "endpoint ratios; onset",
        "checks": checks,
        "summary": f"onset {onset:.1f} keV; sigma(100 keV)={s100:.3f} barn (ref ~0.2); sigma(80 keV)={s80:.1e} barn (ref ~1e-4); checks={checks}",
    }


def run_B10(reg: dict) -> dict:
    from sitewindow.defects.raman import LD_low_defect_regime_nm, id_ig_ratio

    thr = reg["B10"]["thresholds"]
    L = np.array([20.0, 30.0, 50.0])
    r = id_ig_ratio(L, 2.41)
    rel = np.abs(LD_low_defect_regime_nm(r, 2.41) / L - 1.0)
    Lmax = np.linspace(1.0, 20.0, 4000)[
        int(np.argmax(id_ig_ratio(np.linspace(1.0, 20.0, 4000), 2.41)))
    ]
    ok = bool(np.all(rel < thr["relative"]) and 2.0 < Lmax < 4.0)
    return {
        "id": "B10",
        "passed": ok,
        "predicted": {"max_rel_dev_low_defect": float(rel.max()), "L_D_at_max_nm": float(Lmax)},
        "reference": {
            "low_defect_relation": "L_D^2 = 4.3e3 E_L^-4 / (I_D/I_G)",
            "L_D_at_max_nm": "~3",
        },
        "metric": "relative deviation",
        "summary": f"asymptote deviation {rel.max():.3f}; I_D/I_G maximum at L_D={Lmax:.2f} nm",
    }


def run_B1_T2(reg: dict) -> dict:
    """Tier-2 (GFN2-xTB) basal monomer versus B1; reported, not expected to pass without calibration."""
    p = Path("data/processed/verify/atomistic_verify_T2.json")
    if not p.exists():
        return {
            "id": "B1",
            "passed": None,
            "summary": "atomistic verification results not present (run verify-atomistic)",
        }
    import json

    d = json.loads(p.read_text())
    e = d["basal_monomer"]["E_bind_H_atom"]
    return {
        "id": "B1",
        "passed": None,
        "predicted": {"E_bind_H_atom_T2_eV": e},
        "reference": {"DMC_eV": -0.691, "PBE_eV": -0.821},
        "metric": "absolute difference (tier T2, uncalibrated)",
        "summary": f"T2 basal monomer {e:.3f} eV vs DMC -0.691 / PBE -0.821 eV: tier offset {e + 0.691:+.2f} eV (input to C3 calibration)",
    }


def run_all(tier: str = "verify") -> dict:
    reg = _registry()
    results = [run_B8(reg), run_B9_endpoints(reg), run_B10(reg), run_B1_T2(reg)]
    return {"tier": tier, "benchmarks": results}


__all__ = ["run_all", "run_B8", "run_B9_endpoints", "run_B10", "run_B1_T2"]
