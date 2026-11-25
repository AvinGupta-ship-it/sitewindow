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


def run_B1_T1(reg: dict, tier: str = "verify") -> dict:
    """Tier-1 (GPAW, PBE) basal monomer versus the plane-wave PBE value of B1.

    The pass criterion applies to the converged Tier-1 protocol (5x5 cell, 500 eV, 4x4x1
    k-mesh, 15 A vacuum). A verification-tier anchor at reduced settings is reported without a
    pass flag.
    """
    import json

    candidates = [
        Path("data/processed/tier1/tier1_anchor_basal_5x5.json"),
        Path(f"data/processed/{tier}/tier1/tier1_anchor_basal_4x4.json"),
    ]
    thr = reg["B1"]["thresholds"]["T1_PBE_abs_eV"]
    ref = reg["B1"]["reference"]
    for p in candidates:
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        e_pbe = d["binding_eV"]["E_bind_H_atom_PBE"]
        e_d3 = d["binding_eV"]["E_bind_H_atom_PBE_D3"]
        s = d["settings"]
        converged_protocol = s["supercell"] >= 5 and s["ecut_eV"] >= 500 and s["vacuum_A"] >= 15
        dev = e_pbe - ref["PBE"]
        return {
            "id": "B1",
            "passed": bool(abs(dev) < thr) if converged_protocol else None,
            "predicted": {
                "E_bind_H_atom_PBE_eV": e_pbe,
                "E_bind_H_atom_PBE_D3_eV": e_d3,
                "supercell": s["supercell"],
                "ecut_eV": s["ecut_eV"],
                "kmesh": s["kmesh_used"],
                "vacuum_A": s["vacuum_A"],
            },
            "reference": {"PBE_eV": ref["PBE"], "DMC_eV": ref["DMC"], "threshold_eV": thr},
            "metric": "absolute difference of the PBE binding energy (tier T1)",
            "summary": (
                f"T1 PBE basal monomer {e_pbe:.3f} eV ({s['supercell']}x{s['supercell']}, "
                f"{s['ecut_eV']:.0f} eV, k{s['kmesh_used'][0]}x{s['kmesh_used'][1]}) vs PBE ref "
                f"{ref['PBE']} eV: deviation {dev:+.3f} eV"
                + (
                    ""
                    if converged_protocol
                    else " (reduced verification settings; not a pass test)"
                )
            ),
        }
    return {
        "id": "B1-T1",
        "passed": None,
        "summary": "Tier-1 anchor not present (run tier1 anchor)",
    }


def _variant_label(base: dict, settings: dict) -> str:
    """Name the settings a convergence variant changes relative to the base protocol."""
    parts = []
    if settings["supercell"] != base["supercell"]:
        n = settings["supercell"]
        parts.append(f"{n}x{n} cell")
    if list(settings["kmesh_used"]) != list(base["kmesh_used"]):
        k = settings["kmesh_used"]
        parts.append(f"k{k[0]}x{k[1]}x{k[2]}")
    if settings["ecut_eV"] != base["ecut_eV"]:
        parts.append(f"{settings['ecut_eV']:.0f} eV")
    if settings["vacuum_A"] != base["vacuum_A"]:
        parts.append(f"{settings['vacuum_A']:.0f} A vacuum")
    return ", ".join(parts) if parts else "unchanged settings"


def evaluate_convergence(series: dict, threshold_eV: float) -> dict:
    """Acceptance criterion B1 (convergence) applied to one Tier-1 convergence series.

    Every variant must change the PBE binding energy of the basal monomer by less than
    ``threshold_eV`` relative to the base protocol. The pass flag is set only when the base is
    the converged protocol (5x5 cell or larger, 500 eV, 15 A vacuum); a reduced verification
    series is reported without one.
    """
    variants = series["variants"]
    base = variants["base"]
    bs = base["settings"]
    converged_protocol = bs["supercell"] >= 5 and bs["ecut_eV"] >= 500 and bs["vacuum_A"] >= 15
    deltas = {
        name: float(v["E_bind_H_atom_PBE"] - base["E_bind_H_atom_PBE"])
        for name, v in variants.items()
        if name != "base"
    }
    labels = {name: _variant_label(bs, variants[name]["settings"]) for name in deltas}
    checks = {name: bool(abs(d) < threshold_eV) for name, d in deltas.items()}
    all_within = bool(deltas) and all(checks.values())
    k = bs["kmesh_used"]
    base_label = (
        f"{bs['supercell']}x{bs['supercell']}, {bs['ecut_eV']:.0f} eV, k{k[0]}x{k[1]}x{k[2]}, "
        f"{bs['vacuum_A']:.0f} A"
    )
    detail = "; ".join(f"{labels[n]} {deltas[n]:+.4f}" for n in deltas)
    return {
        "id": "B1",
        "passed": all_within if converged_protocol else None,
        "predicted": {
            "E_bind_H_atom_PBE_base_eV": base["E_bind_H_atom_PBE"],
            "delta_PBE_eV": deltas,
            "variant_settings": labels,
            "max_abs_delta_eV": max(abs(d) for d in deltas.values()) if deltas else None,
        },
        "reference": {"threshold_eV": threshold_eV, "base": base_label},
        "checks": checks,
        "metric": "absolute change of the PBE binding energy between the base and each variant (tier T1)",
        "summary": (
            f"T1 convergence of the basal monomer PBE energy from the base ({base_label}): {detail} eV; "
            + (
                f"every change below {threshold_eV:g} eV in magnitude"
                if all_within
                else f"not every change below {threshold_eV:g} eV in magnitude"
            )
            + ("" if converged_protocol else " (reduced verification settings; not a pass test)")
        ),
    }


def run_B1_T1_convergence(reg: dict, tier: str = "verify") -> dict:
    """Acceptance criterion B1 (convergence): the Tier-1 convergence series of the basal monomer."""
    import json

    candidates = [
        Path("data/processed/tier1/tier1_convergence_basal_5x5.json"),
        Path(f"data/processed/{tier}/tier1/tier1_convergence_basal_4x4.json"),
    ]
    thr = reg["B1"]["thresholds"]["T1_convergence_abs_eV"]
    for p in candidates:
        if p.exists():
            return evaluate_convergence(json.loads(p.read_text()), thr)
    return {
        "id": "B1",
        "passed": None,
        "summary": "Tier-1 convergence series not present (run tier1 converge)",
    }


def run_B3_pair_potential(reg: dict) -> dict:
    """Classical pair-potential minimum over coronene versus the CCSD(T)/CBS value (B3).

    The literature parameter sets are the prediction; the B3-calibrated set reproduces the
    reference by construction and is reported, not tested.
    """
    from sitewindow.atomistic import physisorption as ph

    ref = reg["B3"]["reference"]["value"]
    thr = reg["B3"]["thresholds"]["pair_potential_abs_kJ_mol"]
    predicted = {name: ph.coronene_binding_kJ_mol(p)[0] for name, p in ph.PARAMETER_SETS.items()}
    cal = ph.calibrated_potential()
    checks = {name: bool(abs(e - ref) < thr) for name, e in predicted.items()}
    return {
        "id": "B3",
        "passed": bool(all(checks.values())),
        "predicted": {
            **{f"{k}_kJ_mol": v for k, v in predicted.items()},
            "calibrated_epsilon_K": cal.epsilon_K,
        },
        "reference": {"CCSD_T_kJ_mol": ref, "threshold_kJ_mol": thr},
        "metric": "absolute difference of the classical minimum (single-site LJ, no orientation)",
        "checks": checks,
        "summary": "pair potential on coronene: "
        + "; ".join(f"{k} {v:.2f}" for k, v in predicted.items())
        + f" kJ/mol vs {ref} (tol {thr}); calibrated set eps/k={cal.epsilon_K:.1f} K reproduces B3 by construction",
    }


def run_B4_computed(reg: dict) -> dict:
    """Computed basal and slit-pore site energies versus the isosteric-heat interval of B4.

    Evaluated for every parameter set on the basal plane and the verification-tier pore
    widths; the pass flag refers to the B3-calibrated set used for the landscape.
    """
    from sitewindow.atomistic import physisorption as ph
    from sitewindow.atomistic.structures import builders as b
    from sitewindow.landscape import catalogue

    lo, hi = reg["B4"]["reference"]["interval"]
    widths = catalogue()["pore_slit"].recipe("verify")["widths_nm"]
    g = b.graphene_supercell(4)
    pores = {float(w): b.slit_pore(float(w), n=4) for w in widths}
    out: dict[str, dict[str, float]] = {}
    for name, p in ph.all_potentials().items():
        row = {"basal": -ph.sheet_site_energy(g, p)[0]}
        for w, pore in pores.items():
            row[f"slit_{w:.2f}nm"] = -ph.pore_site_energy(pore, p)[0]
        out[name] = row
    cal = ph.calibrated_potential().name
    checks = {name: bool(all(lo <= v <= hi for v in row.values())) for name, row in out.items()}
    return {
        "id": "B4",
        "passed": checks[cal],
        "predicted": out,
        "reference": {
            "interval_kJ_mol": [lo, hi],
            "AX21_average": reg["B4"]["reference"]["AX21_average"],
        },
        "metric": "site energies (minus sign) within the isosteric-heat interval",
        "checks": checks,
        "summary": "B3-calibrated potential: "
        + ", ".join(f"{k} {v:.1f}" for k, v in out[cal].items())
        + f" kJ/mol vs [{lo}, {hi}]; literature sets within interval: "
        + ", ".join(f"{k}={v}" for k, v in checks.items() if k != cal),
    }


def run_all(tier: str = "verify") -> dict:
    reg = _registry()
    results = [
        run_B8(reg),
        run_B9_endpoints(reg),
        run_B10(reg),
        run_B1_T2(reg),
        run_B1_T1(reg, tier),
        run_B1_T1_convergence(reg, tier),
        run_B3_pair_potential(reg),
        run_B4_computed(reg),
    ]
    return {"tier": tier, "benchmarks": results}


__all__ = [
    "run_all",
    "run_B8",
    "run_B9_endpoints",
    "run_B10",
    "run_B1_T2",
    "run_B1_T1",
    "run_B1_T1_convergence",
    "evaluate_convergence",
    "run_B3_pair_potential",
    "run_B4_computed",
]
