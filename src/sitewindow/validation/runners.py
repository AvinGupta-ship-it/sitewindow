"""Benchmark runners. Thresholds are read from ``data/registry/benchmarks.yaml``.

Each runner returns a dict with ``id``, ``passed`` (True/False/None when not executable at the
requested tier), ``predicted``, ``reference``, ``metric`` and ``summary``.
"""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import yaml

from sitewindow.constants import EV_TO_KJ_PER_MOL

if TYPE_CHECKING:
    import pandas as pd


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


def evaluate_B1_T2(
    frame: pd.DataFrame, reg: dict, settings: dict | None = None, tier_label: str = "T2"
) -> dict:
    """Basal monomer of a Tier-2 or Tier-3 landscape table against B1; reported, not expected to
    pass without calibration. The summary names the cell and, for tblite, the electronic
    temperature because the Gamma-only sampling makes the value depend on both."""
    ref = reg["B1"]["reference"]
    rows = frame[
        (frame["site_class"] == "basal")
        & (frame["quantity"] == "E_bind_H_atom")
        & (frame["n_H"] == 1)
    ]
    if rows.empty:
        return {
            "id": "B1",
            "passed": None,
            "predicted": {},
            "reference": {"DMC_eV": ref["DMC"], "PBE_eV": ref["PBE"]},
            "metric": f"absolute difference (tier {tier_label}, uncalibrated)",
            "summary": f"no basal monomer in the {tier_label} landscape table",
        }
    r = rows.iloc[0]
    e = float(r["value_eV"])
    t_el = (settings or {}).get("electronic_temperature")
    energy_kind = (settings or {}).get("energy", "F")
    cell = str(r["supercell"]).replace("graphene_", "")
    return {
        "id": "B1",
        "passed": None,
        "predicted": {
            "E_bind_H_atom_eV": e,
            "tier": tier_label,
            "engine": str(r["engine"]) if "engine" in r else None,
            "supercell": str(r["supercell"]),
            "electronic_temperature_K": t_el,
            "energy": energy_kind,
        },
        "reference": {"DMC_eV": ref["DMC"], "PBE_eV": ref["PBE"]},
        "metric": f"absolute difference (tier {tier_label}, uncalibrated)",
        "summary": (
            f"{tier_label} basal monomer {e:.3f} eV ({cell}"
            + (f", {t_el:.0f} K smearing" if t_el else "")
            + (", zero-temperature extrapolation" if energy_kind == "E0" else "")
            + f") vs DMC {ref['DMC']} / PBE {ref['PBE']} eV: offsets {e - ref['DMC']:+.2f} / "
            f"{e - ref['PBE']:+.2f} eV (input to C3 calibration)"
        ),
    }


def _landscape_table(tier: str, sub: str = "") -> tuple[Path, dict | None] | None:
    """The landscape table of the full tier, else of ``tier``, in the optional subdirectory."""
    import json

    for t in ("full", tier):
        base = Path(f"data/processed/{t}") / sub if sub else Path(f"data/processed/{t}")
        table = base / f"landscape_raw_{t}.parquet"
        if table.exists():
            manifest = base / f"landscape_manifest_{t}.json"
            settings = None
            if manifest.exists():
                m = json.loads(manifest.read_text())
                settings = dict(m.get("calculator", {}).get("settings") or {})
                settings["energy"] = m.get("energy", "F")
            return table, settings
    return None


def run_B1_T3(reg: dict, tier: str = "verify") -> dict:
    """Tier-3 (MACE) basal monomer versus B1 from the ``tier3`` landscape table."""
    import pandas as pd

    found = _landscape_table(tier, "tier3")
    if found is None:
        return {
            "id": "B1",
            "passed": None,
            "summary": "Tier-3 landscape not present (run landscape build --calculator mace)",
        }
    table, settings = found
    return evaluate_B1_T2(pd.read_parquet(table), reg, settings, tier_label="T3")


def run_B1_T2(reg: dict, tier: str = "verify") -> dict:
    """Tier-2 (GFN2-xTB) basal monomer versus B1 from the landscape table of the full tier, else
    of the requested tier; the engine check of ``verify-atomistic`` is the last fallback."""
    import json

    import pandas as pd

    found = _landscape_table(tier)
    if found is not None:
        table, settings = found
        return evaluate_B1_T2(pd.read_parquet(table), reg, settings)
    p = Path("data/processed/verify/atomistic_verify_T2.json")
    if not p.exists():
        return {
            "id": "B1",
            "passed": None,
            "summary": "Tier-2 landscape not present (run landscape build)",
        }
    d = json.loads(p.read_text())
    e = d["basal_monomer"]["E_bind_H_atom"]
    ref = reg["B1"]["reference"]
    return {
        "id": "B1",
        "passed": None,
        "predicted": {"E_bind_H_atom_T2_eV": e},
        "reference": {"DMC_eV": ref["DMC"], "PBE_eV": ref["PBE"]},
        "metric": "absolute difference (tier T2, uncalibrated)",
        "summary": (
            f"T2 basal monomer {e:.3f} eV (engine check, 4x4 cell, 300 K) vs DMC {ref['DMC']} / "
            f"PBE {ref['PBE']} eV"
        ),
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


def evaluate_B2(frame: pd.DataFrame, reg: dict) -> dict:
    """Benchmark B2 applied to a barrier table: the basal monomer sticking barrier per tier and
    method against 0.2 eV. The pass flag is set from the Tier-1 value only; Tier-2 and Tier-3
    values are reported as inputs to the calibration."""
    ref = float(reg["B2"]["reference"]["value"])
    thr = float(reg["B2"]["thresholds"]["T1_abs_eV"])
    rows = frame[
        (frame["quantity"] == "E_barrier_stick")
        & (frame["site_class"] == "basal")
        & (frame["reference_state"] == "H_atom")
    ]
    if rows.empty:
        return {
            "id": "B2",
            "passed": None,
            "predicted": {},
            "reference": {"value_eV": ref, "threshold_T1_eV": thr},
            "metric": "absolute difference of the sticking barrier (tier T1)",
            "summary": "no basal monomer sticking barrier in the barrier table",
        }
    predicted: dict[str, float] = {}
    saddle: dict[str, bool] = {}
    for _, r in rows.iterrows():
        method = str(r["identifier"]).split("@")[-1] if "@" in str(r["identifier"]) else "unknown"
        key = f"{r['tier']}:{method}:{r['supercell']}"
        predicted[key] = float(r["value_eV"])
        saddle[key] = bool(r["converged"]) if "converged" in rows.columns else True
    located = {k: v for k, v in predicted.items() if saddle[k]}
    t1 = {k: v for k, v in located.items() if k.startswith("T1:")}
    passed = None
    if t1:
        passed = bool(all(abs(v - ref) < thr for v in t1.values()))
    parts = "; ".join(
        f"{k} {v:+.3f}" + ("" if saddle[k] else " (path maximum, no saddle located)")
        for k, v in predicted.items()
    )
    return {
        "id": "B2",
        "passed": passed,
        "predicted": predicted,
        "saddle_located": saddle,
        "reference": {"value_eV": ref, "threshold_T1_eV": thr},
        "checks": {k: bool(abs(v - ref) < thr) for k, v in located.items()},
        "metric": "absolute difference of the sticking barrier (tier T1)",
        "summary": f"H sticking barrier on the basal plane, eV from H atom + substrate: {parts} vs {ref} eV"
        + ("" if t1 else " (no Tier-1 saddle value; reported as calibration input)"),
    }


def run_B2(reg: dict, tier: str = "verify") -> dict:
    """Benchmark B2 from the barrier tables of the full and the requested tier, including the
    sticking-cell and NEB image-count convergence runs kept in their own subdirectories."""
    import pandas as pd

    frames = []
    for t in ("full", tier):
        root = Path(f"data/processed/{t}")
        bases = [root, root / "tier3"]
        bases += sorted(root.glob("stick_cell_*")) + sorted((root / "tier3").glob("neb_images_*"))
        for base in bases:
            p = base / f"barriers_{t}.parquet"
            if p.exists():
                frames.append(pd.read_parquet(p))
    for p in sorted(Path("data/processed/tier1").glob("barriers*.parquet")):
        frames.append(pd.read_parquet(p))
    if not frames:
        return {
            "id": "B2",
            "passed": None,
            "summary": "barrier table not present (run landscape barriers)",
        }
    return evaluate_B2(pd.concat(frames, ignore_index=True), reg)


B11_CLASSES = {
    "monovacancy": ("monovacancy_formation", "monovacancy formation"),
    "stone_wales": ("stone_wales_formation", "Stone-Wales formation"),
}


def _b11_deviation(value: float, ref: float | list[float]) -> float:
    """Signed distance of ``value`` from a reference value or from a reference interval (zero
    inside the interval)."""
    if isinstance(ref, list):
        lo, hi = float(ref[0]), float(ref[1])
        return 0.0 if lo <= value <= hi else (value - hi if value > hi else value - lo)
    return value - float(ref)


def evaluate_B11(frames: dict[str, pd.DataFrame], reg: dict) -> dict:
    """Benchmark B11 (defect energetics as engine sanity checks): the monovacancy and Stone-Wales
    formation energies of each tier against the registry values, with the Tier-1 and Tier-2
    thresholds; Tier-3 values are reported without a threshold. ``frames`` maps a tier label to
    a landscape or anchor table holding ``E_formation_defect`` records."""
    ref = reg["B11"]["reference"]
    thr = {k.split("_")[0]: float(v) for k, v in reg["B11"]["thresholds"].items()}
    predicted: dict[str, float] = {}
    checks: dict[str, bool] = {}
    for label, frame in frames.items():
        rows = frame[frame["quantity"] == "E_formation_defect"]
        if "engine" in rows.columns:
            rows = rows[~rows["engine"].astype(str).str.endswith(":PBE")]
        for cls, (ref_key, _) in B11_CLASSES.items():
            sub = rows[rows["site_class"] == cls]
            if sub.empty:
                continue
            r = sub.sort_values("n_atoms", ascending=False).iloc[0]
            key = f"{label}:{cls}:{r['supercell']}"
            value = float(r["value_eV"])
            predicted[key] = value
            if label in thr:
                checks[key] = bool(abs(_b11_deviation(value, ref[ref_key])) < thr[label])
    if not predicted:
        return {
            "id": "B11",
            "passed": None,
            "predicted": {},
            "reference": {k: ref[k] for k in ("monovacancy_formation", "stone_wales_formation")},
            "metric": "absolute deviation of defect formation energies per tier",
            "summary": "no defect formation energies in the landscape or anchor tables",
        }
    missing = [t for t in thr if not any(k.startswith(f"{t}:") for k in predicted)]
    parts = "; ".join(
        f"{k} {v:.2f}" + ("" if k not in checks else (" ok" if checks[k] else " off"))
        for k, v in predicted.items()
    )
    return {
        "id": "B11",
        "passed": bool(all(checks.values())) if checks else None,
        "predicted": predicted,
        "checks": checks,
        "reference": {
            "monovacancy_formation_eV": ref["monovacancy_formation"],
            "stone_wales_formation_eV": ref["stone_wales_formation"],
            "thresholds_eV": thr,
        },
        "metric": "absolute deviation of defect formation energies per tier",
        "summary": f"defect formation energies, eV: {parts} vs monovacancy {ref['monovacancy_formation']} "
        f"and Stone-Wales {ref['stone_wales_formation']} (thresholds T1 {thr.get('T1')}, T2 {thr.get('T2')})"
        + (f"; no {', '.join(missing)} values yet" if missing else ""),
    }


def run_B11(reg: dict, tier: str = "verify") -> dict:
    """Benchmark B11 from the Tier-2 landscape table (full tier, else the requested tier), the
    Tier-3 table and the Tier-1 anchor table."""
    import pandas as pd

    frames: dict[str, pd.DataFrame] = {}
    found = _landscape_table(tier)
    if found is not None:
        frames["T2"] = pd.read_parquet(found[0])
    found3 = _landscape_table(tier, "tier3")
    if found3 is not None:
        frames["T3"] = pd.read_parquet(found3[0])
    anchors = Path("data/processed/tier1/anchors_tier1.parquet")
    if anchors.exists():
        frames["T1"] = pd.read_parquet(anchors)
    return evaluate_B11(frames, reg)


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


def _corpus_summary(path: Path = Path("data/processed/full/isosteric_heats_summary.parquet")):
    import pandas as pd

    return pd.read_parquet(path) if path.exists() else None


def evaluate_B4_corpus(summary, reg: dict) -> dict:
    """Corpus medians of the low-loading isosteric heat by material group against the B4 interval.

    The pass criterion is the registry's: the median over activated carbons lies inside the
    interval. The graphitic group (graphene, graphite, nanotubes) and the subsets with cryogenic
    data (a temperature at or below 100 K) are reported beside it, since Clausius-Clapeyron heats
    from room-temperature pairs of digitized isotherms are poorly conditioned.
    """
    import numpy as np

    lo, hi = reg["B4"]["reference"]["interval"]
    q = summary["q_st_low_kJ_mol"].astype(float)
    cryo = summary["temperatures_K"].apply(lambda t: min(t) <= 100.0)
    groups = {
        "activated_carbon": summary["material_class"] == "activated_carbon",
        "graphitic": summary["graphitic"].astype(bool),
        "all_carbons": np.ones(len(summary), dtype=bool),
    }
    out: dict[str, dict] = {}
    for name, mask in groups.items():
        for label, m in ((name, mask), (name + "_cryogenic", mask & cryo)):
            vals = q[m]
            out[label] = {
                "n": int(len(vals)),
                "median_kJ_mol": float(vals.median()) if len(vals) else None,
                "q25_kJ_mol": float(vals.quantile(0.25)) if len(vals) else None,
                "q75_kJ_mol": float(vals.quantile(0.75)) if len(vals) else None,
                "within_interval": bool(lo <= vals.median() <= hi) if len(vals) else None,
            }
    ac = out["activated_carbon"]
    passed = ac["within_interval"] if ac["n"] >= 3 else None
    return {
        "id": "B4-corpus",
        "passed": passed,
        "predicted": out,
        "reference": {"interval_kJ_mol": [lo, hi], "AX21_average": reg["B4"]["reference"]["AX21_average"]},
        "metric": "corpus median of the low-loading isosteric heat (activated carbons) within the interval",
        "summary": "; ".join(
            f"{k}: n={v['n']}, median {v['median_kJ_mol']:.1f} kJ/mol"
            for k, v in out.items()
            if v["n"]
        )
        + f" vs [{lo}, {hi}]"
        + ("" if ac["n"] >= 3 else "; fewer than three activated carbons, not decided"),
    }


def run_B4_corpus(reg: dict) -> dict:
    """Benchmark B4 from the corpus isosteric heats (``sitewindow corpus isosteric``)."""
    summary = _corpus_summary()
    if summary is None or summary.empty:
        return {
            "id": "B4-corpus",
            "passed": None,
            "summary": "isosteric-heat table not present (run corpus isosteric)",
        }
    return evaluate_B4_corpus(summary, reg)


def evaluate_B7(points, reg: dict, n_boot: int = 2000, seed: int = 0) -> dict:
    """Slope of the maximum 77 K excess uptake (wt%) against the BET area (per 1000 m2/g).

    The Chahine-type rule is a proportionality, so the statistic is the through-origin
    least-squares slope with a bootstrap over materials; the free-intercept fit is reported
    beside it. ``points`` holds one row per material with ``bet_area_m2_g`` and
    ``uptake_wt_percent``.
    """
    import numpy as np

    lo, hi = reg["B7"]["reference"]["interval"]
    x = points["bet_area_m2_g"].to_numpy(dtype=float) / 1000.0
    y = points["uptake_wt_percent"].to_numpy(dtype=float)
    n = len(x)
    if n < 3:
        return {
            "id": "B7",
            "passed": None,
            "predicted": {"n_materials": n},
            "reference": {"interval_wt_percent_per_1000m2_g": [lo, hi]},
            "metric": "through-origin slope with bootstrap 95 % interval overlapping the interval",
            "summary": f"{n} material(s) with a verified BET area and a 77 K isotherm reaching the plateau pressure; three are needed",
        }
    slope = float(np.sum(x * y) / np.sum(x * x))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = x[idx], y[idx]
        boots.append(np.sum(xb * yb) / np.sum(xb * xb))
    ci = [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))]
    a1, a0 = np.polyfit(x, y, 1)
    passed = bool(ci[0] <= hi and ci[1] >= lo)
    return {
        "id": "B7",
        "passed": passed,
        "predicted": {
            "n_materials": n,
            "slope_through_origin_wt_percent_per_1000m2_g": slope,
            "bootstrap_95_interval": ci,
            "free_intercept_fit": {"slope": float(a1), "intercept_wt_percent": float(a0)},
            "materials": points.to_dict(orient="records"),
        },
        "reference": {"interval_wt_percent_per_1000m2_g": [lo, hi]},
        "metric": "through-origin slope with bootstrap 95 % interval overlapping the interval",
        "summary": f"slope {slope:.2f} wt% per 1000 m2/g (95 % bootstrap [{ci[0]:.2f}, {ci[1]:.2f}], n={n}; free-intercept slope {a1:.2f}) vs [{lo}, {hi}]",
    }


def b7_points(plateau_pressure_bar: float = 20.0):
    """Materials with a verified, unambiguous BET area whose 77 K isotherm reaches the plateau
    pressure, with the maximum reported (excess) uptake in wt%."""
    import pandas as pd

    from sitewindow.experimental.corpus import load_corpus, material_id
    from sitewindow.experimental.metadata import load_metadata
    from sitewindow.units import mol_per_kg_to_wt_percent

    path = Path("data/processed/isodb_carbon_h2.parquet")
    meta = load_metadata()
    if not path.exists() or not meta:
        return pd.DataFrame(columns=["material_id", "bet_area_m2_g", "uptake_wt_percent"])
    t = load_corpus(path)
    rows = []
    t = t[(t["temperature_K"] - 77.0).abs() < 1.5]
    for (hk, doi), g in t.groupby(["adsorbent_hashkey", "doi"]):
        mid = material_id(hk, doi)
        m = meta.get(mid)
        if m is None or not m.usable_area:
            continue
        pmax = max(float(p.max()) for p in g["pressure_bar"])
        if pmax < plateau_pressure_bar:
            continue
        nmax = max(float(n.max()) for n in g["loading_mol_kg"])
        rows.append(
            {
                "material_id": mid,
                "adsorbent_name": g["adsorbent_name"].iloc[0],
                "doi": doi,
                "bet_area_m2_g": m.bet_area_m2_g,
                "max_pressure_bar": pmax,
                "uptake_wt_percent": float(mol_per_kg_to_wt_percent(nmax)),
                "where_read": m.where_read,
            }
        )
    return pd.DataFrame(rows)


def run_B7(reg: dict) -> dict:
    """Benchmark B7 from the corpus and the material metadata registry."""
    return evaluate_B7(b7_points(), reg)


def run_all(tier: str = "verify") -> dict:
    reg = _registry()
    results = [
        run_B8(reg),
        run_B9_endpoints(reg),
        run_B10(reg),
        run_B1_T2(reg, tier),
        run_B1_T3(reg, tier),
        run_B1_T1(reg, tier),
        run_B1_T1_convergence(reg, tier),
        run_B2(reg, tier),
        run_B3_pair_potential(reg),
        run_B4_computed(reg),
        run_B4_corpus(reg),
        run_B7(reg),
        run_B11(reg, tier),
    ]
    return {"tier": tier, "benchmarks": results}


__all__ = [
    "run_all",
    "run_B8",
    "run_B9_endpoints",
    "run_B10",
    "run_B1_T2",
    "run_B1_T3",
    "evaluate_B1_T2",
    "run_B1_T1",
    "run_B1_T1_convergence",
    "evaluate_convergence",
    "run_B2",
    "evaluate_B2",
    "run_B3_pair_potential",
    "run_B4_computed",
    "evaluate_B4_corpus",
    "run_B4_corpus",
    "evaluate_B7",
    "b7_points",
    "run_B7",
    "evaluate_B11",
    "run_B11",
]
