"""Electronic-temperature study of the Tier-2 protocol (manual C2, Tier 2).

tblite reports the Mermin free energy F = E - T S of the smeared electronic system, where E is
the internal energy at the Fermi-Dirac occupations and S the electronic entropy of those
occupations (``smearing_terms``). The zero-temperature extrapolation E_0 = F + T S / 2 removes
the smearing error to second order in T when the density of states is smooth across the Fermi
level. This workflow evaluates F, E and E_0 for the quantities whose smearing error matters to
the landscape:

- the far hydrogen atom above the unrelaxed sheet (the separated-reactant reference of
  sticking), whose offset from the separated fragments measures the error directly;
- the basal monomer, whose chemisorption creates a half-filled state at the Fermi level, in
  several cells and at several temperatures;
- one monovacancy and one epoxide configuration, where the carbon-hydrogen bond is localized.

Binding energies are formed from the same quantity (F, E or E_0) of every term, with the
relaxed H2 molecule and the free H atom as references. Every structure is relaxed at the
temperature of the row, as in the landscape workflow, and the decomposition is a single point
on the relaxed structure.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from ase import Atoms

from sitewindow import __version__
from sitewindow.atomistic.calculators.base import CalculatorSpec
from sitewindow.atomistic.calculators.tblite_adapter import TBLiteCalculator, smearing_terms
from sitewindow.atomistic.parallel import is_master
from sitewindow.atomistic.structures.builders import (
    add_hydrogen_top,
    add_hydrogens,
    central_atom_index,
    graphene_supercell,
    hydrogen_atom,
    hydrogen_molecule,
    set_unpaired_electrons,
    structure_hash,
)
from sitewindow.atomistic.workflows.adsorption import relax, unpaired_electrons
from sitewindow.landscape.classes import Configuration, catalogue

log = logging.getLogger(__name__)

ENERGY_KINDS = ("F", "E_internal", "E0")
_KEYS = {"F": "F_eV", "E_internal": "E_internal_eV", "E0": "E0_eV"}

STUDY_CONFIGURATIONS: dict[str, str] = {
    "basal": "monomer",
    "monovacancy": "seq1",
    "epoxide": "monomer_adjacent",
}


@dataclass(frozen=True)
class SmearingSpec:
    cells: tuple[int, ...] = (9, 10, 12)
    temperatures_K: tuple[float, ...] = (1500.0, 3000.0, 6000.0)
    far_temperatures_K: tuple[float, ...] = (300.0, 1500.0, 3000.0, 6000.0)
    classes: tuple[str, ...] = ("basal", "monovacancy", "epoxide")
    far_height_A: float = 6.0
    fmax: float = 0.05
    steps: int = 400


def tier2_calculator(temperature_K: float) -> TBLiteCalculator:
    """The Tier-2 calculator of the landscape protocol at another electronic temperature."""
    return TBLiteCalculator(electronic_temperature=float(temperature_K))


def smearing_single_point(
    atoms: Atoms, calc: TBLiteCalculator, cache_dir: Path | None = None
) -> dict[str, Any]:
    """Single point with the free-energy decomposition, cached like the other single points."""
    spec = CalculatorSpec(calc.name, calc.tier, calc.settings())
    key = f"{structure_hash(atoms)}_{spec.settings_hash}_smear_u{unpaired_electrons(atoms)}"
    meta_p = None
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        meta_p = cache_dir / f"{key}.json"
        if meta_p.exists():
            return {**json.loads(meta_p.read_text()), "cached": True}
    a = atoms.copy()
    calc.attach(a)
    t0 = time.perf_counter()
    a.get_potential_energy()
    terms = smearing_terms(a.calc, calc.electronic_temperature)
    charges = np.asarray(a.calc.results.get("charges", np.zeros(len(a))))
    h_mask = a.numbers == 1
    out: dict[str, Any] = {
        **terms,
        "electronic_temperature_K": calc.electronic_temperature,
        "n_atoms": len(a),
        "H_mulliken_charge": float(charges[h_mask].sum()) if h_mask.any() else None,
        "wall_s": time.perf_counter() - t0,
    }
    a.calc = None
    if meta_p is not None and is_master():
        meta_p.write_text(json.dumps(out))
    return {**out, "cached": False}


def _study_configuration(cls: str, n: int) -> Configuration:
    name = STUDY_CONFIGURATIONS[cls]
    for cfg in catalogue()[cls].configurations("full", {"supercell": n}):
        if cfg.name == name:
            return cfg
    raise KeyError(f"{cls}/{name} is not a catalogue configuration")


def _parity(atoms: Atoms) -> int:
    return int(np.sum(atoms.numbers)) % 2


def _references(calc: TBLiteCalculator, cache_dir: Path, fmax: float) -> dict[str, dict[str, Any]]:
    h2, _ = relax(hydrogen_molecule(), calc, cache_dir, fmax=0.01)
    h = hydrogen_atom()
    return {
        "H2": {
            **smearing_single_point(h2, calc, cache_dir),
            "bond_A": float(h2.get_distance(0, 1)),
        },
        "H": smearing_single_point(h, calc, cache_dir),
    }


def _binding(sub: dict[str, Any], ads: dict[str, Any], refs: dict[str, dict[str, Any]]) -> dict:
    out: dict[str, float] = {}
    for kind in ENERGY_KINDS:
        k = _KEYS[kind]
        out[f"E_bind_H_atom_{kind}"] = ads[k] - sub[k] - refs["H"][k]
        out[f"E_bind_H_half_H2_{kind}"] = ads[k] - sub[k] - 0.5 * refs["H2"][k]
    return out


def far_hydrogen_series(spec: SmearingSpec, cache_dir: Path) -> list[dict[str, Any]]:
    """E(sheet + H at ``far_height_A``) - E(sheet) - E(H) on the unrelaxed sheet for every cell
    and temperature, as F, E and E_0; zero up to physisorption when the reference is consistent."""
    rows: list[dict[str, Any]] = []
    for n in spec.cells:
        g = graphene_supercell(n)
        set_unpaired_electrons(g, 0)
        gh = add_hydrogen_top(g, central_atom_index(g), height=spec.far_height_A)
        set_unpaired_electrons(gh, 1)
        for T in spec.far_temperatures_K:
            calc = tier2_calculator(T)
            s = smearing_single_point(g, calc, cache_dir)
            a = smearing_single_point(gh, calc, cache_dir)
            h = smearing_single_point(hydrogen_atom(), calc, cache_dir)
            row: dict[str, Any] = {
                "cell": f"{n}x{n}",
                "n": n,
                "T_K": T,
                "height_A": spec.far_height_A,
                "TS_sheet_eV": s["TS_eV"],
                "TS_sheet_H_eV": a["TS_eV"],
                "H_mulliken_charge": a["H_mulliken_charge"],
                "n_fractional_orbitals": a["n_fractional_orbitals"],
            }
            for kind in ENERGY_KINDS:
                k = _KEYS[kind]
                row[f"offset_{kind}_eV"] = a[k] - s[k] - h[k]
            rows.append(row)
            log.info(
                "far H %s T=%.0f K offsets F %+.3f E %+.3f E0 %+.3f eV",
                row["cell"],
                T,
                row["offset_F_eV"],
                row["offset_E_internal_eV"],
                row["offset_E0_eV"],
            )
    return rows


def binding_series(spec: SmearingSpec, cache_dir: Path) -> list[dict[str, Any]]:
    """Binding energies of the study configurations in every cell and at every temperature,
    with each structure relaxed at its own temperature, as F, E and E_0."""
    rows: list[dict[str, Any]] = []
    for cls in spec.classes:
        for n in spec.cells:
            cfg = _study_configuration(cls, n)
            for T in spec.temperatures_K:
                calc = tier2_calculator(T)
                refs = _references(calc, cache_dir, spec.fmax)
                sub0 = cfg.substrate.copy()
                set_unpaired_electrons(sub0, _parity(sub0))
                sub, sub_res = relax(sub0, calc, cache_dir, fmax=spec.fmax, steps=spec.steps)
                sub.info.update(cfg.substrate.info)
                start = add_hydrogens(sub, cfg.placements)
                ads, ads_res = relax(start, calc, cache_dir, fmax=spec.fmax, steps=spec.steps)
                s = smearing_single_point(sub, calc, cache_dir)
                a = smearing_single_point(ads, calc, cache_dir)
                site_c = int(cfg.placements[-1][0])
                h_index = len(ads) - 1
                row: dict[str, Any] = {
                    "site_class": cls,
                    "configuration": cfg.name,
                    "cell": f"{n}x{n}",
                    "n": n,
                    "T_K": T,
                    "n_atoms": len(ads),
                    "converged": bool(sub_res.converged and ads_res.converged),
                    "TS_substrate_eV": s["TS_eV"],
                    "TS_adsorbed_eV": a["TS_eV"],
                    "TS_H_eV": refs["H"]["TS_eV"],
                    "TS_H2_eV": refs["H2"]["TS_eV"],
                    "delta_TS_eV": a["TS_eV"] - s["TS_eV"] - refs["H"]["TS_eV"],
                    "n_fractional_orbitals_substrate": s["n_fractional_orbitals"],
                    "n_fractional_orbitals_adsorbed": a["n_fractional_orbitals"],
                    "H_mulliken_charge": a["H_mulliken_charge"],
                    "C_H_bond_A": float(ads.get_distance(site_c, h_index, mic=True)),
                    **_binding(s, a, refs),
                }
                rows.append(row)
                log.info(
                    "%s/%s %s T=%.0f K E_b(H atom): F %.3f E %.3f E0 %.3f eV (dTS %.3f)",
                    cls,
                    cfg.name,
                    row["cell"],
                    T,
                    row["E_bind_H_atom_F"],
                    row["E_bind_H_atom_E_internal"],
                    row["E_bind_H_atom_E0"],
                    row["delta_TS_eV"],
                )
    return rows


def build_smearing_study(spec: SmearingSpec, cache_dir: Path, out_dir: Path) -> dict[str, Any]:
    """Run both series and write ``smearing_extrapolation.json`` and the binding table."""
    t0 = time.perf_counter()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    far = far_hydrogen_series(spec, cache_dir)
    binding = binding_series(spec, cache_dir)
    result: dict[str, Any] = {
        "tier": "T2",
        "engine": "tblite:GFN2-xTB",
        "package_version": __version__,
        "reported_energy": (
            "Mermin free energy F = E - T S (tblite adds the entropy term inside the "
            "self-consistent loop); E_internal = F + T S; E0 = F + T S / 2"
        ),
        "spec": asdict(spec),
        "far_hydrogen": far,
        "binding": binding,
        "wall_s": time.perf_counter() - t0,
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "smearing_extrapolation.json").write_text(json.dumps(result, indent=2) + "\n")
        pd.DataFrame(binding).to_csv(out_dir / "smearing_binding.csv", index=False)
        pd.DataFrame(far).to_csv(out_dir / "smearing_far_hydrogen.csv", index=False)
    return result


__all__ = [
    "ENERGY_KINDS",
    "STUDY_CONFIGURATIONS",
    "SmearingSpec",
    "tier2_calculator",
    "smearing_single_point",
    "far_hydrogen_series",
    "binding_series",
    "build_smearing_study",
]
