"""Tier-1 density-functional anchors and convergence study (manual C2; benchmark B1).

The basal-monomer anchor relaxes graphene and graphene+H at Tier 2 (GFN2-xTB, when available)
for a starting geometry, re-relaxes them with GPAW in LCAO dzp mode (PBE-D3(BJ), spin-polarized),
and evaluates plane-wave single points for the energies. The PBE and D3 parts are reported
separately, so that the PBE binding energy can be compared with the plane-wave PBE value of
benchmark B1 (-0.821 eV, Dumi et al. 2022) at the registry threshold. The convergence series
evaluates the plane-wave binding energy for the supercell, k-mesh, cutoff and vacuum variants of
acceptance criterion B1 on the LCAO-relaxed geometries.

Run under MPI with ``mpirun -np N gpaw python --module sitewindow.workflows.cli tier1 ...``.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms

from sitewindow import __version__
from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
from sitewindow.atomistic.calculators.gpaw_adapter import GPAWCalculator, dispersion_energy_eV
from sitewindow.atomistic.parallel import is_master
from sitewindow.atomistic.structures.builders import (
    add_hydrogens,
    central_atom_index,
    graphene_supercell,
    hydrogen_atom,
    hydrogen_molecule,
    set_unpaired_electrons,
    structure_hash,
)
from sitewindow.atomistic.workflows.adsorption import relax, single_point
from sitewindow.errors import CalculatorUnavailable, ConvergenceError

log = logging.getLogger(__name__)

B1_PBE_EV = -0.821
B1_DMC_EV = -0.691
B1_THRESHOLD_EV = 0.10


@dataclass(frozen=True)
class Tier1Settings:
    supercell: int = 5
    vacuum_A: float = 15.0
    ecut_eV: float = 500.0
    kmesh: tuple[int, int, int] | None = None
    fmax_eV_A: float = 0.02
    relax_steps: int = 150
    prerelax_tier2: bool = True
    dispersion: str | None = "d3bj"
    box_A: float = 15.0

    def lcao(self, txt: str | None = None) -> GPAWCalculator:
        return GPAWCalculator(mode="lcao", kmesh=self.kmesh, dispersion=self.dispersion, txt=txt)

    def pw(self, txt: str | None = None) -> GPAWCalculator:
        return GPAWCalculator(
            mode="pw", ecut_eV=self.ecut_eV, kmesh=self.kmesh, dispersion=self.dispersion, txt=txt
        )

    def molecule_pw(self, txt: str | None = None) -> GPAWCalculator:
        return GPAWCalculator(
            mode="pw",
            ecut_eV=self.ecut_eV,
            kmesh=(1, 1, 1),
            dispersion=self.dispersion,
            dipole_correction=False,
            txt=txt,
        )


def _log_path(log_dir: Path | None, *parts: str) -> str | None:
    return None if log_dir is None else str(Path(log_dir) / ("_".join(parts) + ".txt"))


def check_h2_bond(bond_A: float, window_A: tuple[float, float] = (0.65, 0.85)) -> None:
    """Reject an H2 reference whose relaxed bond length is outside ``window_A`` (dissociation)."""
    lo, hi = window_A
    if not lo < bond_A < hi:
        raise ConvergenceError(
            f"H2 reference relaxed to a bond length of {bond_A:.3f} A; the molecule dissociated "
            "(the initial magnetic moments must start H2 spin-paired)"
        )


def _tier2_calculator() -> Calculator | None:
    try:
        import tblite.ase  # noqa: F401

        from sitewindow.atomistic.calculators.tblite_adapter import TBLiteCalculator

        return TBLiteCalculator()
    except Exception:
        return None


def _require_gpaw() -> None:
    try:
        import gpaw  # noqa: F401
    except ImportError as e:
        raise CalculatorUnavailable(
            "gpaw is not installed; use the conda-forge environment in environment.yml"
        ) from e


def _structures(settings: Tier1Settings) -> tuple[Atoms, Atoms]:
    g = graphene_supercell(settings.supercell, vacuum=settings.vacuum_A)
    set_unpaired_electrons(g, 0)
    gh = add_hydrogens(g, [(central_atom_index(g), "up")])
    return g, gh


def _relaxed_geometries(
    settings: Tier1Settings, cache: Path, log_dir: Path | None
) -> dict[str, Any]:
    """Tier-2 pre-relaxation followed by Tier-1 LCAO relaxation of graphene and graphene+H."""
    g, gh = _structures(settings)
    t2 = _tier2_calculator() if settings.prerelax_tier2 else None
    if t2 is not None:
        g, _ = relax(g, t2, cache, fmax=0.05, steps=400)
        gh, _ = relax(gh, t2, cache, fmax=0.05, steps=400)
    tag = f"{settings.supercell}x{settings.supercell}"
    lcao_g, res_g = relax(
        g,
        settings.lcao(_log_path(log_dir, tag, "lcao_graphene")),
        cache,
        settings.fmax_eV_A,
        settings.relax_steps,
    )
    lcao_gh, res_gh = relax(
        gh,
        settings.lcao(_log_path(log_dir, tag, "lcao_graphene_H")),
        cache,
        settings.fmax_eV_A,
        settings.relax_steps,
    )
    return {"graphene": (lcao_g, res_g), "graphene_H": (lcao_gh, res_gh)}


def _references(settings: Tier1Settings, cache: Path, log_dir: Path | None) -> dict[str, Any]:
    h2 = hydrogen_molecule(settings.box_A, pbc=True)
    set_unpaired_electrons(h2, 0)
    lcao_mol = replace(
        settings.lcao(), kmesh=(1, 1, 1), dipole_correction=False, txt=_log_path(log_dir, "lcao_H2")
    )
    h2_relaxed, _ = relax(h2, lcao_mol, cache, fmax=0.01, steps=100)
    check_h2_bond(float(h2_relaxed.get_distance(0, 1)))
    h = hydrogen_atom(settings.box_A, pbc=True)
    tag = f"{settings.ecut_eV:.0f}eV"
    sp_h2 = single_point(h2_relaxed, settings.molecule_pw(_log_path(log_dir, tag, "pw_H2")), cache)
    sp_h = single_point(h, settings.molecule_pw(_log_path(log_dir, tag, "pw_H")), cache)
    d3 = {} if settings.dispersion is None else {"H2": dispersion_energy_eV(h2_relaxed), "H": 0.0}
    return {
        "E_H2": sp_h2.energy_eV,
        "E_H": sp_h.energy_eV,
        "H2_bond_A": float(h2_relaxed.get_distance(0, 1)),
        "D3": d3,
    }


def anchor_basal_monomer(
    settings: Tier1Settings, cache_dir: Path, out_dir: Path, log_dir: Path | None = None
) -> dict[str, Any]:
    """Tier-1 basal-monomer anchor: energies, binding energies, geometry and the B1 comparison."""
    _require_gpaw()
    t0 = time.perf_counter()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    if log_dir is not None and is_master():
        Path(log_dir).mkdir(parents=True, exist_ok=True)
    geoms = _relaxed_geometries(settings, cache_dir, log_dir)
    refs = _references(settings, cache_dir, log_dir)
    pw = settings.pw
    g, res_g = geoms["graphene"]
    gh, res_gh = geoms["graphene_H"]
    tag = f"{settings.supercell}x{settings.supercell}_{settings.ecut_eV:.0f}eV"
    sp_g = single_point(g, pw(_log_path(log_dir, tag, "pw_graphene")), cache_dir)
    sp_gh = single_point(gh, pw(_log_path(log_dir, tag, "pw_graphene_H")), cache_dir)
    d3_g = 0.0 if settings.dispersion is None else dispersion_energy_eV(g)
    d3_gh = 0.0 if settings.dispersion is None else dispersion_energy_eV(gh)
    d3_h2 = refs["D3"].get("H2", 0.0)
    e_b_atom = sp_gh.energy_eV - sp_g.energy_eV - refs["E_H"]
    e_b_half = sp_gh.energy_eV - sp_g.energy_eV - 0.5 * refs["E_H2"]
    d3_part_atom = d3_gh - d3_g
    d3_part_half = d3_gh - d3_g - 0.5 * d3_h2
    e_b_atom_pbe = e_b_atom - d3_part_atom
    e_b_half_pbe = e_b_half - d3_part_half
    e_b_atom_lcao = res_gh.energy_eV - res_g.energy_eV - refs["E_H"]
    c_index = central_atom_index(g)
    h_index = len(gh) - 1
    ch_bond = float(gh.get_distance(c_index, h_index))
    pucker = float(gh.positions[c_index, 2] - np.median(gh.positions[: len(gh) - 1, 2]))
    kmesh = settings.pw().kpoints(g)
    result: dict[str, Any] = {
        "tier": "T1",
        "engine": "gpaw",
        "package_version": __version__,
        "settings": {**asdict(settings), "kmesh_used": list(kmesh)},
        "settings_hash_pw": CalculatorSpec("gpaw", "T1", pw().settings()).settings_hash,
        "structure_hashes": {"graphene": structure_hash(g), "graphene_H": structure_hash(gh)},
        "energies_eV": {
            "pw": {
                "graphene": sp_g.energy_eV,
                "graphene_H": sp_gh.energy_eV,
                **{k: refs[k] for k in ("E_H2", "E_H")},
            },
            "lcao": {"graphene": res_g.energy_eV, "graphene_H": res_gh.energy_eV},
            "d3": {"graphene": d3_g, "graphene_H": d3_gh, "H2": d3_h2},
        },
        "binding_eV": {
            "E_bind_H_atom_PBE_D3": e_b_atom,
            "E_bind_H_atom_PBE": e_b_atom_pbe,
            "E_bind_H_half_H2_PBE_D3": e_b_half,
            "E_bind_H_half_H2_PBE": e_b_half_pbe,
            "D3_contribution_atom": d3_part_atom,
            "E_bind_H_atom_PBE_D3_lcao": e_b_atom_lcao,
        },
        "geometry": {"C_H_bond_A": ch_bond, "puckering_A": pucker, "H2_bond_A": refs["H2_bond_A"]},
        "relaxation": {
            "graphene": {
                "converged": res_g.converged,
                "max_force_eV_A": res_g.max_force_eV_A,
                "n_steps": res_g.n_steps,
            },
            "graphene_H": {
                "converged": res_gh.converged,
                "max_force_eV_A": res_gh.max_force_eV_A,
                "n_steps": res_gh.n_steps,
            },
        },
        "benchmark_B1": {
            "reference_PBE_eV": B1_PBE_EV,
            "reference_DMC_eV": B1_DMC_EV,
            "threshold_eV": B1_THRESHOLD_EV,
            "deviation_from_PBE_eV": e_b_atom_pbe - B1_PBE_EV,
            "within_threshold": bool(abs(e_b_atom_pbe - B1_PBE_EV) < B1_THRESHOLD_EV),
            "PBE_minus_DMC_offset_eV": e_b_atom_pbe - B1_DMC_EV,
        },
        "wall_s": time.perf_counter() - t0,
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        n = settings.supercell
        (out_dir / f"tier1_anchor_basal_{n}x{n}.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
    return result


def convergence_variants(base: Tier1Settings) -> dict[str, Tier1Settings]:
    """The variants of acceptance criterion B1: cell, k-mesh, cutoff and vacuum."""
    g = graphene_supercell(base.supercell, vacuum=base.vacuum_A)
    k_base = base.pw().kpoints(g)
    k_dense = (k_base[0] + 2, k_base[1] + 2, 1)
    return {
        "base": base,
        "cell_plus_one": replace(base, supercell=base.supercell + 1),
        "kmesh_dense": replace(base, kmesh=k_dense),
        "cutoff_600": replace(base, ecut_eV=600.0),
        "vacuum_20": replace(base, vacuum_A=20.0),
    }


def _reembed(atoms: Atoms, vacuum_A: float) -> Atoms:
    a = atoms.copy()
    cell = np.array(a.cell)
    cell[2, 2] = vacuum_A + float(np.ptp(a.positions[:, 2]))
    a.set_cell(cell, scale_atoms=False)
    a.center(axis=2)
    return a


def convergence_series(
    base: Tier1Settings,
    cache_dir: Path,
    out_dir: Path,
    variants: list[str] | None = None,
    log_dir: Path | None = None,
) -> dict[str, Any]:
    """Plane-wave binding energies for the B1 convergence variants on LCAO-relaxed geometries.

    The k-mesh, cutoff and vacuum variants reuse the geometry of the base cell; the larger cell is
    relaxed on its own. Energy differences are reported relative to the base settings.
    """
    _require_gpaw()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    all_variants = convergence_variants(base)
    names = variants or list(all_variants)
    base_geoms = _relaxed_geometries(base, cache_dir, log_dir)
    rows: dict[str, Any] = {}
    for name in names:
        s = all_variants[name]
        if s.supercell != base.supercell:
            geoms = _relaxed_geometries(s, cache_dir, log_dir)
        else:
            geoms = base_geoms
        g, _ = geoms["graphene"]
        gh, _ = geoms["graphene_H"]
        if s.vacuum_A != base.vacuum_A:
            g, gh = _reembed(g, s.vacuum_A), _reembed(gh, s.vacuum_A)
        refs = _references(s, cache_dir, log_dir)
        sp_g = single_point(g, s.pw(_log_path(log_dir, name, "pw_graphene")), cache_dir)
        sp_gh = single_point(gh, s.pw(_log_path(log_dir, name, "pw_graphene_H")), cache_dir)
        e_b = sp_gh.energy_eV - sp_g.energy_eV - refs["E_H"]
        d3 = 0.0 if s.dispersion is None else dispersion_energy_eV(gh) - dispersion_energy_eV(g)
        rows[name] = {
            "settings": {**asdict(s), "kmesh_used": list(s.pw().kpoints(g))},
            "E_bind_H_atom_PBE_D3": e_b,
            "E_bind_H_atom_PBE": e_b - d3,
            "n_atoms": len(gh),
        }
        log.info("convergence %s E_b(H atom)=%.4f eV", name, e_b)
    base_e = rows.get("base", {}).get("E_bind_H_atom_PBE")
    for row in rows.values():
        row["delta_vs_base_PBE_eV"] = None if base_e is None else row["E_bind_H_atom_PBE"] - base_e
    result = {
        "tier": "T1",
        "engine": "gpaw",
        "package_version": __version__,
        "criterion": "B1: |delta| < 0.02 eV for cell, k-mesh, cutoff and vacuum variants",
        "variants": rows,
        "all_within_0.02_eV": bool(
            all(
                r["delta_vs_base_PBE_eV"] is not None and abs(r["delta_vs_base_PBE_eV"]) < 0.02
                for n, r in rows.items()
                if n != "base"
            )
        )
        if len(rows) > 1
        else None,
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / f"tier1_convergence_basal_{base.supercell}x{base.supercell}.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
    return result


__all__ = [
    "Tier1Settings",
    "check_h2_bond",
    "anchor_basal_monomer",
    "convergence_variants",
    "convergence_series",
    "B1_PBE_EV",
    "B1_DMC_EV",
    "B1_THRESHOLD_EV",
]
