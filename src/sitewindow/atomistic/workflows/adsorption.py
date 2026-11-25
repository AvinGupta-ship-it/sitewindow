"""Relaxation, single-point and binding-energy workflows with content-addressed caching.

Binding energies (DERIVATIONS and manual C2):
    E_b(H; 1/2 H2) = E(sub+H) - E(sub) - E(H2)/2
    E_b(H; atom)   = E(sub+H) - E(sub) - E(H)
    E_b(H2)        = E(sub+H2) - E(sub) - E(H2)

Cache keys combine the structure hash, the calculator settings hash, the force tolerance and the
number of unpaired electrons, so that spin states of the same geometry are cached separately.
Under MPI every rank computes and only rank 0 writes.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import CalculationFailed
from ase.io import read, write
from ase.optimize import BFGS, FIRE

from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
from sitewindow.atomistic.parallel import barrier, broadcast_array, is_master
from sitewindow.atomistic.structures.builders import (
    hydrogen_atom,
    hydrogen_molecule,
    structure_hash,
)

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class RelaxResult:
    energy_eV: float
    max_force_eV_A: float
    n_steps: int
    converged: bool
    structure_hash: str
    settings_hash: str
    wall_s: float
    cached: bool
    scf_retried: bool = False


@dataclass(frozen=True)
class SinglePointResult:
    energy_eV: float
    structure_hash: str
    settings_hash: str
    wall_s: float
    cached: bool
    scf_retried: bool = False


def _robust(calc: Calculator) -> Calculator | None:
    """A variant of ``calc`` for retrying a failed self-consistent field, when the adapter has one."""
    robust = getattr(calc, "robust", None)
    return robust() if callable(robust) else None


def unpaired_electrons(atoms: Atoms) -> int:
    return int(round(float(np.sum(atoms.get_initial_magnetic_moments()))))


def constraint_tag(constraints: Sequence[Any] | None) -> str:
    """Short deterministic tag of a list of ASE constraints (from their ``todict`` forms)."""
    if not constraints:
        return ""
    payload = json.dumps(
        [c.todict() for c in constraints], sort_keys=True, default=lambda o: o.tolist()
    )
    return "_c" + hashlib.sha256(payload.encode()).hexdigest()[:8]


def _cache_key(
    atoms: Atoms,
    spec: CalculatorSpec,
    fmax: float,
    fixed_cell: bool,
    constraints: Sequence[Any] | None = None,
) -> str:
    return (
        f"{structure_hash(atoms)}_{spec.settings_hash}_f{fmax:g}_{int(fixed_cell)}"
        f"_u{unpaired_electrons(atoms)}{constraint_tag(constraints)}"
    )


def _cached(meta_p: Path) -> bool:
    """Whether a cache entry exists, decided by rank 0 so that every rank takes the same path."""
    return bool(broadcast_array([float(meta_p.exists())])[0])


def _store(cache_dir: Path, key: str, atoms: Atoms, meta: dict) -> None:
    """Write the cache entry. ``ase.io.write`` and ``ase.io.read`` are collective under an MPI
    world (rank 0 performs the I/O and broadcasts), so every rank must call them; only the plain
    JSON write is restricted to rank 0."""
    out = atoms.copy()
    out.info = {}
    out.calc = None
    write(cache_dir / f"{key}.extxyz", out)
    if is_master():
        (cache_dir / f"{key}.json").write_text(json.dumps(meta))
    barrier()


def relax(
    atoms: Atoms,
    calc: Calculator,
    cache_dir: Path | None = None,
    fmax: float = 0.05,
    steps: int = 400,
    optimizer: str = "BFGS",
    fixed_cell: bool = True,
    constraints: Sequence[Any] | None = None,
) -> tuple[Atoms, RelaxResult]:
    """Relax ``atoms`` with ``calc``; results are cached by structure, settings and spin state.

    The returned structure keeps ``atoms.info`` (recipe metadata) and its initial magnetic
    moments; only positions change. ``constraints`` (ASE constraint objects, for example a fixed
    bond length or a fixed plane along a reaction coordinate) replace any constraint carried by
    ``atoms``; they are part of the cache key and the reported maximum force is the constrained
    one.
    """
    spec = CalculatorSpec(calc.name, calc.tier, calc.settings())
    key = _cache_key(atoms, spec, fmax, fixed_cell, constraints)
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        if is_master():
            cache_dir.mkdir(parents=True, exist_ok=True)
        barrier()
        meta_p = cache_dir / f"{key}.json"
        if _cached(meta_p):
            meta = json.loads(meta_p.read_text())
            relaxed = atoms.copy()
            relaxed.positions = read(cache_dir / f"{key}.extxyz").positions
            return relaxed, RelaxResult(**{**meta, "cached": True})
    a = atoms.copy()
    if constraints is not None:
        a.set_constraint(list(constraints))
    calc.attach(a)
    t0 = time.perf_counter()
    opt = (BFGS if optimizer == "BFGS" else FIRE)(a, logfile=None)
    retried = False
    try:
        converged = bool(opt.run(fmax=fmax, steps=steps))
    except CalculationFailed:
        robust = _robust(calc)
        if robust is None:
            raise
        log.warning("relax %s: SCF failed, retrying with the damped mixer", key)
        retried = True
        robust.attach(a)
        done = int(opt.get_number_of_steps())
        opt = (BFGS if optimizer == "BFGS" else FIRE)(a, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=max(steps - done, 1)))
    e = float(a.get_potential_energy())
    fmax_final = float(np.linalg.norm(a.get_forces(), axis=1).max()) if len(a) > 1 else 0.0
    a.calc = None
    # Under MPI every rank holds rank 0's result, whatever the calculator's own parallelism.
    a.positions = broadcast_array(a.positions)
    scalars = broadcast_array([e, fmax_final, float(opt.get_number_of_steps()), float(converged)])
    res = RelaxResult(
        float(scalars[0]),
        float(scalars[1]),
        int(scalars[2]),
        bool(scalars[3]),
        structure_hash(atoms),
        spec.settings_hash,
        time.perf_counter() - t0,
        False,
        retried,
    )
    log.info(
        "relax %s n=%d steps=%d fmax=%.3f converged=%s wall=%.1fs",
        key,
        len(a),
        res.n_steps,
        res.max_force_eV_A,
        res.converged,
        res.wall_s,
    )
    if cache_dir is not None:
        _store(cache_dir, key, a, {k: v for k, v in asdict(res).items() if k != "cached"})
    return a, res


def single_point(
    atoms: Atoms, calc: Calculator, cache_dir: Path | None = None
) -> SinglePointResult:
    """Total energy of ``atoms`` with ``calc`` without relaxation, cached like :func:`relax`."""
    spec = CalculatorSpec(calc.name, calc.tier, calc.settings())
    key = f"{structure_hash(atoms)}_{spec.settings_hash}_sp_u{unpaired_electrons(atoms)}"
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        if is_master():
            cache_dir.mkdir(parents=True, exist_ok=True)
        barrier()
        meta_p = cache_dir / f"{key}.json"
        if _cached(meta_p):
            meta = json.loads(meta_p.read_text())
            return SinglePointResult(**{**meta, "cached": True})
    a = atoms.copy()
    calc.attach(a)
    t0 = time.perf_counter()
    retried = False
    try:
        e_local = a.get_potential_energy()
    except CalculationFailed:
        robust = _robust(calc)
        if robust is None:
            raise
        log.warning("single point %s: SCF failed, retrying with the damped mixer", key)
        retried = True
        robust.attach(a)
        e_local = a.get_potential_energy()
    e = float(broadcast_array([e_local])[0])
    res = SinglePointResult(
        e, structure_hash(atoms), spec.settings_hash, time.perf_counter() - t0, False, retried
    )
    log.info("single point %s n=%d wall=%.1fs", key, len(a), res.wall_s)
    if cache_dir is not None:
        _store(cache_dir, key, a, {k: v for k, v in asdict(res).items() if k != "cached"})
    a.calc = None
    return res


def reference_energies(
    calc: Calculator, cache_dir: Path | None = None, pbc: bool = False
) -> dict[str, float]:
    """Relaxed energies of H2 and the H atom with the same calculator."""
    h2, r2 = relax(hydrogen_molecule(pbc=pbc), calc, cache_dir, fmax=0.01)
    h, r1 = relax(hydrogen_atom(pbc=pbc), calc, cache_dir, fmax=0.01, steps=0)
    return {"E_H2": r2.energy_eV, "E_H": r1.energy_eV, "H2_bond_A": float(h2.get_distance(0, 1))}


def binding_energies(
    E_sub: float, E_sub_ads: float, refs: dict[str, float], adsorbate: str = "H", n: int = 1
) -> dict[str, float]:
    """Binding energy per adsorbed entity for ``n`` H atoms or one H2 molecule."""
    if adsorbate == "H":
        return {
            "E_bind_H_half_H2": (E_sub_ads - E_sub - 0.5 * n * refs["E_H2"]) / n,
            "E_bind_H_atom": (E_sub_ads - E_sub - n * refs["E_H"]) / n,
        }
    if adsorbate == "H2":
        return {"E_bind_H2": E_sub_ads - E_sub - refs["E_H2"]}
    raise ValueError(adsorbate)


__all__ = [
    "constraint_tag",
    "relax",
    "single_point",
    "RelaxResult",
    "SinglePointResult",
    "reference_energies",
    "binding_energies",
    "unpaired_electrons",
]
