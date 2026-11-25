"""Relaxation and binding-energy workflows with content-addressed caching.

Binding energies (DERIVATIONS and manual C2):
    E_b(H; 1/2 H2) = E(sub+H) - E(sub) - E(H2)/2
    E_b(H; atom)   = E(sub+H) - E(sub) - E(H)
    E_b(H2)        = E(sub+H2) - E(sub) - E(H2)
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, FIRE

from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
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


def _cache_key(atoms: Atoms, spec: CalculatorSpec, fmax: float, fixed_cell: bool) -> str:
    return f"{structure_hash(atoms)}_{spec.settings_hash}_f{fmax:g}_{int(fixed_cell)}"


def relax(
    atoms: Atoms,
    calc: Calculator,
    cache_dir: Path | None = None,
    fmax: float = 0.05,
    steps: int = 400,
    optimizer: str = "BFGS",
    fixed_cell: bool = True,
) -> tuple[Atoms, RelaxResult]:
    """Relax ``atoms`` with ``calc``; results are cached by structure and settings hash."""
    spec = CalculatorSpec(calc.name, calc.tier, calc.settings())
    key = _cache_key(atoms, spec, fmax, fixed_cell)
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        meta_p = cache_dir / f"{key}.json"
        if meta_p.exists():
            meta = json.loads(meta_p.read_text())
            relaxed = read(cache_dir / f"{key}.extxyz")
            return relaxed, RelaxResult(**{**meta, "cached": True})
    a = atoms.copy()
    calc.attach(a)
    t0 = time.perf_counter()
    opt = (BFGS if optimizer == "BFGS" else FIRE)(a, logfile=None)
    converged = bool(opt.run(fmax=fmax, steps=steps))
    e = float(a.get_potential_energy())
    fmax_final = float(np.linalg.norm(a.get_forces(), axis=1).max()) if len(a) > 1 else 0.0
    res = RelaxResult(
        e,
        fmax_final,
        int(opt.get_number_of_steps()),
        converged,
        structure_hash(atoms),
        spec.settings_hash,
        time.perf_counter() - t0,
        False,
    )
    if cache_dir is not None:
        write(cache_dir / f"{key}.extxyz", a)
        (cache_dir / f"{key}.json").write_text(
            json.dumps({k: v for k, v in asdict(res).items() if k != "cached"})
        )
    a.calc = None
    return a, res


def reference_energies(calc: Calculator, cache_dir: Path | None = None) -> dict[str, float]:
    """Relaxed energies of H2 and the H atom with the same calculator."""
    h2, r2 = relax(hydrogen_molecule(), calc, cache_dir, fmax=0.01)
    h, r1 = relax(hydrogen_atom(), calc, cache_dir, fmax=0.01, steps=1)
    return {"E_H2": r2.energy_eV, "E_H": r1.energy_eV, "H2_bond_A": float(h2.get_distance(0, 1))}


def binding_energies(
    E_sub: float, E_sub_ads: float, refs: dict[str, float], adsorbate: str = "H"
) -> dict[str, float]:
    if adsorbate == "H":
        return {
            "E_bind_H_half_H2": E_sub_ads - E_sub - 0.5 * refs["E_H2"],
            "E_bind_H_atom": E_sub_ads - E_sub - refs["E_H"],
        }
    if adsorbate == "H2":
        return {"E_bind_H2": E_sub_ads - E_sub - refs["E_H2"]}
    raise ValueError(adsorbate)


__all__ = ["relax", "RelaxResult", "reference_energies", "binding_energies"]
