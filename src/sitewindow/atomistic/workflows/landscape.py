"""Landscape construction: relax every catalogue configuration at one fidelity tier and write
``SiteEnergyRecord`` tables (manual C2, experiment E1).

For each configuration the substrate is relaxed once, hydrogen atoms are placed on the relaxed
substrate, the hydrogenated structure is relaxed in the lowest spin multiplicity (and, when
requested, in higher ones, keeping the lowest energy), and the records written are

- ``E_bind_H_half_H2`` and ``E_bind_H_atom``: binding energy per H with the 1/2 H2 and the free
  H atom references;
- ``E_sequential``: the binding energy of the last hydrogen relative to the parent configuration,
  with both references;
- ``E_formation_defect``: the defect formation energy of the substrate relative to the pristine
  cell, when the class defines one.

Configurations whose relaxation does not converge are excluded from the table and logged.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from ase import Atoms
from ase.calculators.calculator import CalculationFailed
from ase.io import write as ase_write

from sitewindow import __version__
from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
from sitewindow.atomistic.descriptors import site_descriptors
from sitewindow.atomistic.parallel import is_master
from sitewindow.atomistic.structures.builders import add_hydrogens, set_unpaired_electrons
from sitewindow.atomistic.workflows.adsorption import (
    RelaxResult,
    binding_energies,
    reference_energies,
    relax,
    relax_tube_radius,
)
from sitewindow.io import sha256_file, write_table
from sitewindow.landscape.classes import Configuration, catalogue, chemisorption_classes
from sitewindow.schema import (
    Extraction,
    Provenance,
    ProvenanceKind,
    Quantity,
    ReferenceState,
    SiteDescriptors,
    SiteEnergyRecord,
    Tier,
    UncertainScalar,
    records_to_frame,
)

log = logging.getLogger(__name__)


@dataclass
class LandscapeSpec:
    """What to compute: tier, classes, recipe overrides and relaxation settings."""

    tier: str = "verify"
    classes: list[str] | None = None
    recipes: dict[str, dict[str, Any]] = field(default_factory=dict)
    fmax: float = 0.05
    steps: int = 400
    spin_states: tuple[int, ...] = (0,)
    numerical_sigma_eV: float = 0.05


@dataclass
class ConfigurationResult:
    config: Configuration
    config_hash: str
    energy_eV: float
    relax: RelaxResult
    atoms: Atoms
    n_unpaired: int


def _parity(atoms: Atoms) -> int:
    return int(np.sum(atoms.numbers)) % 2


def _config_hash(start: Atoms, settings_hash: str, n_unpaired: int) -> str:
    from sitewindow.atomistic.structures.builders import structure_hash

    return hashlib.sha256(
        f"{structure_hash(start)}:{settings_hash}:{n_unpaired}".encode()
    ).hexdigest()[:16]


def _relax_spin_states(
    start: Atoms, calc: Calculator, cache: Path, spec: LandscapeSpec
) -> tuple[Atoms, RelaxResult, int]:
    """Relax in each requested spin state and return the lowest converged one (or the lowest)."""
    base = _parity(start)
    best: tuple[Atoms, RelaxResult, int] | None = None
    for extra in spec.spin_states:
        n_unp = base + extra
        a0 = start.copy()
        set_unpaired_electrons(a0, n_unp)
        a, res = relax(a0, calc, cache, fmax=spec.fmax, steps=spec.steps)
        if best is None:
            best = (a, res, n_unp)
            continue
        if (res.converged and not best[1].converged) or (
            res.converged == best[1].converged and res.energy_eV < best[1].energy_eV
        ):
            best = (a, res, n_unp)
    assert best is not None
    return best


def _extxyz(atoms: Atoms) -> str:
    buf = io.StringIO()
    out = atoms.copy()
    out.info = {}
    ase_write(buf, out, format="extxyz")
    return buf.getvalue()


def _record(
    cfg: Configuration,
    config_hash: str,
    quantity: Quantity,
    value: float,
    reference: ReferenceState,
    settings_hash: str,
    tier: Tier,
    engine: str,
    converged: bool,
    descriptors: dict[str, Any] | None,
    n_atoms: int,
    sigma: float,
    n_H: int,
) -> SiteEnergyRecord:
    return SiteEnergyRecord(
        config_hash=config_hash,
        tier=tier,
        engine=engine,
        quantity=quantity,
        value=UncertainScalar(value=value, sigma=sigma, components={"numerical": sigma}),
        reference_state=reference,
        converged=converged,
        provenance=Provenance(
            kind=ProvenanceKind.computed_here,
            source_id=config_hash,
            identifier=cfg.key,
            extraction=Extraction.computed,
            settings_hash=settings_hash,
        ),
        site_class=cfg.site_class,
        site_id=cfg.site.site_id,
        n_H=n_H,
        descriptors=None if descriptors is None else SiteDescriptors(**descriptors),
        n_atoms=n_atoms,
        supercell=cfg.substrate_name,
    )


def build_landscape(
    spec: LandscapeSpec, calc: Calculator, cache_dir: Path, out_dir: Path
) -> dict[str, Any]:
    """Compute all configurations of the requested classes and write the landscape tables."""
    t_start = time.perf_counter()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    cat = catalogue()
    classes = spec.classes or chemisorption_classes()
    tier = Tier(calc.tier)
    engine = f"{calc.name}:{calc.settings().get('method', calc.settings().get('xc', ''))}".rstrip(
        ":"
    )
    settings_hash = CalculatorSpec(calc.name, calc.tier, calc.settings()).settings_hash
    sigma = spec.numerical_sigma_eV
    refs = reference_energies(calc, cache_dir)
    log.info("references E_H2=%.6f E_H=%.6f settings=%s", refs["E_H2"], refs["E_H"], settings_hash)

    records: list[SiteEnergyRecord] = []
    exclusions: list[dict[str, Any]] = []
    config_rows: list[dict[str, Any]] = []
    substrates: dict[str, tuple[Atoms, RelaxResult]] = {}
    results: dict[str, ConfigurationResult] = {}

    radius_opt: dict[str, dict[str, Any]] = {}

    def relaxed_substrate(name: str, atoms: Atoms) -> tuple[Atoms, RelaxResult]:
        if name not in substrates:
            a0 = atoms.copy()
            set_unpaired_electrons(a0, _parity(a0))
            if a0.info.get("geometry") == "tube":
                a0, rad = relax_tube_radius(a0, calc, cache_dir)
                radius_opt[name] = rad
                log.info(
                    "substrate %s: radial scale %.5f, dE=%.4f eV (%d single points)",
                    name,
                    rad["scale"],
                    rad["dE_eV"],
                    rad["n_single_points"],
                )
            a, res = relax(a0, calc, cache_dir, fmax=spec.fmax, steps=spec.steps)
            a.info.update(atoms.info)
            substrates[name] = (a, res)
        return substrates[name]

    for cls in classes:
        if cls not in cat:
            raise KeyError(f"unknown site class {cls!r}")
        defn = cat[cls]
        if "T2" not in defn.tiers:
            continue
        configs = defn.configurations(spec.tier, spec.recipes.get(cls))
        formation_done: set[str] = set()
        for cfg in configs:
            try:
                sub, sub_res = relaxed_substrate(cfg.substrate_name, cfg.substrate)
            except (CalculationFailed, RuntimeError) as exc:
                exclusions.append(
                    {
                        "key": cfg.key,
                        "reason": f"substrate calculation failed: {str(exc).splitlines()[0]}",
                        "max_force_eV_A": None,
                        "n_steps": None,
                    }
                )
                log.warning("%s excluded: %s", cfg.key, str(exc).splitlines()[0])
                continue
            if not sub_res.converged:
                exclusions.append(
                    {
                        "key": cfg.key,
                        "reason": "substrate relaxation did not converge",
                        "max_force_eV_A": sub_res.max_force_eV_A,
                        "n_steps": sub_res.n_steps,
                    }
                )
                continue
            if cfg.pristine is not None and cfg.substrate_name not in formation_done:
                pris, pris_res = relaxed_substrate(cfg.pristine_name or "pristine", cfg.pristine)
                n_c_sub = int(np.sum(sub.numbers == 6))
                n_c_pris = int(np.sum(pris.numbers == 6))
                e_form = sub_res.energy_eV - (n_c_sub / n_c_pris) * pris_res.energy_eV
                h_sub = _config_hash(sub, settings_hash, _parity(sub))
                centre = cfg.substrate.info.get("defect_center")
                site_c = int(cfg.site.carbon_indices[0])
                records.append(
                    _record(
                        cfg,
                        h_sub,
                        Quantity.E_formation_defect,
                        e_form,
                        ReferenceState.pristine_cell,
                        settings_hash,
                        tier,
                        engine,
                        sub_res.converged and pris_res.converged,
                        site_descriptors(
                            sub, site_c, None if centre is None else np.asarray(centre)
                        ),
                        len(sub),
                        sigma,
                        0,
                    )
                )
                formation_done.add(cfg.substrate_name)
            start = add_hydrogens(sub, cfg.placements)
            try:
                atoms, res, n_unp = _relax_spin_states(start, calc, cache_dir, spec)
            except (CalculationFailed, RuntimeError) as exc:
                exclusions.append(
                    {
                        "key": cfg.key,
                        "reason": f"calculation failed: {str(exc).splitlines()[0]}",
                        "max_force_eV_A": None,
                        "n_steps": None,
                    }
                )
                log.warning("%s excluded: %s", cfg.key, str(exc).splitlines()[0])
                continue
            atoms.info.update(cfg.substrate.info)
            chash = _config_hash(start, settings_hash, n_unp)
            results[cfg.key] = ConfigurationResult(cfg, chash, res.energy_eV, res, atoms, n_unp)
            if not res.converged:
                exclusions.append(
                    {
                        "key": cfg.key,
                        "reason": f"relaxation did not converge (fmax {spec.fmax} eV/A, {spec.steps} steps)",
                        "max_force_eV_A": res.max_force_eV_A,
                        "n_steps": res.n_steps,
                    }
                )
                continue
            n = cfg.n_H
            be = binding_energies(sub_res.energy_eV, res.energy_eV, refs, "H", n)
            site_c = int(cfg.placements[-1][0])
            centre = cfg.substrate.info.get("defect_center")
            desc = site_descriptors(atoms, site_c, None if centre is None else np.asarray(centre))
            common = (settings_hash, tier, engine, True, desc, len(atoms), sigma, n)
            records.append(
                _record(
                    cfg,
                    chash,
                    Quantity.E_bind_H_half_H2,
                    be["E_bind_H_half_H2"],
                    ReferenceState.half_H2,
                    *common,
                )
            )
            records.append(
                _record(
                    cfg,
                    chash,
                    Quantity.E_bind_H_atom,
                    be["E_bind_H_atom"],
                    ReferenceState.H_atom,
                    *common,
                )
            )
            if cfg.parent is not None:
                parent_key = f"{cfg.site_class.value}/{cfg.substrate_name}/{cfg.parent}"
                parent = results.get(parent_key)
                if parent is not None and parent.relax.converged:
                    d_e = res.energy_eV - parent.energy_eV
                    records.append(
                        _record(
                            cfg,
                            chash,
                            Quantity.E_sequential,
                            d_e - 0.5 * refs["E_H2"],
                            ReferenceState.half_H2,
                            *common,
                        )
                    )
                    records.append(
                        _record(
                            cfg,
                            chash,
                            Quantity.E_sequential,
                            d_e - refs["E_H"],
                            ReferenceState.H_atom,
                            *common,
                        )
                    )
            config_rows.append(
                {
                    "config_hash": chash,
                    "site_class": cfg.site_class.value,
                    "site_id": cfg.site.site_id,
                    "configuration": cfg.name,
                    "substrate": cfg.substrate_name,
                    "n_H": n,
                    "n_atoms": len(atoms),
                    "spin_multiplicity": n_unp + 1,
                    "energy_eV": res.energy_eV,
                    "max_force_eV_A": res.max_force_eV_A,
                    "n_steps": res.n_steps,
                    "wall_s": res.wall_s,
                    "cached": res.cached,
                    "structure": _extxyz(atoms),
                    "cell": json.dumps(np.asarray(atoms.cell).round(6).tolist()),
                    "pbc": json.dumps([bool(x) for x in atoms.pbc]),
                }
            )
            log.info(
                "%s n_H=%d E_b(H atom)=%.3f eV spin=%d wall=%.1fs",
                cfg.key,
                n,
                be["E_bind_H_atom"],
                n_unp + 1,
                res.wall_s,
            )

    frame = records_to_frame(records)
    out: dict[str, Any] = {
        "records": frame,
        "configurations": pd.DataFrame(config_rows),
        "excluded": pd.DataFrame(exclusions),
        "references": refs,
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        label = spec.tier
        p_rec = write_table(
            frame,
            out_dir / f"landscape_raw_{label}.parquet",
            {"tier": label, "calculator_settings_hash": settings_hash},
        )
        p_cfg = write_table(
            out["configurations"],
            out_dir / f"landscape_configurations_{label}.parquet",
            {"tier": label},
        )
        p_exc = out_dir / f"landscape_raw_{label}_excluded.csv"
        out["excluded"].to_csv(p_exc, index=False)
        manifest = {
            "tier": label,
            "package_version": __version__,
            "calculator": {
                "name": calc.name,
                "tier": calc.tier,
                "settings": calc.settings(),
                "settings_hash": settings_hash,
            },
            "classes": classes,
            "relaxation": {
                "fmax_eV_A": spec.fmax,
                "max_steps": spec.steps,
                "spin_states": list(spec.spin_states),
                "tube_radius_optimization": radius_opt,
            },
            "references_eV": refs,
            "n_records": int(len(frame)),
            "n_configurations": int(len(config_rows)),
            "n_excluded": int(len(exclusions)),
            "wall_s": time.perf_counter() - t_start,
            "files": {p.name: sha256_file(p) for p in (p_rec, p_cfg, p_exc)},
        }
        (out_dir / f"landscape_manifest_{label}.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        out["manifest"] = manifest
    return out


__all__ = ["LandscapeSpec", "ConfigurationResult", "build_landscape"]
