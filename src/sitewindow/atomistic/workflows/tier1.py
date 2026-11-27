"""Tier-1 density-functional anchors and convergence study (manual C2; benchmark B1).

The basal-monomer anchor relaxes graphene and graphene+H at Tier 2 (GFN2-xTB, when available)
for a starting geometry, re-relaxes them with GPAW in LCAO dzp mode (PBE-D3(BJ), spin-polarized),
and evaluates plane-wave single points for the energies. The PBE and D3 parts are reported
separately, so that the PBE binding energy can be compared with the plane-wave PBE value of
benchmark B1 (-0.821 eV, Dumi et al. 2022) at the registry threshold. The convergence series
evaluates the plane-wave binding energy for the supercell, k-mesh, cutoff and vacuum variants of
acceptance criterion B1 on the LCAO-relaxed geometries.

The anchor set (``run_anchor_set``) applies the same protocol to every configuration of the
registry ``data/registry/tier1_anchors.yaml``: Tier-2 pre-relaxation of the substrate (with the
radial breathing minimization for tubes) and of the hydrogenated structure under the landscape
protocol, LCAO relaxation of both, plane-wave single points, the D3 energy stored separately, the
defect formation energy where the class defines a pristine reference, and the site descriptors of
the relaxed structure. Every structure is checkpointed through the calculation cache and every
finished anchor is written as its own JSON record, so an interrupted set resumes where it stopped;
the anchor table ``anchors_tier1.parquet`` is rebuilt from the records after each anchor.

Spin state: the nominal number of unpaired electrons is the electron-count parity; every GPAW
relaxation starts from the symmetry-broken guess of ``broken_symmetry_moments`` (one Bohr magneton
on each carbon with fewer than three carbon neighbours, the sign alternating with the sublattice) with the total moment
free, the plane-wave single point starts from the moments the LCAO relaxation converged to, and
the converged moments of both stages are part of the record.

Run under MPI with ``mpirun -np N gpaw python --module sitewindow.workflows.cli tier1 ...``.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import yaml
from ase import Atoms
from ase.io import write as ase_write

from sitewindow import __version__
from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
from sitewindow.atomistic.calculators.gpaw_adapter import GPAWCalculator, dispersion_energy_eV
from sitewindow.atomistic.descriptors import site_descriptors
from sitewindow.atomistic.parallel import is_master
from sitewindow.atomistic.structures.builders import (
    BROKEN_SYMMETRY_MOMENT,
    add_hydrogens,
    broken_symmetry_moments,
    central_atom_index,
    graphene_supercell,
    hydrogen_atom,
    hydrogen_molecule,
    set_unpaired_electrons,
    structure_hash,
)
from sitewindow.atomistic.workflows.adsorption import (
    RelaxResult,
    relax,
    relax_tube_radius,
    single_point,
)
from sitewindow.errors import CalculatorUnavailable, ConvergenceError
from sitewindow.io import sha256_file, write_table
from sitewindow.landscape.classes import Configuration, catalogue
from sitewindow.schema import (
    Extraction,
    Provenance,
    ProvenanceKind,
    Quantity,
    ReferenceState,
    SiteClassName,
    SiteDescriptors,
    SiteEnergyRecord,
    Tier,
    UncertainScalar,
    records_to_frame,
)

log = logging.getLogger(__name__)

B1_PBE_EV = -0.821
B1_DMC_EV = -0.691
B1_THRESHOLD_EV = 0.10
T1_NUMERICAL_SIGMA_EV = 0.02
ANCHOR_REGISTRY = Path("data/registry/tier1_anchors.yaml")
PRERELAX_ELECTRONIC_TEMPERATURE_K = 3000.0


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


def _tier2_calculator(electronic_temperature_K: float | None = None) -> Calculator | None:
    try:
        import tblite.ase  # noqa: F401

        from sitewindow.atomistic.calculators.tblite_adapter import TBLiteCalculator

        if electronic_temperature_K is None:
            return TBLiteCalculator()
        return TBLiteCalculator(electronic_temperature=float(electronic_temperature_K))
    except Exception:
        return None


class Tier1Engine(Protocol):
    """The calculators of one Tier-1 protocol: pre-relaxation, relaxation, energies, dispersion."""

    def prerelax(self) -> Calculator | None: ...

    def lcao(self, txt: str | None = None) -> Calculator: ...

    def pw(self, txt: str | None = None) -> Calculator: ...

    def molecule_pw(self, txt: str | None = None) -> Calculator: ...

    def dispersion(self, atoms: Atoms) -> float: ...

    def name(self) -> str: ...


@dataclass(frozen=True)
class GPAWEngine:
    """The production engine: GPAW LCAO/plane-wave with D3(BJ) through ``dftd3-python`` and a
    Tier-2 (GFN2-xTB, landscape protocol) pre-relaxation when tblite is installed."""

    settings: Tier1Settings
    prerelax_temperature_K: float = PRERELAX_ELECTRONIC_TEMPERATURE_K
    first_attempt_maxiter: int = 150

    def prerelax(self) -> Calculator | None:
        if not self.settings.prerelax_tier2:
            return None
        return _tier2_calculator(self.prerelax_temperature_K)

    def lcao(self, txt: str | None = None) -> Calculator:
        return replace(self.settings.lcao(txt), maxiter=self.first_attempt_maxiter)

    def pw(self, txt: str | None = None) -> Calculator:
        return replace(self.settings.pw(txt), maxiter=self.first_attempt_maxiter)

    def molecule_pw(self, txt: str | None = None) -> Calculator:
        return replace(self.settings.molecule_pw(txt), maxiter=self.first_attempt_maxiter)

    def dispersion(self, atoms: Atoms) -> float:
        return 0.0 if self.settings.dispersion is None else dispersion_energy_eV(atoms)

    def name(self) -> str:
        return "gpaw"


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


def _molecule_lcao(engine: Tier1Engine, txt: str | None) -> Calculator:
    """The LCAO calculator of ``engine`` for an isolated molecule (Gamma point, no dipole layer)."""
    calc = engine.lcao(txt)
    if isinstance(calc, GPAWCalculator):
        return replace(calc, kmesh=(1, 1, 1), dipole_correction=False)
    return calc


def _references(
    settings: Tier1Settings,
    cache: Path,
    log_dir: Path | None,
    engine: Tier1Engine | None = None,
) -> dict[str, Any]:
    engine = engine or GPAWEngine(settings)
    h2 = hydrogen_molecule(settings.box_A, pbc=True)
    set_unpaired_electrons(h2, 0)
    lcao_mol = _molecule_lcao(engine, _log_path(log_dir, "lcao_H2"))
    h2_relaxed, _ = relax(h2, lcao_mol, cache, fmax=0.01, steps=100)
    check_h2_bond(float(h2_relaxed.get_distance(0, 1)))
    h = hydrogen_atom(settings.box_A, pbc=True)
    tag = f"{settings.ecut_eV:.0f}eV"
    sp_h2 = single_point(h2_relaxed, engine.molecule_pw(_log_path(log_dir, tag, "pw_H2")), cache)
    sp_h = single_point(h, engine.molecule_pw(_log_path(log_dir, tag, "pw_H")), cache)
    d3 = {} if settings.dispersion is None else {"H2": engine.dispersion(h2_relaxed), "H": 0.0}
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


# --------------------------------------------------------------------------- anchor set


def _parity(atoms: Atoms) -> int:
    return int(np.sum(atoms.numbers)) % 2


def _total_moment(moments: list[float] | None) -> float | None:
    return None if moments is None else float(np.sum(moments))


def _extxyz(atoms: Atoms) -> str:
    buf = io.StringIO()
    out = atoms.copy()
    out.info = {}
    ase_write(buf, out, format="extxyz")
    return buf.getvalue()


def anchor_tag(cfg: Configuration) -> str:
    """File stem of an anchor record: class, substrate and configuration name."""
    return f"{cfg.site_class.value}__{cfg.substrate_name}__{cfg.name}"


def load_anchor_registry(path: Path = ANCHOR_REGISTRY) -> dict[str, Any]:
    """The anchor registry: protocol summary and the list of (round, class, configuration, recipe)."""
    reg = yaml.safe_load(Path(path).read_text()) or {}
    if "anchors" not in reg:
        raise ValueError(f"{path} has no 'anchors' list")
    return reg


def resolve_anchors(registry: dict[str, Any]) -> list[tuple[int, Configuration]]:
    """Resolve every registry entry to its catalogue configuration, ordered by round and size."""
    cat = catalogue()
    out: list[tuple[int, Configuration]] = []
    for entry in registry["anchors"]:
        cls = entry["site_class"]
        if cls not in cat:
            raise KeyError(f"unknown site class {cls!r} in the anchor registry")
        matches = [
            c
            for c in cat[cls].configurations("full", entry.get("recipe"))
            if c.name == entry["configuration"]
        ]
        if len(matches) != 1:
            raise KeyError(
                f"anchor {cls}/{entry['configuration']} with recipe {entry.get('recipe')} "
                f"resolves to {len(matches)} configurations"
            )
        out.append((int(entry.get("round", 1)), matches[0]))
    out.sort(key=lambda rc: (rc[0], len(rc[1].substrate) + rc[1].n_H, rc[1].key))
    return out


def _relaxation_record(res: RelaxResult) -> dict[str, Any]:
    return {
        "converged": res.converged,
        "max_force_eV_A": res.max_force_eV_A,
        "n_steps": res.n_steps,
        "scf_retried": res.scf_retried,
        "cached": res.cached,
    }


def with_converged_moments(atoms: Atoms, res: RelaxResult) -> Atoms:
    """The relaxed structure with the magnetic moments its relaxation converged to as initial
    moments, so that the plane-wave single point starts in the magnetic state the LCAO
    relaxation found (ribbon edges and defects can hold several states, and a self-consistent
    field started from a one-atom guess sloshes between them)."""
    out = atoms.copy()
    if res.magmoms is not None and len(res.magmoms) == len(out):
        out.set_initial_magnetic_moments(np.asarray(res.magmoms, dtype=float))
    return out


def _prerelaxed(
    atoms: Atoms, t2: Calculator | None, cache: Path, info: dict[str, Any]
) -> tuple[Atoms, dict[str, Any] | None]:
    """Tier-2 pre-relaxation (radial breathing minimization first for tubes); the structure keeps
    the builder's metadata so that hydrogen placement and descriptors see the defect centre."""
    if t2 is None:
        return atoms, None
    a0 = atoms.copy()
    record: dict[str, Any] = {"settings": t2.settings()}
    if a0.info.get("geometry") == "tube":
        a0, rad = relax_tube_radius(a0, t2, cache)
        record["tube_radius"] = rad
    a, res = relax(a0, t2, cache, fmax=0.05, steps=400)
    a.info.update(info)
    record["relaxation"] = _relaxation_record(res)
    return a, record


def anchor_configuration(
    cfg: Configuration,
    settings: Tier1Settings,
    cache_dir: Path,
    out_dir: Path,
    log_dir: Path | None = None,
    engine: Tier1Engine | None = None,
) -> dict[str, Any]:
    """Tier-1 anchor of one catalogue configuration; writes ``anchors/<tag>.json`` in ``out_dir``.

    The binding energies use the plane-wave energies of the relaxed substrate and hydrogenated
    structure with the plane-wave references of ``_references``; PBE and PBE-D3(BJ) values are both
    recorded. When the class defines a pristine reference the defect formation energy of the
    substrate is computed the same way as in the landscape workflow.
    """
    if engine is None:
        _require_gpaw()
        engine = GPAWEngine(settings)
    t0 = time.perf_counter()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    tag = anchor_tag(cfg)
    if log_dir is not None and is_master():
        Path(log_dir).mkdir(parents=True, exist_ok=True)
    t2 = engine.prerelax()

    sub0 = cfg.substrate.copy()
    set_unpaired_electrons(sub0, _parity(sub0))
    sub0, pre_sub = _prerelaxed(sub0, t2, cache_dir, cfg.substrate.info)
    guess_sub = broken_symmetry_moments(sub0, _parity(sub0))
    sub, res_sub = relax(
        sub0,
        engine.lcao(_log_path(log_dir, tag, "lcao_substrate")),
        cache_dir,
        settings.fmax_eV_A,
        settings.relax_steps,
    )
    sub.info.update(cfg.substrate.info)

    start = add_hydrogens(sub, cfg.placements)
    n_unp = _parity(start)
    set_unpaired_electrons(start, n_unp)
    start_hash = structure_hash(start)
    start_pre, pre_ads = _prerelaxed(start, t2, cache_dir, cfg.substrate.info)
    guess_ads = broken_symmetry_moments(start_pre, n_unp)
    ads, res_ads = relax(
        start_pre,
        engine.lcao(_log_path(log_dir, tag, "lcao_configuration")),
        cache_dir,
        settings.fmax_eV_A,
        settings.relax_steps,
    )
    ads.info.update(cfg.substrate.info)

    refs = _references(settings, cache_dir, log_dir, engine)
    pw_tag = f"{tag}_{settings.ecut_eV:.0f}eV"
    sp_sub = single_point(
        with_converged_moments(sub, res_sub),
        engine.pw(_log_path(log_dir, pw_tag, "pw_substrate")),
        cache_dir,
    )
    sp_ads = single_point(
        with_converged_moments(ads, res_ads),
        engine.pw(_log_path(log_dir, pw_tag, "pw_configuration")),
        cache_dir,
    )
    d3_sub, d3_ads = engine.dispersion(sub), engine.dispersion(ads)
    d3_h2 = refs["D3"].get("H2", 0.0)
    n = cfg.n_H
    e_b_atom = (sp_ads.energy_eV - sp_sub.energy_eV - n * refs["E_H"]) / n
    e_b_half = (sp_ads.energy_eV - sp_sub.energy_eV - 0.5 * n * refs["E_H2"]) / n
    d3_atom = (d3_ads - d3_sub) / n
    d3_half = (d3_ads - d3_sub - 0.5 * n * d3_h2) / n
    e_b_atom_lcao = (res_ads.energy_eV - res_sub.energy_eV - n * refs["E_H"]) / n

    formation: dict[str, Any] | None = None
    pristine_hash: str | None = None
    if cfg.pristine is not None:
        pris0 = cfg.pristine.copy()
        set_unpaired_electrons(pris0, _parity(pris0))
        pris0, _ = _prerelaxed(pris0, t2, cache_dir, cfg.pristine.info)
        broken_symmetry_moments(pris0, _parity(pris0))
        pris, res_pris = relax(
            pris0,
            engine.lcao(_log_path(log_dir, cfg.pristine_name or "pristine", "lcao_pristine")),
            cache_dir,
            settings.fmax_eV_A,
            settings.relax_steps,
        )
        sp_pris = single_point(
            with_converged_moments(pris, res_pris),
            engine.pw(
                _log_path(
                    log_dir, f"{cfg.pristine_name or 'pristine'}_{settings.ecut_eV:.0f}eV", "pw"
                )
            ),
            cache_dir,
        )
        d3_pris = engine.dispersion(pris)
        ratio = int(np.sum(sub.numbers == 6)) / int(np.sum(pris.numbers == 6))
        e_form = sp_sub.energy_eV - ratio * sp_pris.energy_eV
        d3_form = d3_sub - ratio * d3_pris
        pristine_hash = structure_hash(pris)
        formation = {
            "pristine": cfg.pristine_name,
            "E_formation_PBE_D3": e_form,
            "E_formation_PBE": e_form - d3_form,
            "E_pw_pristine": sp_pris.energy_eV,
            "d3_pristine": d3_pris,
            "relaxation": _relaxation_record(res_pris),
            "converged": res_pris.converged,
        }

    site_c = int(cfg.placements[-1][0])
    h_index = len(ads) - 1
    centre = cfg.substrate.info.get("defect_center")
    desc = site_descriptors(ads, site_c, None if centre is None else np.asarray(centre))
    geometry: dict[str, Any] = {
        "C_H_bond_A": float(ads.get_distance(site_c, h_index, mic=True)),
        "H2_bond_A": refs["H2_bond_A"],
    }
    if cfg.substrate.info.get("geometry") == "sheet" and bool(ads.pbc[0] and ads.pbc[1]):
        carbons = ads.numbers == 6
        geometry["puckering_A"] = float(
            ads.positions[site_c, 2] - np.median(ads.positions[carbons, 2])
        )
    pw_settings = engine.pw().settings()
    settings_hash_pw = CalculatorSpec(engine.name(), "T1", pw_settings).settings_hash
    config_hash = hashlib.sha256(f"{start_hash}:{settings_hash_pw}:{n_unp}".encode()).hexdigest()[
        :16
    ]
    kmesh_used = getattr(engine.pw(), "kpoints", lambda a: (1, 1, 1))(sub)
    result: dict[str, Any] = {
        "key": cfg.key,
        "tag": tag,
        "site_class": cfg.site_class.value,
        "site_id": cfg.site.site_id,
        "configuration": cfg.name,
        "substrate": cfg.substrate_name,
        "parent": cfg.parent,
        "n_H": n,
        "n_atoms": len(ads),
        "spin_multiplicity": n_unp + 1,
        "tier": "T1",
        "engine": engine.name(),
        "package_version": __version__,
        "settings": {**asdict(settings), "kmesh_used": list(kmesh_used)},
        "settings_hash_pw": settings_hash_pw,
        "config_hash": config_hash,
        "structure_hashes": {
            "substrate": structure_hash(sub),
            "configuration": structure_hash(ads),
            "start": start_hash,
            "pristine": pristine_hash,
        },
        "energies_eV": {
            "pw": {
                "substrate": sp_sub.energy_eV,
                "configuration": sp_ads.energy_eV,
                "E_H2": refs["E_H2"],
                "E_H": refs["E_H"],
            },
            "lcao": {"substrate": res_sub.energy_eV, "configuration": res_ads.energy_eV},
            "d3": {"substrate": d3_sub, "configuration": d3_ads, "H2": d3_h2},
        },
        "binding_eV": {
            "E_bind_H_atom_PBE_D3": e_b_atom,
            "E_bind_H_atom_PBE": e_b_atom - d3_atom,
            "E_bind_H_half_H2_PBE_D3": e_b_half,
            "E_bind_H_half_H2_PBE": e_b_half - d3_half,
            "D3_contribution_atom": d3_atom,
            "E_bind_H_atom_PBE_D3_lcao": e_b_atom_lcao,
        },
        "formation_eV": formation,
        "descriptors": desc,
        "geometry": geometry,
        "relaxation": {
            "substrate": _relaxation_record(res_sub),
            "configuration": _relaxation_record(res_ads),
        },
        "magnetic_moments_lcao": {
            "substrate": _total_moment(res_sub.magmoms),
            "configuration": _total_moment(res_ads.magmoms),
        },
        "magnetic_moments_pw": {
            "substrate": _total_moment(sp_sub.magmoms),
            "configuration": _total_moment(sp_ads.magmoms),
        },
        "spin_guess": {
            "substrate_sites": guess_sub,
            "configuration_sites": guess_ads,
            "moment_muB": BROKEN_SYMMETRY_MOMENT,
        },
        "converged": bool(
            res_sub.converged
            and res_ads.converged
            and (formation is None or bool(formation["converged"]))
        ),
        "tier2_prerelaxation": {"substrate": pre_sub, "configuration": pre_ads},
        "structure": _extxyz(ads),
        "wall_s": time.perf_counter() - t0,
    }
    if is_master():
        (out_dir / "anchors").mkdir(parents=True, exist_ok=True)
        (out_dir / "anchors" / f"{tag}.json").write_text(json.dumps(result, indent=2) + "\n")
    log.info(
        "anchor %s E_b(H atom) PBE-D3 %.4f eV (PBE %.4f) converged=%s wall=%.0fs",
        cfg.key,
        e_b_atom,
        e_b_atom - d3_atom,
        result["converged"],
        result["wall_s"],
    )
    return result


def _anchor_record(
    d: dict[str, Any],
    quantity: Quantity,
    value: float,
    reference: ReferenceState,
    engine: str,
    n_H: int,
    converged: bool,
) -> SiteEnergyRecord:
    desc = d.get("descriptors")
    return SiteEnergyRecord(
        config_hash=d["config_hash"],
        tier=Tier.T1,
        engine=engine,
        quantity=quantity,
        value=UncertainScalar(
            value=value,
            sigma=T1_NUMERICAL_SIGMA_EV,
            components={"numerical": T1_NUMERICAL_SIGMA_EV},
        ),
        reference_state=reference,
        converged=converged,
        provenance=Provenance(
            kind=ProvenanceKind.computed_here,
            source_id=d["config_hash"],
            identifier=d["key"],
            extraction=Extraction.computed,
            settings_hash=d["settings_hash_pw"],
        ),
        site_class=SiteClassName(d["site_class"]),
        site_id=d["site_id"],
        n_H=n_H,
        descriptors=None if desc is None else SiteDescriptors(**desc),
        n_atoms=d["n_atoms"],
        supercell=d["substrate"],
    )


ENGINE_LABELS = {"PBE_D3": "gpaw:PBE-D3(BJ)", "PBE": "gpaw:PBE"}


def collect_anchor_table(out_dir: Path) -> dict[str, Any]:
    """Rebuild ``anchors_tier1.parquet`` and its manifest from every anchor record in ``out_dir``.

    Records are written for both the PBE-D3(BJ) (the Tier-1 protocol) and the PBE-only energies,
    distinguished by the engine label, so that the PBE values can be compared with published PBE
    anchors while the calibration of C3 uses the dispersion-corrected tier. Sequential binding
    energies are formed when the parent configuration's anchor exists.
    """
    out_dir = Path(out_dir)
    files = sorted((out_dir / "anchors").glob("*.json")) if (out_dir / "anchors").exists() else []
    anchors = [json.loads(p.read_text()) for p in files]
    by_key = {d["key"]: d for d in anchors}
    records: list[SiteEnergyRecord] = []
    for d in anchors:
        n = int(d["n_H"])
        conv = bool(d["converged"])
        b = d["binding_eV"]
        for suffix, engine in ENGINE_LABELS.items():
            records.append(
                _anchor_record(
                    d,
                    Quantity.E_bind_H_half_H2,
                    b[f"E_bind_H_half_H2_{suffix}"],
                    ReferenceState.half_H2,
                    engine,
                    n,
                    conv,
                )
            )
            records.append(
                _anchor_record(
                    d,
                    Quantity.E_bind_H_atom,
                    b[f"E_bind_H_atom_{suffix}"],
                    ReferenceState.H_atom,
                    engine,
                    n,
                    conv,
                )
            )
            form = d.get("formation_eV")
            if form is not None:
                records.append(
                    _anchor_record(
                        d,
                        Quantity.E_formation_defect,
                        form[f"E_formation_{suffix}"],
                        ReferenceState.pristine_cell,
                        engine,
                        0,
                        conv and bool(form["converged"]),
                    )
                )
            parent = d.get("parent")
            if parent is not None:
                pkey = f"{d['site_class']}/{d['substrate']}/{parent}"
                pd_ = by_key.get(pkey)
                if pd_ is not None:
                    e_cfg = d["energies_eV"]["pw"]["configuration"]
                    e_par = pd_["energies_eV"]["pw"]["configuration"]
                    d_e = e_cfg - e_par
                    if suffix == "PBE":
                        d_e -= d["energies_eV"]["d3"]["configuration"]
                        d_e += pd_["energies_eV"]["d3"]["configuration"]
                    d3_h2 = d["energies_eV"]["d3"]["H2"] if suffix == "PBE" else 0.0
                    records.append(
                        _anchor_record(
                            d,
                            Quantity.E_sequential,
                            d_e - 0.5 * (d["energies_eV"]["pw"]["E_H2"] - d3_h2),
                            ReferenceState.half_H2,
                            engine,
                            n,
                            conv and bool(pd_["converged"]),
                        )
                    )
                    records.append(
                        _anchor_record(
                            d,
                            Quantity.E_sequential,
                            d_e - d["energies_eV"]["pw"]["E_H"],
                            ReferenceState.H_atom,
                            engine,
                            n,
                            conv and bool(pd_["converged"]),
                        )
                    )
    frame = records_to_frame(records)
    manifest: dict[str, Any] = {
        "tier": "T1",
        "package_version": __version__,
        "n_anchors": len(anchors),
        "n_records": int(len(frame)),
        "classes": sorted({d["site_class"] for d in anchors}),
        "anchors": [
            {
                "key": d["key"],
                "n_atoms": d["n_atoms"],
                "converged": d["converged"],
                "E_bind_H_atom_PBE_D3": d["binding_eV"]["E_bind_H_atom_PBE_D3"],
                "wall_s": d["wall_s"],
            }
            for d in anchors
        ],
        "wall_s_total": float(sum(d["wall_s"] for d in anchors)),
        "failed": [
            json.loads(p.read_text()) for p in sorted((out_dir / "anchors_failed").glob("*.json"))
        ]
        if (out_dir / "anchors_failed").exists()
        else [],
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        p = write_table(frame, out_dir / "anchors_tier1.parquet", {"tier": "T1"})
        manifest["files"] = {p.name: sha256_file(p)}
        (out_dir / "anchors_manifest_tier1.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"records": frame, "manifest": manifest}


def run_anchor_set(
    settings: Tier1Settings,
    cache_dir: Path,
    out_dir: Path,
    registry_path: Path = ANCHOR_REGISTRY,
    rounds: list[int] | None = None,
    keys: list[str] | None = None,
    log_dir: Path | None = None,
    engine: Tier1Engine | None = None,
    dry_run: bool = False,
    retry_failed: bool = False,
) -> dict[str, Any]:
    """Run the anchor registry in order, skipping anchors whose record exists, rebuilding the
    anchor table after each one; returns the plan with the status of every anchor.

    An anchor whose self-consistent field or relaxation fails is recorded under
    ``anchors_failed/`` with the reason and the set continues; failed anchors are skipped on
    later runs unless ``retry_failed`` is set.
    """
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    plan: list[dict[str, Any]] = []
    for rnd, cfg in resolve_anchors(load_anchor_registry(registry_path)):
        selected = (rounds is None or rnd in rounds) and (keys is None or cfg.key in keys)
        record = out_dir / "anchors" / f"{anchor_tag(cfg)}.json"
        failed = out_dir / "anchors_failed" / f"{anchor_tag(cfg)}.json"
        status = "done" if record.exists() else ("failed" if failed.exists() else "pending")
        plan.append(
            {
                "round": rnd,
                "key": cfg.key,
                "n_atoms": len(cfg.substrate) + cfg.n_H,
                "selected": selected,
                "status": status,
                "record": str(record),
            }
        )
    if dry_run:
        return {"plan": plan}
    if engine is None:
        _require_gpaw()
        engine = GPAWEngine(settings)
    for item, (_, cfg) in zip(plan, resolve_anchors(load_anchor_registry(registry_path))):
        if not item["selected"] or item["status"] == "done":
            continue
        if item["status"] == "failed" and not retry_failed:
            continue
        log.info("anchor %s (round %d, %d atoms) starting", cfg.key, item["round"], item["n_atoms"])
        t0 = time.perf_counter()
        try:
            res = anchor_configuration(cfg, settings, cache_dir, out_dir, log_dir, engine)
        except (ConvergenceError, RuntimeError) as exc:
            reason = str(exc).splitlines()[0]
            log.error("anchor %s failed: %s", cfg.key, reason)
            item["status"] = "failed"
            item["reason"] = reason
            if is_master():
                (out_dir / "anchors_failed").mkdir(parents=True, exist_ok=True)
                (out_dir / "anchors_failed" / f"{anchor_tag(cfg)}.json").write_text(
                    json.dumps(
                        {"key": cfg.key, "reason": reason, "wall_s": time.perf_counter() - t0},
                        indent=2,
                    )
                    + "\n"
                )
            continue
        item["status"] = "done"
        item["wall_s"] = res["wall_s"]
        item["E_bind_H_atom_PBE_D3"] = res["binding_eV"]["E_bind_H_atom_PBE_D3"]
        collect_anchor_table(out_dir)
    return {"plan": plan, **collect_anchor_table(out_dir)}


__all__ = [
    "Tier1Settings",
    "Tier1Engine",
    "GPAWEngine",
    "check_h2_bond",
    "anchor_basal_monomer",
    "anchor_configuration",
    "anchor_tag",
    "load_anchor_registry",
    "resolve_anchors",
    "collect_anchor_table",
    "run_anchor_set",
    "with_converged_moments",
    "convergence_variants",
    "convergence_series",
    "B1_PBE_EV",
    "B1_DMC_EV",
    "B1_THRESHOLD_EV",
    "T1_NUMERICAL_SIGMA_EV",
    "ANCHOR_REGISTRY",
]
