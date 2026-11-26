"""Charging and release barriers along designated pathways (manual C2; benchmark B2).

Two methods are implemented.

- ``constrained_scan``: a reaction coordinate is held fixed while every other degree of freedom
  relaxes, on a 0.1 A grid refined to 0.02 A around the maximum (manual: Tier 2, and the Tier-1
  anchors under MPI). The coordinate is a carbon-hydrogen distance (fixed bond length): of the
  adsorbed hydrogen for sticking, and of the first hydrogen of the pair for recombination, so
  that the leaving hydrogen can abstract its partner while every other degree of freedom relaxes.
- ``ci_neb``: climbing-image nudged elastic band with seven images and a spring constant of
  0.1 eV/A between the relaxed initial state and the relaxed product state: the hydrogen atom
  held far from the surface for sticking; for recombination the relaxed end of the scan when the
  scan formed H2, and otherwise an H2 molecule relaxed above the pair, with its centre of mass
  held when a free relaxation lets it fall back onto the surface (the vacancy hole); a band whose
  final state is not the product is not run. The band starts from the
  geometries of the constrained scan when that scan reached the product (an IDPP interpolation
  otherwise), is first relaxed without the climbing image to a loose tolerance and then with it
  to the tolerance of the tier, and is resumable from its trajectory (manual: Tier 3; run at
  Tier 2 as well, where it is cheap).

Every method result is checked against the product. The sticking product is the hydrogen atom
separated from the surface (final carbon-hydrogen distance within 1 A of the far end of the grid);
the recombination product is H2 (final hydrogen-hydrogen distance within a factor
``PRODUCT_H_H_FACTOR`` of the engine's own H2 bond length). A method has located the saddle only
when every point or image converged, the path ends in the product, the maximum lies in the
interior of the path, the constrained path is continuous (no energy step steeper than
``SCAN_SLOPE_LIMIT_EV_A``, about twice the largest restoring force of a carbon-hydrogen bond,
between adjacent points: a steeper step means the constrained minimum jumped between valleys and
the maximum before the jump is not a saddle) and, for sticking, the far end of the path lies
within ``FAR_POINT_TOLERANCE_EV`` of the separated reactants (the far point reproduces them up to
physisorption, so a larger offset means the reference is inconsistent with the path, as happens
when Fermi smearing shares electrons between a distant hydrogen atom and the sheet); otherwise
``saddle_located`` is false and ``saddle_caveat`` names the reason. A single carbon-hydrogen
coordinate does not reach H2 from every initial state (a para pair, or the two dangling-bond
hydrogens of a monovacancy, only detach one hydrogen), and a scan that ends in the wrong product
is not a bound on the recombination barrier in either direction.

The table ``barriers_<tier>.parquet`` holds one record per pathway from the primary method: the
constrained scan for sticking, whose single coordinate follows the minimum-energy path (the NEB
confirms it), and the NEB for recombination, where a constrained scan can miss the saddle in
either direction (the constrained minimum can jump between valleys, or follow a higher valley
than the true path). When the preferred method has not located the saddle the other method is
used if it has; when neither has, the scan maximum is written and the record carries
``converged = False``, so that a record with ``converged = True`` is a saddle-point energy and
one without is the maximum of a path that stopped early, ended in the wrong product, or rose to
the end of the grid (manual: failure modes). At Tier 2 the GFN2-xTB self-consistent field does
not converge at the hydrogen-abstraction saddle, so recombination records there are of the
second kind. The profile file keeps every method with its energies, flags and structures, and
for sticking the energy of the far end of the grid relative to the separated reactants (the
``far_point_offset_eV``), which is zero up to physisorption when the reference is consistent.

Pathways (manual C2, governing principles) are taken from the site-class catalogue so that the
initial states are the landscape configurations themselves: monomer sticking on the basal plane,
ortho- and para-dimer recombination, recombination of the two dangling-bond hydrogens of a
monovacancy, and ortho-dimer recombination on the exterior of armchair nanotubes at two diameters.

Conventions (manual Section 4 notation): the charging barrier E_a of sticking is measured from
the separated reactants, the relaxed substrate plus a free H atom (reference state ``H_atom``);
the release barrier E_d of recombination is measured from the relaxed adsorbed initial state
(reference state ``adsorbed_state``). A sticking barrier below zero means that no barrier exists
along the coordinate; the value is kept and flagged in the profile.
"""

from __future__ import annotations

import functools
import hashlib
import io
import json
import logging
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, is_dataclass, replace
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from ase import Atoms
from ase.calculators.calculator import CalculationFailed
from ase.constraints import FixBondLength, FixSubsetCom
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.mep import NEB
from ase.optimize import FIRE

from sitewindow import __version__
from sitewindow.atomistic.calculators.base import Calculator, CalculatorSpec
from sitewindow.atomistic.descriptors import site_descriptors
from sitewindow.atomistic.parallel import is_master, size
from sitewindow.atomistic.structures.builders import (
    add_hydrogens,
    set_unpaired_electrons,
    structure_hash,
    surface_normal,
)
from sitewindow.atomistic.workflows.adsorption import (
    RelaxResult,
    reference_energies,
    relax,
    relax_tube_radius,
)
from sitewindow.io import sha256_file, write_table
from sitewindow.landscape.classes import Configuration, catalogue
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

STICK = "stick"
RECOMB = "recomb"
CONSTRAINED_SCAN = "constrained_scan"
CI_NEB = "ci_neb"
METHODS = (CONSTRAINED_SCAN, CI_NEB)
H2_BOND_A = 0.75
PRODUCT_H_H_FACTOR = 1.3
SCAN_SLOPE_LIMIT_EV_A = 10.0
FAR_POINT_TOLERANCE_EV = 0.1


@dataclass(frozen=True)
class Pathway:
    """A barrier pathway: the catalogue configuration of its initial state and its coordinate.

    ``coordinate`` holds the carbon and hydrogen indices of the constrained bond in the
    hydrogenated structure (the adsorbed hydrogen for sticking, the first hydrogen of the pair
    for recombination); ``hydrogens`` lists the hydrogen indices of the pathway.
    """

    name: str
    kind: str
    config: Configuration
    coordinate: tuple[int, int]
    tiers: tuple[str, ...]

    @property
    def key(self) -> str:
        return f"{self.config.site_class.value}/{self.config.substrate_name}/{self.name}"

    @property
    def quantity(self) -> Quantity:
        return Quantity.E_barrier_stick if self.kind == STICK else Quantity.E_barrier_recomb

    @property
    def reference_state(self) -> ReferenceState:
        return ReferenceState.H_atom if self.kind == STICK else ReferenceState.adsorbed_state

    @property
    def site_carbons(self) -> tuple[int, ...]:
        return tuple(int(c) for c, _ in self.config.placements)

    @property
    def hydrogens(self) -> tuple[int, ...]:
        n_sub = len(self.config.substrate)
        return tuple(range(n_sub, n_sub + self.config.n_H))

    @property
    def side(self) -> str:
        return str(self.config.placements[0][1])


@dataclass
class BarrierSpec:
    """What to compute: tier, pathways, methods, grids and relaxation settings."""

    tier: str = "verify"
    pathways: list[str] | None = None
    recipes: dict[str, dict[str, Any]] = field(default_factory=dict)
    methods: tuple[str, ...] = (CONSTRAINED_SCAN, CI_NEB)
    fmax: float = 0.05
    steps: int = 400
    grid_step_A: float = 0.1
    refine_step_A: float = 0.02
    refine_halfwidth_A: float = 0.1
    far_distance_A: float = 3.6
    recomb_distance_A: float = 3.0
    h2_height_A: float = 3.0
    h2_far_height_A: float = 6.0
    neb_images: int = 7
    neb_spring_eV_A: float = 0.1
    neb_steps: int = 600
    neb_prerelax_fmax: float = 0.3
    numerical_sigma_eV: float = 0.05
    tube_diameters_nm: tuple[float, float] = (0.68, 1.36)
    stick_recipes: dict[str, dict[str, Any]] = field(default_factory=dict)
    stick_calculator_settings: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ScanPoint:
    coordinate_A: float
    energy_eV: float
    converged: bool
    atoms: Atoms | None
    excluded: str | None = None

    @property
    def accepted(self) -> bool:
        """A finite point that is not excluded from the constrained path."""
        return bool(np.isfinite(self.energy_eV)) and self.excluded is None


@dataclass(frozen=True)
class NEBResult:
    energies_eV: tuple[float, ...]
    forward_barrier_eV: float
    reverse_barrier_eV: float
    converged: bool
    n_steps: int
    images: tuple[Atoms, ...]
    wall_s: float
    cached: bool = False


# --------------------------------------------------------------------------- pathway catalogue


def _nearest(values: Sequence[float], target: float) -> float:
    return min(values, key=lambda v: abs(float(v) - target))


def pathways(spec: BarrierSpec) -> list[Pathway]:
    """The designated pathways of the manual for ``spec.tier``, built from the catalogue."""
    cat = catalogue()
    out: list[Pathway] = []

    def find(cls: str, name: str, override: dict[str, Any] | None = None) -> Configuration | None:
        recipe = dict(spec.recipes.get(cls, {}))
        if override:
            recipe.update(override)
        for cfg in cat[cls].configurations(spec.tier, recipe or None):
            if cfg.name == name:
                return cfg
        return None

    def add(name: str, kind: str, cfg: Configuration | None) -> None:
        if cfg is None:
            return
        coord = (int(cfg.placements[0][0]), len(cfg.substrate))
        out.append(Pathway(name, kind, cfg, coord, cfg.tiers))

    add("monomer_stick", STICK, find("basal", "monomer", spec.stick_recipes.get("basal")))
    add("ortho_recomb", RECOMB, find("basal", "ortho"))
    add("para_recomb", RECOMB, find("basal", "para"))
    add("dangling_pair_recomb", RECOMB, find("monovacancy", "seq2"))
    cnt_recipe = cat["cnt_armchair"].recipe(spec.tier, spec.recipes.get("cnt_armchair"))
    diameters = sorted({_nearest(cnt_recipe["diameters_nm"], d) for d in spec.tube_diameters_nm})
    for d in diameters:
        for cfg in cat["cnt_armchair"].configurations(
            spec.tier, {**cnt_recipe, "diameters_nm": [d]}
        ):
            if cfg.name == "ortho_exterior":
                add("ortho_exterior_recomb", RECOMB, cfg)
    if spec.pathways:
        wanted = set(spec.pathways)
        out = [p for p in out if p.key in wanted or p.name in wanted]
    return out


# --------------------------------------------------------------------------- coordinates


def set_coordinate(
    atoms: Atoms, coordinate: tuple[int, int], distance_A: float, kind: str = STICK
) -> Atoms:
    """Return a copy with the coordinate pair at ``distance_A``.

    For sticking only the second atom (the hydrogen) moves along the pair direction; otherwise
    both atoms move symmetrically about their midpoint.
    """
    a = atoms.copy()
    i, j = coordinate
    v = a.get_distance(i, j, mic=True, vector=True)
    u = v / np.linalg.norm(v)
    if kind == STICK:
        a.positions[j] = a.positions[i] + u * distance_A
    else:
        mid = a.positions[i] + 0.5 * v
        a.positions[i] = mid - 0.5 * u * distance_A
        a.positions[j] = mid + 0.5 * u * distance_A
    return a


def h2_above(
    atoms: Atoms,
    hydrogens: Sequence[int],
    carbons: Sequence[int],
    normal: np.ndarray,
    height_A: float,
    bond_A: float = H2_BOND_A,
) -> Atoms:
    """Return a copy in which the two hydrogens form an H2 molecule ``height_A`` above the
    midpoint of ``carbons`` along ``normal``, oriented parallel to the surface."""
    a = atoms.copy()
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    c = list(carbons)
    if len(c) >= 2:
        axis = a.get_distance(c[0], c[1], mic=True, vector=True)
        mid = a.positions[c[0]] + 0.5 * axis
    else:
        axis = np.array([1.0, 0.0, 0.0])
        mid = a.positions[c[0]].copy()
    axis = axis - np.dot(axis, n) * n
    if np.linalg.norm(axis) < 1e-6:
        axis = np.cross(n, [0.0, 0.0, 1.0])
        if np.linalg.norm(axis) < 1e-6:
            axis = np.cross(n, [1.0, 0.0, 0.0])
    axis /= np.linalg.norm(axis)
    centre = mid + height_A * n
    h1, h2 = hydrogens
    a.positions[h1] = centre - 0.5 * bond_A * axis
    a.positions[h2] = centre + 0.5 * bond_A * axis
    return a


def resample_path(geometries: Sequence[Atoms], n_images: int) -> list[Atoms]:
    """``n_images`` geometries at equal fractions of the cumulative displacement along an ordered
    sequence of structures, interpolating positions linearly between neighbours; the first and
    last geometries are kept."""
    if len(geometries) < 2:
        raise ValueError("a path needs at least two geometries")
    pos = [g.positions.copy() for g in geometries]
    steps = [float(np.sqrt(((b - a) ** 2).sum())) for a, b in zip(pos, pos[1:])]
    cum = np.concatenate([[0.0], np.cumsum(steps)])
    if cum[-1] <= 0.0:
        return [geometries[0].copy() for _ in range(n_images - 1)] + [geometries[-1].copy()]
    targets = np.linspace(0.0, cum[-1], n_images)
    out: list[Atoms] = []
    for t in targets:
        k = int(np.clip(np.searchsorted(cum, t, side="right") - 1, 0, len(steps) - 1))
        frac = 0.0 if steps[k] == 0.0 else (t - cum[k]) / steps[k]
        img = geometries[0].copy()
        img.positions = (1.0 - frac) * pos[k] + frac * pos[k + 1]
        out.append(img)
    out[0].positions = pos[0]
    out[-1].positions = pos[-1]
    return out


def scan_profile(
    evaluate: Callable[[float, ScanPoint | None], ScanPoint],
    grid: Sequence[float],
    refine_step: float = 0.02,
    refine_halfwidth: float = 0.1,
) -> list[ScanPoint]:
    """Evaluate the coarse grid in order, each point starting from the previous one, then refine
    around the highest converged point; returns all points sorted by coordinate."""
    points: list[ScanPoint] = []
    prev: ScanPoint | None = None
    for d in grid:
        point = evaluate(float(d), prev)
        points.append(point)
        if point.atoms is not None:
            prev = point
    points = mark_islands(points)
    accepted = [p for p in points if p.accepted]
    if not accepted:
        return points
    converged = [p for p in accepted if p.converged] or accepted
    top = max(converged, key=lambda p: p.energy_eV)
    lo, hi = min(grid), max(grid)
    fine = np.arange(
        top.coordinate_A - refine_halfwidth,
        top.coordinate_A + refine_halfwidth + 1e-9,
        refine_step,
    )
    for d in fine:
        if d < lo - 1e-9 or d > hi + 1e-9:
            continue
        if any(abs(d - p.coordinate_A) < 1e-6 for p in points):
            continue
        points.append(evaluate(float(round(d, 4)), top))
    return mark_islands(points)


ISLAND_REASON = "self-consistent field solution above the ground state of the neighbouring points"


def mark_islands(
    points: Sequence[ScanPoint], slope_limit_eV_A: float = SCAN_SLOPE_LIMIT_EV_A
) -> list[ScanPoint]:
    """Mark the islands of a scan: maximal runs of finite points that are separated from both
    of their neighbours by energy steps steeper than ``slope_limit_eV_A`` and lie entirely above
    both neighbours. Such a run is a self-consistent field solution that is not the ground state
    of the constrained system, since lower states exist at the adjacent coordinates; its points
    are kept with their energies but excluded from the path. Earlier marks are recomputed, and
    the result is sorted by coordinate."""
    ordered = sorted((replace(p, excluded=None) for p in points), key=lambda p: p.coordinate_A)
    finite = [k for k, p in enumerate(ordered) if np.isfinite(p.energy_eV)]
    if len(finite) < 3:
        return ordered
    segments: list[list[int]] = [[finite[0]]]
    for a, b in zip(finite, finite[1:]):
        pa, pb = ordered[a], ordered[b]
        gap = abs(pb.coordinate_A - pa.coordinate_A)
        if gap > 0 and abs(pb.energy_eV - pa.energy_eV) / gap > slope_limit_eV_A:
            segments.append([b])
        else:
            segments[-1].append(b)
    for n in range(1, len(segments) - 1):
        seg = segments[n]
        left = ordered[segments[n - 1][-1]].energy_eV
        right = ordered[segments[n + 1][0]].energy_eV
        if min(ordered[k].energy_eV for k in seg) > max(left, right):
            for k in seg:
                ordered[k] = replace(ordered[k], excluded=ISLAND_REASON)
    return ordered


def _grid(start: float, end: float, step: float) -> list[float]:
    if end >= start:
        first = np.ceil(start / step) * step
        vals = np.arange(first, end + 1e-9, step)
    else:
        first = np.floor(start / step) * step
        vals = np.arange(first, end - 1e-9, -step)
    return [float(round(v, 4)) for v in vals if abs(v - start) > 1e-6]


# --------------------------------------------------------------------------- climbing-image NEB


def neb_barrier(
    initial: Atoms,
    final: Atoms,
    calc: Calculator,
    n_images: int = 7,
    spring_eV_A: float = 0.1,
    fmax: float = 0.05,
    steps: int = 300,
    trajectory: Path | None = None,
    climb: bool = True,
    prerelax_fmax: float = 0.3,
    seed: Sequence[Atoms] | None = None,
) -> NEBResult:
    """Climbing-image NEB between two relaxed states.

    With ``trajectory`` the band is checkpointed: a converged band recorded in the sidecar JSON
    next to the trajectory is returned without recomputation; an unconverged band resumes from
    its last images with a fresh step budget until it has taken twice the budget, after which it
    is returned as it stands, and the recorded step count is cumulative. A ``seed`` (an ordered
    sequence of geometries from ``initial`` to
    ``final``, for example the constrained scan) replaces the IDPP interpolation of a new band so
    that the band starts in the valley the scan found. Not available under MPI (the Tier-1
    anchors use the constrained scan).
    """
    if size() > 1:
        raise RuntimeError("NEB runs in a single process; use the constrained scan under MPI")
    images = [initial.copy() for _ in range(n_images - 1)] + [final.copy()]
    for img in images:
        img.set_constraint()
    sidecar = None if trajectory is None else Path(trajectory).with_suffix(".json")
    saved: list[Atoms] = []
    if trajectory is not None and Path(trajectory).exists():
        try:
            saved = list(ase_read(f"{trajectory}@-{n_images}:"))
        except Exception:  # a truncated file restarts the band
            saved = []
        if len(saved) != n_images or not (
            np.allclose(saved[0].positions, initial.positions, atol=1e-3)
            and np.allclose(saved[-1].positions, final.positions, atol=1e-3)
        ):
            saved = []  # a band between other end points is not resumed
    prior_steps = 0
    if saved and sidecar is not None and sidecar.exists():
        meta = json.loads(sidecar.read_text())
        prior_steps = int(meta.get("n_steps", 0))
        complete = len(meta.get("energies_eV", [])) == n_images
        if complete and (meta.get("converged") or prior_steps >= 2 * steps):
            for img, s in zip(images, saved):
                img.positions = s.positions
            energies = tuple(float("nan") if e is None else float(e) for e in meta["energies_eV"])
            top = float(np.nanmax(energies)) if np.isfinite(energies).any() else float("nan")
            return NEBResult(
                energies,
                top - energies[0],
                top - energies[-1],
                bool(meta.get("converged")),
                prior_steps,
                tuple(images),
                0.0,
                True,
            )
    for img in images:
        calc.attach(img)
    neb = NEB(
        images,
        k=spring_eV_A,
        climb=climb,
        method="improvedtangent",
        allow_shared_calculator=False,
    )
    if saved:
        for img, s in zip(images, saved):
            img.positions = s.positions
    elif seed is not None and len(seed) >= 2:
        path = resample_path([initial, *seed, final], n_images)
        for img, s in zip(images[1:-1], path[1:-1]):
            img.positions = s.positions
    else:
        neb.interpolate(method="idpp", mic=bool(initial.pbc.any()))
    t0 = time.perf_counter()
    traj_arg = None if trajectory is None else str(trajectory)
    n_steps = 0
    if climb and not saved:
        neb.climb = False
        pre = FIRE(neb, trajectory=traj_arg, logfile=None)
        try:
            pre.run(fmax=prerelax_fmax, steps=steps // 2)
        except (CalculationFailed, RuntimeError) as exc:
            log.warning("NEB pre-relaxation stopped: %s", str(exc).splitlines()[0])
        n_steps += int(pre.get_number_of_steps())
        neb.climb = True
    opt = FIRE(neb, trajectory=traj_arg, logfile=None)
    try:
        converged = bool(opt.run(fmax=fmax, steps=steps - n_steps))
    except (CalculationFailed, RuntimeError) as exc:
        log.warning("NEB stopped: %s", str(exc).splitlines()[0])
        converged = False
    n_steps += int(opt.get_number_of_steps()) + prior_steps
    energies_list: list[float] = []
    for img in images:
        try:
            energies_list.append(float(img.get_potential_energy()))
        except (CalculationFailed, RuntimeError):
            energies_list.append(float("nan"))
            converged = False
    energies = tuple(energies_list)
    for img in images:
        img.calc = None
    if sidecar is not None and is_master():
        sidecar.write_text(
            json.dumps(
                {
                    "energies_eV": list(energies),
                    "converged": converged,
                    "n_steps": n_steps,
                    "fmax_eV_A": fmax,
                    "spring_eV_A": spring_eV_A,
                }
            )
        )
    top = float(np.nanmax(energies)) if np.isfinite(energies).any() else float("nan")
    return NEBResult(
        energies,
        top - energies[0],
        top - energies[-1],
        converged,
        n_steps,
        tuple(images),
        time.perf_counter() - t0,
    )


# --------------------------------------------------------------------------- records


def _extxyz(atoms: Atoms) -> str:
    buf = io.StringIO()
    out = atoms.copy()
    out.info = {}
    out.set_constraint()
    ase_write(buf, out, format="extxyz")
    return buf.getvalue()


def _hash(*parts: str) -> str:
    return hashlib.sha256(":".join(parts).encode()).hexdigest()[:16]


def _record(
    pw: Pathway,
    config_hash: str,
    method: str,
    value: float,
    settings_hash: str,
    tier: Tier,
    engine: str,
    converged: bool,
    descriptors: dict[str, Any],
    n_atoms: int,
    sigma: float,
) -> SiteEnergyRecord:
    cfg = pw.config
    return SiteEnergyRecord(
        config_hash=config_hash,
        tier=tier,
        engine=engine,
        quantity=pw.quantity,
        value=UncertainScalar(value=value, sigma=sigma, components={"numerical": sigma}),
        reference_state=pw.reference_state,
        converged=converged,
        provenance=Provenance(
            kind=ProvenanceKind.computed_here,
            source_id=config_hash,
            identifier=f"{pw.key}@{method}",
            extraction=Extraction.computed,
            settings_hash=settings_hash,
        ),
        site_class=cfg.site_class,
        site_id=cfg.site.site_id,
        n_H=cfg.n_H,
        descriptors=SiteDescriptors(**descriptors),
        n_atoms=n_atoms,
        supercell=cfg.substrate_name,
    )


def _located(result: dict[str, Any]) -> bool:
    return bool(result.get("saddle_located", result.get("converged", False)))


def primary_method(kind: str, results: dict[str, dict[str, Any]]) -> str | None:
    """The method whose barrier enters the table: the constrained scan for sticking, the NEB for
    recombination, falling back to the other method when the preferred one is absent or has not
    located the saddle. When neither method has located it the scan is reported before the NEB,
    because the constrained minima are meaningful points of the surface while an unconverged
    band is not."""
    order = [CONSTRAINED_SCAN, CI_NEB] if kind == STICK else [CI_NEB, CONSTRAINED_SCAN]
    usable = {
        m: r
        for m, r in results.items()
        if r.get("barrier_eV") is not None and np.isfinite(r["barrier_eV"])
    }
    for m in order:
        if m in usable and _located(usable[m]):
            return m
    for m in (CONSTRAINED_SCAN, CI_NEB):
        if m in usable:
            return m
    return None


def product_distance(
    kind: str,
    atoms: Atoms,
    coordinate: tuple[int, int],
    hydrogens: Sequence[int],
    far_distance_A: float,
    h2_bond_A: float,
) -> tuple[float, bool]:
    """The distance that defines the product of a pathway and whether ``atoms`` is that product:
    the carbon-hydrogen distance of the scanned bond for sticking (product: within 1 A of the far
    end of the grid), the hydrogen-hydrogen distance for recombination (product: H2, within
    ``PRODUCT_H_H_FACTOR`` of the engine's bond length)."""
    if kind == STICK:
        d = float(atoms.get_distance(*coordinate, mic=True))
        return d, bool(d >= far_distance_A - 1.0)
    h1, h2 = hydrogens
    d = float(atoms.get_distance(h1, h2, mic=True))
    return d, bool(d <= PRODUCT_H_H_FACTOR * h2_bond_A)


def saddle_caveat(
    kind: str,
    converged: bool,
    product_reached: bool,
    maximum_at_end: bool,
    final_distance_A: float | None,
    method: str = CONSTRAINED_SCAN,
    discontinuity: tuple[float, float, float] | None = None,
    far_offset_eV: float | None = None,
) -> str | None:
    """Why a method result is not a saddle-point energy, or ``None`` when it is.
    ``discontinuity`` is the largest energy step of a scan with the coordinates it joins, and
    ``far_offset_eV`` the energy of the far end of a sticking path relative to the separated
    reactants."""
    reasons = []
    if not converged:
        reasons.append(
            "not every point converged"
            if method == CONSTRAINED_SCAN
            else "the band did not converge"
        )
    if not product_reached:
        label = "C-H" if kind == STICK else "H-H"
        where = "" if final_distance_A is None else f" ({label} distance {final_distance_A:.2f} A)"
        reasons.append(f"the path does not end in the product{where}")
    if maximum_at_end:
        reasons.append(
            "the maximum is at the end of the grid"
            if method == CONSTRAINED_SCAN
            else "the highest image is an end point"
        )
    if discontinuity is not None:
        jump, d1, d2 = discontinuity
        if abs(d2 - d1) > 0 and jump / abs(d2 - d1) > SCAN_SLOPE_LIMIT_EV_A:
            reasons.append(
                f"the constrained path is discontinuous ({jump:.2f} eV between {d1:.2f} and {d2:.2f} A)"
            )
    if far_offset_eV is not None and abs(far_offset_eV) > FAR_POINT_TOLERANCE_EV:
        reasons.append(
            f"the far end of the path lies {far_offset_eV:+.2f} eV from the separated reactants"
        )
    return "; ".join(reasons) or None


def largest_step(points: Sequence[ScanPoint]) -> tuple[float, float, float] | None:
    """The largest energy difference between adjacent finite points of a scan, with the two
    coordinates it joins; ``None`` with fewer than two finite points."""
    finite = sorted((p for p in points if np.isfinite(p.energy_eV)), key=lambda p: p.coordinate_A)
    if len(finite) < 2:
        return None
    best = max(zip(finite, finite[1:]), key=lambda pq: abs(pq[1].energy_eV - pq[0].energy_eV))
    return (
        float(abs(best[1].energy_eV - best[0].energy_eV)),
        float(best[0].coordinate_A),
        float(best[1].coordinate_A),
    )


# --------------------------------------------------------------------------- workflow


def _bond_evaluator(
    init: Atoms,
    coordinate: tuple[int, int],
    relax_with: Callable[[Atoms, list[Any]], tuple[Atoms, RelaxResult]],
) -> Callable[[float, ScanPoint | None], ScanPoint]:
    """Constrained relaxation at a fixed carbon-hydrogen distance, chained from the previous point."""
    i, j = coordinate

    def evaluate(d: float, prev: ScanPoint | None) -> ScanPoint:
        base = init if prev is None or prev.atoms is None else prev.atoms
        try:
            a, res = relax_with(set_coordinate(base, (i, j), d, STICK), [FixBondLength(i, j)])
        except (CalculationFailed, RuntimeError) as exc:
            log.warning("scan point %.2f A failed: %s", d, str(exc).splitlines()[0])
            return ScanPoint(d, float("nan"), False, None)
        return ScanPoint(d, res.energy_eV, res.converged, a)

    return evaluate


def _points(pts: Sequence[ScanPoint]) -> list[dict[str, Any]]:
    return [
        {
            "coordinate_A": p.coordinate_A,
            "energy_eV": p.energy_eV if np.isfinite(p.energy_eV) else None,
            "converged": p.converged,
            "excluded": p.excluded,
        }
        for p in pts
    ]


def _engine(calc: Calculator, cache_dir: Path) -> dict[str, Any]:
    """A calculator with its tier, engine label, settings hash and reference energies."""
    settings = calc.settings()
    return {
        "calc": calc,
        "tier": Tier(calc.tier),
        "settings": settings,
        "engine": f"{calc.name}:{settings.get('method', settings.get('xc', ''))}".rstrip(":"),
        "settings_hash": CalculatorSpec(calc.name, calc.tier, settings).settings_hash,
        "refs": reference_energies(calc, cache_dir),
    }


def build_barriers(
    spec: BarrierSpec, calc: Calculator, cache_dir: Path, out_dir: Path
) -> dict[str, Any]:
    """Compute the barriers of every pathway of the tier and write the barrier tables.

    Sticking pathways use ``spec.stick_recipes`` for their substrate and the calculator with
    ``spec.stick_calculator_settings`` applied (the same calculator when empty), because a
    barrier measured from the separated reactants needs integer occupations in a cell whose
    far-hydrogen state reproduces them (manual: Tier-2 protocol)."""
    t_start = time.perf_counter()
    cache_dir, out_dir = Path(cache_dir), Path(out_dir)
    sigma = spec.numerical_sigma_eV
    for m in spec.methods:
        if m not in METHODS:
            raise ValueError(f"unknown barrier method {m!r}")
    methods = tuple(m for m in spec.methods if not (m == CI_NEB and size() > 1))
    if methods != tuple(spec.methods):
        log.info("NEB skipped under MPI; the constrained scan is the only method")
    main = _engine(calc, cache_dir)
    engines = {RECOMB: main, STICK: main}
    if spec.stick_calculator_settings:
        if not is_dataclass(calc):
            raise ValueError("stick calculator settings require a dataclass calculator")
        calc_stick = cast(Calculator, replace(cast(Any, calc), **spec.stick_calculator_settings))
        engines[STICK] = _engine(calc_stick, cache_dir)

    records: list[SiteEnergyRecord] = []
    profiles: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    substrates: dict[tuple[str, str], tuple[Atoms, RelaxResult]] = {}

    radius_opt: dict[str, dict[str, Any]] = {}

    def relaxed_substrate(
        name: str, atoms: Atoms, c: Calculator, key: str
    ) -> tuple[Atoms, RelaxResult]:
        if (name, key) not in substrates:
            a0 = atoms.copy()
            set_unpaired_electrons(a0, int(np.sum(a0.numbers)) % 2)
            if a0.info.get("geometry") == "tube":
                a0, rad = relax_tube_radius(a0, c, cache_dir)
                radius_opt[f"{name}@{key}"] = rad
                log.info(
                    "substrate %s: radial scale %.5f, dE=%.4f eV (%d single points)",
                    name,
                    rad["scale"],
                    rad["dE_eV"],
                    rad["n_single_points"],
                )
            a, res = relax(a0, c, cache_dir, fmax=spec.fmax, steps=spec.steps)
            a.info.update(atoms.info)
            substrates[(name, key)] = (a, res)
        return substrates[(name, key)]

    def relax_free(a0: Atoms, c: Calculator) -> tuple[Atoms, RelaxResult]:
        return relax(a0, c, cache_dir, fmax=spec.fmax, steps=spec.steps)

    def relax_with(a0: Atoms, constraints: list[Any], c: Calculator) -> tuple[Atoms, RelaxResult]:
        return relax(a0, c, cache_dir, fmax=spec.fmax, steps=spec.steps, constraints=constraints)

    for pw in pathways(spec):
        cfg = pw.config
        eng = engines[pw.kind]
        c: Calculator = eng["calc"]
        refs: dict[str, float] = eng["refs"]
        settings_hash: str = eng["settings_hash"]
        tier: Tier = eng["tier"]
        engine: str = eng["engine"]
        try:
            sub, sub_res = relaxed_substrate(cfg.substrate_name, cfg.substrate, c, settings_hash)
            if not sub_res.converged:
                exclusions.append(
                    {"key": pw.key, "reason": "substrate relaxation did not converge"}
                )
                continue
            start = add_hydrogens(sub, cfg.placements)
            n_unp = int(np.sum(start.numbers)) % 2
            init, init_res = relax_free(start, c)
        except (CalculationFailed, RuntimeError) as exc:
            exclusions.append(
                {"key": pw.key, "reason": f"calculation failed: {str(exc).splitlines()[0]}"}
            )
            log.warning("pathway %s excluded: %s", pw.key, str(exc).splitlines()[0])
            continue
        if not init_res.converged:
            exclusions.append({"key": pw.key, "reason": "initial state did not converge"})
            continue
        chash = _hash(structure_hash(start), settings_hash, str(n_unp))
        i, j = pw.coordinate
        hydrogens = list(pw.hydrogens)
        normal = surface_normal(init, pw.site_carbons[0], pw.side)  # type: ignore[arg-type]
        centre = cfg.substrate.info.get("defect_center")
        desc = site_descriptors(
            init, int(cfg.placements[-1][0]), None if centre is None else np.asarray(centre)
        )
        e_ref = sub_res.energy_eV + refs["E_H"] if pw.kind == STICK else init_res.energy_eV
        e_gas = sub_res.energy_eV + (refs["E_H"] if pw.kind == STICK else refs["E_H2"])
        d0 = float(init.get_distance(i, j, mic=True))
        profile: dict[str, Any] = {
            "pathway": pw.key,
            "kind": pw.kind,
            "quantity": pw.quantity.value,
            "reference_state": pw.reference_state.value,
            "coordinate_indices": [i, j],
            "coordinate": "C-H distance",
            "hydrogen_indices": hydrogens,
            "initial_coordinate_A": d0,
            "initial_H_H_distance_A": (
                float(init.get_distance(hydrogens[0], hydrogens[1], mic=True))
                if len(hydrogens) == 2
                else None
            ),
            "E_initial_eV": init_res.energy_eV,
            "E_reference_eV": e_ref,
            "E_substrate_eV": sub_res.energy_eV,
            "E_H_eV": refs["E_H"],
            "E_H2_eV": refs["E_H2"],
            "E_desorbed_gas_from_initial_eV": e_gas - init_res.energy_eV,
            "spin_multiplicity": n_unp + 1,
            "n_atoms": len(init),
            "config_hash": chash,
            "initial_structure": _extxyz(init),
            "methods": {},
        }

        scan_seed: list[Atoms] = []
        scan_final: tuple[Atoms, RelaxResult] | None = None
        if CONSTRAINED_SCAN in methods:
            t0 = time.perf_counter()
            evaluate = _bond_evaluator(init, (i, j), functools.partial(relax_with, c=c))
            end = spec.far_distance_A if pw.kind == STICK else spec.recomb_distance_A
            grid = _grid(d0, end, spec.grid_step_A)
            first = ScanPoint(d0, init_res.energy_eV, True, init)
            points = scan_profile(evaluate, grid, spec.refine_step_A, spec.refine_halfwidth_A)
            all_points = sorted([first, *points], key=lambda p: p.coordinate_A)
            finite = [p for p in all_points if p.accepted]
            conv = [p for p in finite if p.converged] or finite
            top = max(conv, key=lambda p: p.energy_eV)
            last = next(p for p in reversed(all_points) if p.atoms is not None and p.accepted)
            final_distance: float | None = None
            product = False
            try:
                final_free, final_res = relax_free(
                    last.atoms if last.atoms is not None else init, c
                )
                final_energy: float | None = final_res.energy_eV
                final_converged = final_res.converged
                final_distance, product = product_distance(
                    pw.kind, final_free, (i, j), hydrogens, spec.far_distance_A, refs["H2_bond_A"]
                )
                if final_res.converged and product:
                    scan_final = (final_free, final_res)
            except (CalculationFailed, RuntimeError) as exc:
                log.warning("final state of %s failed: %s", pw.key, str(exc).splitlines()[0])
                final_free, final_energy, final_converged = init, None, False
            scan_seed = [p.atoms for p in all_points if p.atoms is not None]
            all_converged = (
                all(p.converged for p in points if p.excluded is None) and final_converged
            )
            barrier = top.energy_eV - e_ref
            at_end = top is last
            far_offset = (
                float(last.energy_eV - e_ref)
                if pw.kind == STICK and np.isfinite(last.energy_eV)
                else None
            )
            step = largest_step(finite)
            caveat = saddle_caveat(
                pw.kind,
                all_converged,
                product,
                at_end,
                final_distance,
                discontinuity=step,
                far_offset_eV=far_offset,
            )
            k_top = next(k for k, q in enumerate(all_points) if q is top)
            if any(q.excluded for q in all_points[max(0, k_top - 1) : k_top + 2]):
                caveat = "; ".join(
                    x for x in (caveat, "the maximum neighbours an excluded point") if x
                )
            profile["methods"][CONSTRAINED_SCAN] = {
                "barrier_eV": barrier,
                "barrier_from_initial_eV": top.energy_eV - init_res.energy_eV,
                "barrierless": bool(barrier <= 0.0),
                "transition_coordinate_A": top.coordinate_A,
                "maximum_at_grid_end": bool(at_end),
                "far_point_offset_eV": far_offset,
                "largest_step_eV": None if step is None else step[0],
                "largest_step_between_A": None if step is None else [step[1], step[2]],
                "E_final_eV": final_energy,
                "reaction_energy_eV": (
                    None if final_energy is None else final_energy - init_res.energy_eV
                ),
                "final_distance_A": final_distance,
                "product_reached": product,
                "saddle_located": caveat is None,
                "saddle_caveat": caveat,
                "points": _points(all_points),
                "n_points": len(points),
                "n_excluded_points": sum(q.excluded is not None for q in all_points),
                "converged": all_converged,
                "transition_structure": _extxyz(top.atoms if top.atoms is not None else init),
                "final_structure": _extxyz(final_free),
                "wall_s": time.perf_counter() - t0,
            }
            log.info(
                "barrier %s %s E=%.3f eV at %.2f A (%d points, saddle=%s%s) wall=%.1fs",
                pw.key,
                CONSTRAINED_SCAN,
                barrier,
                top.coordinate_A,
                len(points),
                caveat is None,
                "" if caveat is None else ": " + caveat,
                time.perf_counter() - t0,
            )

        neb_final: tuple[Atoms, RelaxResult] | None = None
        final_label = ""
        if CI_NEB in methods:

            def intact(
                a: Atoms,
                pair: tuple[int, int] = (i, j),
                hs: Sequence[int] = hydrogens,
                bond: float = refs["H2_bond_A"],
            ) -> bool:
                return product_distance(RECOMB, a, pair, hs, spec.far_distance_A, bond)[1]

            try:
                if pw.kind == STICK:
                    final0 = set_coordinate(init, (i, j), spec.far_distance_A, STICK)
                    neb_final = relax_with(final0, [FixBondLength(i, j)], c)
                    final_label = "H atom at the far distance"
                elif scan_final is not None:
                    neb_final = scan_final
                    final_label = "relaxed end point of the scan"
                else:
                    held = [FixSubsetCom(list(hydrogens))]
                    final0 = h2_above(init, hydrogens, pw.site_carbons, normal, spec.h2_height_A)
                    neb_final = relax_free(final0, c)
                    final_label = "H2 above the pair"
                    if not intact(neb_final[0]):
                        log.info(
                            "H2 above %s does not stay intact when relaxed freely; "
                            "its centre of mass is held",
                            pw.key,
                        )
                        neb_final = relax_with(final0, held, c)
                        final_label = "H2 above the pair, centre of mass held"
                    if not intact(neb_final[0]):
                        log.info(
                            "H2 dissociates %.1f A above %s with its centre of mass held; "
                            "the desorbed molecule %.1f A above the pair is the product",
                            spec.h2_height_A,
                            pw.key,
                            spec.h2_far_height_A,
                        )
                        far0 = h2_above(
                            init, hydrogens, pw.site_carbons, normal, spec.h2_far_height_A
                        )
                        neb_final = relax_with(far0, held, c)
                        final_label = "desorbed H2, centre of mass held"
                    scan_seed = []
            except (CalculationFailed, RuntimeError) as exc:
                log.warning("NEB final state of %s failed: %s", pw.key, str(exc).splitlines()[0])
                profile["methods"][CI_NEB] = {
                    "barrier_eV": None,
                    "converged": False,
                    "failure": f"final state: {str(exc).splitlines()[0]}",
                    "final_state": final_label,
                }
        if neb_final is not None:
            end_distance, end_is_product = product_distance(
                pw.kind, neb_final[0], (i, j), hydrogens, spec.far_distance_A, refs["H2_bond_A"]
            )
            if not end_is_product:
                label = "C-H" if pw.kind == STICK else "H-H"
                log.warning(
                    "NEB final state of %s is not the product (%s distance %.2f A); band skipped",
                    pw.key,
                    label,
                    end_distance,
                )
                profile["methods"][CI_NEB] = {
                    "barrier_eV": None,
                    "converged": False,
                    "failure": f"final state is not the product ({label} distance {end_distance:.2f} A)",
                    "final_state": final_label,
                }
                neb_final = None
        if neb_final is not None:
            final, final_res = neb_final
            traj_dir = cache_dir / "neb"
            if is_master():
                traj_dir.mkdir(parents=True, exist_ok=True)
            traj = traj_dir / f"{chash}_{structure_hash(final)[:8]}_{spec.neb_images}.traj"
            res = neb_barrier(
                init,
                final,
                c,
                spec.neb_images,
                spec.neb_spring_eV_A,
                spec.fmax,
                spec.neb_steps,
                traj,
                prerelax_fmax=spec.neb_prerelax_fmax,
                seed=scan_seed or None,
            )
            if not np.isfinite(res.energies_eV).all():
                failed = [k for k, e in enumerate(res.energies_eV) if not np.isfinite(e)]
                log.warning(
                    "barrier %s %s: images %s have no converged energy after %d steps; "
                    "the path maximum is undefined",
                    pw.key,
                    CI_NEB,
                    failed,
                    res.n_steps,
                )
                profile["methods"][CI_NEB] = {
                    "barrier_eV": None,
                    "converged": False,
                    "failure": f"images {failed} have no converged energy",
                    "final_state": final_label,
                    "image_energies_eV": [e if np.isfinite(e) else None for e in res.energies_eV],
                    "n_images": spec.neb_images,
                    "seeded_from_scan": bool(scan_seed),
                    "n_steps": res.n_steps,
                    "wall_s": res.wall_s,
                }
                neb_final = None
        if neb_final is not None:
            finite_e = np.where(np.isfinite(res.energies_eV), res.energies_eV, -np.inf)
            top_i = int(np.argmax(finite_e))
            barrier = float(np.nanmax(res.energies_eV)) - e_ref
            neb_distance, neb_product = product_distance(
                pw.kind, res.images[-1], (i, j), hydrogens, spec.far_distance_A, refs["H2_bond_A"]
            )
            neb_converged = res.converged and final_res.converged
            neb_far_offset = (
                float(res.energies_eV[-1] - e_ref)
                if pw.kind == STICK and np.isfinite(res.energies_eV[-1])
                else None
            )
            neb_caveat = saddle_caveat(
                pw.kind,
                neb_converged,
                neb_product,
                top_i in (0, spec.neb_images - 1),
                neb_distance,
                CI_NEB,
                far_offset_eV=neb_far_offset,
            )
            profile["methods"][CI_NEB] = {
                "barrier_eV": barrier,
                "barrier_from_initial_eV": res.forward_barrier_eV,
                "reverse_barrier_eV": res.reverse_barrier_eV,
                "barrierless": bool(barrier <= 0.0),
                "image_energies_eV": [e if np.isfinite(e) else None for e in res.energies_eV],
                "n_images": spec.neb_images,
                "spring_eV_A": spec.neb_spring_eV_A,
                "seeded_from_scan": bool(scan_seed),
                "final_state": final_label,
                "converged": neb_converged,
                "far_point_offset_eV": neb_far_offset,
                "final_distance_A": neb_distance,
                "product_reached": neb_product,
                "saddle_located": neb_caveat is None,
                "saddle_caveat": neb_caveat,
                "n_steps": res.n_steps,
                "cached": res.cached,
                "E_final_eV": res.energies_eV[-1],
                "reaction_energy_eV": res.energies_eV[-1] - init_res.energy_eV,
                "transition_structure": _extxyz(res.images[top_i]),
                "final_structure": _extxyz(res.images[-1]),
                "wall_s": res.wall_s,
            }
            log.info(
                "barrier %s %s E=%.3f eV (image %d of %d, saddle=%s%s) wall=%.1fs",
                pw.key,
                CI_NEB,
                barrier,
                top_i,
                spec.neb_images,
                neb_caveat is None,
                "" if neb_caveat is None else ": " + neb_caveat,
                res.wall_s,
            )

        primary = primary_method(pw.kind, profile["methods"])
        profile["primary_method"] = primary
        if primary is not None:
            r = profile["methods"][primary]
            records.append(
                _record(
                    pw,
                    chash,
                    primary,
                    r["barrier_eV"],
                    settings_hash,
                    tier,
                    engine,
                    bool(r.get("saddle_located", r["converged"])),
                    desc,
                    len(init),
                    sigma,
                )
            )
        profiles.append(profile)

    frame = records_to_frame(records)
    manifest: dict[str, Any] = {
        "tier": spec.tier,
        "calculator": {
            "name": calc.name,
            "tier": calc.tier,
            "settings": main["settings"],
            "settings_hash": main["settings_hash"],
        },
        "stick_calculator_settings": spec.stick_calculator_settings,
        "stick_recipes": spec.stick_recipes,
        "stick_settings_hash": engines[STICK]["settings_hash"],
        "methods": list(methods),
        "primary_method": {STICK: CONSTRAINED_SCAN, RECOMB: CI_NEB},
        "grid_step_A": spec.grid_step_A,
        "refine_step_A": spec.refine_step_A,
        "fmax_eV_A": spec.fmax,
        "tube_radius_optimization": radius_opt,
        "neb_images": spec.neb_images,
        "neb_spring_eV_A": spec.neb_spring_eV_A,
        "reference_energies_eV": refs,
        "n_records": len(frame),
        "n_pathways": len(profiles),
        "n_excluded": len(exclusions),
        "package_version": __version__,
        "wall_s": time.perf_counter() - t_start,
    }
    if is_master():
        out_dir.mkdir(parents=True, exist_ok=True)
        table = out_dir / f"barriers_{spec.tier}.parquet"
        write_table(frame, table, {"tier": spec.tier, "settings_hash": main["settings_hash"]})
        prof = out_dir / f"barrier_profiles_{spec.tier}.json"
        prof.write_text(json.dumps(profiles, indent=1))
        excl = out_dir / f"barriers_{spec.tier}_excluded.csv"
        pd.DataFrame(exclusions, columns=["key", "reason"]).to_csv(excl, index=False)
        manifest["files"] = {p.name: sha256_file(p) for p in (table, prof, excl)}
        (out_dir / f"barriers_manifest_{spec.tier}.json").write_text(json.dumps(manifest, indent=2))
    return {"records": frame, "profiles": profiles, "exclusions": exclusions, "manifest": manifest}


__all__ = [
    "BarrierSpec",
    "NEBResult",
    "Pathway",
    "ScanPoint",
    "build_barriers",
    "h2_above",
    "largest_step",
    "mark_islands",
    "neb_barrier",
    "pathways",
    "primary_method",
    "product_distance",
    "resample_path",
    "saddle_caveat",
    "scan_profile",
    "set_coordinate",
]
