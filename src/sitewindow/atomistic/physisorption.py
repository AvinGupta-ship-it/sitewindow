"""Pairwise H2-carbon potentials with the Feynman-Hibbs correction and physisorption site
energies on carbon structures (manual C2, physisorption set; DERIVATIONS 2.10; benchmark B3).

H2 is a single Lennard-Jones site. Two parameter sets, both Lorentz-Berthelot combinations of
Steele's graphite-carbon parameters with a published H2 parameter set, serve as the model-form
alternatives of the manual: the Buch set and the Darkrim-Levesque set (values and sources in
``data/registry/physisorption_potentials.yaml``). Site energies are the minimum of the summed
pair potential over probe positions, with periodic images of in-plane periodic structures.
Energies are reported in kJ mol^-1 per H2 (negative for binding).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from ase import Atoms
from numpy.typing import NDArray
from scipy.optimize import minimize

from sitewindow.constants import AMU_KG, HBAR, K_B_J, M_H2_AMU, R_GAS

K_TO_KJ_PER_MOL = R_GAS / 1000.0
CUTOFF_A = 15.0


@dataclass(frozen=True)
class PairPotential:
    """Lennard-Jones H2-C pair potential: well depth epsilon/k_B (K) and diameter sigma (A)."""

    name: str
    epsilon_K: float
    sigma_A: float
    source: str

    @property
    def settings_hash(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:16]


STEELE_CARBON = {"epsilon_K": 28.0, "sigma_A": 3.40, "source": "Steele, Surf. Sci. 36, 317 (1973)"}
H2_BUCH = {"epsilon_K": 34.2, "sigma_A": 2.96, "source": "Buch, J. Chem. Phys. 100, 7610 (1994)"}
H2_DARKRIM_LEVESQUE = {
    "epsilon_K": 36.7,
    "sigma_A": 2.958,
    "source": "Darkrim and Levesque, J. Chem. Phys. 109, 4981 (1998)",
}


def lorentz_berthelot(name: str, carbon: dict[str, Any], hydrogen: dict[str, Any]) -> PairPotential:
    return PairPotential(
        name,
        float(np.sqrt(carbon["epsilon_K"] * hydrogen["epsilon_K"])),
        0.5 * (carbon["sigma_A"] + hydrogen["sigma_A"]),
        f"Lorentz-Berthelot: {carbon['source']}; {hydrogen['source']}",
    )


STEELE_BUCH = lorentz_berthelot("steele_buch", STEELE_CARBON, H2_BUCH)
STEELE_DARKRIM_LEVESQUE = lorentz_berthelot(
    "steele_darkrim_levesque", STEELE_CARBON, H2_DARKRIM_LEVESQUE
)
PARAMETER_SETS: dict[str, PairPotential] = {
    p.name: p for p in (STEELE_BUCH, STEELE_DARKRIM_LEVESQUE)
}
B3_CORONENE_KJ_MOL = -4.7


def calibrate_to_coronene(
    p: PairPotential, target_kJ_mol: float = B3_CORONENE_KJ_MOL
) -> PairPotential:
    """Scale the well depth so that the classical coronene minimum equals the B3 reference.

    The summed potential is linear in epsilon at fixed sigma, so the minimum position is
    unchanged and the scale factor is exact: epsilon' = epsilon * target / E_min(p).
    """
    e_min, _ = coronene_binding_kJ_mol(p)
    scale = target_kJ_mol / e_min
    return PairPotential(
        f"{p.name}_b3",
        p.epsilon_K * scale,
        p.sigma_A,
        f"{p.source}; epsilon scaled by {scale:.4f} to the B3 coronene value {target_kJ_mol} kJ/mol",
    )


_CALIBRATED: dict[str, PairPotential] = {}


def calibrated_potential(base: PairPotential = STEELE_BUCH) -> PairPotential:
    """The B3-calibrated potential (computed once per process)."""
    if base.name not in _CALIBRATED:
        _CALIBRATED[base.name] = calibrate_to_coronene(base)
    return _CALIBRATED[base.name]


def all_potentials() -> dict[str, PairPotential]:
    """Literature sets plus the B3-calibrated set."""
    cal = calibrated_potential()
    return {**PARAMETER_SETS, cal.name: cal}


def lj_energy_K(r_A: NDArray[np.float64], p: PairPotential) -> NDArray[np.float64]:
    x = (p.sigma_A / r_A) ** 6
    return 4.0 * p.epsilon_K * (x * x - x)


def lj_laplacian_K_per_A2(r_A: NDArray[np.float64], p: PairPotential) -> NDArray[np.float64]:
    """Laplacian of the pair potential, 4 eps [132 sigma^12/r^14 - 30 sigma^6/r^8] (DERIVATIONS 2.10)."""
    s6 = p.sigma_A**6
    return 4.0 * p.epsilon_K * (132.0 * s6 * s6 / r_A**14 - 30.0 * s6 / r_A**8)


def feynman_hibbs_coefficient_A2(T_K: float, mass_amu: float = M_H2_AMU) -> float:
    """hbar^2 / (24 mu k_B T) in A^2 for a molecule of mass ``mass_amu`` against a rigid surface."""
    mu = mass_amu * AMU_KG
    return HBAR**2 / (24.0 * mu * K_B_J * T_K) * 1.0e20


def pair_energy_K(
    r_A: NDArray[np.float64], p: PairPotential, T_K: float | None = None
) -> NDArray[np.float64]:
    """Classical (``T_K=None``) or Feynman-Hibbs-corrected pair energy (K)."""
    u = lj_energy_K(r_A, p)
    if T_K is None:
        return u
    return u + feynman_hibbs_coefficient_A2(T_K) * lj_laplacian_K_per_A2(r_A, p)


def carbon_images(atoms: Atoms, cutoff_A: float = CUTOFF_A) -> NDArray[np.float64]:
    """Carbon positions including the periodic images needed within ``cutoff_A`` of the cell."""
    pos = atoms.positions[atoms.numbers == 6]
    cell = np.array(atoms.cell)
    reps = []
    for axis in range(3):
        if atoms.pbc[axis]:
            n = int(np.ceil(cutoff_A / np.linalg.norm(cell[axis]))) + 1
            reps.append(range(-n, n + 1))
        else:
            reps.append(range(0, 1))
    shifts = np.array(
        [i * cell[0] + j * cell[1] + k * cell[2] for i in reps[0] for j in reps[1] for k in reps[2]]
    )
    return (pos[None, :, :] + shifts[:, None, :]).reshape(-1, 3)


def probe_energy_kJ_mol(
    points: NDArray[np.float64],
    carbons: NDArray[np.float64],
    p: PairPotential,
    T_K: float | None = None,
    cutoff_A: float = CUTOFF_A,
) -> NDArray[np.float64]:
    """Summed H2-C energy (kJ mol^-1) at each probe point."""
    pts = np.atleast_2d(points)
    out = np.empty(len(pts))
    for i, q in enumerate(pts):
        d = np.linalg.norm(carbons - q, axis=1)
        d = d[(d < cutoff_A) & (d > 0.5)]
        out[i] = float(np.sum(pair_energy_K(d, p, T_K))) * K_TO_KJ_PER_MOL
    return out


def minimize_site_energy(
    atoms: Atoms,
    start_points: NDArray[np.float64],
    p: PairPotential,
    T_K: float | None = None,
    bounds: list[tuple[float, float]] | None = None,
) -> tuple[float, NDArray[np.float64]]:
    """Minimum of the summed potential starting from the best of ``start_points``."""
    carbons = carbon_images(atoms)
    e = probe_energy_kJ_mol(start_points, carbons, p, T_K)
    best = int(np.argmin(e))
    res = minimize(
        lambda q: probe_energy_kJ_mol(q[None, :], carbons, p, T_K)[0],
        start_points[best],
        method="L-BFGS-B",
        bounds=bounds,
    )
    return float(res.fun), np.asarray(res.x, dtype=float)


def _lateral_grid(atoms: Atoms, spacing_A: float) -> NDArray[np.float64]:
    cell = np.array(atoms.cell)
    n0 = max(2, int(np.ceil(np.linalg.norm(cell[0]) / spacing_A)))
    n1 = max(2, int(np.ceil(np.linalg.norm(cell[1]) / spacing_A)))
    f0, f1 = np.meshgrid(np.arange(n0) / n0, np.arange(n1) / n1, indexing="ij")
    return f0.reshape(-1, 1) * cell[0][None, :2] + f1.reshape(-1, 1) * cell[1][None, :2]


def sheet_site_energy(
    atoms: Atoms,
    p: PairPotential,
    T_K: float | None = None,
    heights_A: tuple[float, ...] = (2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0),
    spacing_A: float = 0.3,
) -> tuple[float, NDArray[np.float64]]:
    """Strongest H2 site above a sheet (graphene, defective graphene, coronene, ribbon)."""
    carbons = atoms.positions[atoms.numbers == 6]
    z0 = float(np.max(carbons[:, 2]))
    if any(atoms.pbc[:2]):
        lateral = _lateral_grid(atoms, spacing_A)
    else:
        xy = carbons[:, :2]
        xs = np.arange(xy[:, 0].min(), xy[:, 0].max() + 1e-9, spacing_A)
        ys = np.arange(xy[:, 1].min(), xy[:, 1].max() + 1e-9, spacing_A)
        lateral = np.array([(x, y) for x in xs for y in ys])
    pts = np.array([(x, y, z0 + h) for x, y in lateral for h in heights_A])
    lo, hi = z0 + 2.0, z0 + 6.0
    return minimize_site_energy(atoms, pts, p, T_K, bounds=[(None, None), (None, None), (lo, hi)])


def pore_site_energy(
    pore: Atoms, p: PairPotential, T_K: float | None = None, spacing_A: float = 0.3
) -> tuple[float, NDArray[np.float64]]:
    """Strongest H2 site inside a slit pore (two sheets)."""
    z = pore.positions[pore.numbers == 6, 2]
    z_lo, z_hi = float(z.min()), float(z.max())
    lateral = _lateral_grid(pore, spacing_A)
    heights = (
        np.linspace(z_lo + 2.5, z_hi - 2.5, 9)
        if z_hi - z_lo > 5.0
        else np.array([0.5 * (z_lo + z_hi)])
    )
    pts = np.array([(x, y, h) for x, y in lateral for h in heights])
    return minimize_site_energy(
        pore, pts, p, T_K, bounds=[(None, None), (None, None), (z_lo + 1.5, z_hi - 1.5)]
    )


def tube_site_energy(
    tube: Atoms, p: PairPotential, interior: bool, T_K: float | None = None, n_angles: int = 24
) -> tuple[float, NDArray[np.float64]]:
    """Strongest H2 site inside (``interior``) or outside a nanotube."""
    carbons = tube.positions[tube.numbers == 6]
    axis_xy = carbons[:, :2].mean(axis=0)
    radius = float(np.mean(np.linalg.norm(carbons[:, :2] - axis_xy, axis=1)))
    zs = np.linspace(0.0, float(tube.cell[2, 2]), 5, endpoint=False)
    if interior:
        radii = np.linspace(0.0, max(radius - 2.6, 0.0), 6)
    else:
        radii = radius + np.array([2.8, 3.1, 3.4, 3.8])
    angles = np.linspace(0.0, 2 * np.pi, n_angles, endpoint=False)
    pts = np.array(
        [
            (axis_xy[0] + r * np.cos(a), axis_xy[1] + r * np.sin(a), z)
            for r in radii
            for a in angles
            for z in zs
        ]
    )
    if interior:
        bounds = None
    else:
        bounds = None
    e, q = minimize_site_energy(tube, pts, p, T_K, bounds=bounds)
    r_final = float(np.linalg.norm(q[:2] - axis_xy))
    if interior and r_final > radius:
        e, q = (
            float(
                probe_energy_kJ_mol(
                    pts[np.argmin(probe_energy_kJ_mol(pts, carbon_images(tube), p, T_K))][None, :],
                    carbon_images(tube),
                    p,
                    T_K,
                )[0]
            ),
            pts[0],
        )
    return e, q


def bundle(tube: Atoms, gap_A: float = 3.2) -> Atoms:
    """Three tubes on a triangular lattice (van der Waals gap ``gap_A``) for groove and
    interstitial sites, periodic along the tube axis."""
    carbons = tube.positions[tube.numbers == 6]
    axis_xy = carbons[:, :2].mean(axis=0)
    radius = float(np.mean(np.linalg.norm(carbons[:, :2] - axis_xy, axis=1)))
    d = 2.0 * radius + gap_A
    out = tube.copy()
    out.info = dict(tube.info)
    for shift in ((d, 0.0), (0.5 * d, 0.5 * np.sqrt(3.0) * d)):
        t = tube.copy()
        t.positions[:, 0] += shift[0]
        t.positions[:, 1] += shift[1]
        out += t
    cell = np.array(out.cell)
    cell[0, 0] = max(cell[0, 0], 3 * d + 2 * radius)
    cell[1, 1] = max(cell[1, 1], 3 * d + 2 * radius)
    out.set_cell(cell)
    out.info.update(
        geometry="bundle",
        tube_centres=[axis_xy, axis_xy + (d, 0.0), axis_xy + (0.5 * d, 0.5 * np.sqrt(3.0) * d)],
        radius_A=radius,
    )
    return out


def bundle_site_energy(
    b: Atoms, p: PairPotential, site: str, T_K: float | None = None
) -> tuple[float, NDArray[np.float64]]:
    """Groove (between two tubes, outside the bundle) or interstitial (between three tubes) site."""
    c = np.array(b.info["tube_centres"])
    zs = np.linspace(0.0, float(b.cell[2, 2]), 5, endpoint=False)
    if site == "interstitial":
        centre = c.mean(axis=0)
        pts = np.array(
            [
                (centre[0] + dx, centre[1] + dy, z)
                for dx in (-0.3, 0.0, 0.3)
                for dy in (-0.3, 0.0, 0.3)
                for z in zs
            ]
        )
    elif site == "groove":
        mid = 0.5 * (c[0] + c[1])
        outward = np.array([0.0, -1.0])
        pts = np.array(
            [(mid[0], mid[1] + outward[1] * h, z) for h in np.linspace(1.5, 4.5, 7) for z in zs]
        )
    else:
        raise ValueError("site must be 'groove' or 'interstitial'")
    return minimize_site_energy(b, pts, p, T_K)


def coronene_binding_kJ_mol(p: PairPotential = STEELE_BUCH) -> tuple[float, NDArray[np.float64]]:
    """Classical minimum of the pair potential over coronene (benchmark B3, no quantum correction)."""
    from sitewindow.atomistic.structures.builders import coronene

    return sheet_site_energy(coronene(), p, None, spacing_A=0.25)


__all__ = [
    "PairPotential",
    "STEELE_BUCH",
    "STEELE_DARKRIM_LEVESQUE",
    "PARAMETER_SETS",
    "B3_CORONENE_KJ_MOL",
    "calibrate_to_coronene",
    "calibrated_potential",
    "all_potentials",
    "lorentz_berthelot",
    "lj_energy_K",
    "lj_laplacian_K_per_A2",
    "feynman_hibbs_coefficient_A2",
    "pair_energy_K",
    "carbon_images",
    "probe_energy_kJ_mol",
    "minimize_site_energy",
    "sheet_site_energy",
    "pore_site_energy",
    "tube_site_energy",
    "bundle",
    "bundle_site_energy",
    "coronene_binding_kJ_mol",
    "K_TO_KJ_PER_MOL",
]
