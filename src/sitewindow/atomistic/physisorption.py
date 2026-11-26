"""Pairwise H2-carbon potentials with the Feynman-Hibbs correction and physisorption site
energies on carbon structures (manual C2, physisorption set; DERIVATIONS 2.10; benchmark B3).

H2 is a single Lennard-Jones site. Two parameter sets, both Lorentz-Berthelot combinations of
Steele's graphite-carbon parameters with a published H2 parameter set, serve as the model-form
alternatives of the manual: the Buch set and the Darkrim-Levesque set (values and sources in
``data/registry/physisorption_potentials.yaml``). Oxygen atoms of surface groups interact with
H2 through Lorentz-Berthelot combinations of the same H2 set with the OPLS-AA oxygen parameters
of the group type (ether oxygen for epoxides, alcohol oxygen for hydroxyls and the acid hydroxyl,
carbonyl oxygen for carbonyls), assigned from the oxygen's bonding environment; group hydrogens
carry no Lennard-Jones centre, as in OPLS-AA. Site energies are the minimum of the summed pair
potential over probe positions, with periodic images of in-plane periodic structures. Energies
are reported in kJ mol^-1 per H2 (negative for binding).
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
    """Lennard-Jones H2-C pair potential: well depth epsilon/k_B (K) and diameter sigma (A).

    ``h2_epsilon_K`` and ``h2_sigma_A`` are the H2 parameters the carbon pair was combined from;
    they combine with the oxygen parameters of surface groups (``oxygen_pair``). The B3
    calibration scales the carbon pair only, since B3 is a carbon-only benchmark.
    """

    name: str
    epsilon_K: float
    sigma_A: float
    source: str
    h2_epsilon_K: float | None = None
    h2_sigma_A: float | None = None

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
OPLS_SOURCE = "Jorgensen, Maxwell and Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"
KCAL_PER_MOL_TO_K = 503.2197
OXYGEN_TYPES: dict[str, dict[str, Any]] = {
    "ether": {"epsilon_K": 0.140 * KCAL_PER_MOL_TO_K, "sigma_A": 2.90, "source": OPLS_SOURCE},
    "alcohol": {"epsilon_K": 0.170 * KCAL_PER_MOL_TO_K, "sigma_A": 3.12, "source": OPLS_SOURCE},
    "carbonyl": {"epsilon_K": 0.210 * KCAL_PER_MOL_TO_K, "sigma_A": 2.96, "source": OPLS_SOURCE},
}
O_C_BOND_A = 1.55
O_C_DOUBLE_BOND_A = 1.30
O_H_BOND_A = 1.2


def lorentz_berthelot(name: str, carbon: dict[str, Any], hydrogen: dict[str, Any]) -> PairPotential:
    return PairPotential(
        name,
        float(np.sqrt(carbon["epsilon_K"] * hydrogen["epsilon_K"])),
        0.5 * (carbon["sigma_A"] + hydrogen["sigma_A"]),
        f"Lorentz-Berthelot: {carbon['source']}; {hydrogen['source']}",
        float(hydrogen["epsilon_K"]),
        float(hydrogen["sigma_A"]),
    )


def oxygen_type(atoms: Atoms, index: int) -> str:
    """Group type of oxygen ``index`` from its bonding: a hydrogen within ``O_H_BOND_A`` ->
    alcohol (hydroxyl, acid hydroxyl); a carbon within the double-bond length ``O_C_DOUBLE_BOND_A``
    -> carbonyl; two carbons within ``O_C_BOND_A`` -> ether (epoxide); otherwise carbonyl."""
    d = atoms.get_distances(index, range(len(atoms)), mic=bool(atoms.pbc.any()))
    carbon = (atoms.numbers == 6) & (d > 0.0)
    if int(np.sum((atoms.numbers == 1) & (d < O_H_BOND_A))) >= 1:
        return "alcohol"
    if int(np.sum(carbon & (d < O_C_DOUBLE_BOND_A))) >= 1:
        return "carbonyl"
    if int(np.sum(carbon & (d < O_C_BOND_A))) >= 2:
        return "ether"
    return "carbonyl"


def oxygen_pair(p: PairPotential, kind: str) -> tuple[float, float]:
    """(epsilon_K, sigma_A) of the H2-O pair for oxygen ``kind`` with the H2 set of ``p``."""
    if p.h2_epsilon_K is None or p.h2_sigma_A is None:
        raise ValueError(f"potential {p.name} carries no H2 parameters to combine with oxygen")
    o = OXYGEN_TYPES[kind]
    return float(np.sqrt(o["epsilon_K"] * p.h2_epsilon_K)), 0.5 * (o["sigma_A"] + p.h2_sigma_A)


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
        p.h2_epsilon_K,
        p.h2_sigma_A,
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


def lj_energy_K(
    r_A: NDArray[np.float64],
    p: PairPotential,
    epsilon_K: NDArray[np.float64] | float | None = None,
    sigma_A: NDArray[np.float64] | float | None = None,
) -> NDArray[np.float64]:
    eps = p.epsilon_K if epsilon_K is None else epsilon_K
    sig = p.sigma_A if sigma_A is None else sigma_A
    x = (sig / r_A) ** 6
    return 4.0 * eps * (x * x - x)


def lj_laplacian_K_per_A2(
    r_A: NDArray[np.float64],
    p: PairPotential,
    epsilon_K: NDArray[np.float64] | float | None = None,
    sigma_A: NDArray[np.float64] | float | None = None,
) -> NDArray[np.float64]:
    """Laplacian of the pair potential, 4 eps [132 sigma^12/r^14 - 30 sigma^6/r^8] (DERIVATIONS 2.10)."""
    eps = p.epsilon_K if epsilon_K is None else epsilon_K
    sig = p.sigma_A if sigma_A is None else sigma_A
    s6 = sig**6
    return 4.0 * eps * (132.0 * s6 * s6 / r_A**14 - 30.0 * s6 / r_A**8)


def feynman_hibbs_coefficient_A2(T_K: float, mass_amu: float = M_H2_AMU) -> float:
    """hbar^2 / (24 mu k_B T) in A^2 for a molecule of mass ``mass_amu`` against a rigid surface."""
    mu = mass_amu * AMU_KG
    return HBAR**2 / (24.0 * mu * K_B_J * T_K) * 1.0e20


def pair_energy_K(
    r_A: NDArray[np.float64],
    p: PairPotential,
    T_K: float | None = None,
    epsilon_K: NDArray[np.float64] | float | None = None,
    sigma_A: NDArray[np.float64] | float | None = None,
) -> NDArray[np.float64]:
    """Classical (``T_K=None``) or Feynman-Hibbs-corrected pair energy (K)."""
    u = lj_energy_K(r_A, p, epsilon_K, sigma_A)
    if T_K is None:
        return u
    return u + feynman_hibbs_coefficient_A2(T_K) * lj_laplacian_K_per_A2(r_A, p, epsilon_K, sigma_A)


def _image_shifts(atoms: Atoms, cutoff_A: float) -> NDArray[np.float64]:
    cell = np.array(atoms.cell)
    reps = []
    for axis in range(3):
        if atoms.pbc[axis]:
            n = int(np.ceil(cutoff_A / np.linalg.norm(cell[axis]))) + 1
            reps.append(range(-n, n + 1))
        else:
            reps.append(range(0, 1))
    return np.array(
        [i * cell[0] + j * cell[1] + k * cell[2] for i in reps[0] for j in reps[1] for k in reps[2]]
    )


def carbon_images(atoms: Atoms, cutoff_A: float = CUTOFF_A) -> NDArray[np.float64]:
    """Carbon positions including the periodic images needed within ``cutoff_A`` of the cell."""
    pos = atoms.positions[atoms.numbers == 6]
    shifts = _image_shifts(atoms, cutoff_A)
    return (pos[None, :, :] + shifts[:, None, :]).reshape(-1, 3)


@dataclass(frozen=True)
class SiteImages:
    """Interaction centres of a structure with their periodic images and pair parameters: the
    carbons with the H2-C pair of the potential and the oxygens of surface groups with the H2-O
    pair of their group type. Group hydrogens carry no centre."""

    positions: NDArray[np.float64]
    epsilon_K: NDArray[np.float64]
    sigma_A: NDArray[np.float64]
    n_oxygen: int

    def __len__(self) -> int:
        return len(self.positions)


def site_images(atoms: Atoms, p: PairPotential, cutoff_A: float = CUTOFF_A) -> SiteImages:
    """Carbon and oxygen centres with periodic images and their pair parameters under ``p``."""
    shifts = _image_shifts(atoms, cutoff_A)
    carbons = atoms.positions[atoms.numbers == 6]
    pos = [carbons]
    eps = [np.full(len(carbons), p.epsilon_K)]
    sig = [np.full(len(carbons), p.sigma_A)]
    oxygens = np.flatnonzero(atoms.numbers == 8)
    for k in oxygens:
        e_o, s_o = oxygen_pair(p, oxygen_type(atoms, int(k)))
        pos.append(atoms.positions[k][None, :])
        eps.append(np.array([e_o]))
        sig.append(np.array([s_o]))
    centres = np.vstack(pos)
    eps_a = np.concatenate(eps)
    sig_a = np.concatenate(sig)
    n_shift = len(shifts)
    return SiteImages(
        (centres[None, :, :] + shifts[:, None, :]).reshape(-1, 3),
        np.tile(eps_a, n_shift),
        np.tile(sig_a, n_shift),
        int(len(oxygens)),
    )


def probe_energy_kJ_mol(
    points: NDArray[np.float64],
    carbons: NDArray[np.float64] | SiteImages,
    p: PairPotential,
    T_K: float | None = None,
    cutoff_A: float = CUTOFF_A,
) -> NDArray[np.float64]:
    """Summed pair energy (kJ mol^-1) at each probe point over carbon positions (the H2-C pair of
    ``p``) or over :class:`SiteImages` (per-centre parameters)."""
    pts = np.atleast_2d(points)
    out = np.empty(len(pts))
    if isinstance(carbons, SiteImages):
        for i, q in enumerate(pts):
            d = np.linalg.norm(carbons.positions - q, axis=1)
            mask = (d < cutoff_A) & (d > 0.5)
            out[i] = (
                float(
                    np.sum(
                        pair_energy_K(
                            d[mask], p, T_K, carbons.epsilon_K[mask], carbons.sigma_A[mask]
                        )
                    )
                )
                * K_TO_KJ_PER_MOL
            )
        return out
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
    """Minimum of the summed potential starting from the best of ``start_points``; oxygen atoms of
    surface groups contribute through their group-type pair parameters."""
    carbons: NDArray[np.float64] | SiteImages = (
        site_images(atoms, p) if np.any(atoms.numbers == 8) else carbon_images(atoms)
    )
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
    near: NDArray[np.float64] | None = None,
    radius_A: float = 4.0,
) -> tuple[float, NDArray[np.float64]]:
    """Strongest H2 site above a sheet (graphene, defective graphene, coronene, ribbon).

    With ``near`` the probe grid is restricted to the lateral disc of ``radius_A`` around that
    point (the site adjacent to a surface group), and the heights are measured from the sheet's
    carbon plane and extended above the group's atoms.
    """
    carbons = atoms.positions[atoms.numbers == 6]
    z0 = float(np.median(carbons[:, 2])) if near is not None else float(np.max(carbons[:, 2]))
    if any(atoms.pbc[:2]):
        lateral = _lateral_grid(atoms, spacing_A)
    else:
        xy = carbons[:, :2]
        xs = np.arange(xy[:, 0].min(), xy[:, 0].max() + 1e-9, spacing_A)
        ys = np.arange(xy[:, 1].min(), xy[:, 1].max() + 1e-9, spacing_A)
        lateral = np.array([(x, y) for x in xs for y in ys])
    heights = list(heights_A)
    if near is not None:
        c = np.asarray(near, dtype=float)[:2]
        cell = np.array(atoms.cell)
        d = lateral - c[None, :]
        for axis in range(2):
            if atoms.pbc[axis]:
                frac = np.rint(d @ cell[axis, :2] / float(cell[axis, :2] @ cell[axis, :2]))
                d = d - frac[:, None] * cell[axis, :2][None, :]
        lateral = lateral[np.linalg.norm(d, axis=1) <= radius_A]
        top = float(np.max(atoms.positions[:, 2]))
        heights += [top - z0 + h for h in (2.8, 3.2, 3.6)]
    pts = np.array([(x, y, z0 + h) for x, y in lateral for h in heights])
    lo, hi = z0 + 2.0, float(np.max(atoms.positions[:, 2])) + 6.0
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
