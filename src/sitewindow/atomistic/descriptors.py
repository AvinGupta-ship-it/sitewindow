"""Local geometric descriptors of hydrogen binding sites (manual C2; DATA_SCHEMA SiteSpec).

Descriptors are evaluated on the relaxed structure that carries the site:

- pyramidalization angle theta_p = theta_sigma-pi - 90 deg of the site carbon from the
  pi-orbital axis vector (POAV1) construction: the vector making equal angles with the three
  C-C bonds; 0 for planar sp2 carbon, 11.6 deg in C60, 19.47 deg for sp3;
- coordination: number of carbon neighbours within the C-C cutoff;
- n_neighbor_H: hydrogen atoms bonded to other carbons within a radius of the site carbon;
- distance to the defect centre (nm, minimum image);
- curvature 1/d (nm^-1) for nanotubes, 0 for sheets.
"""

from __future__ import annotations

import math

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.neighborlist import neighbor_list
from numpy.typing import NDArray

C_C_CUTOFF_A = 1.75
C_H_CUTOFF_A = 1.30
NEIGHBOR_H_RADIUS_A = 3.0


def carbon_neighbors(
    atoms: Atoms, index: int, cutoff: float = C_C_CUTOFF_A
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Indices of carbon neighbours of ``index`` and the minimum-image vectors to them."""
    i, j, D = neighbor_list("ijD", atoms, cutoff)
    sel = (i == index) & (atoms.numbers[j] == 6)
    return j[sel], D[sel]


def coordination(atoms: Atoms, index: int, cutoff: float = C_C_CUTOFF_A) -> int:
    return int(len(carbon_neighbors(atoms, index, cutoff)[0]))


def pyramidalization_angle_deg(atoms: Atoms, index: int, cutoff: float = C_C_CUTOFF_A) -> float:
    """POAV1 pyramidalization angle (degrees); NaN unless the carbon has three carbon neighbours."""
    _, D = carbon_neighbors(atoms, index, cutoff)
    if len(D) != 3:
        return math.nan
    u = D / np.linalg.norm(D, axis=1)[:, None]
    n = np.cross(u[1] - u[0], u[2] - u[0])
    nn = float(np.linalg.norm(n))
    if nn < 1e-10:
        return math.nan
    n /= nn
    c = float(np.mean(u @ n))
    if c > 0:
        c = -c
    theta_sigma_pi = math.degrees(math.acos(max(-1.0, min(1.0, c))))
    return theta_sigma_pi - 90.0


def n_neighbor_hydrogens(
    atoms: Atoms,
    index: int,
    radius: float = NEIGHBOR_H_RADIUS_A,
    ch_cutoff: float = C_H_CUTOFF_A,
) -> int:
    """Hydrogen atoms bonded to carbons other than ``index`` within ``radius`` of ``index``."""
    i, j, d = neighbor_list("ijd", atoms, radius)
    near_carbons = {int(k) for k in j[(i == index) & (atoms.numbers[j] == 6)]}
    near_carbons.discard(index)
    ci, hj, dd = neighbor_list("ijd", atoms, ch_cutoff)
    bonded = (atoms.numbers[ci] == 6) & (atoms.numbers[hj] == 1)
    hydrogens = {
        int(h) for c, h in zip(ci[bonded], hj[bonded], strict=True) if int(c) in near_carbons
    }
    return len(hydrogens)


def distance_to_point_nm(atoms: Atoms, index: int, point: NDArray[np.float64]) -> float:
    """Minimum-image distance (nm) from atom ``index`` to a point in the cell."""
    D = np.asarray(point, dtype=float)[None, :] - atoms.positions[index][None, :]
    _, dist = find_mic(D, atoms.cell, atoms.pbc)
    return float(dist[0]) / 10.0


def curvature_inv_nm(atoms: Atoms) -> float:
    d = atoms.info.get("diameter_nm")
    return 0.0 if not d else 1.0 / float(d)


def site_descriptors(
    atoms: Atoms, index: int, defect_center: NDArray[np.float64] | None = None
) -> dict[str, float | int | None]:
    """All descriptors of the carbon ``index`` in one dictionary (DATA_SCHEMA SiteDescriptors)."""
    pyr = pyramidalization_angle_deg(atoms, index)
    dist = None if defect_center is None else distance_to_point_nm(atoms, index, defect_center)
    return {
        "pyramidalization_deg": None if math.isnan(pyr) else pyr,
        "coordination": coordination(atoms, index),
        "n_neighbor_H": n_neighbor_hydrogens(atoms, index),
        "distance_to_defect_nm": dist,
        "curvature_inv_nm": curvature_inv_nm(atoms),
    }


__all__ = [
    "carbon_neighbors",
    "coordination",
    "pyramidalization_angle_deg",
    "n_neighbor_hydrogens",
    "distance_to_point_nm",
    "curvature_inv_nm",
    "site_descriptors",
]
