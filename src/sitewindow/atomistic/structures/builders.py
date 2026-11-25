"""Deterministic construction of defect-class structures (graphene supercells, vacancies,
Stone-Wales rotation, nanotubes) and hydrogen placement.

Every builder is a pure function of its parameters so that a (class, recipe) pair maps to one
structure; the structure hash is the cache key for all calculations.
"""

from __future__ import annotations

import hashlib

import numpy as np
from ase import Atoms
from ase.build import graphene, nanotube
from ase.neighborlist import natural_cutoffs, neighbor_list

from sitewindow.constants import GRAPHENE_LATTICE_A

VACUUM_A = 15.0


def graphene_supercell(
    n: int, m: int | None = None, vacuum: float = VACUUM_A, a: float = GRAPHENE_LATTICE_A
) -> Atoms:
    """n x m graphene supercell (orthogonal to z), periodic in-plane with vacuum along z."""
    m = n if m is None else m
    g = graphene(a=a, size=(n, m, 1), vacuum=vacuum / 2.0)
    g.pbc = (True, True, False)
    g.center(axis=2)
    return g


def _neighbors(atoms: Atoms, cutoff: float = 1.7) -> tuple[np.ndarray, np.ndarray]:
    i, j = neighbor_list("ij", atoms, cutoff)
    return i, j


def central_atom_index(atoms: Atoms) -> int:
    """Index of the carbon closest to the in-plane cell centre."""
    cell = atoms.cell
    centre = 0.5 * (cell[0] + cell[1])
    centre[2] = atoms.positions[:, 2].mean()
    d = np.linalg.norm(atoms.positions - centre, axis=1)
    return int(np.argmin(d))


def monovacancy(n: int = 6, vacuum: float = VACUUM_A) -> tuple[Atoms, int]:
    """Graphene with one carbon removed near the centre; returns (atoms, removed_index)."""
    g = graphene_supercell(n, vacuum=vacuum)
    idx = central_atom_index(g)
    removed_position = g.positions[idx].copy()
    del g[idx]
    g.info["removed_position"] = removed_position
    return g, idx


def divacancy_585(n: int = 7, vacuum: float = VACUUM_A) -> Atoms:
    """Graphene with two adjacent carbons removed (5-8-5 after relaxation)."""
    g = graphene_supercell(n, vacuum=vacuum)
    i0 = central_atom_index(g)
    i, j = _neighbors(g)
    partner = int(j[i == i0][0])
    for k in sorted([i0, partner], reverse=True):
        del g[k]
    return g


def stone_wales(n: int = 6, vacuum: float = VACUUM_A) -> Atoms:
    """Graphene with the central C-C bond rotated by 90 degrees about its midpoint."""
    g = graphene_supercell(n, vacuum=vacuum)
    i0 = central_atom_index(g)
    i, j = _neighbors(g)
    i1 = int(j[i == i0][0])
    p0, p1 = g.positions[i0].copy(), g.positions[i1].copy()
    mid = 0.5 * (p0 + p1)
    v = p1 - p0
    rot = np.array([-v[1], v[0], 0.0])
    g.positions[i0] = mid - 0.5 * rot
    g.positions[i1] = mid + 0.5 * rot
    g.info["sw_pair"] = (i0, i1)
    return g


def swcnt(n: int, m: int, length: int = 3, vacuum: float = VACUUM_A) -> Atoms:
    """Single-walled nanotube (n, m) with ``length`` unit cells, periodic along z."""
    t = nanotube(n, m, length=length, bond=1.42, vacuum=vacuum / 2.0)
    t.pbc = (False, False, True)
    t.center()
    return t


def tube_diameter_nm(n: int, m: int, a: float = GRAPHENE_LATTICE_A) -> float:
    return a * np.sqrt(n * n + n * m + m * m) / np.pi / 10.0


def add_hydrogen_top(
    atoms: Atoms, carbon_index: int, height: float = 1.1, outward: bool = True
) -> Atoms:
    """Place one H atom above carbon ``carbon_index`` along the local surface normal.

    For planar sheets the normal is +z; for tubes the normal is radial (outward or inward).
    """
    a = atoms.copy()
    p = a.positions[carbon_index]
    if a.pbc[2] and not a.pbc[0]:  # nanotube along z
        axis_centre = a.positions[:, :2].mean(axis=0)
        r = p[:2] - axis_centre
        nrm = np.array([r[0], r[1], 0.0]) / np.linalg.norm(r)
        if not outward:
            nrm = -nrm
    else:
        nrm = np.array([0.0, 0.0, 1.0])
    a.append("H")
    a.positions[-1] = p + height * nrm
    return a


def hydrogen_molecule(box: float = 15.0) -> Atoms:
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[box, box, box], pbc=False)
    h2.center()
    return h2


def hydrogen_atom(box: float = 15.0) -> Atoms:
    h = Atoms("H", positions=[[0, 0, 0]], cell=[box, box, box], pbc=False)
    h.center()
    h.set_initial_magnetic_moments([1.0])
    return h


def ortho_para_meta_partners(atoms: Atoms, carbon_index: int) -> dict[str, int]:
    """Indices of ortho (nearest), meta (second), and para (third) neighbours on a hexagon."""
    i, j = _neighbors(atoms)
    first = set(j[i == carbon_index].tolist())
    second = set()
    for k in first:
        second |= set(j[i == k].tolist())
    second -= first | {carbon_index}
    third = set()
    for k in second:
        third |= set(j[i == k].tolist())
    third -= first | second | {carbon_index}

    # para partner across the hexagon: among third neighbours, the one at ~2.84 A
    def dist(k: int) -> float:
        return float(np.linalg.norm(atoms.get_distance(carbon_index, k, mic=True)))

    para = min(third, key=dist) if third else -1
    return {"ortho": int(sorted(first)[0]), "meta": int(sorted(second)[0]), "para": int(para)}


def structure_hash(atoms: Atoms) -> str:
    """SHA-256 of symbols, rounded positions, cell and pbc."""
    h = hashlib.sha256()
    h.update("".join(atoms.get_chemical_symbols()).encode())
    h.update(np.round(atoms.positions, 4).tobytes())
    h.update(np.round(np.array(atoms.cell), 4).tobytes())
    h.update(np.array(atoms.pbc).tobytes())
    return h.hexdigest()[:16]


__all__ = [
    "graphene_supercell",
    "monovacancy",
    "divacancy_585",
    "stone_wales",
    "swcnt",
    "tube_diameter_nm",
    "add_hydrogen_top",
    "hydrogen_molecule",
    "hydrogen_atom",
    "ortho_para_meta_partners",
    "central_atom_index",
    "structure_hash",
    "natural_cutoffs",
]
