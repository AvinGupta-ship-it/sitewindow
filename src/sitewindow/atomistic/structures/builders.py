"""Deterministic construction of defect-class structures and hydrogen placement.

Graphene supercells, mono- and divacancies (5-8-5 and 555-777 starting geometries), the
Stone-Wales rotation, nanoribbons with zigzag and armchair edges, single-walled nanotubes with
and without defects, oxygen functionalities, the coronene molecule used by benchmark B3, and slit
pores. Every builder is a pure function of its parameters so that a (class, recipe) pair maps to
one structure; the structure hash is the cache key for all calculations.

Each structure records its geometry type in ``atoms.info['geometry']`` (``sheet``, ``tube``,
``molecule``, ``pore``) and, for defects, the defect centre in ``atoms.info['defect_center']``.
"""

from __future__ import annotations

import hashlib
from collections import deque
from collections.abc import Sequence

import numpy as np
from ase import Atoms
from ase.build import graphene, graphene_nanoribbon, nanotube
from ase.geometry import find_mic
from ase.neighborlist import natural_cutoffs, neighbor_list

from sitewindow.constants import GRAPHENE_CC_BOND_A, GRAPHENE_LATTICE_A

VACUUM_A = 15.0
C_C_CUTOFF_A = 1.75
C_H_BOND_A = 1.09
H_HEIGHT_A = 1.1

Side = str  # "up" (outward for tubes) or "down" (inward for tubes)


# --------------------------------------------------------------------------- graph helpers


def _neighbors(atoms: Atoms, cutoff: float = C_C_CUTOFF_A) -> tuple[np.ndarray, np.ndarray]:
    i, j = neighbor_list("ij", atoms, cutoff)
    return i, j


def neighbor_map(atoms: Atoms, cutoff: float = C_C_CUTOFF_A) -> dict[int, list[int]]:
    """Carbon-carbon neighbour lists (minimum image), sorted for determinism."""
    i, j = _neighbors(atoms, cutoff)
    carbon = atoms.numbers == 6
    out: dict[int, list[int]] = {int(k): [] for k in np.flatnonzero(carbon)}
    for a, b in zip(i, j, strict=True):
        if carbon[a] and carbon[b]:
            out[int(a)].append(int(b))
    return {k: sorted(v) for k, v in out.items()}


def mic_vector(atoms: Atoms, i: int, j: int) -> np.ndarray:
    """Minimum-image vector from atom ``i`` to atom ``j``."""
    D = (atoms.positions[j] - atoms.positions[i])[None, :]
    v, _ = find_mic(D, atoms.cell, atoms.pbc)
    return np.asarray(v[0], dtype=float)


def carbon_indices(atoms: Atoms) -> np.ndarray:
    return np.flatnonzero(atoms.numbers == 6)


def central_atom_index(atoms: Atoms) -> int:
    """Index of the carbon closest to the geometric centre of the structure."""
    carbons = carbon_indices(atoms)
    centre = 0.5 * (atoms.cell[0] + atoms.cell[1] + atoms.cell[2])
    for axis in range(3):
        if not atoms.pbc[axis]:
            centre[axis] = atoms.positions[carbons, axis].mean()
    d = np.linalg.norm(atoms.positions[carbons] - centre, axis=1)
    return int(carbons[int(np.argmin(d))])


def sublattices(atoms: Atoms) -> np.ndarray:
    """Two-colouring of the carbon network (0/1 per atom; -1 for non-carbon atoms)."""
    nb = neighbor_map(atoms)
    colour = -np.ones(len(atoms), dtype=int)
    for start in sorted(nb):
        if colour[start] >= 0:
            continue
        colour[start] = 0
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for b in nb[a]:
                if colour[b] < 0:
                    colour[b] = 1 - colour[a]
                    queue.append(b)
    return colour


def hexagon(atoms: Atoms, c0: int, nb: dict[int, list[int]] | None = None) -> list[int]:
    """Six carbons of a hexagon containing ``c0`` in ring order [c0, ortho, meta, para, meta', ortho']."""
    nb = nb or neighbor_map(atoms)
    first = nb[c0]
    for n1 in first:
        for n5 in first:
            if n5 <= n1:
                continue
            for n2 in nb[n1]:
                if n2 == c0:
                    continue
                for n4 in nb[n5]:
                    if n4 == c0 or n4 == n2:
                        continue
                    common = set(nb[n2]) & set(nb[n4]) - {n1, n5, c0}
                    if common:
                        n3 = min(common)
                        return [c0, n1, n2, n3, n4, n5]
    raise ValueError(f"no hexagon found around atom {c0}")


def ortho_para_meta_partners(atoms: Atoms, carbon_index: int) -> dict[str, int]:
    """Indices of ortho (nearest), meta (second) and para (opposite) partners on one hexagon."""
    h = hexagon(atoms, carbon_index)
    return {"ortho": h[1], "meta": h[2], "para": h[3]}


# --------------------------------------------------------------------------- sheets


def graphene_supercell(
    n: int, m: int | None = None, vacuum: float = VACUUM_A, a: float = GRAPHENE_LATTICE_A
) -> Atoms:
    """n x m graphene supercell, periodic in-plane with ``vacuum`` along z."""
    m = n if m is None else m
    g = graphene(a=a, size=(n, m, 1), vacuum=vacuum / 2.0)
    g.pbc = (True, True, False)
    g.center(axis=2)
    g.info.update(geometry="sheet", supercell=(n, m))
    return g


def monovacancy(n: int = 6, vacuum: float = VACUUM_A) -> tuple[Atoms, int]:
    """Graphene with one carbon removed near the centre; returns (atoms, removed_index)."""
    g = graphene_supercell(n, vacuum=vacuum)
    idx = central_atom_index(g)
    removed_position = g.positions[idx].copy()
    del g[idx]
    g.info.update(removed_position=removed_position, defect_center=removed_position.copy())
    return g, idx


def _dimer_sides(atoms: Atoms) -> tuple[int, int, list[int], list[int], dict[int, list[int]]]:
    nb = neighbor_map(atoms)
    A = central_atom_index(atoms)
    B = nb[A][0]
    a_side = [k for k in nb[A] if k != B]
    b_side = [k for k in nb[B] if k != A]
    return A, B, a_side, b_side, nb


def _shifted(indices: list[int], deleted: list[int]) -> list[int]:
    """Indices after deleting ``deleted`` from the structure."""
    return [k - sum(1 for d in deleted if d < k) for k in indices]


def _pinch(atoms: Atoms, p: int, q: int, amount: float) -> None:
    v = mic_vector(atoms, p, q)
    u = v / np.linalg.norm(v)
    atoms.positions[p] += amount * u
    atoms.positions[q] -= amount * u


def divacancy_585(n: int = 7, vacuum: float = VACUUM_A, pinch: float = 0.35) -> Atoms:
    """Two adjacent carbons removed; the dangling pairs are moved towards each other so that
    relaxation forms the two pentagons of the 5-8-5 reconstruction."""
    g = graphene_supercell(n, vacuum=vacuum)
    A, B, a_side, b_side, _ = _dimer_sides(g)
    centre = g.positions[A] + 0.5 * mic_vector(g, A, B)
    removed = [g.positions[A].copy(), g.positions[B].copy()]
    if pinch:
        _pinch(g, a_side[0], a_side[1], pinch)
        _pinch(g, b_side[0], b_side[1], pinch)
    recon = _shifted(a_side + b_side, [A, B])
    del g[[A, B]]
    g.info.update(removed_positions=removed, defect_center=centre, reconstruction_carbons=recon)
    return g


def _rotate_bond_in_plane(atoms: Atoms, x: int, y: int) -> None:
    v = mic_vector(atoms, x, y)
    mid = atoms.positions[x] + 0.5 * v
    rot = np.array([-v[1], v[0], 0.0])
    atoms.positions[x] = mid - 0.5 * rot
    atoms.positions[y] = mid + 0.5 * rot


def divacancy_555777(n: int = 7, vacuum: float = VACUUM_A, pinch: float = 0.35) -> Atoms:
    """Divacancy with one octagon bond rotated by 90 degrees: the starting geometry of the
    555-777 reconstruction (three pentagons, three heptagons)."""
    g = graphene_supercell(n, vacuum=vacuum)
    A, B, a_side, b_side, nb = _dimer_sides(g)
    centre = g.positions[A] + 0.5 * mic_vector(g, A, B)
    removed = [g.positions[A].copy(), g.positions[B].copy()]
    pair: tuple[int, int] | None = None
    for a in a_side:
        for y in nb[a]:
            if y == A:
                continue
            for b in b_side:
                for x in nb[b]:
                    if x == B or x == y:
                        continue
                    if y in nb[x]:
                        pair = (x, y)
                        break
                if pair:
                    break
            if pair:
                break
        if pair:
            break
    if pair is None:  # pragma: no cover - impossible on a graphene lattice
        raise ValueError("octagon bond for the 555-777 rotation not found")
    if pinch:
        _pinch(g, a_side[0], a_side[1], pinch)
        _pinch(g, b_side[0], b_side[1], pinch)
    _rotate_bond_in_plane(g, pair[0], pair[1])
    rotated = [g.positions[pair[0]].copy(), g.positions[pair[1]].copy()]
    recon = _shifted(a_side + b_side, [A, B])
    rotated_pair = _shifted([pair[0], pair[1]], [A, B])
    del g[[A, B]]
    g.info.update(
        removed_positions=removed,
        defect_center=centre,
        rotated_positions=rotated,
        reconstruction_carbons=recon,
        rotated_pair=rotated_pair,
    )
    return g


def stone_wales(n: int = 6, vacuum: float = VACUUM_A) -> Atoms:
    """Graphene with the central C-C bond rotated by 90 degrees about its midpoint."""
    g = graphene_supercell(n, vacuum=vacuum)
    nb = neighbor_map(g)
    i0 = central_atom_index(g)
    i1 = nb[i0][0]
    mid = g.positions[i0] + 0.5 * mic_vector(g, i0, i1)
    _rotate_bond_in_plane(g, i0, i1)
    g.info.update(sw_pair=(i0, i1), defect_center=mid)
    return g


# --------------------------------------------------------------------------- edges


def nanoribbon(
    edge: str = "zigzag",
    width: int = 6,
    length: int = 4,
    vacuum: float = VACUUM_A,
    saturated: bool = True,
) -> Atoms:
    """Graphene nanoribbon in the xy plane, periodic along x, with H-terminated edges.

    ``width`` counts zigzag chains (zigzag) or dimer lines (armchair); ``length`` counts unit
    cells along the periodic direction.
    """
    if edge not in ("zigzag", "armchair"):
        raise ValueError("edge must be 'zigzag' or 'armchair'")
    r = graphene_nanoribbon(
        width,
        length,
        type=edge,
        saturated=saturated,
        C_H=C_H_BOND_A,
        C_C=GRAPHENE_CC_BOND_A,
        vacuum=vacuum / 2.0,
    )
    lengths = r.cell.lengths()
    pos = r.positions[:, [2, 0, 1]]
    out = Atoms(
        symbols=r.get_chemical_symbols(),
        positions=pos,
        cell=[lengths[2], lengths[0], lengths[1]],
        pbc=(True, False, False),
    )
    out.center(axis=(1, 2))
    out.info.update(geometry="sheet", edge=edge, supercell=(width, length))
    return out


def edge_carbons(atoms: Atoms, ch_cutoff: float = 1.25) -> list[int]:
    """Carbons bonded to a terminating hydrogen."""
    i, j = neighbor_list("ij", atoms, ch_cutoff)
    return sorted(
        {
            int(a)
            for a, b in zip(i, j, strict=True)
            if atoms.numbers[a] == 6 and atoms.numbers[b] == 1
        }
    )


def subedge_carbons(atoms: Atoms) -> list[int]:
    """Carbons bonded to an edge carbon that are not edge carbons themselves."""
    nb = neighbor_map(atoms)
    edge = set(edge_carbons(atoms))
    return sorted({b for a in edge for b in nb[a] if b not in edge})


def remove_edge_hydrogen(
    atoms: Atoms, carbon_index: int, ch_cutoff: float = 1.25
) -> tuple[Atoms, int]:
    """Remove the hydrogen bonded to ``carbon_index``; returns (atoms, new carbon index)."""
    i, j, d = neighbor_list("ijd", atoms, ch_cutoff)
    hs = [int(b) for a, b in zip(i, j, strict=True) if a == carbon_index and atoms.numbers[b] == 1]
    if not hs:
        raise ValueError(f"carbon {carbon_index} has no terminating hydrogen")
    h = hs[0]
    a = atoms.copy()
    del a[h]
    new_index = carbon_index - 1 if h < carbon_index else carbon_index
    a.info.update(bare_edge_carbon=new_index, defect_center=a.positions[new_index].copy())
    return a, new_index


# --------------------------------------------------------------------------- tubes


def tube_diameter_nm(n: int, m: int, a: float = GRAPHENE_LATTICE_A) -> float:
    return a * np.sqrt(n * n + n * m + m * m) / np.pi / 10.0


def swcnt(n: int, m: int, length: int = 3, vacuum: float = VACUUM_A) -> Atoms:
    """Single-walled nanotube (n, m) with ``length`` unit cells, periodic along z."""
    t = nanotube(n, m, length=length, bond=GRAPHENE_CC_BOND_A, vacuum=vacuum / 2.0)
    t.pbc = (False, False, True)
    t.center()
    t.info.update(
        geometry="tube", chirality=(n, m), diameter_nm=tube_diameter_nm(n, m), length=length
    )
    return t


def scale_tube_radius(atoms: Atoms, scale: float) -> Atoms:
    """Scale the distance of every atom from the tube axis (the z axis through the cell centre)
    by ``scale``; the uniform radial breathing coordinate of a tube."""
    out = atoms.copy()
    centre = np.array([out.cell[0, 0], out.cell[1, 1]]) / 2.0
    out.positions[:, :2] = centre + scale * (out.positions[:, :2] - centre)
    return out


def chirality_for_diameter(diameter_nm: float, chirality: str = "armchair") -> tuple[int, int]:
    """(n, m) of the armchair or zigzag tube whose diameter is closest to ``diameter_nm``."""
    a = GRAPHENE_LATTICE_A
    if chirality == "armchair":
        n = max(3, round(diameter_nm * 10.0 * np.pi / (a * np.sqrt(3.0))))
        return n, n
    if chirality == "zigzag":
        n = max(5, round(diameter_nm * 10.0 * np.pi / a))
        return n, 0
    raise ValueError("chirality must be 'armchair' or 'zigzag'")


def swcnt_for_diameter(diameter_nm: float, chirality: str = "armchair", length: int = 3) -> Atoms:
    n, m = chirality_for_diameter(diameter_nm, chirality)
    return swcnt(n, m, length=length)


def tube_middle_atom(atoms: Atoms) -> int:
    """Carbon nearest the middle of the periodic direction (ties broken by angular position)."""
    z = atoms.positions[:, 2]
    zmid = 0.5 * atoms.cell[2, 2]
    order = np.lexsort(
        (
            np.arctan2(
                atoms.positions[:, 1] - atoms.cell[1, 1] / 2,
                atoms.positions[:, 0] - atoms.cell[0, 0] / 2,
            ),
            np.abs(z - zmid),
        )
    )
    for k in order:
        if atoms.numbers[k] == 6:
            return int(k)
    raise ValueError("no carbon atoms")


def _radial_unit(atoms: Atoms, point: np.ndarray) -> np.ndarray:
    axis_xy = atoms.positions[:, :2].mean(axis=0)
    r = np.array([point[0] - axis_xy[0], point[1] - axis_xy[1], 0.0])
    return r / np.linalg.norm(r)


def swcnt_monovacancy(n: int, m: int, length: int = 3) -> Atoms:
    t = swcnt(n, m, length=length)
    idx = tube_middle_atom(t)
    removed = t.positions[idx].copy()
    del t[idx]
    t.info.update(removed_position=removed, defect_center=removed.copy())
    return t


def swcnt_stone_wales(n: int, m: int, length: int = 3) -> Atoms:
    """Nanotube with one bond at mid-length rotated by 90 degrees about the local radial axis."""
    t = swcnt(n, m, length=length)
    nb = neighbor_map(t)
    i0 = tube_middle_atom(t)
    i1 = nb[i0][0]
    v = mic_vector(t, i0, i1)
    mid = t.positions[i0] + 0.5 * v
    r = _radial_unit(t, mid)
    v_rot = np.cross(r, v) + r * np.dot(r, v)
    t.positions[i0] = mid - 0.5 * v_rot
    t.positions[i1] = mid + 0.5 * v_rot
    t.info.update(sw_pair=(i0, i1), defect_center=mid)
    return t


# --------------------------------------------------------------------------- oxygen groups


def epoxide(n: int = 6, vacuum: float = VACUUM_A, height: float = 1.25) -> Atoms:
    g = graphene_supercell(n, vacuum=vacuum)
    nb = neighbor_map(g)
    c0 = central_atom_index(g)
    c1 = nb[c0][0]
    mid = g.positions[c0] + 0.5 * mic_vector(g, c0, c1)
    g.append("O")
    g.positions[-1] = mid + [0.0, 0.0, height]
    g.info.update(defect_center=mid, functionalized_carbons=(c0, c1), group_atoms=(len(g) - 1,))
    return g


def hydroxyl(n: int = 6, vacuum: float = VACUUM_A) -> Atoms:
    g = graphene_supercell(n, vacuum=vacuum)
    c0 = central_atom_index(g)
    p = g.positions[c0].copy()
    g.append("O")
    g.positions[-1] = p + [0.0, 0.0, 1.49]
    g.append("H")
    g.positions[-1] = p + [0.60, 0.0, 1.49 + 0.76]
    g.info.update(
        defect_center=p, functionalized_carbons=(c0,), group_atoms=(len(g) - 2, len(g) - 1)
    )
    return g


def carbonyl_pair(n: int = 7, vacuum: float = VACUUM_A, co_bond: float = 1.23) -> Atoms:
    """Divacancy hole with two carbonyl oxygens on diagonal dangling carbons."""
    g = graphene_supercell(n, vacuum=vacuum)
    A, B, a_side, b_side, _ = _dimer_sides(g)
    centre = g.positions[A] + 0.5 * mic_vector(g, A, B)
    oxygens = []
    for c, removed in ((a_side[0], A), (b_side[1], B)):
        u = mic_vector(g, c, removed)
        u /= np.linalg.norm(u)
        oxygens.append(g.positions[c] + co_bond * u)
    removed_pos = [g.positions[A].copy(), g.positions[B].copy()]
    carbonyl_carbons = _shifted([a_side[0], b_side[1]], [A, B])
    bare_carbons = _shifted([a_side[1], b_side[0]], [A, B])
    del g[[A, B]]
    for o in oxygens:
        g.append("O")
        g.positions[-1] = o
    g.info.update(
        removed_positions=removed_pos,
        defect_center=centre,
        group_atoms=(len(g) - 2, len(g) - 1),
        carbonyl_carbons=carbonyl_carbons,
        bare_carbons=bare_carbons,
    )
    return g


def carboxyl(n: int = 6, vacuum: float = VACUUM_A) -> Atoms:
    g = graphene_supercell(n, vacuum=vacuum)
    c0 = central_atom_index(g)
    p = g.positions[c0].copy()
    cprime = p + [0.0, 0.0, 1.55]
    o1 = cprime + [0.0, 1.05, 0.62]
    o2 = cprime + [0.0, -1.14, 0.72]
    h = o2 + [0.0, -0.60, 0.76]
    for sym, pos in (("C", cprime), ("O", o1), ("O", o2), ("H", h)):
        g.append(sym)
        g.positions[-1] = pos
    g.info.update(
        defect_center=p, functionalized_carbons=(c0,), group_atoms=tuple(range(len(g) - 4, len(g)))
    )
    return g


# --------------------------------------------------------------------------- molecules and pores


def coronene(box: float = 25.0) -> Atoms:
    """Coronene C24H12 cut from the graphene lattice with 1.09 A C-H bonds (benchmark B3)."""
    g = graphene(a=GRAPHENE_LATTICE_A, size=(8, 8, 1), vacuum=10.0)
    g.pbc = (True, True, False)
    nb = neighbor_map(g)
    c0 = central_atom_index(g)
    centre = g.positions[c0] + mic_vector(g, c0, nb[c0][0]) + mic_vector(g, c0, nb[c0][1])
    D, dist = find_mic(g.positions - centre, g.cell, g.pbc)
    keep = np.flatnonzero(dist < 3.9)
    if len(keep) != 24:  # pragma: no cover - lattice geometry is fixed
        raise ValueError(f"coronene cut selected {len(keep)} carbons")
    pos = D[keep]
    mol = Atoms("C" * 24, positions=pos)
    nbm = neighbor_map(mol)
    for k in range(24):
        if len(nbm[k]) == 2:
            u = pos[k].copy()
            u[2] = 0.0
            u /= np.linalg.norm(u)
            mol.append("H")
            mol.positions[-1] = pos[k] + C_H_BOND_A * u
    mol.cell = [box, box, box]
    mol.pbc = False
    mol.center()
    mol.info.update(geometry="molecule", name="coronene")
    return mol


def slit_pore(width_nm: float, n: int = 5, vacuum: float = VACUUM_A) -> Atoms:
    """Two graphene sheets separated by ``width_nm`` (centre-to-centre), periodic in-plane."""
    g = graphene_supercell(n, vacuum=0.0)
    w = width_nm * 10.0
    top = g.copy()
    top.positions[:, 2] += w
    pore = g + top
    pore.cell[2, 2] = w + vacuum
    pore.pbc = (True, True, False)
    pore.center(axis=2)
    pore.info.update(geometry="pore", width_nm=width_nm, supercell=(n, n))
    return pore


def hydrogen_molecule(box: float = 15.0, pbc: bool = False) -> Atoms:
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[box, box, box], pbc=pbc)
    h2.center()
    h2.info.update(geometry="molecule", name="H2")
    return h2


def hydrogen_atom(box: float = 15.0, pbc: bool = False) -> Atoms:
    h = Atoms("H", positions=[[0, 0, 0]], cell=[box, box, box], pbc=pbc)
    h.center()
    h.set_initial_magnetic_moments([1.0])
    h.info.update(geometry="molecule", name="H")
    return h


# --------------------------------------------------------------------------- hydrogen placement


def surface_normal(atoms: Atoms, index: int, side: Side = "up") -> np.ndarray:
    """Direction of a C-H bond at carbon ``index``.

    ``up``/``down``: the local surface normal (+z or -z for sheets and molecules, radial outward
    or inward for tubes). ``hole``: within the surface towards ``atoms.info['defect_center']``
    (dangling bonds at vacancies). ``edge``: within the sheet away from a ribbon's axis.
    """
    geometry = atoms.info.get("geometry")
    if geometry is None:
        geometry = "tube" if (atoms.pbc[2] and not atoms.pbc[0]) else "sheet"
    p = atoms.positions[index]
    nrm = _radial_unit(atoms, p) if geometry == "tube" else np.array([0.0, 0.0, 1.0])
    if side == "up":
        return nrm
    if side == "down":
        return -nrm
    if side == "hole":
        centre = np.asarray(atoms.info["defect_center"], dtype=float)
        d, _ = find_mic((centre - p)[None, :], atoms.cell, atoms.pbc)
        v = d[0] - np.dot(d[0], nrm) * nrm
        return v / np.linalg.norm(v)
    if side == "edge":
        carbons = carbon_indices(atoms)
        v = p - atoms.positions[carbons].mean(axis=0)
        for axis in range(3):
            if atoms.pbc[axis]:
                v[axis] = 0.0
        v -= np.dot(v, nrm) * nrm
        return v / np.linalg.norm(v)
    raise ValueError(f"unknown side {side!r}")


def add_hydrogen_top(
    atoms: Atoms,
    carbon_index: int,
    height: float = H_HEIGHT_A,
    outward: bool = True,
    side: Side | None = None,
) -> Atoms:
    """Place one H atom above carbon ``carbon_index`` along the local surface normal.

    ``side='up'`` (default) is the outward normal for tubes and +z for sheets; ``side='down'``
    is the opposite face. ``outward=False`` is an alias of ``side='down'``.
    """
    s: Side = side if side is not None else ("up" if outward else "down")
    a = atoms.copy()
    a.append("H")
    a.positions[-1] = atoms.positions[carbon_index] + height * surface_normal(
        atoms, carbon_index, s
    )
    return a


def add_hydrogens(
    atoms: Atoms, placements: Sequence[tuple[int, Side]], height: float = H_HEIGHT_A
) -> Atoms:
    """Add several H atoms; ``placements`` are (carbon index, side) pairs."""
    a = atoms.copy()
    for idx, side in placements:
        a.append("H")
        a.positions[-1] = atoms.positions[idx] + height * surface_normal(atoms, idx, side)
    set_unpaired_electrons(a, len(placements) % 2, hydrogen_first=True)
    return a


def set_unpaired_electrons(atoms: Atoms, n_unpaired: int, hydrogen_first: bool = True) -> None:
    """Set initial magnetic moments so that their sum equals ``n_unpaired``.

    tblite derives the number of unpaired electrons from the sum of the initial magnetic
    moments; GPAW uses them as the starting magnetization. Moments are placed on hydrogen
    atoms first (the last-added atoms) and otherwise on the first carbon atoms.
    """
    magmoms = np.zeros(len(atoms))
    order = list(range(len(atoms)))
    if hydrogen_first:
        order = [k for k in reversed(order) if atoms.numbers[k] == 1] + [
            k for k in order if atoms.numbers[k] != 1
        ]
    for k in order[: max(0, n_unpaired)]:
        magmoms[k] = 1.0
    atoms.set_initial_magnetic_moments(magmoms)
    atoms.info["n_unpaired"] = int(n_unpaired)


BROKEN_SYMMETRY_MOMENT = 1.0


def undercoordinated_carbons(atoms: Atoms, cutoff: float = C_C_CUTOFF_A) -> list[int]:
    """Carbons with fewer than three carbon neighbours within ``cutoff``: the dangling bonds of
    vacancies and the atoms of edges, hydrogen-terminated or bare (the π edge states of a zigzag
    edge carry its magnetism whether or not the σ bond is saturated)."""
    return [k for k, nb in neighbor_map(atoms, cutoff).items() if len(nb) < 3]


def broken_symmetry_moments(
    atoms: Atoms,
    n_unpaired: int,
    moment: float = BROKEN_SYMMETRY_MOMENT,
    hydrogen_first: bool = True,
) -> list[int]:
    """Initial magnetic moments for a spin-polarized density-functional start that can break
    spin symmetry: the nominal ``n_unpaired`` electrons placed as :func:`set_unpaired_electrons`,
    plus ``±moment`` on every carbon with fewer than three carbon neighbours, the sign
    alternating with the sublattice (Lieb's theorem for the bipartite network: opposite spins on
    the two sublattices, the antiferromagnetic coupling of the two zigzag edges, parallel moments
    on the three same-sublattice neighbours of a monovacancy). The total moment is left free in
    the calculation; ``atoms.info['n_unpaired']`` keeps the nominal count. Returns the indices
    that received a guess (empty for a fully three-coordinated network, whose moments stay the
    canonical placement). A self-consistent field started from zero moments cannot leave the
    non-magnetic state of a system whose ground state is magnetic, and one started from a
    one-atom guess finds the state nearest to that guess rather than the lowest.
    """
    set_unpaired_electrons(atoms, n_unpaired, hydrogen_first)
    under = undercoordinated_carbons(atoms)
    if not under:
        return []
    colour = sublattices(atoms)
    magmoms = np.array(atoms.get_initial_magnetic_moments(), dtype=float)
    for k in under:
        magmoms[k] += moment if colour[k] == 0 else -moment
    atoms.set_initial_magnetic_moments(magmoms)
    return under


def graphane_placements(atoms: Atoms) -> list[tuple[int, Side]]:
    """Chair graphane: every carbon hydrogenated, sides alternating with the sublattice."""
    sub = sublattices(atoms)
    return [(int(k), "up" if sub[k] == 0 else "down") for k in carbon_indices(atoms)]


def structure_hash(atoms: Atoms) -> str:
    """SHA-256 of symbols, rounded positions, cell and pbc."""
    h = hashlib.sha256()
    h.update("".join(atoms.get_chemical_symbols()).encode())
    h.update(np.round(atoms.positions, 4).tobytes())
    h.update(np.round(np.array(atoms.cell), 4).tobytes())
    h.update(np.array(atoms.pbc).tobytes())
    return h.hexdigest()[:16]


__all__ = [
    "scale_tube_radius",
    "VACUUM_A",
    "GRAPHENE_LATTICE_A",
    "neighbor_map",
    "mic_vector",
    "carbon_indices",
    "central_atom_index",
    "sublattices",
    "hexagon",
    "ortho_para_meta_partners",
    "graphene_supercell",
    "monovacancy",
    "divacancy_585",
    "divacancy_555777",
    "stone_wales",
    "nanoribbon",
    "edge_carbons",
    "subedge_carbons",
    "remove_edge_hydrogen",
    "tube_diameter_nm",
    "swcnt",
    "chirality_for_diameter",
    "swcnt_for_diameter",
    "tube_middle_atom",
    "swcnt_monovacancy",
    "swcnt_stone_wales",
    "epoxide",
    "hydroxyl",
    "carbonyl_pair",
    "carboxyl",
    "coronene",
    "slit_pore",
    "hydrogen_molecule",
    "hydrogen_atom",
    "surface_normal",
    "add_hydrogen_top",
    "add_hydrogens",
    "set_unpaired_electrons",
    "undercoordinated_carbons",
    "broken_symmetry_moments",
    "BROKEN_SYMMETRY_MOMENT",
    "graphane_placements",
    "structure_hash",
    "natural_cutoffs",
]
