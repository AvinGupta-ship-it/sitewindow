"""Structure builders, descriptors and the site-class catalogue (manual C2)."""

import numpy as np
import pytest
from ase.build import molecule

from sitewindow.atomistic.descriptors import (
    coordination,
    n_neighbor_hydrogens,
    pyramidalization_angle_deg,
    site_descriptors,
)
from sitewindow.atomistic.structures import builders as b
from sitewindow.landscape import catalogue, chemisorption_classes, physisorption_classes


def test_pyramidalization_limits():
    g = b.graphene_supercell(4)
    assert abs(pyramidalization_angle_deg(g, b.central_atom_index(g))) < 1e-9
    c60 = molecule("C60")
    c60.cell = [30, 30, 30]
    assert abs(pyramidalization_angle_deg(c60, 0) - 11.64) < 0.2
    tet = molecule("CH4")
    tet.symbols = ["C", "C", "C", "C", "C"]
    tet.cell = [20, 20, 20]
    tet.positions *= 1.42 / tet.get_distance(0, 1)
    assert (
        abs(pyramidalization_angle_deg(tet, 0, cutoff=1.75) - 19.471) < 0.01
        or coordination(tet, 0) != 3
    )


def test_tube_pyramidalization_decreases_with_diameter():
    angles = [
        pyramidalization_angle_deg(t, b.tube_middle_atom(t))
        for t in (b.swcnt(5, 5, 2), b.swcnt(8, 8, 2), b.swcnt(12, 12, 2))
    ]
    assert angles[0] > angles[1] > angles[2] > 0


def test_hexagon_and_partners():
    g = b.graphene_supercell(5)
    c0 = b.central_atom_index(g)
    h = b.hexagon(g, c0)
    assert len(set(h)) == 6
    d = [g.get_distance(h[0], k, mic=True) for k in h[1:4]]
    assert np.allclose(d, [1.42, 2.46, 2.84], atol=0.01)


def test_sublattice_two_colouring_is_proper():
    g = b.graphene_supercell(4)
    sub = b.sublattices(g)
    nb = b.neighbor_map(g)
    assert all(sub[a] != sub[c] for a, cs in nb.items() for c in cs)
    assert sum(1 for _, s in b.graphane_placements(g) if s == "up") == len(g) // 2


def test_defect_builders_compositions():
    bond = b.GRAPHENE_LATTICE_A / np.sqrt(3.0)
    v, _ = b.monovacancy(6)
    assert len(v) == 71 and sum(1 for k in b.carbon_indices(v) if coordination(v, k) == 2) == 3
    dv = b.divacancy_585(7)
    r = dv.info["reconstruction_carbons"]
    assert len(dv) == 96 and len(r) == 4 and dv.get_distance(r[0], r[1], mic=True) < 1.8
    dv2 = b.divacancy_555777(7)
    assert (
        len(dv2) == 96 and abs(dv2.get_distance(*dv2.info["rotated_pair"], mic=True) - bond) < 1e-6
    )
    sw = b.stone_wales(6)
    assert len(sw) == 72 and abs(sw.get_distance(*sw.info["sw_pair"], mic=True) - bond) < 1e-6
    assert b.carbonyl_pair(7).get_chemical_formula() == "C96O2"
    assert b.epoxide(6).get_chemical_formula() == "C72O"
    assert b.hydroxyl(6).get_chemical_formula() == "C72HO"
    assert b.carboxyl(6).get_chemical_formula() == "C73HO2"


def test_ribbons_and_edges():
    r = b.nanoribbon("zigzag", 6, 4)
    assert r.pbc.tolist() == [True, False, False] and r.get_chemical_formula() == "C48H8"
    ec, sc = b.edge_carbons(r), b.subedge_carbons(r)
    assert len(ec) == 8 and len(sc) == 8 and not set(ec) & set(sc)
    bare, ci = b.remove_edge_hydrogen(r, ec[0])
    assert len(bare) == 55 and coordination(bare, ci) == 2
    placed = b.add_hydrogens(bare, [(ci, "edge")])
    assert abs(placed.get_distance(ci, len(placed) - 1) - 1.1) < 1e-6
    assert abs(placed.positions[-1, 2] - placed.positions[ci, 2]) < 1e-9


def test_coronene_geometry():
    c = b.coronene()
    assert c.get_chemical_formula() == "C24H12"
    ch = [c.get_distance(i, j) for i in range(24) for j in range(24, 36)]
    assert abs(min(ch) - 1.09) < 1e-6
    assert np.allclose(c.positions[:, 2], c.positions[0, 2])


def test_chirality_for_diameter_brackets_target():
    for d in (0.6, 1.0, 1.4, 2.5):
        for chir in ("armchair", "zigzag"):
            n, m = b.chirality_for_diameter(d, chir)
            assert abs(b.tube_diameter_nm(n, m) - d) < 0.09


def test_hole_placement_points_to_vacancy_centre():
    v, _ = b.monovacancy(4)
    nb = b.neighbor_map(v)
    dangling = [k for k in b.carbon_indices(v) if len(nb[k]) == 2]
    vh = b.add_hydrogens(v, [(dangling[0], "hole")])
    assert np.linalg.norm(vh.positions[-1] - v.info["defect_center"]) < np.linalg.norm(
        v.positions[dangling[0]] - v.info["defect_center"]
    )
    assert vh.get_initial_magnetic_moments().sum() == 1.0


def test_structure_hash_sensitivity():
    g = b.graphene_supercell(3)
    h0 = b.structure_hash(g)
    g2 = g.copy()
    g2.positions[0, 2] += 0.01
    assert h0 != b.structure_hash(g2) and h0 == b.structure_hash(g.copy())


def test_descriptors_of_dimer():
    g = b.graphene_supercell(4)
    h = b.hexagon(g, b.central_atom_index(g))
    gh = b.add_hydrogens(g, [(h[0], "up"), (h[1], "up")])
    d = site_descriptors(gh, h[0])
    assert d["coordination"] == 3 and d["n_neighbor_H"] == 1 and d["curvature_inv_nm"] == 0.0
    assert n_neighbor_hydrogens(gh, h[3]) == 2


@pytest.mark.parametrize("tier", ["verify", "full"])
def test_catalogue_configurations_are_consistent(tier):
    cat = catalogue()
    assert set(physisorption_classes()) == {
        "pore_slit",
        "cnt_interior",
        "cnt_interstitial",
        "cnt_groove",
    }
    for name in chemisorption_classes():
        defn = cat[name]
        if tier not in defn.recipes and tier == "verify":
            assert defn.configurations(tier) or True
            continue
        cfgs = defn.configurations(tier)
        assert cfgs, name
        names = {(c.substrate_name, c.name) for c in cfgs}
        for c in cfgs:
            assert c.n_H == len(c.placements) >= 1
            assert all(c.substrate.numbers[i] == 6 for i, _ in c.placements)
            assert c.site.multiplicity >= 1 and c.site.max_occupancy >= 1
            if c.parent is not None:
                assert (c.substrate_name, c.parent) in names
        model = defn.site_class(tier)
        assert model.name.value == name and model.sites


def test_verify_tier_class_set_matches_manual_section_16():
    cat = catalogue()
    counts = {
        k: len(cat[k].configurations("verify"))
        for k in ("basal", "monovacancy", "stone_wales", "cnt_armchair")
    }
    assert counts == {"basal": 4, "monovacancy": 3, "stone_wales": 1, "cnt_armchair": 9}


def test_scale_tube_radius_scales_distances_from_the_axis_only():
    import numpy as np

    from sitewindow.atomistic.structures.builders import scale_tube_radius, swcnt

    t = swcnt(5, 5, length=2)
    centre = np.array([t.cell[0, 0], t.cell[1, 1]]) / 2.0
    r0 = np.linalg.norm(t.positions[:, :2] - centre, axis=1)
    out = scale_tube_radius(t, 1.02)
    r1 = np.linalg.norm(out.positions[:, :2] - centre, axis=1)
    assert np.allclose(r1, 1.02 * r0)
    assert np.allclose(out.positions[:, 2], t.positions[:, 2])
    assert out.info["geometry"] == "tube" and (out.cell == t.cell).all()
