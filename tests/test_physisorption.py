"""Pairwise physisorption potential: identities, Feynman-Hibbs correction, B3 calibration."""

import numpy as np

from sitewindow.atomistic import physisorption as ph
from sitewindow.atomistic.structures import builders as b
from sitewindow.constants import AMU_KG, HBAR, K_B_J, M_H2_AMU


def test_lennard_jones_minimum_and_laplacian():
    p = ph.STEELE_BUCH
    r_min = 2 ** (1 / 6) * p.sigma_A
    assert abs(float(ph.lj_energy_K(np.array([r_min]), p)[0]) + p.epsilon_K) < 1e-9
    r = np.linspace(2.8, 6.0, 50)
    h = 1e-4
    u = lambda x: ph.lj_energy_K(x, p)  # noqa: E731
    lap_num = (u(r + h) - 2 * u(r) + u(r - h)) / h**2 + 2.0 / r * (u(r + h) - u(r - h)) / (2 * h)
    assert np.allclose(ph.lj_laplacian_K_per_A2(r, p), lap_num, rtol=1e-4, atol=1e-6)


def test_feynman_hibbs_coefficient_and_sign():
    c = ph.feynman_hibbs_coefficient_A2(77.0)
    expected = HBAR**2 / (24 * M_H2_AMU * AMU_KG * K_B_J * 77.0) * 1e20
    assert abs(c - expected) < 1e-12
    g = b.graphene_supercell(4)
    e_cl, _ = ph.sheet_site_energy(g, ph.STEELE_BUCH, None)
    e_77, _ = ph.sheet_site_energy(g, ph.STEELE_BUCH, 77.0)
    e_298, _ = ph.sheet_site_energy(g, ph.STEELE_BUCH, 298.15)
    assert e_cl < e_298 < e_77 < 0.0


def test_lorentz_berthelot_values():
    assert abs(ph.STEELE_BUCH.epsilon_K - np.sqrt(28.0 * 34.2)) < 1e-9
    assert abs(ph.STEELE_BUCH.sigma_A - 3.18) < 1e-9
    assert ph.STEELE_BUCH.settings_hash != ph.STEELE_DARKRIM_LEVESQUE.settings_hash


def test_calibration_reproduces_b3_exactly():
    cal = ph.calibrated_potential()
    e, _ = ph.coronene_binding_kJ_mol(cal)
    assert abs(e - ph.B3_CORONENE_KJ_MOL) < 1e-6
    assert cal.sigma_A == ph.STEELE_BUCH.sigma_A and cal.epsilon_K > ph.STEELE_BUCH.epsilon_K


def test_pore_energy_stronger_than_sheet_and_monotone_for_wide_pores():
    p = ph.STEELE_BUCH
    e_sheet, _ = ph.sheet_site_energy(b.graphene_supercell(4), p)
    e_pores = [ph.pore_site_energy(b.slit_pore(w, n=4), p)[0] for w in (0.7, 0.9, 1.2)]
    assert e_pores[0] < e_pores[1] < e_pores[2] < e_sheet < 0.0
    assert e_pores[2] > 2.0 * e_sheet


def test_tube_interior_stronger_than_exterior():
    t = b.swcnt(10, 10, 2)
    e_in, _ = ph.tube_site_energy(t, ph.STEELE_BUCH, True)
    e_out, _ = ph.tube_site_energy(t, ph.STEELE_BUCH, False)
    assert e_in < e_out < 0.0


def test_oxygen_types_follow_the_bonding_environment():
    for build, expected in (
        (b.epoxide, {"ether"}),
        (b.hydroxyl, {"alcohol"}),
        (b.carbonyl_pair, {"carbonyl"}),
    ):
        a = build(4) if build is not b.carbonyl_pair else build(5)
        kinds = {ph.oxygen_type(a, int(k)) for k in np.flatnonzero(a.numbers == 8)}
        assert kinds == expected
    acid = b.carboxyl(4)
    kinds = sorted(ph.oxygen_type(acid, int(k)) for k in np.flatnonzero(acid.numbers == 8))
    assert kinds == ["alcohol", "carbonyl"]
    eps, sig = ph.oxygen_pair(ph.STEELE_BUCH, "ether")
    assert (
        abs(eps - np.sqrt(0.140 * 503.2197 * 34.2)) < 1e-6 and abs(sig - 0.5 * (2.90 + 2.96)) < 1e-9
    )
    cal = ph.calibrated_potential()
    assert ph.oxygen_pair(cal, "ether") == ph.oxygen_pair(ph.STEELE_BUCH, "ether")


def test_oxygen_adjacent_site_differs_from_the_clean_sheet_and_carbon_only_paths_agree():
    g = b.graphene_supercell(4)
    e_clean, _ = ph.sheet_site_energy(g, ph.STEELE_BUCH, None)
    imgs = ph.site_images(g, ph.STEELE_BUCH)
    assert imgs.n_oxygen == 0
    pts = np.array([[3.0, 3.0, g.positions[:, 2].max() + 3.1]])
    e_a = ph.probe_energy_kJ_mol(pts, ph.carbon_images(g), ph.STEELE_BUCH)
    e_b = ph.probe_energy_kJ_mol(pts, imgs, ph.STEELE_BUCH)
    assert abs(float(e_a[0] - e_b[0])) < 1e-9
    epo = b.epoxide(4)
    e_adj, q = ph.sheet_site_energy(epo, ph.STEELE_BUCH, None, near=epo.info["defect_center"])
    assert e_adj < 0.0 and abs(e_adj - e_clean) > 0.05
    centre = np.asarray(epo.info["defect_center"])
    assert np.linalg.norm(q[:2] - centre[:2]) < 6.0
    assert q[2] > epo.positions[:, 2].max() - 1.0
