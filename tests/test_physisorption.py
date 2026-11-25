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
