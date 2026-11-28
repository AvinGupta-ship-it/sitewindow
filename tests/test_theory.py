"""Tests of the window theory against analytic identities (DERIVATIONS 2, 5; acceptance A1-A5)."""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from sitewindow.constants import EV_TO_KJ_PER_MOL, R_GAS
from sitewindow.theory.deliverable import (
    SiteTable,
    bound_wt_percent,
    deliverable,
    random_landscape_bound_check,
)
from sitewindow.theory.entropy import ConstantEntropy, HarmonicOscillatorEntropy, TwoDGasEntropy
from sitewindow.theory.window import (
    BHATIA_MYERS,
    W1,
    W2,
    W3,
    OperatingWindow,
    heterogeneity_threshold,
    occupancy,
    optimum_energy,
    optimum_enthalpy_isothermal_kJ,
    second_virial_lj,
    swing,
    swing_max_isothermal,
)


@pytest.mark.parametrize("r", [2.0, 5.0, 20.0, 100.0])
def test_isothermal_swing_maximum_formula(r):
    """A1: numerical maximum of the swing equals (sqrt r - 1)/(sqrt r + 1) to 1e-8."""
    w = OperatingWindow("t", 298.15, r * 1.0, 298.15, 1.0, fugacity=False)
    _, eta = optimum_energy(w, ConstantEntropy())
    assert abs(eta - swing_max_isothermal(r)) < 1e-8


def test_bhatia_myers_optimum_enthalpy():
    """A2 / B8: -15.1 kJ/mol at 298 K between 30 and 1.5 bar with dS = -8R."""
    dH = optimum_enthalpy_isothermal_kJ(298.15, 30.0, 1.5, -8.0 * R_GAS)
    assert abs(dH - (-15.1)) < 0.1
    eps_star, _ = optimum_energy(BHATIA_MYERS, ConstantEntropy())
    assert abs(eps_star * EV_TO_KJ_PER_MOL - dH) < 0.05


def test_bhatia_myers_carbon_optimum_temperature():
    """A2 / B8: for dH = -5.8 kJ/mol the swing between 30 and 1.5 bar is maximal near 115 K."""
    eps = -5.8 / EV_TO_KJ_PER_MOL
    Ts = np.linspace(60.0, 300.0, 2401)
    vals = [
        float(swing(eps, OperatingWindow("t", T, 30.0, T, 1.5, fugacity=False), ConstantEntropy()))
        for T in Ts
    ]
    T_opt = Ts[int(np.argmax(vals))]
    assert abs(T_opt - 115.0) < 5.0


def test_limiting_cases():
    """A3: r -> 1 gives zero swing; eps -> -inf gives theta -> 1 and zero swing."""
    w = OperatingWindow("t", 298.15, 10.0, 298.15, 10.0, fugacity=False)
    assert abs(float(swing(-0.1, w))) < 1e-12
    assert abs(float(occupancy(-5.0, 298.15, 1.0)) - 1.0) < 1e-9
    assert abs(float(swing(-5.0, W1))) < 1e-6


@given(
    st.floats(min_value=-1.0, max_value=-0.001),
    st.floats(min_value=100.0, max_value=400.0),
    st.floats(min_value=0.1, max_value=200.0),
)
@settings(max_examples=200, deadline=None)
def test_occupancy_bounds_and_monotonicity(eps, T, P):
    th = float(occupancy(eps, T, P, use_fugacity=False))
    assert 0.0 <= th <= 1.0
    th2 = float(occupancy(eps, T, 2.0 * P, use_fugacity=False))
    assert th2 >= th - 1e-12


def test_swing_bounded_by_isothermal_maximum():
    from sitewindow.theory.window import effective_pressure

    eps = np.linspace(-1.0, -0.001, 500)
    s = swing(eps, W1, ConstantEntropy())
    r_eff = effective_pressure(W1.T_c, W1.P_c) / effective_pressure(W1.T_r, W1.P_r)
    assert np.all(s <= swing_max_isothermal(r_eff) + 1e-9)


def test_bound_never_violated_random_landscapes():
    """A4: 10^3 random molecular landscapes never exceed D_max(W1) (full 10^4 in the slow tier)."""
    assert random_landscape_bound_check(W1, np.random.default_rng(12345), n_trials=1000)


def test_heterogeneity_threshold_exists_and_is_unique():
    """A5: for eps1 weaker than eps*, eps_c exists on the strong side and swing(eps_c) = swing(eps1)."""
    rng = np.random.default_rng(7)
    for _ in range(100):
        Pc = rng.uniform(20.0, 150.0)
        Pr = rng.uniform(1.0, 10.0)
        w = OperatingWindow("t", 298.15, Pc, 298.15, Pr, fugacity=False)
        eps_star, _ = optimum_energy(w)
        eps1 = eps_star * 0.5  # weaker than optimum
        eps_c = heterogeneity_threshold(eps1, w)
        assert eps_c is not None and eps_c < eps_star
        assert abs(float(swing(eps_c, w)) - float(swing(eps1, w))) < 1e-6
    assert heterogeneity_threshold(-1.4, W1) is None


def test_entropy_models_finite_and_ordered():
    for m in (ConstantEntropy(), HarmonicOscillatorEntropy(), TwoDGasEntropy()):
        ds = m.delta_s(298.15, 1.0)
        assert np.isfinite(ds) and -150.0 < ds < 0.0


def test_virial_fugacity_fallback_reasonable():
    """A6 (fallback branch): B(298 K) within the 10-16 cm3/mol range of tabulated hydrogen data; B(77 K) negative."""
    assert 10.0 < second_virial_lj(298.15) < 16.0
    assert second_virial_lj(77.0) < 0.0


def test_deliverable_homogeneous_site_at_optimum_recovers_bound():
    """G4 negative control: a homogeneous landscape at the mass-accounted optimum eps_D with the
    monolayer site count hits D_max, and one at the swing optimum eps* stays within 1e-3 of it."""
    from sitewindow.constants import M_C_AMU
    from sitewindow.theory.deliverable import monolayer_capacity_mol_per_kg

    b = bound_wt_percent(W1, "molecular")
    assert abs(b["eps_D_eV"] - b["eps_star_eV"]) < 0.01
    n_sites_per_C = monolayer_capacity_mol_per_kg() / 1000.0 * M_C_AMU
    t_star = SiteTable(
        np.array([b["eps_star_eV"]]),
        np.array([np.nan]),
        np.array([1e13]),
        np.array([1]),
        np.array([n_sites_per_C]),
        np.array([1.0]),
        np.array([0.0]),
        np.array([0]),
    )
    d_star = deliverable(t_star, W1).D_wt_percent
    assert 0.0 <= b["D_max_wt_percent"] - d_star < 1e-3 * b["D_max_wt_percent"]
    t = SiteTable(
        np.array([b["eps_D_eV"]]),
        np.array([np.nan]),
        np.array([1e13]),
        np.array([1]),
        np.array([n_sites_per_C]),
        np.array([1.0]),
        np.array([0.0]),
        np.array([0]),
    )
    res = deliverable(t, W1)
    assert abs(res.D_wt_percent - b["D_max_wt_percent"]) < 1e-6 * max(1.0, b["D_max_wt_percent"])


def test_pristine_graphene_physisorption_small_at_W1():
    """G4 negative control: a -5 kJ/mol basal site delivers a small fraction of its capacity at W1."""
    s = float(swing(-5.0 / EV_TO_KJ_PER_MOL, W1))
    assert 0.0 < s < 0.3


def test_temperature_swing_increases_efficiency():
    _, eta1 = optimum_energy(W1)
    _, eta2 = optimum_energy(W2)
    _, eta3 = optimum_energy(W3)
    assert eta2 > eta1 > 0.5 and 0.5 < eta3 < 0.7


def test_isothermal_efficiency_is_entropy_model_independent():
    """For an isothermal swing K is pressure independent, so eta = (sqrt r - 1)/(sqrt r + 1) for every model."""
    w = OperatingWindow("t", 298.15, 30.0, 298.15, 1.5, fugacity=False)
    for m in (ConstantEntropy(), HarmonicOscillatorEntropy(), TwoDGasEntropy()):
        _, eta = optimum_energy(w, m)
        assert abs(eta - swing_max_isothermal(20.0)) < 1e-7
