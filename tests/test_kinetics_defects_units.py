"""Tests for kinetics (DERIVATIONS 4), knock-on physics (DERIVATIONS 6-7), Raman relations (8), units."""

import numpy as np
import pytest

from sitewindow.constants import K_B_EV
from sitewindow.defects.knockon import (
    debye_mean_square_velocity,
    mckinley_feshbach_barn,
    t_max_eV,
    thermally_smeared_cross_section_barn,
)
from sitewindow.defects.raman import (
    LD_low_defect_regime_nm,
    RamanParameters,
    defects_per_carbon_from_LD,
    id_idprime_mixture,
    id_ig_ratio,
)
from sitewindow.defects.rates import integrate_population
from sitewindow.kinetics.polanyi_wigner import (
    redhead_peak_temperature,
    release_fraction_isothermal,
    t90_first_order,
    t_release,
    tpd_linear_ramp,
)
from sitewindow.units import (
    UnitError,
    loading_to_mol_per_kg,
    mol_per_kg_to_wt_percent,
    pressure_to_bar,
)


def test_release_fraction_closed_forms():
    f1 = float(release_fraction_isothermal(1.0, 1e13, 1, 500.0, 1000.0))
    k = 1e13 * np.exp(-1.0 / (K_B_EV * 500.0))
    assert abs(f1 - (1 - np.exp(-k * 1000.0))) < 1e-12
    f2 = float(release_fraction_isothermal(1.0, 1e13, 2, 500.0, 1000.0, 0.5))
    assert abs(f2 - (k * 0.5 * 1000.0) / (1 + k * 0.5 * 1000.0)) < 1e-12


def test_t90_consistency():
    T = float(t90_first_order(1.2, 1e13, 1000.0))
    assert abs(float(release_fraction_isothermal(1.2, 1e13, 1, T, 1000.0)) - 0.9) < 1e-9
    assert abs(t_release(1.2, 1e13, 1, 1000.0) - T) < 1e-4


def test_redhead_peak_matches_ode_first_order():
    """D1: ODE peak within 1 K of the Redhead relation."""
    E, nu, beta = 1.3, 1e13, 1.0
    T, theta, rate = tpd_linear_ramp(E, nu, 1, beta, T0=200.0, T1=900.0, n_points=7001)
    Tp_ode = T[int(np.argmax(rate))]
    Tp_red = redhead_peak_temperature(E, nu, beta, 1)
    assert abs(Tp_ode - Tp_red) < 1.0
    assert abs(theta[0] - 1.0) < 1e-6 and theta[-1] < 1e-3


def test_second_order_peak_shifts_with_initial_coverage():
    E, nu, beta = 1.3, 1e13, 1.0
    T, _, r_full = tpd_linear_ramp(E, nu, 2, beta, theta0=1.0, T0=200.0, T1=1000.0, n_points=4001)
    _, _, r_half = tpd_linear_ramp(E, nu, 2, beta, theta0=0.25, T0=200.0, T1=1000.0, n_points=4001)
    assert T[int(np.argmax(r_half))] > T[int(np.argmax(r_full))]


def test_mckinley_feshbach_static_onset_and_scaling():
    """E1/B9: T_max(E) = 22 eV at ~108.7 keV; zero below; positive and increasing above."""
    from scipy.optimize import brentq

    E_on = brentq(lambda E: float(t_max_eV(E)) - 22.0, 50.0, 200.0)
    assert abs(E_on - 108.7) < 0.5
    assert float(mckinley_feshbach_barn(100.0, 22.0)) == 0.0
    s110, s120 = (
        float(mckinley_feshbach_barn(110.0, 22.0)),
        float(mckinley_feshbach_barn(120.0, 22.0)),
    )
    assert 0.0 < s110 < s120


def test_mckinley_feshbach_independent_evaluation():
    """E1: compare against a direct re-implementation of the closed form at several energies."""
    from sitewindow.constants import ALPHA_FS, R_E_M

    for E in (110.0, 130.0, 160.0, 200.0, 300.0):
        Ek = E * 1e3
        gamma = 1 + Ek / 510998.95
        b2 = 1 - 1 / gamma**2
        b = np.sqrt(b2)
        Mc2 = 12.0 * 931494102.42
        tmax = 2 * Ek * (Ek + 2 * 510998.95) / Mc2
        x = tmax / 22.0
        a = 6 * ALPHA_FS
        ref = (
            np.pi
            * 36
            * R_E_M**2
            * (1 - b2)
            / b2**2
            * ((x - 1) - b2 * np.log(x) + np.pi * a * b * (2 * (np.sqrt(x) - 1) - np.log(x)))
            / 1e-28
        )
        assert abs(float(mckinley_feshbach_barn(E, 22.0)) / ref - 1.0) < 1e-6


def test_debye_velocity_zero_point_limit():
    v2 = debye_mean_square_velocity(1287.0, 0.0)
    from sitewindow.constants import AMU_KG, K_B_J

    assert abs(v2 - 3 * K_B_J * 1287.0 / (8 * 12.0 * AMU_KG)) / v2 < 1e-6
    assert debye_mean_square_velocity(1287.0, 300.0) > v2


def test_smeared_cross_section_finite_below_onset():
    v2 = debye_mean_square_velocity(1287.0, 300.0)
    s100 = thermally_smeared_cross_section_barn(100.0, 22.0, v2)
    s90 = thermally_smeared_cross_section_barn(90.0, 22.0, v2)
    assert 0.0 < s90 < s100 < 1.0


def test_rate_equations_conserve_carbon():
    """E3: conservation residual below 1e-10."""
    out = integrate_population(
        1e4, 100.0, 300.0, sigma_V_barn=0.2, sigma_DV_barn=0.01, sigma_SW_barn=0.001
    )
    assert np.max(np.abs(out["conservation_residual"])) < 1e-10
    assert np.all(np.diff(out["x_total"]) >= -1e-12)


def test_raman_low_defect_limit_matches_published_factor():
    """E4/B10: beyond 20 nm the full model reduces to L_D^2 = 4.3e3 E_L^-4 / (I_D/I_G) within 5 %."""
    p = RamanParameters()
    for L in (20.0, 30.0, 50.0):
        r = float(id_ig_ratio(L, 2.41, p))
        assert abs(float(LD_low_defect_regime_nm(r, 2.41)) / L - 1.0) < 0.05


def test_raman_maximum_near_3nm():
    L = np.linspace(1.0, 20.0, 2000)
    r = id_ig_ratio(L, 2.41)
    assert 2.0 < L[int(np.argmax(r))] < 4.0


def test_defect_density_conversion():
    x = defects_per_carbon_from_LD(1.0)
    assert abs(x - 1 / (np.pi * 38.2)) < 1e-12
    assert abs(id_idprime_mixture({"sp3": 1.0}) - 13.0) < 1e-12


@pytest.mark.parametrize(
    "unit,val,expected",
    [
        ("mmol/g", 10.0, 10.0),
        ("mol/kg", 10.0, 10.0),
        ("wt%", 2.0, (0.02 / 0.98) / 2.01588 * 1000.0),
        ("cm3(STP)/g", 224.1397, 10.0),
        ("mg/g", 20.1588, 10.0),
        ("mmol/kg", 1000.0, 1.0),
    ],
)
def test_loading_units(unit, val, expected):
    assert abs(float(loading_to_mol_per_kg(val, unit)) - expected) < 1e-6 * max(1.0, expected)


def test_wt_percent_round_trip():
    n = np.array([0.5, 5.0, 40.0])
    w = mol_per_kg_to_wt_percent(n)
    back = loading_to_mol_per_kg(w, "wt%")
    assert np.allclose(back, n)


@pytest.mark.parametrize(
    "unit,factor", [("bar", 1.0), ("kPa", 0.01), ("MPa", 10.0), ("atm", 1.01325), ("Pa", 1e-5)]
)
def test_pressure_units(unit, factor):
    assert abs(float(pressure_to_bar(2.0, unit)) - 2.0 * factor) < 1e-12


def test_unit_errors():
    with pytest.raises(UnitError):
        loading_to_mol_per_kg(1.0, "g/l")
    with pytest.raises(UnitError):
        pressure_to_bar(1.0, "furlongs")
