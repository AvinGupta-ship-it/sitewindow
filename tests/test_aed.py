"""Adsorption-energy-distribution inversion (C6): kernel, windows, synthetic data, pooling."""

import numpy as np
import pytest

from sitewindow.constants import EV_TO_KJ_PER_MOL, R_GAS
from sitewindow.experimental import aed
from sitewindow.theory.entropy import ConstantEntropy
from sitewindow.theory.window import W1, W3, occupancy, optimum_energy


def test_energy_grid_matches_the_manual():
    g = aed.energy_grid()
    assert g[0] == 2.0 and g[-1] == 40.0 and len(g) == 77
    assert np.allclose(np.diff(g), 0.5)
    assert len(aed.energy_grid(1.0)) == 39


def test_kernel_reproduces_the_langmuir_occupancy_of_the_theory_module():
    g = aed.energy_grid()
    for T, P in [(77.0, 0.5), (298.15, 30.0)]:
        th = aed.occupancy_matrix(g, np.array([T]), np.array([P]))[0]
        ref = occupancy(-g / EV_TO_KJ_PER_MOL, T, P, ConstantEntropy(), use_fugacity=False)
        assert np.allclose(th, ref, atol=1e-12)
    # An entropy offset of +R raises ln K by one.
    a = aed.log_kernel(g, np.array([100.0]), np.array([1.0]), 0.0)
    b = aed.log_kernel(g, np.array([100.0]), np.array([1.0]), R_GAS)
    assert np.allclose(b - a, 1.0)


def test_window_band_contains_the_optimum_energy_and_is_contiguous():
    g = aed.energy_grid()
    for w in (W1, W3):
        mask = aed.window_mask(g, w)
        eps_star = -optimum_energy(w, ConstantEntropy())[0] * EV_TO_KJ_PER_MOL
        inside = g[mask]
        assert inside.min() - 0.5 <= eps_star <= inside.max() + 0.5
        idx = np.flatnonzero(mask)
        assert np.all(np.diff(idx) == 1)
    # W3 (77 K) sits at weak energies, W1 (298 K) at 7-17 kJ/mol.
    assert g[aed.window_mask(g, W3)].max() < g[aed.window_mask(g, W1)].min()


def test_synthetic_material_follows_the_planted_distribution_and_error_model():
    g = aed.energy_grid()
    w = aed.planted_weights(g, [(5.0, 0.7), (12.0, 0.3)])
    assert abs(w.sum() - 1.0) < 1e-12 and w[g == 5.0] == 0.7 and w[g == 12.0] == 0.3
    ms = aed.synthetic_material(w, g, n_m=10.0, pore_volume_cm3_g=0.0, temperatures=(77.0, 298.0), seed=0)
    assert ms.n_temperatures == 2 and all(i.n_points == 20 for i in ms.isotherms)
    iso77 = ms.isotherms[0]
    # Without noise the 77 K loading at 50 bar approaches n_m (both modes nearly saturated).
    theta = aed.occupancy_matrix(g, np.array([77.0]), np.array([50.0]))[0]
    assert abs(iso77.loading_mol_kg[-1] - 10.0 * theta @ w) < 1.0
    data = aed.material_data(ms)
    assert data.multi_temperature and data.y.shape == (40,)
    assert data.dlnP.max() <= aed.ErrorScales().dig_pressure_cap + 1e-12
    assert data.dlnP.min() >= aed.ErrorScales().dig_pressure_rel - 1e-12


def test_mode_energies_split_the_grid_at_the_given_energies():
    g = aed.energy_grid()
    w = aed.planted_weights(g, [(5.0, 0.7), (12.0, 0.3)])
    e = aed.mode_energies(w[None, :], g, [8.5])
    assert np.allclose(e[0], [5.0, 12.0])
    assert np.allclose(aed.mode_energies(w[None, :], g, [])[0], [0.7 * 5 + 0.3 * 12])


def test_pooled_scales_weight_precise_estimates_and_leave_one_source_out():
    def rec(doi, a, b, spread):
        return {
            "metadata": {"doi": doi},
            "error_scales": {
                "estimated": True,
                "a_dig": {"q50": a, "q05": a / spread, "q95": a * spread},
                "b_dig": {"q50": b, "q05": b / spread, "q95": b * spread},
            },
        }

    stage1 = [rec("x", 0.04, 0.05, 1.2), rec("y", 0.08, 0.10, 3.0)]
    pooled = aed.pooled_scales(stage1)
    # The precise estimate dominates the precision-weighted mean of the logarithms.
    assert 0.04 < pooled.a_dig < 0.05 and 0.05 < pooled.b_dig < 0.06
    loo = aed.pooled_scales(stage1, exclude_doi="x")
    assert abs(loo.a_dig - 0.08) < 1e-9 and "excluding x" in loo.source
    assert aed.pooled_scales([]).source.startswith("manual defaults")


def test_outlier_classification_uses_metadata_features():
    base = {
        "strong_fraction": {"q50": 0.3},
        "label": "single_T_entropy_fixed",
        "temperatures_K": [298.0],
        "ppc": {"chi2_p_value": 0.5},
        "metadata": {"digitized_only": True, "generic_name": True, "max_pressure_bar": 1.0,
                     "quality_flags": ["nonmonotonic"], "material_class": "activated_carbon"},
    }
    flags = aed.classify_outlier(base)
    assert {"single_T_entropy_fixed", "digitized_only", "generic_adsorbent_name", "low_coverage_only"} <= set(flags)
    assert any(f.startswith("quality_flags") for f in flags)
    assert aed.classify_outlier({**base, "strong_fraction": {"q50": 0.05}}) == []


@pytest.mark.slow
def test_single_site_recovery_with_pymc():
    pytest.importorskip("pymc")
    g = aed.energy_grid()
    w0 = aed.planted_weights(g, [(6.0, 1.0)])
    ms = aed.synthetic_material(w0, g, temperatures=(77.0, 87.0), seed=1)
    res = aed.fit_material(
        ms, draws=300, tune=400, chains=2, cores=2, seed=1, entropy_free=False
    )
    e = aed.mode_energies(res.draws, g, [])[:, 0]
    assert np.mean(np.abs(e - 6.0) <= 1.0) >= 0.9
    assert res.strong_fraction["q50"] < 0.05 and abs(res.n_m["q50"] - 10.0) < 1.5
    assert res.diagnostics["divergences"] <= 5 and 0.0 <= res.ppc["chi2_p_value"] <= 1.0
    assert set(res.window_occupancy) == {"W1", "W2", "W3"}


def test_resolvability_of_planted_energies_by_design():
    P = np.logspace(-2, np.log10(50.0), 20)
    # A 12 kJ/mol site is saturated at every cryogenic point but steps within reach at 298 K.
    assert not aed.resolvable(12.0, (77.0, 87.0), P)
    assert aed.resolvable(12.0, (77.0, 298.0), P)
    assert aed.resolvable(6.0, (77.0,), P)
    # A 40 kJ/mol site is saturated at every point of the 77 and 298 K design.
    assert not aed.resolvable(40.0, (77.0, 298.0), P)
