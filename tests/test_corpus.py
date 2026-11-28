"""Corpus inference (C6): equation of state, material sets, isosteric heats, benchmark runners."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sitewindow.constants import R_GAS
from sitewindow.experimental import eos
from sitewindow.experimental.corpus import (
    MaterialSet,
    load_corpus,
    material_id,
    material_sets,
    multi_temperature,
)
from sitewindow.experimental.isosteric import corpus_isosteric_heats, isosteric_heat
from sitewindow.units import loading_to_mol_per_kg, pressure_to_bar

CORPUS = Path("data/processed/isodb_carbon_h2.parquet")


def _langmuir_table(q_kJ: float, temps=(77.0, 87.0), n_m=10.0, ds_R=-8.0, n_pts=14, noise=0.0):
    """Langmuir isotherms with site energy ``q_kJ`` and constant entropy ``ds_R`` R."""
    rows = []
    rng = np.random.default_rng(3)
    for k, T in enumerate(temps):
        P = np.logspace(-2.5, 1.5, n_pts)
        K = np.exp(q_kJ * 1000.0 / (R_GAS * T) + ds_R)
        n = n_m * K * P / (1.0 + K * P)
        n = n * (1.0 + noise * rng.standard_normal(n_pts))
        rows.append(
            {
                "isodb_filename": f"synthetic.Isotherm{k}",
                "doi": "10.0000/synthetic",
                "adsorbent_name": "Synthetic carbon",
                "adsorbent_hashkey": "NIST-MATDB-0000000000000000",
                "temperature_K": T,
                "pressure_bar": P,
                "loading_mol_kg": n,
                "digitized": True,
                "tabular": False,
                "quality_flags": "",
                "category": "",
                "material_class": "activated_carbon",
            }
        )
    return pd.DataFrame(rows)


def test_h2_density_ideal_gas_limit_and_fallback():
    rho = eos.h2_density_mol_m3(298.15, [0.01])[0]
    assert abs(rho / (1.0e3 / (R_GAS * 298.15)) - 1.0) < 1e-3
    assert eos.h2_density_mol_m3(77.0, [0.0, -1.0]).tolist() == [0.0, 0.0]
    # Pore-gas term: 1 cm^3/g at the ideal density of 1 bar and 77 K is 0.156 mol/kg.
    g = eos.gas_in_pores_mol_kg(77.0, [1.0], 1.0)[0]
    assert abs(g - 1.0e5 / (R_GAS * 77.0) * 1.0e-3) < 0.01
    # Saturation pressure brackets the normal boiling point and vanishes above T_c.
    assert 0.9 < eos.h2_saturation_pressure_bar(20.4) < 1.15
    assert eos.h2_saturation_pressure_bar(40.0) == float("inf")
    # Fugacity is used only above 50 bar.
    f = eos.effective_pressures_bar(77.0, [10.0, 100.0])
    assert f[0] == 10.0 and 80.0 < f[1] < 110.0 and f[1] != 100.0


def test_material_sets_group_by_hashkey_and_doi_and_log_exclusions():
    t = _langmuir_table(6.0)
    bad = t.iloc[[0]].copy()
    bad["isodb_filename"] = "bad.Isotherm"
    bad["pressure_bar"] = [np.array([1e-20, 1.0, 2.0e5])]
    bad["loading_mol_kg"] = [np.array([0.1, 0.2, 0.3])]
    mod = t.iloc[[1]].copy()
    mod["isodb_filename"] = "mod.Isotherm"
    mod["category"] = "mod"
    sub = t.iloc[[1]].copy()
    sub["isodb_filename"] = "sub.Isotherm"
    sub["temperature_K"] = 20.0
    other = t.iloc[[0]].copy()
    other["isodb_filename"] = "other.Isotherm"
    other["doi"] = "10.0000/other"
    table = pd.concat([t, bad, mod, sub, other], ignore_index=True)
    sets, excluded = material_sets(table)
    ids = {s.material_id for s in sets}
    assert len(sets) == 2 and material_id("NIST-MATDB-0000000000000000", "10.0000/other") in ids
    reasons = " ".join(excluded["reason"])
    assert len(excluded) == 3
    assert "pressure axis" in reasons and "modelled" in reasons and "subcritical" in reasons
    assert len(multi_temperature(sets)) == 1


def test_isosteric_heat_recovers_planted_langmuir_heat():
    sets, _ = material_sets(_langmuir_table(6.0))
    res = isosteric_heat(sets[0], n_grid=12, n_boot=30, seed=1)
    assert res.usable
    assert np.allclose(res.table["q_st_kJ_mol"], 6.0, atol=0.05)
    # Noise-free data: the bootstrap spread comes from the interpolation alone.
    assert res.table["q_st_sd_kJ_mol"].median() < 0.3
    # Three temperatures: the regression over 1/T reproduces the same heat.
    sets3, _ = material_sets(_langmuir_table(5.0, temps=(77.0, 87.0, 100.0)))
    res3 = isosteric_heat(sets3[0], n_grid=8, n_boot=10)
    assert np.allclose(res3.table["q_st_kJ_mol"], 5.0, atol=0.05)


def test_isosteric_heat_requires_overlapping_loadings_and_two_temperatures():
    t = _langmuir_table(6.0, temps=(77.0, 298.0))
    # At 298 K the Langmuir loadings stay below the 77 K range: no common loading.
    t.loc[1, "pressure_bar"] = np.logspace(-6, -4, 14)
    t.loc[1, "loading_mol_kg"] = np.linspace(1e-6, 1e-4, 14)
    sets, _ = material_sets(t)
    res = isosteric_heat(sets[0], n_grid=6, n_boot=5)
    assert not res.usable and "overlap" in res.reason
    single, _ = material_sets(_langmuir_table(6.0, temps=(77.0,)))
    assert "two temperatures" in isosteric_heat(single[0]).reason


def test_corpus_isosteric_heats_returns_curves_summary_and_log():
    t = pd.concat([_langmuir_table(6.0), _langmuir_table(4.0, temps=(298.0,))], ignore_index=True)
    t.loc[2, "doi"] = "10.0000/single"
    sets, _ = material_sets(t)
    curves, summary, log = corpus_isosteric_heats(sets, n_grid=8, n_boot=5)
    assert len(summary) == 1 and abs(summary["q_st_low_kJ_mol"].iloc[0] - 6.0) < 0.1
    assert set(curves["material_id"]) == {summary["material_id"].iloc[0]}
    assert summary["graphitic"].iloc[0] is np.False_ or not summary["graphitic"].iloc[0]


@pytest.mark.skipif(not CORPUS.exists(), reason="harmonized corpus not present")
def test_every_unit_string_of_the_carbon_subset_converts():
    t = load_corpus(CORPUS)
    for unit in set(t["adsorption_units"]):
        assert np.isfinite(loading_to_mol_per_kg(1.0, unit))
    for unit in set(t["pressure_units"]):
        assert np.isfinite(pressure_to_bar(1.0, unit))
    sets, excluded = material_sets(t)
    assert len(sets) >= 50 and set(excluded.columns) >= {"isodb_filename", "reason"}
    assert all(isinstance(s, MaterialSet) and s.n_temperatures >= 1 for s in sets)


def test_metadata_registry_loads_and_filters(tmp_path):
    from sitewindow.experimental.metadata import load_metadata, pore_volumes, simulated_ids

    reg = tmp_path / "meta.yaml"
    reg.write_text(
        "materials:\n"
        "  - {doi: 10.1/a, adsorbent_hashkey: NIST-MATDB-aaaaaaaaaaaa, adsorbent_name: A,\n"
        "     bet_area_m2_g: 2000, pore_volume_cm3_g: 1.1, where_read: 'Table 1', access: open}\n"
        "  - {doi: 10.1/b, adsorbent_hashkey: NIST-MATDB-bbbbbbbbbbbb, adsorbent_name: B,\n"
        "     bet_area_m2_g: 1500, sample_identity: ambiguous, where_read: abstract}\n"
        "  - {doi: 10.1/c, adsorbent_hashkey: NIST-MATDB-cccccccccccc, adsorbent_name: C, simulated: true}\n"
        "  - {doi: 10.1/d, adsorbent_hashkey: NIST-MATDB-dddddddddddd, adsorbent_name: D,\n"
        "     bet_area_m2_g: 1500, value_source: secondary, where_read: 'other paper'}\n"
    )
    meta = load_metadata(reg)
    a = meta[material_id("NIST-MATDB-aaaaaaaaaaaa", "10.1/a")]
    b = meta[material_id("NIST-MATDB-bbbbbbbbbbbb", "10.1/b")]
    d = meta[material_id("NIST-MATDB-dddddddddddd", "10.1/d")]
    assert a.usable_area and not b.usable_area and not d.usable_area
    assert pore_volumes(meta) == {a.material_id: 1.1}
    assert simulated_ids(meta) == {material_id("NIST-MATDB-cccccccccccc", "10.1/c")}
    assert load_metadata(tmp_path / "absent.yaml") == {}
    # The checked-in registry parses and marks the model-calculation records as simulated.
    checked_in = load_metadata()
    assert any(m.simulated for m in checked_in.values())
    assert all(m.where_read for m in checked_in.values())


def test_b7_evaluator_recovers_a_planted_proportionality():
    from sitewindow.validation.runners import _registry, evaluate_B7

    reg = _registry()
    rng = np.random.default_rng(0)
    areas = np.array([800.0, 1500.0, 2200.0, 2800.0, 3300.0])
    uptake = 2.0 * areas / 1000.0 * (1.0 + 0.05 * rng.standard_normal(len(areas)))
    pts = pd.DataFrame({"material_id": list("abcde"), "bet_area_m2_g": areas, "uptake_wt_percent": uptake})
    r = evaluate_B7(pts, reg, n_boot=500)
    assert r["id"] == "B7" and r["passed"] is True
    s = r["predicted"]["slope_through_origin_wt_percent_per_1000m2_g"]
    assert abs(s - 2.0) < 0.2
    lo, hi = r["predicted"]["bootstrap_95_interval"]
    assert lo <= s <= hi
    assert evaluate_B7(pts.iloc[:2], reg)["passed"] is None


def test_b4_corpus_evaluator_uses_the_activated_carbon_median():
    from sitewindow.validation.runners import _registry, evaluate_B4_corpus

    reg = _registry()
    summary = pd.DataFrame(
        {
            "material_class": ["activated_carbon"] * 3 + ["SWCNT", "graphene"],
            "graphitic": [False] * 3 + [True, True],
            "q_st_low_kJ_mol": [5.0, 6.0, 12.0, 4.0, 4.5],
            "temperatures_K": [[77.0, 87.0], [77.0, 298.0], [298.0, 318.0], [77.0, 87.0], [77.0, 100.0]],
        }
    )
    r = evaluate_B4_corpus(summary, reg)
    assert r["id"] == "B4-corpus" and r["passed"] is True
    assert r["predicted"]["activated_carbon"]["median_kJ_mol"] == 6.0
    assert r["predicted"]["activated_carbon_cryogenic"]["n"] == 2
    assert r["predicted"]["graphitic"]["median_kJ_mol"] == 4.25
    summary.loc[:2, "q_st_low_kJ_mol"] = [9.0, 10.0, 12.0]
    assert evaluate_B4_corpus(summary, reg)["passed"] is False
