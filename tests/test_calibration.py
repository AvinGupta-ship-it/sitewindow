"""Published-value registry, multi-fidelity calibration and calibrated-landscape assembly (C3)."""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from sitewindow.landscape import calibration as cal
from sitewindow.landscape.assemble import assemble_calibrated, atomization_constant, calibrate_table
from sitewindow.landscape.registry import (
    cell_atoms,
    check_registry,
    published_values,
    registry_frame,
)

REGISTRY = Path(__file__).resolve().parents[1] / "data" / "registry" / "published_values.yaml"
DESC = list(cal.DESCRIPTORS)


# --------------------------------------------------------------------------- registry


def test_registry_loads_site_energies_in_ev_with_sources():
    values = published_values(REGISTRY)
    assert values, "the registry holds site energies"
    for v in values:
        assert v.doi or v.arxiv
        assert v.unit == "eV"
        assert v.extraction.value in {
            "tabular_text",
            "abstract_statement",
            "digitized_figure",
            "database_record",
        }
    frame = registry_frame(values)
    assert not frame["compilation"].any()
    assert (frame["sigma_eV"] > 0).all()
    dmc = frame[frame["level"] == "DMC"]
    assert len(dmc) >= 1 and dmc["tier0"].all()
    assert check_registry(REGISTRY)["problems"] == []


def test_registry_rejects_a_value_without_identifier(tmp_path):
    bad = [
        {
            "id": "PV-9999",
            "system": "graphene",
            "site_class": "basal",
            "quantity": "E_bind_H_atom",
            "value": -0.8,
            "unit": "eV",
            "reference_state": "H_atom",
            "method": {"functional": "PBE"},
            "extraction": "tabular_text",
            "source_id": "S4",
        }
    ]
    p = tmp_path / "pv.yaml"
    p.write_text(yaml.safe_dump(bad))
    assert any("no DOI or arXiv" in s for s in check_registry(p)["problems"])
    with pytest.raises(ValueError):
        published_values(p)


def test_cell_atoms_and_unit_conversion(tmp_path):
    assert cell_atoms("5x5") == 50 and cell_atoms("graphene 12×12 supercell") == 288
    assert math.isnan(cell_atoms(None)) and math.isnan(cell_atoms("(5,5) tube"))
    entry = [
        {
            "id": "PV-0100",
            "doi": "10.1000/test",
            "system": "coronene",
            "site_class": "basal",
            "quantity": "E_bind_H2",
            "value": -4.7,
            "unit": "kJ/mol",
            "reference_state": "H2_molecule",
            "method": {"level": "CCSD(T)/CBS"},
            "extraction": "abstract_statement",
            "source_id": "S8",
        },
        {
            "id": "PV-0101",
            "doi": "10.1000/test2",
            "system": "graphene",
            "site_class": "basal",
            "quantity": "E_bind_H_atom",
            "value": [-0.98, -0.79],
            "unit": "eV",
            "reference_state": "H_atom",
            "method": {"functional": "PBE"},
            "extraction": "tabular_text",
            "source_id": "S4",
        },
    ]
    p = tmp_path / "pv.yaml"
    p.write_text(yaml.safe_dump(entry))
    values = published_values(p)
    assert values[0].value.value == pytest.approx(-4.7 / 96.48533212)
    frame = registry_frame(values)
    assert len(frame) == 1 and frame["tier0"].iloc[0]
    assert len(registry_frame(values, include_compilations=True)) == 2


# --------------------------------------------------------------------------- transfer and GP


def _pairs(rng, n_per_class, a_b, noise=0.02):
    rows = []
    for k, (n, (a, b)) in enumerate(zip(n_per_class, a_b, strict=True)):
        e_t = rng.normal(-1.5, 0.6, n)
        x = rng.normal(0, 1, (n, len(DESC)))
        for i in range(n):
            row = dict(zip(DESC, x[i], strict=True))
            row.update(
                site_class=f"c{k}",
                configuration=f"cfg{i}",
                n_H=1,
                E_t=e_t[i],
                E_T1=a + b * e_t[i] + rng.normal(0, noise),
                tier="T2",
            )
            rows.append(row)
    return pd.DataFrame(rows)


def test_transfer_recovers_planted_slope_and_fixes_slope_with_few_anchors():
    rng = np.random.default_rng(1)
    pairs = _pairs(rng, [12], [(0.3, 0.9)], noise=0.01)
    tr = cal.fit_transfer(pairs["E_t"], pairs["E_T1"], np.full(12, 0.02), "c0", "T2")
    assert tr.a == pytest.approx(0.3, abs=0.05) and tr.b == pytest.approx(0.9, abs=0.05)
    assert not tr.b_fixed
    two = cal.fit_transfer(
        np.array([-1.0, -2.0]), np.array([-0.8, -1.8]), np.array([0.02, 0.02]), "c", "T2"
    )
    assert two.b == 1.0 and two.b_fixed and two.a == pytest.approx(0.2)
    one = cal.fit_transfer(
        np.array([-1.0]), np.array([-0.8]), np.array([0.02]), "c", "T2", pooled_sd=0.3
    )
    assert one.b_fixed and one.a == pytest.approx(0.2) and math.sqrt(one.cov[0, 0]) >= 0.3


def test_matern_kernels_and_gp_interpolate_training_points():
    assert cal.matern52(np.array([0.0]))[0] == 1.0 and cal.matern32(np.array([0.0]))[0] == 1.0
    assert cal.matern52(np.array([3.0]))[0] < 0.1
    hyper = cal.GPHyper(0.5, np.ones(2), 1e-3)
    z = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    r = np.array([0.3, -0.2, 0.1])
    mean, cov = cal.gp_predict(z, r, np.full(3, 1e-3), z, hyper)
    assert np.allclose(mean, r, atol=1e-2)
    assert np.all(np.diag(cov) < 0.01)
    far_mean, far_cov = cal.gp_predict(z, r, np.full(3, 1e-3), np.array([[30.0, 30.0]]), hyper)
    assert abs(far_mean[0]) < 1e-6 and far_cov[0, 0] == pytest.approx(0.25, rel=1e-3)


def test_calibration_statuses_and_prediction_shapes():
    rng = np.random.default_rng(2)
    pairs = _pairs(rng, [5, 2, 1], [(0.2, 1.0), (0.4, 1.0), (-0.3, 1.0)])
    c = cal.Calibration("T2").fit(pairs)
    assert c.status("c0") == "calibrated"
    assert c.status("c1") == "offset_calibrated" and c.n_anchors("c1") == 2
    assert c.status("c2") == "offset_calibrated"
    assert c.status("unknown") == "awaiting_anchors"
    frame = pd.DataFrame(
        {
            "site_class": ["c0", "c1", "c2", "unknown"],
            "value_eV": [-1.5, -1.0, -2.0, -1.2],
            **{d: [0.0, 0.5, -0.5, 0.1] for d in DESC},
        }
    )
    pred = c.predict(frame)
    assert list(pred["calibration_status"]) == [
        "calibrated",
        "offset_calibrated",
        "offset_calibrated",
        "awaiting_anchors",
    ]
    assert (pred["calibrated_sd_eV"] > 0).all() and np.isfinite(pred["calibrated_sd_eV"]).all()
    assert pred["calibrated_mean_eV"].iloc[0] == pytest.approx(-1.5 + 0.2, abs=0.1)
    rep = c.report()
    assert rep["gp"]["kernel"] == "matern52" and rep["classes"]["c0"]["n_anchors"] == 5


def test_loo_coverage_and_held_out_groups():
    rng = np.random.default_rng(3)
    pairs = _pairs(rng, [6, 6], [(0.2, 1.0), (0.1, 1.05)], noise=0.03)
    pairs.loc[pairs.site_class == "c1", "site_class"] = "epoxide"
    loo = cal.loo_coverage(pairs, "T2")
    assert loo["n"] == 12 and 0.0 <= loo["coverage"] <= 1.0 and len(loo["anchors"]) == 12
    held = cal.held_out_class_errors(
        pairs, "T2", {"oxygen": ("epoxide",), "curvature": ("cnt_armchair",)}
    )
    assert held["oxygen"]["n_held_out"] == 6 and held["curvature"]["n_held_out"] == 0
    assert held["oxygen"]["rmse_eV"] < 0.5


def test_synthetic_recovery_meets_the_acceptance_criterion():
    rec = cal.synthetic_recovery(n_replicates=40, seed=5)
    assert rec["n_replicates"] == 40
    assert rec["coverage_of_planted_class_means"] >= 0.85


def test_tier0_offset_from_the_basal_anchor():
    anchors = pd.DataFrame(
        {
            "site_class": ["basal", "basal"],
            "quantity": ["E_bind_H_atom", "E_bind_H_half_H2"],
            "reference_state": ["H_atom", "half_H2"],
            "engine": [cal.T1_ENGINE] * 2,
            "n_H": [1, 1],
            "identifier": ["basal/graphene_5x5/monomer"] * 2,
            "value_eV": [-0.82, 1.44],
        }
    )
    registry = pd.DataFrame(
        {
            "id": ["PV-0001"],
            "site_class": ["basal"],
            "site_id": ["top_monomer"],
            "quantity": ["E_bind_H_atom"],
            "reference_state": ["H_atom"],
            "value_eV": [-0.691],
            "sigma_eV": [0.019],
            "level": ["DMC"],
            "tier0": [True],
        }
    )
    off = cal.tier0_offset(anchors, registry)
    assert off["n_pairs"] == 1 and off["offset_eV"] == pytest.approx(-0.129)
    assert off["sigma_carried_eV"] == pytest.approx(0.05)
    assert off["sigma_statistical_eV"] == pytest.approx(math.hypot(0.019, 0.02))
    none = cal.tier0_offset(anchors.iloc[0:0], registry)
    assert none["n_pairs"] == 0 and none["offset_eV"] == 0.0


# --------------------------------------------------------------------------- assembly


def _landscape_frame(rows):
    base = {
        "config_hash": "h",
        "tier": "T2",
        "engine": "tblite:GFN2-xTB",
        "sigma_eV": 0.03,
        "distribution": "normal",
        "temperature_K": None,
        "converged": True,
        "n_atoms": 100,
        "settings_hash": "s",
        "provenance_kind": "computed_here",
        "source_id": "x",
        "extraction": "computed",
        "package_version": "0",
        "seed": None,
        "created": "now",
        "pyramidalization_deg": 10.0,
        "coordination": 3,
        "n_neighbor_H": 0,
        "distance_to_defect_nm": float("nan"),
        "curvature_inv_nm": 0.0,
    }
    return pd.DataFrame([{**base, **r} for r in rows])


def test_assemble_calibrated_landscape(tmp_path):
    processed = tmp_path / "full"
    processed.mkdir()
    t2 = _landscape_frame(
        [
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 1,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -0.82,
                "identifier": "basal/graphene_12x12/monomer",
                "supercell": "graphene_12x12",
            },
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 1,
                "quantity": "E_bind_H_half_H2",
                "reference_state": "half_H2",
                "value_eV": 0.87,
                "identifier": "basal/graphene_12x12/monomer",
                "supercell": "graphene_12x12",
            },
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 2,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -1.1,
                "identifier": "basal/graphene_12x12/ortho",
                "supercell": "graphene_12x12",
            },
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 2,
                "quantity": "E_sequential",
                "reference_state": "H_atom",
                "value_eV": -1.4,
                "identifier": "basal/graphene_12x12/ortho",
                "supercell": "graphene_12x12",
            },
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 2,
                "quantity": "E_sequential",
                "reference_state": "half_H2",
                "value_eV": 0.3,
                "identifier": "basal/graphene_12x12/ortho",
                "supercell": "graphene_12x12",
            },
            {
                "site_class": "monovacancy",
                "site_id": "dangling",
                "n_H": 1,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -3.29,
                "identifier": "monovacancy/monovacancy_12x12/seq1",
                "supercell": "monovacancy_12x12",
            },
            {
                "site_class": "monovacancy",
                "site_id": "dangling",
                "n_H": 0,
                "quantity": "E_formation_defect",
                "reference_state": "pristine_cell",
                "value_eV": 7.69,
                "identifier": "monovacancy/monovacancy_12x12/seq1",
                "supercell": "monovacancy_12x12",
            },
            {
                "site_class": "epoxide",
                "site_id": "adjacent",
                "n_H": 1,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -1.24,
                "identifier": "epoxide/epoxide_12x12/monomer_adjacent",
                "supercell": "epoxide_12x12",
            },
        ]
    )
    t2.to_parquet(processed / "landscape_raw_full.parquet")
    anchors = _landscape_frame(
        [
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 1,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -0.820,
                "identifier": "basal/graphene_5x5/monomer",
                "supercell": "graphene_5x5",
                "tier": "T1",
                "engine": cal.T1_ENGINE,
            },
            {
                "site_class": "basal",
                "site_id": "top",
                "n_H": 1,
                "quantity": "E_bind_H_half_H2",
                "reference_state": "half_H2",
                "value_eV": 1.440,
                "identifier": "basal/graphene_5x5/monomer",
                "supercell": "graphene_5x5",
                "tier": "T1",
                "engine": cal.T1_ENGINE,
            },
            {
                "site_class": "monovacancy",
                "site_id": "dangling",
                "n_H": 1,
                "quantity": "E_bind_H_atom",
                "reference_state": "H_atom",
                "value_eV": -3.797,
                "identifier": "monovacancy/monovacancy_6x6/seq1",
                "supercell": "monovacancy_6x6",
                "tier": "T1",
                "engine": cal.T1_ENGINE,
            },
            {
                "site_class": "monovacancy",
                "site_id": "dangling",
                "n_H": 0,
                "quantity": "E_formation_defect",
                "reference_state": "pristine_cell",
                "value_eV": 7.689,
                "identifier": "monovacancy/monovacancy_6x6/seq1",
                "supercell": "monovacancy_6x6",
                "tier": "T1",
                "engine": cal.T1_ENGINE,
            },
        ]
    )
    anchors_path = tmp_path / "anchors_tier1.parquet"
    anchors.to_parquet(anchors_path)
    result = assemble_calibrated(
        processed, "full", anchors_path, REGISTRY, meta=False, recovery_replicates=5
    )
    table = result["table"]
    status = result["status"].set_index("site_class")
    assert (
        status.loc["basal", "status"] == "offset_calibrated"
        and status.loc["basal", "n_anchors"] == 1
    )
    assert status.loc["monovacancy", "status"] == "offset_calibrated"
    assert status.loc["epoxide", "status"] == "awaiting_anchors"
    mono = table[
        (table.identifier == "basal/graphene_12x12/monomer") & (table.quantity == "E_bind_H_atom")
    ].iloc[0]
    assert abs(mono["calibrated_mean_T1_eV"] - (-0.820)) < mono["calibrated_sd_T1_eV"]
    half = table[
        (table.identifier == "basal/graphene_12x12/monomer")
        & (table.quantity == "E_bind_H_half_H2")
    ].iloc[0]
    assert half["calibrated_mean_T1_eV"] == pytest.approx(
        mono["calibrated_mean_T1_eV"] + atomization_constant(anchors), abs=1e-6
    )
    seq = table[(table.quantity == "E_sequential") & (table.reference_state == "half_H2")].iloc[0]
    assert (
        np.isfinite(seq["calibrated_mean_T1_eV"])
        and seq["calibration_status"] == "offset_calibrated"
    )
    form = table[table.quantity == "E_formation_defect"].iloc[0]
    assert form["calibrated_mean_T1_eV"] == pytest.approx(7.689) and form[
        "calibrated_mean_T0_eV"
    ] == pytest.approx(7.689)
    offset = result["report"]["tier0_offset"]
    assert offset["n_pairs"] == 1 and mono["calibrated_mean_T0_eV"] == pytest.approx(
        mono["calibrated_mean_T1_eV"] - offset["offset_eV"]
    )
    assert (processed / "landscape_calibrated.parquet").exists()
    report = json.loads((processed / "calibration_report.json").read_text())
    assert report["synthetic_recovery_calibration"]["n_replicates"] == 5
    assert {r["site_class"] for r in report["class_status"]} == {"basal", "monovacancy", "epoxide"}
    assert set(table["calibration_status"]) <= {
        "calibrated",
        "offset_calibrated",
        "awaiting_anchors",
        "not_calibrated",
    }
    out, c, pairs = calibrate_table(t2, anchors, "T2", offset)
    assert len(pairs) == 2 and len(out) == len(t2)


@pytest.mark.slow
def test_meta_analysis_fits_the_registry_with_pymc():
    pytest.importorskip("pymc")
    from sitewindow.landscape.meta_analysis import fit_meta_analysis

    frame = registry_frame(published_values(REGISTRY))
    group = frame[
        (frame.quantity == "E_bind_H_atom")
        & (frame.reference_state == "H_atom")
        & frame.site_class.notna()
    ]
    if len(group) < 3:
        pytest.skip("fewer than three published values")
    res = fit_meta_analysis(group, draws=150, tune=150, chains=2, seed=1)
    assert res.n_values == len(group)
    assert set(res.class_means["site_class"]) == set(group["site_class"])
    assert "rhat_max" in res.diagnostics and res.posterior_predictive["coverage_90"] >= 0.0
    assert any(
        e["effect"] == "gamma_cell_size" for e in res.method_effects.to_dict(orient="records")
    )
