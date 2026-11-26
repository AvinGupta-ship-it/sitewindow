"""Validation runners that need no external engines, and the Markdown summary (C8)."""

from sitewindow.validation.report import render_markdown
from sitewindow.validation.runners import (
    _registry,
    evaluate_B1_T2,
    evaluate_B2,
    evaluate_B11,
    evaluate_convergence,
    run_B3_pair_potential,
    run_B8,
    run_B10,
)


def test_runners_return_the_report_contract():
    reg = _registry()
    for r in (run_B8(reg), run_B10(reg), run_B3_pair_potential(reg)):
        assert {"id", "passed", "predicted", "reference", "metric", "summary"} <= set(r)
        assert isinstance(r["passed"], bool)


def test_b3_pair_potential_reports_both_literature_sets():
    r = run_B3_pair_potential(_registry())
    assert set(r["checks"]) == {"steele_buch", "steele_darkrim_levesque"}
    assert r["predicted"]["calibrated_epsilon_K"] > 40.0


def test_render_markdown_has_one_row_per_benchmark():
    report = {
        "tier": "verify",
        "benchmarks": [
            {"id": "B8", "passed": True, "summary": "ok"},
            {"id": "B9", "passed": False, "summary": "a | b"},
            {"id": "B1", "passed": None, "summary": "reported"},
        ],
    }
    md = render_markdown(report, "Benchmarks B2 are not listed.")
    rows = [line for line in md.splitlines() if line.startswith("| B")]
    assert len(rows) == 3
    assert "| B9 | fail | a / b |" in md and "| B1 | reported |" in md
    assert md.rstrip().endswith("not listed.")


def _series(deltas: dict[str, float], supercell: int = 5) -> dict:
    def settings(**over):
        base = {"supercell": supercell, "vacuum_A": 15.0, "ecut_eV": 500.0, "kmesh_used": [4, 4, 1]}
        base.update(over)
        return base

    e0 = -0.78
    variants = {"base": {"settings": settings(), "E_bind_H_atom_PBE": e0}}
    changes = {
        "cell_plus_one": settings(supercell=supercell + 1, kmesh_used=[3, 3, 1]),
        "kmesh_dense": settings(kmesh_used=[6, 6, 1]),
        "cutoff_600": settings(ecut_eV=600.0),
        "vacuum_20": settings(vacuum_A=20.0),
    }
    for name, d in deltas.items():
        variants[name] = {"settings": changes[name], "E_bind_H_atom_PBE": e0 + d}
    return {"variants": variants}


def test_convergence_criterion_passes_only_when_every_variant_is_within_threshold():
    thr = _registry()["B1"]["thresholds"]["T1_convergence_abs_eV"]
    assert thr == 0.02
    ok = evaluate_convergence(
        _series(
            {
                "cell_plus_one": -0.0199,
                "kmesh_dense": -0.0136,
                "cutoff_600": -0.0023,
                "vacuum_20": -0.0011,
            }
        ),
        thr,
    )
    assert ok["passed"] is True
    assert ok["predicted"]["variant_settings"]["cell_plus_one"] == "6x6 cell, k3x3x1"
    assert ok["predicted"]["variant_settings"]["vacuum_20"] == "20 A vacuum"
    assert abs(ok["predicted"]["max_abs_delta_eV"] - 0.0199) < 1e-9
    bad = evaluate_convergence(_series({"cell_plus_one": -0.021, "cutoff_600": 0.001}), thr)
    assert bad["passed"] is False and bad["checks"] == {"cell_plus_one": False, "cutoff_600": True}
    reduced = evaluate_convergence(_series({"cutoff_600": 0.001}, supercell=4), thr)
    assert reduced["passed"] is None and "not a pass test" in reduced["summary"]


def test_b2_evaluation_reports_tier2_and_passes_only_on_tier1():
    import pandas as pd

    reg = _registry()
    cols = [
        "quantity",
        "site_class",
        "reference_state",
        "identifier",
        "tier",
        "supercell",
        "value_eV",
    ]
    t2 = pd.DataFrame(
        [
            [
                "E_barrier_stick",
                "basal",
                "H_atom",
                "basal/graphene_4x4/monomer_stick@constrained_scan",
                "T2",
                "graphene_4x4",
                0.31,
            ]
        ],
        columns=cols,
    )
    r = evaluate_B2(t2, reg)
    assert r["passed"] is None and "T2:constrained_scan:graphene_4x4" in r["predicted"]
    t1 = pd.DataFrame(
        [
            [
                "E_barrier_stick",
                "basal",
                "H_atom",
                "basal/graphene_5x5/monomer_stick@constrained_scan",
                "T1",
                "graphene_5x5",
                0.23,
            ]
        ],
        columns=cols,
    )
    assert evaluate_B2(pd.concat([t2, t1]), reg)["passed"] is True
    t1.loc[0, "value_eV"] = 0.40
    assert evaluate_B2(pd.concat([t2, t1]), reg)["passed"] is False
    assert evaluate_B2(t2.iloc[:0], reg)["passed"] is None
    flagged = pd.concat([t2, t1]).assign(converged=[True, False])
    r = evaluate_B2(flagged, reg)
    assert r["passed"] is None and "no saddle located" in r["summary"]
    assert r["saddle_located"] == {
        "T2:constrained_scan:graphene_4x4": True,
        "T1:constrained_scan:graphene_5x5": False,
    }
    assert set(r["checks"]) == {"T2:constrained_scan:graphene_4x4"}


def test_b1_tier2_evaluation_reports_cell_and_smearing():
    import pandas as pd

    reg = _registry()
    frame = pd.DataFrame(
        [
            ["basal", "E_bind_H_atom", 1, "graphene_9x9", -0.818],
            ["basal", "E_bind_H_atom", 2, "graphene_9x9", -1.1],
        ],
        columns=["site_class", "quantity", "n_H", "supercell", "value_eV"],
    )
    r = evaluate_B1_T2(frame, reg, {"electronic_temperature": 3000.0})
    assert r["passed"] is None and r["predicted"]["E_bind_H_atom_eV"] == -0.818
    assert "9x9, 3000 K smearing" in r["summary"] and "-0.13" in r["summary"]
    assert evaluate_B1_T2(frame.iloc[:0], reg)["passed"] is None
    t3 = evaluate_B1_T2(frame, reg, {"model": "medium"}, tier_label="T3")
    assert t3["summary"].startswith("T3 basal monomer") and "smearing" not in t3["summary"]


def test_b11_evaluation_applies_tier_thresholds_and_reports_tier3():
    import pandas as pd

    reg = _registry()

    def frame(rows):
        return pd.DataFrame(
            [
                {
                    "quantity": "E_formation_defect",
                    "site_class": cls,
                    "supercell": cell,
                    "n_atoms": 100,
                    "engine": eng,
                    "value_eV": v,
                }
                for cls, cell, eng, v in rows
            ]
        )

    t2 = frame(
        [
            ("monovacancy", "monovacancy_12x12", "tblite:GFN2-xTB", 7.59),
            ("stone_wales", "stone_wales_12x12", "tblite:GFN2-xTB", 4.67),
        ]
    )
    t3 = frame([("monovacancy", "monovacancy_9x9", "mace", 5.69)])
    r = evaluate_B11({"T2": t2, "T3": t3}, reg)
    assert r["passed"] is True
    assert r["checks"] == {
        "T2:monovacancy:monovacancy_12x12": True,
        "T2:stone_wales:stone_wales_12x12": True,
    }
    assert "T3:monovacancy:monovacancy_9x9" in r["predicted"] and "no T1 values yet" in r["summary"]
    t1 = frame(
        [
            ("monovacancy", "monovacancy_6x6", "gpaw:PBE-D3(BJ)", 8.2),
            ("monovacancy", "monovacancy_6x6", "gpaw:PBE", 8.1),
            ("stone_wales", "stone_wales_6x6", "gpaw:PBE-D3(BJ)", 5.1),
        ]
    )
    r1 = evaluate_B11({"T1": t1, "T2": t2}, reg)
    assert r1["passed"] is False
    assert r1["checks"]["T1:monovacancy:monovacancy_6x6"] is False
    assert r1["checks"]["T1:stone_wales:stone_wales_6x6"] is True
    assert r1["predicted"]["T1:monovacancy:monovacancy_6x6"] == 8.2
    assert evaluate_B11({}, reg)["passed"] is None
