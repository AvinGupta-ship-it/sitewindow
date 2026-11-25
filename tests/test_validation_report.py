"""Validation runners that need no external engines, and the Markdown summary (C8)."""

from sitewindow.validation.report import render_markdown
from sitewindow.validation.runners import (
    _registry,
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
