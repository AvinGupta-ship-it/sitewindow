"""Validation runners that need no external engines, and the Markdown summary (C8)."""

from sitewindow.validation.report import render_markdown
from sitewindow.validation.runners import _registry, run_B3_pair_potential, run_B8, run_B10


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
