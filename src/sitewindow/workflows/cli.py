"""Command-line interface for SiteWindow workflows."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import typer

from sitewindow import __version__

app = typer.Typer(
    add_completion=False, help="SiteWindow: site-energy landscapes and deliverable-capacity limits."
)


def _data_root() -> Path:
    import os

    return Path(os.environ.get("SITEWINDOW_DATA", "data"))


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", help="Print version and exit."),
) -> None:
    if version:
        typer.echo(f"sitewindow {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@app.command("fetch-isodb")
def fetch_isodb(
    dest: Path = typer.Option(
        None, help="Clone directory (default: $SITEWINDOW_DATA/raw/isodb-library)"
    ),
) -> None:
    """Clone or update the NIST ISODB GitHub mirror and harmonize the carbon H2 isotherms."""
    from sitewindow.experimental.isodb import clone_mirror, load_carbon_h2, write_parquet

    dest = dest or _data_root() / "raw" / "isodb-library"
    commit = clone_mirror(dest)
    res = load_carbon_h2(dest)
    out = write_parquet(res, _data_root() / "processed", commit)
    (_data_root() / "sources.lock").write_text(json.dumps({"isodb_commit": commit}, indent=2))
    typer.echo(json.dumps({"commit": commit, "output": str(out), **res.notes}, indent=2))


@app.command("theory-tables")
def theory_tables(
    out: Path = typer.Option(None, help="Output directory (default: $SITEWINDOW_DATA/processed)"),
) -> None:
    """Window efficiencies, optimum and threshold energies, and bounds for the standard windows (E3)."""
    import pandas as pd

    from sitewindow.constants import EV_TO_KJ_PER_MOL
    from sitewindow.theory.deliverable import bound_wt_percent
    from sitewindow.theory.entropy import DEFAULT_MODEL_SET
    from sitewindow.theory.window import (
        BHATIA_MYERS,
        W1,
        W2,
        W3,
        heterogeneity_threshold,
        optimum_energy,
    )

    rows = []
    for w in (BHATIA_MYERS, W1, W2, W3):
        for m in DEFAULT_MODEL_SET:
            eps_star, eta = optimum_energy(w, m)
            b = bound_wt_percent(w, "molecular", m)
            eps_c = heterogeneity_threshold(-5.8 / EV_TO_KJ_PER_MOL, w, m)
            rows.append(
                {
                    "window": w.name,
                    "entropy_model": m.name,
                    "eps_star_kJ_mol": eps_star * EV_TO_KJ_PER_MOL,
                    "eta": eta,
                    "D_max_molecular_wt_percent": b["D_max_wt_percent"],
                    "eps_c_for_-5.8kJ_site_kJ_mol": None
                    if eps_c is None
                    else eps_c * EV_TO_KJ_PER_MOL,
                }
            )
    df = pd.DataFrame(rows)
    out = out or _data_root() / "processed"
    out.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out / "theory_window_tables.parquet", index=False)
    df.to_csv(out / "theory_window_tables.csv", index=False)
    typer.echo(df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))


@app.command("knockon")
def knockon(
    energies_keV: str = typer.Option("80,85,90,95,100,105,110,120"),
    E_d_eV: float = 22.0,
    theta_D_K: float = 1287.0,
    T_K: float = 300.0,
) -> None:
    """Knock-on cross sections (static and Debye-smeared) for graphene (C5, benchmark B9)."""
    from sitewindow.defects.knockon import (
        debye_mean_square_velocity,
        mckinley_feshbach_barn,
        thermally_smeared_cross_section_barn,
    )

    v2 = debye_mean_square_velocity(theta_D_K, T_K)
    for E in [float(x) for x in energies_keV.split(",")]:
        typer.echo(
            f"{E:7.1f} keV  static {float(mckinley_feshbach_barn(E, E_d_eV)):.4e} barn  "
            f"smeared {thermally_smeared_cross_section_barn(E, E_d_eV, v2):.4e} barn"
        )


@app.command("validate")
def validate(
    tier: str = typer.Option("verify", help="verify or full"), out: Path = typer.Option(None)
) -> None:
    """Run the benchmark runners available at the requested tier and write validation_report.json."""
    from sitewindow.validation.runners import run_all

    report = run_all(tier=tier)
    out = out or _data_root() / "processed" / tier
    out.mkdir(parents=True, exist_ok=True)
    (out / "validation_report.json").write_text(json.dumps(report, indent=2, default=float))
    for r in report["benchmarks"]:
        status = "PASS" if r["passed"] else ("FAIL" if r["passed"] is False else "NOT RUN")
        typer.echo(f"{r['id']:4s} {status:8s} {r['summary']}")


@app.command("verify-atomistic")
def verify_atomistic(cache: Path = typer.Option(None), supercell: int = 4) -> None:
    """Verification-tier GFN2-xTB calculations: basal monomer and vacancy dangling-bond H (C2, tier T2)."""
    from sitewindow.atomistic.calculators.tblite_adapter import TBLiteCalculator
    from sitewindow.atomistic.structures.builders import (
        add_hydrogen_top,
        central_atom_index,
        graphene_supercell,
        monovacancy,
    )
    from sitewindow.atomistic.workflows.adsorption import (
        binding_energies,
        reference_energies,
        relax,
    )

    cache = cache or _data_root() / "raw" / "cache_verify_tblite"
    calc = TBLiteCalculator()
    refs = reference_energies(calc, cache)
    g, rg = relax(graphene_supercell(supercell), calc, cache)
    ci = central_atom_index(g)
    gh, rgh = relax(add_hydrogen_top(g, ci), calc, cache)
    be = binding_energies(rg.energy_eV, rgh.energy_eV, refs)
    v, _ = monovacancy(supercell)
    v, rv = relax(v, calc, cache)
    rp = v.info["removed_position"]
    db = int(np.argmin(np.linalg.norm(v.positions - rp, axis=1)))
    vh = add_hydrogen_top(v, db)
    vh.positions[-1] = v.positions[db] + 1.1 * (rp - v.positions[db]) / np.linalg.norm(
        rp - v.positions[db]
    )
    vh, rvh = relax(vh, calc, cache)
    bev = binding_energies(rv.energy_eV, rvh.energy_eV, refs)
    result = {
        "tier": "T2",
        "engine": "GFN2-xTB (tblite)",
        "supercell": supercell,
        "basal_monomer": be,
        "monovacancy_formation_eV_approx": rv.energy_eV - rg.energy_eV * (len(v) / len(g)),
        "vacancy_dangling_bond_H": bev,
        "settings_hash": rg.settings_hash,
        "converged": all(r.converged for r in (rg, rgh, rv, rvh)),
    }
    out = _data_root() / "processed" / "verify"
    out.mkdir(parents=True, exist_ok=True)
    (out / "atomistic_verify_T2.json").write_text(json.dumps(result, indent=2))
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":  # pragma: no cover
    app()
