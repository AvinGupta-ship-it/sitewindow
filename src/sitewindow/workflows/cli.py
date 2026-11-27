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
landscape_app = typer.Typer(help="Site-energy landscape construction (C2).")
tier1_app = typer.Typer(help="Tier-1 density-functional anchors and convergence (C2, B1).")
tier2_app = typer.Typer(help="Tier-2 (GFN2-xTB) protocol studies (C2).")
schema_app = typer.Typer(help="Canonical data schema.")
app.add_typer(landscape_app, name="landscape")
app.add_typer(tier1_app, name="tier1")
app.add_typer(tier2_app, name="tier2")
app.add_typer(schema_app, name="schema")


def _data_root() -> Path:
    import os

    return Path(os.environ.get("SITEWINDOW_DATA", "data"))


def _tier_config(tier: str, config: Path | None) -> dict:
    import yaml

    p = config or Path("configs") / "tiers" / f"{tier}.yaml"
    return (yaml.safe_load(p.read_text()) or {}) if p.exists() else {}


def _calculator_settings(cfg: dict, calculator: str) -> dict:
    """Calculator keyword settings from the tier configuration: ``calculator_settings`` applies to
    every engine, ``calculator_settings_<name>`` to one engine."""
    out = dict(cfg.get("calculator_settings", {}) or {})
    out.update(cfg.get(f"calculator_settings_{calculator}", {}) or {})
    return out


def _is_master() -> bool:
    from sitewindow.atomistic.parallel import is_master

    return is_master()


def _configure_logging(verbose: bool) -> None:
    """Timing and settings logging on rank 0 only (manual Section 10.5)."""
    import logging

    logging.basicConfig(
        level=logging.INFO if (verbose and _is_master()) else logging.WARNING,
        format="%(asctime)s %(name)s %(message)s",
    )


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
    update: bool = typer.Option(
        False,
        "--update",
        help="Use the mirror's current HEAD instead of the commit pinned in sources.lock.",
    ),
) -> None:
    """Clone the NIST ISODB GitHub mirror at the pinned commit and harmonize the carbon H2 isotherms."""
    from sitewindow.experimental.isodb import (
        ISODB_GIT_URL,
        clone_mirror,
        load_carbon_h2,
        write_parquet,
    )

    dest = dest or _data_root() / "raw" / "isodb-library"
    lock_path = _data_root() / "sources.lock"
    pinned = None
    if lock_path.exists() and not update:
        pinned = json.loads(lock_path.read_text()).get("isodb_commit")
    commit = clone_mirror(dest, commit=pinned)
    res = load_carbon_h2(dest)
    out = write_parquet(res, _data_root() / "processed", commit)
    lock_path.write_text(
        json.dumps({"isodb_commit": commit, "isodb_url": ISODB_GIT_URL}, indent=2) + "\n"
    )
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
    tier: str = typer.Option("verify", help="verify or full"),
    out: Path = typer.Option(None),
    markdown: Path = typer.Option(None, help="Also write the human-readable summary (Markdown)"),
) -> None:
    """Run the benchmark runners available at the requested tier and write validation_report.json."""
    from sitewindow.validation.report import render_markdown
    from sitewindow.validation.runners import run_all

    report = run_all(tier=tier)
    out = out or _data_root() / "processed" / tier
    out.mkdir(parents=True, exist_ok=True)
    (out / "validation_report.json").write_text(json.dumps(report, indent=2, default=float))
    if markdown is not None:
        ids = {r["id"].split("-")[0] for r in report["benchmarks"]}
        missing = [f"B{k}" for k in range(1, 13) if f"B{k}" not in ids]
        note = (
            "Benchmarks "
            + ", ".join(missing)
            + " require components not yet executed and are not listed."
            if missing
            else None
        )
        markdown.write_text(render_markdown(report, note))
    for r in report["benchmarks"]:
        status = "PASS" if r["passed"] else ("FAIL" if r["passed"] is False else "REPORTED")
        typer.echo(f"{r['id']:5s} {status:8s} {r['summary']}")


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


@landscape_app.command("build")
def landscape_build(
    tier: str = typer.Option("verify", help="verify or full (configs/tiers/<tier>.yaml)"),
    calculator: str = typer.Option("tblite", help="tblite, gpaw, mace or a registered plugin"),
    classes: str = typer.Option(None, help="Comma-separated site classes (default: tier config)"),
    config: Path = typer.Option(None, help="Tier configuration file"),
    cache: Path = typer.Option(
        None, help="Cache directory (default: $SITEWINDOW_DATA/raw/cache_<tier>_<calculator>)"
    ),
    out: Path = typer.Option(
        None, help="Output directory (default: $SITEWINDOW_DATA/processed/<tier>)"
    ),
    fmax: float = typer.Option(None, help="Force tolerance (eV/A); default from the tier config"),
    spin_search: bool = typer.Option(False, help="Also relax in the next higher spin multiplicity"),
    energy: str = typer.Option(
        None,
        help="F (reported energy) or E0 (zero-temperature extrapolation, tblite only); default "
        "from the tier config key energy_<calculator>, else F",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Relax every catalogue configuration of the tier and write landscape_raw_<tier>.parquet (E1)."""
    from sitewindow.atomistic.calculators import make_calculator
    from sitewindow.atomistic.workflows.landscape import LandscapeSpec, build_landscape
    from sitewindow.landscape import catalogue

    _configure_logging(verbose)
    cfg = _tier_config(tier, config).get("atomistic", {})
    cls_list = [c.strip() for c in classes.split(",")] if classes else list(cfg.get("classes", []))
    cls_list = cls_list or None
    recipes: dict[str, dict] = {}
    if "supercell" in cfg:
        for name, d in catalogue().items():
            if "supercell" in d.recipe(tier if tier in d.recipes else "full"):
                recipes.setdefault(name, {})["supercell"] = int(cfg["supercell"])
    cnt_cfg = dict(cfg.get("cnt", {}))
    if "cnt_diameters_nm" in cfg:
        cnt_cfg.setdefault("diameters_nm", list(cfg["cnt_diameters_nm"]))
    if cnt_cfg:
        for name in ("cnt_armchair", "cnt_zigzag"):
            recipes.setdefault(name, {}).update(
                {k: v for k, v in cnt_cfg.items() if k in ("diameters_nm", "length")}
            )
        lengths = {
            "cnt_armchair": cnt_cfg.get("length_armchair"),
            "cnt_zigzag": cnt_cfg.get("length_zigzag"),
            "mixed_vacancy_cnt": cnt_cfg.get("length_armchair", cnt_cfg.get("length")),
            "mixed_sw_cnt": cnt_cfg.get("length_armchair", cnt_cfg.get("length")),
        }
        for name, length in lengths.items():
            if length is not None:
                recipes.setdefault(name, {})["length"] = int(length)
    spec = LandscapeSpec(
        tier=tier,
        classes=cls_list,
        recipes=recipes,
        fmax=float(fmax if fmax is not None else cfg.get("fmax_eV_A", 0.05)),
        spin_states=(0, 2) if spin_search else (0,),
        energy=str(energy or cfg.get(f"energy_{calculator}", "F")),
    )
    calc = make_calculator(calculator, **_calculator_settings(cfg, calculator))
    cache = cache or _data_root() / "raw" / f"cache_{tier}_{calculator}"
    out = out or _data_root() / "processed" / tier
    result = build_landscape(spec, calc, cache, out)
    m = result.get("manifest", {})
    typer.echo(
        json.dumps(
            {
                k: m.get(k)
                for k in ("tier", "energy", "n_records", "n_configurations", "n_excluded", "wall_s")
            },
            indent=2,
        )
    )


@landscape_app.command("barriers")
def landscape_barriers(
    tier: str = typer.Option("verify", help="verify or full (configs/tiers/<tier>.yaml)"),
    calculator: str = typer.Option("tblite", help="tblite, gpaw, mace or a registered plugin"),
    pathways: str = typer.Option(
        None, help="Comma-separated pathway names or keys (default: all designated pathways)"
    ),
    method: str = typer.Option(
        "constrained_scan", help="Comma-separated: constrained_scan, ci_neb"
    ),
    config: Path = typer.Option(None, help="Tier configuration file"),
    cache: Path = typer.Option(None, help="Cache directory (default: landscape cache of the tier)"),
    out: Path = typer.Option(
        None, help="Output directory (default: $SITEWINDOW_DATA/processed/<tier>)"
    ),
    fmax: float = typer.Option(None, help="Force tolerance (eV/A); default from the tier config"),
    neb_images: int = typer.Option(None, help="Number of NEB images (default: tier config, 7)"),
    seed_profiles: Path = typer.Option(
        None,
        help="Profile file of another run whose transition structures seed the bands of "
        "--seed-pathways (default: every pathway it holds)",
    ),
    seed_pathways: str = typer.Option(None, help="Comma-separated pathway names or keys to seed"),
    seed_label: str = typer.Option("seeded", help="Label of the seeded band trajectories"),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Sticking and recombination barriers along the designated pathways (C2; benchmark B2)."""
    from sitewindow.atomistic.calculators import make_calculator
    from sitewindow.atomistic.workflows.barriers import BarrierSpec, build_barriers
    from sitewindow.landscape import catalogue

    _configure_logging(verbose)
    cfg = _tier_config(tier, config).get("atomistic", {})
    recipes: dict[str, dict] = {}
    if "supercell" in cfg:
        for name, d in catalogue().items():
            if "supercell" in d.recipe(tier if tier in d.recipes else "full"):
                recipes.setdefault(name, {})["supercell"] = int(cfg["supercell"])
    cnt_cfg = dict(cfg.get("cnt", {}))
    if "cnt_diameters_nm" in cfg:
        cnt_cfg.setdefault("diameters_nm", list(cfg["cnt_diameters_nm"]))
    if cnt_cfg:
        recipes.setdefault("cnt_armchair", {}).update(
            {k: v for k, v in cnt_cfg.items() if k in ("diameters_nm", "length")}
        )
        if cnt_cfg.get("length_armchair") is not None:
            recipes["cnt_armchair"]["length"] = int(cnt_cfg["length_armchair"])
    stick_cfg = dict((cfg.get("barriers", {}) or {}).get("stick", {}) or {})
    stick_recipes: dict[str, dict] = {}
    if "supercell" in stick_cfg:
        stick_recipes["basal"] = {"supercell": int(stick_cfg["supercell"])}
    spec = BarrierSpec(
        tier=tier,
        pathways=[x.strip() for x in pathways.split(",")] if pathways else None,
        recipes=recipes,
        methods=tuple(m.strip() for m in method.split(",")),
        fmax=float(fmax if fmax is not None else cfg.get("fmax_eV_A", 0.05)),
        stick_recipes=stick_recipes,
        stick_calculator_settings=_calculator_settings(stick_cfg, calculator),
        neb_images=int(
            neb_images
            if neb_images is not None
            else (cfg.get("tier3", {}) or {}).get("neb_images", 7)
        ),
        seed_profiles=seed_profiles,
        seed_pathways=(
            tuple(x.strip() for x in seed_pathways.split(",")) if seed_pathways else None
        ),
        seed_label=seed_label,
    )
    barrier_cfg = {**cfg, **{k: v for k, v in (cfg.get("barriers") or {}).items() if k != "stick"}}
    calc = make_calculator(calculator, **_calculator_settings(barrier_cfg, calculator))
    cache = cache or _data_root() / "raw" / f"cache_{tier}_{calculator}"
    out = out or _data_root() / "processed" / tier
    result = build_barriers(spec, calc, cache, out)
    if not _is_master():
        return
    rows = []
    for prof in result["profiles"]:
        for m, r in prof["methods"].items():
            rows.append(
                {
                    "pathway": prof["pathway"],
                    "method": m,
                    "barrier_eV": None if r["barrier_eV"] is None else round(r["barrier_eV"], 3),
                    "reference": prof["reference_state"],
                    "saddle_located": r.get("saddle_located", False),
                }
            )
    if rows:
        import pandas as pd

        typer.echo(pd.DataFrame(rows).to_string(index=False))
    m = result["manifest"]
    typer.echo(
        json.dumps({k: m.get(k) for k in ("tier", "n_records", "n_excluded", "wall_s")}, indent=2)
    )


@landscape_app.command("physisorption")
def landscape_physisorption(
    tier: str = typer.Option("verify", help="verify or full"),
    out: Path = typer.Option(
        None, help="Output directory (default: $SITEWINDOW_DATA/processed/<tier>)"
    ),
    supercell: int = typer.Option(4),
) -> None:
    """Pair-potential H2 site energies for the physisorption set with Feynman-Hibbs corrections (E5)."""
    from sitewindow.atomistic.workflows.physisorption import PhysisorptionSpec, build_physisorption

    out = out or _data_root() / "processed" / tier
    res = build_physisorption(PhysisorptionSpec(tier=tier, supercell=supercell), out)
    t = res["table"]
    classical = t[t["T_K"].isna()].pivot(
        index=["geometry", "site_id"], columns="potential", values="E_kJ_mol"
    )
    typer.echo(classical.round(2).to_string())
    typer.echo(
        json.dumps({k: res["manifest"][k] for k in ("tier", "n_records", "wall_s")}, indent=2)
    )


@tier1_app.command("anchor")
def tier1_anchor(
    supercell: int = typer.Option(5),
    ecut: float = typer.Option(500.0, help="Plane-wave cutoff (eV)"),
    kmesh: str = typer.Option(
        None, help="Explicit k-mesh 'k1,k2,k3' (default: scaled from 4x4x1 on 5x5)"
    ),
    vacuum: float = typer.Option(15.0),
    fmax: float = typer.Option(0.02),
    no_prerelax: bool = typer.Option(False, help="Skip the Tier-2 pre-relaxation"),
    cache: Path = typer.Option(None),
    out: Path = typer.Option(None),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Tier-1 basal-monomer anchor with GPAW (PBE-D3(BJ)); compares the PBE value with B1."""
    from sitewindow.atomistic.workflows.tier1 import Tier1Settings, anchor_basal_monomer

    _configure_logging(verbose)
    km = tuple(int(x) for x in kmesh.split(",")) if kmesh else None
    settings = Tier1Settings(
        supercell=supercell,
        vacuum_A=vacuum,
        ecut_eV=ecut,
        kmesh=km,
        fmax_eV_A=fmax,
        prerelax_tier2=not no_prerelax,
    )
    cache = cache or _data_root() / "raw" / "cache_tier1_gpaw"
    out = out or _data_root() / "processed" / "tier1"
    res = anchor_basal_monomer(settings, cache, out, log_dir=cache / "logs")
    if not _is_master():
        return
    typer.echo(
        json.dumps(
            {
                "binding_eV": res["binding_eV"],
                "benchmark_B1": res["benchmark_B1"],
                "wall_s": res["wall_s"],
            },
            indent=2,
        )
    )


@tier1_app.command("anchors")
def tier1_anchors(
    registry: Path = typer.Option(
        None, help="Anchor registry (default: $SITEWINDOW_DATA/registry/tier1_anchors.yaml)"
    ),
    rounds: str = typer.Option(None, help="Comma-separated rounds to run (default: all)"),
    keys: str = typer.Option(None, help="Comma-separated anchor keys to run (default: all)"),
    ecut: float = typer.Option(500.0, help="Plane-wave cutoff (eV)"),
    vacuum: float = typer.Option(15.0),
    fmax: float = typer.Option(0.02),
    relax_steps: int = typer.Option(200),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print the ordered plan and exit"),
    cache: Path = typer.Option(None),
    out: Path = typer.Option(None),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Tier-1 anchor set (PBE-D3(BJ) at the converged protocol) over the registry, checkpointed
    per structure and per anchor; rebuilds anchors_tier1.parquet after every anchor (C2, C3)."""
    from sitewindow.atomistic.workflows.tier1 import Tier1Settings, run_anchor_set

    _configure_logging(verbose)
    settings = Tier1Settings(vacuum_A=vacuum, ecut_eV=ecut, fmax_eV_A=fmax, relax_steps=relax_steps)
    registry = registry or _data_root() / "registry" / "tier1_anchors.yaml"
    cache = cache or _data_root() / "raw" / "cache_tier1_gpaw"
    out = out or _data_root() / "processed" / "tier1"
    res = run_anchor_set(
        settings,
        cache,
        out,
        registry_path=registry,
        rounds=[int(r) for r in rounds.split(",")] if rounds else None,
        keys=[k.strip() for k in keys.split(",")] if keys else None,
        log_dir=cache / "logs",
        dry_run=dry_run,
    )
    if not _is_master():
        return
    import pandas as pd

    plan = pd.DataFrame(res["plan"])
    typer.echo(plan.drop(columns=["record"]).to_string(index=False))
    if "manifest" in res:
        m = res["manifest"]
        typer.echo(
            json.dumps(
                {k: m.get(k) for k in ("n_anchors", "n_records", "classes", "wall_s_total")},
                indent=2,
            )
        )


@tier1_app.command("converge")
def tier1_converge(
    supercell: int = typer.Option(5),
    ecut: float = typer.Option(500.0),
    vacuum: float = typer.Option(15.0),
    variants: str = typer.Option(
        None, help="Comma-separated subset of base,cell_plus_one,kmesh_dense,cutoff_600,vacuum_20"
    ),
    cache: Path = typer.Option(None),
    out: Path = typer.Option(None),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Convergence of the Tier-1 basal-monomer binding energy (acceptance criterion B1)."""
    from sitewindow.atomistic.workflows.tier1 import Tier1Settings, convergence_series

    _configure_logging(verbose)
    base = Tier1Settings(supercell=supercell, vacuum_A=vacuum, ecut_eV=ecut)
    cache = cache or _data_root() / "raw" / "cache_tier1_gpaw"
    out = out or _data_root() / "processed" / "tier1"
    res = convergence_series(
        base, cache, out, variants.split(",") if variants else None, log_dir=cache / "logs"
    )
    if not _is_master():
        return
    typer.echo(
        json.dumps(
            {
                k: {kk: v[kk] for kk in ("E_bind_H_atom_PBE", "delta_vs_base_PBE_eV")}
                for k, v in res["variants"].items()
            },
            indent=2,
        )
    )


@tier2_app.command("smearing")
def tier2_smearing(
    cells: str = typer.Option("9,10,12", help="Comma-separated graphene cell sizes"),
    temperatures: str = typer.Option("1500,3000,6000", help="Electronic temperatures (K)"),
    classes: str = typer.Option("basal,monovacancy,epoxide"),
    cache: Path = typer.Option(None),
    out: Path = typer.Option(
        None, help="Output directory (default: $SITEWINDOW_DATA/processed/tier2)"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Free energy, internal energy and zero-temperature extrapolation of Tier-2 binding energies
    under Fermi smearing: far hydrogen, basal monomer, monovacancy and epoxide (C2, Tier 2)."""
    from sitewindow.atomistic.workflows.smearing import SmearingSpec, build_smearing_study

    _configure_logging(verbose)
    spec = SmearingSpec(
        cells=tuple(int(x) for x in cells.split(",")),
        temperatures_K=tuple(float(x) for x in temperatures.split(",")),
        classes=tuple(c.strip() for c in classes.split(",")),
    )
    cache = cache or _data_root() / "raw" / "cache_tier2_smearing"
    out = out or _data_root() / "processed" / "tier2"
    res = build_smearing_study(spec, cache, out)
    if not _is_master():
        return
    import pandas as pd

    far = pd.DataFrame(res["far_hydrogen"])
    typer.echo(
        far[["cell", "T_K", "offset_F_eV", "offset_E_internal_eV", "offset_E0_eV"]]
        .round(3)
        .to_string(index=False)
    )
    b = pd.DataFrame(res["binding"])
    cols = [
        "site_class",
        "cell",
        "T_K",
        "E_bind_H_atom_F",
        "E_bind_H_atom_E_internal",
        "E_bind_H_atom_E0",
        "delta_TS_eV",
    ]
    typer.echo(b[cols].round(3).to_string(index=False))


@schema_app.command("export")
def schema_export(
    out: Path = typer.Option(Path("docs/schema"), help="Directory for JSON Schema files"),
) -> None:
    """Write the JSON Schema of every canonical model."""
    from sitewindow import schema as sch

    out.mkdir(parents=True, exist_ok=True)
    names = [
        "Provenance",
        "UncertainScalar",
        "SiteClass",
        "SiteConfiguration",
        "SiteEnergyRecord",
        "PublishedValue",
        "Isotherm",
        "Material",
        "DefectPopulation",
        "OperatingWindow",
        "CapacityResult",
        "BenchmarkRecord",
    ]
    for name in names:
        (out / f"{name}.json").write_text(
            json.dumps(sch.json_schema(getattr(sch, name)), indent=2, sort_keys=True) + "\n"
        )
    typer.echo(f"wrote {len(names)} schema files to {out}")


if __name__ == "__main__":  # pragma: no cover
    app()
