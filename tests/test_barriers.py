"""Barrier workflow (C2, B2): coordinate scan, climbing-image NEB, pathway catalogue, tables."""

import json
from dataclasses import dataclass, field

import numpy as np
import pytest
from ase import Atoms
from ase.build import add_adsorbate, fcc100
from ase.constraints import FixAtoms, FixBondLength

from sitewindow.atomistic.structures.builders import hydrogen_molecule
from sitewindow.atomistic.workflows.adsorption import relax
from sitewindow.atomistic.workflows.barriers import (
    CI_NEB,
    CONSTRAINED_SCAN,
    BarrierSpec,
    ScanPoint,
    build_barriers,
    h2_above,
    largest_step,
    mark_islands,
    neb_barrier,
    pathways,
    primary_method,
    product_distance,
    resample_path,
    saddle_caveat,
    scan_profile,
    set_coordinate,
)
from sitewindow.io import read_table
from sitewindow.schema import frame_to_records


@dataclass
class EMTCalculator:
    """Effective-medium test calculator (crude C/H parameters); tier label T3 for the pipeline."""

    name: str = "emt"
    tier: str = "T3"
    extra: dict = field(default_factory=dict)

    def attach(self, atoms: Atoms) -> None:
        from ase.calculators.emt import EMT

        atoms.calc = EMT()

    def settings(self) -> dict:
        return {"model": "emt", **self.extra}


def test_scan_profile_refines_the_maximum_to_the_fine_grid():
    calls: list[float] = []

    def evaluate(d, prev):
        calls.append(d)
        return ScanPoint(d, float(np.exp(-((d - 1.57) ** 2) / 0.02)), True, None)

    pts = scan_profile(evaluate, [1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8], 0.02, 0.1)
    top = max(pts, key=lambda p: p.energy_eV)
    assert abs(top.coordinate_A - 1.57) <= 0.011
    assert [p.coordinate_A for p in pts] == sorted(p.coordinate_A for p in pts)
    assert len(set(calls)) == len(calls)
    assert calls[:7] == [1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]


def test_mark_islands_excludes_upper_scf_solutions_but_not_valley_jumps():
    # The 10x10 sticking scan pattern: a physical barrier, then isolated points that converged
    # to a solution 1.6 eV above the ground state of their neighbours.
    energies = {1.5: 0.12, 1.6: 0.23, 1.7: 0.25, 1.8: 0.21, 2.9: -0.04, 3.0: 1.53, 3.1: -0.04}
    energies.update({3.2: 1.54, 3.3: 1.54, 3.32: 1.66, 3.34: 1.66, 3.36: 1.55, 3.44: -0.03})
    energies.update({3.5: -0.03, 3.6: -0.02})
    pts = [ScanPoint(d, e, True, None) for d, e in energies.items()]
    marked = mark_islands(pts)
    excluded = {p.coordinate_A for p in marked if p.excluded}
    assert excluded == {3.0, 3.2, 3.3, 3.32, 3.34, 3.36}
    assert all(p.excluded is None for p in marked if p.coordinate_A < 2.95)
    assert [p.coordinate_A for p in marked] == sorted(energies)
    # A jump between two valleys keeps both branches: the lower branch holds the end point.
    valley = [ScanPoint(d, 0.8 * (d - 1.1), True, None) for d in (1.1, 1.5, 1.9, 2.1, 2.16)]
    valley += [ScanPoint(d, -1.1, True, None) for d in (2.18, 2.5, 3.0)]
    assert all(p.excluded is None for p in mark_islands(valley))
    # A run above its left neighbour but below its right one is not an island.
    middle = [ScanPoint(1.0, 0.0, True, None), ScanPoint(1.1, 1.0, True, None)]
    middle += [ScanPoint(1.12, 2.0, True, None), ScanPoint(1.14, 3.0, True, None)]
    assert all(p.excluded is None for p in mark_islands(middle))
    assert mark_islands([]) == []


def test_scan_profile_refines_around_the_ground_state_maximum_not_an_island():
    def evaluate(d, prev):
        spike = 1.6 if 2.95 < d < 3.45 and round(d * 100) % 20 == 0 else 0.0
        return ScanPoint(d, spike + float(np.exp(-((d - 1.73) ** 2) / 0.02)), True, None)

    grid = [1.2 + 0.1 * k for k in range(25)]
    pts = scan_profile(evaluate, grid, 0.02, 0.1)
    accepted = [p for p in pts if p.accepted]
    top = max(accepted, key=lambda p: p.energy_eV)
    assert abs(top.coordinate_A - 1.73) <= 0.011
    assert any(p.excluded for p in pts)
    assert all(abs(p.coordinate_A - 1.73) < 0.12 for p in pts if p.coordinate_A * 10 % 1 > 1e-6)


def test_set_coordinate_moves_only_the_hydrogen_for_sticking_and_both_for_recombination():
    h2 = hydrogen_molecule()
    stick = set_coordinate(h2, (0, 1), 1.5, "stick")
    assert np.allclose(stick.positions[0], h2.positions[0])
    assert abs(stick.get_distance(0, 1) - 1.5) < 1e-9
    rec = set_coordinate(h2, (0, 1), 1.5, "recomb")
    assert abs(rec.get_distance(0, 1) - 1.5) < 1e-9
    assert np.allclose(
        0.5 * (rec.positions[0] + rec.positions[1]), 0.5 * (h2.positions[0] + h2.positions[1])
    )


def test_relax_with_fixed_bond_keeps_the_distance_and_caches_separately(tmp_path):
    h2 = set_coordinate(hydrogen_molecule(), (0, 1), 1.0, "recomb")
    calc = EMTCalculator()
    fixed, res = relax(h2, calc, tmp_path, fmax=0.05, constraints=[FixBondLength(0, 1)])
    assert abs(fixed.get_distance(0, 1) - 1.0) < 1e-6
    free, res_free = relax(h2, calc, tmp_path, fmax=0.05)
    assert abs(free.get_distance(0, 1) - 1.0) > 1e-3
    assert res_free.energy_eV < res.energy_eV
    again, res2 = relax(h2, calc, tmp_path, fmax=0.05, constraints=[FixBondLength(0, 1)])
    assert res2.cached and abs(res2.energy_eV - res.energy_eV) < 1e-12


def _adatom_states():
    slab = fcc100("Al", size=(2, 2, 2), vacuum=6.0)
    add_adsorbate(slab, "Au", 1.7, "hollow")
    slab.set_constraint(FixAtoms(mask=[a.tag > 1 for a in slab]))
    slab.pbc = (True, True, False)
    final = slab.copy()
    final.positions[-1, 0] += slab.cell[0, 0] / 2.0
    return slab, final


def test_neb_barrier_on_the_emt_adatom_hop_is_positive_and_resumes(tmp_path):
    calc = EMTCalculator()
    initial, final = _adatom_states()
    initial, _ = relax(initial, calc, None, fmax=0.05)
    final, _ = relax(final, calc, None, fmax=0.05)
    traj = tmp_path / "band.traj"
    res = neb_barrier(initial, final, calc, n_images=3, fmax=0.3, steps=60, trajectory=traj)
    assert len(res.energies_eV) == 3 and len(res.images) == 3
    assert res.forward_barrier_eV > 0.0 and res.reverse_barrier_eV > 0.0
    assert traj.exists()
    again = neb_barrier(initial, final, calc, n_images=3, fmax=0.3, steps=60, trajectory=traj)
    if res.converged:
        assert again.cached and again.n_steps == res.n_steps
        assert again.energies_eV == res.energies_eV
    else:
        assert not again.cached and again.n_steps > res.n_steps
    # a band that has taken twice the budget without converging is returned as it stands
    sidecar = traj.with_suffix(".json")
    meta = json.loads(sidecar.read_text())
    meta.update(converged=False, n_steps=120)
    sidecar.write_text(json.dumps(meta))
    stale = neb_barrier(initial, final, calc, n_images=3, fmax=0.3, steps=60, trajectory=traj)
    assert stale.cached and not stale.converged and stale.n_steps == 120
    assert stale.energies_eV == tuple(meta["energies_eV"])


def test_pathways_follow_the_catalogue_configurations():
    spec = BarrierSpec(
        tier="verify", recipes={"basal": {"supercell": 3}, "monovacancy": {"supercell": 3}}
    )
    pws = pathways(spec)
    names = {p.key for p in pws}
    assert "basal/graphene_3x3/monomer_stick" in names
    assert "basal/graphene_3x3/ortho_recomb" in names
    assert "basal/graphene_3x3/para_recomb" in names
    assert "monovacancy/monovacancy_3x3/dangling_pair_recomb" in names
    tubes = [p for p in pws if p.name == "ortho_exterior_recomb"]
    assert len(tubes) == 2 and len({p.config.substrate_name for p in tubes}) == 2
    for p in pws:
        n_sub = len(p.config.substrate)
        assert p.coordinate == (p.config.placements[0][0], n_sub)
        if p.kind == "stick":
            assert p.hydrogens == (n_sub,)
            assert p.quantity.value == "E_barrier_stick" and p.reference_state.value == "H_atom"
        else:
            assert p.hydrogens == (n_sub, n_sub + 1) and p.config.n_H == 2
            assert p.quantity.value == "E_barrier_recomb"
            assert p.reference_state.value == "adsorbed_state"
    only = pathways(BarrierSpec(tier="verify", pathways=["monomer_stick"], recipes=spec.recipes))
    bigger = pathways(
        BarrierSpec(tier="verify", recipes=spec.recipes, stick_recipes={"basal": {"supercell": 4}})
    )
    by_name = {p.name: p for p in bigger}
    assert by_name["monomer_stick"].config.substrate_name == "graphene_4x4"
    assert by_name["ortho_recomb"].config.substrate_name == "graphene_3x3"
    assert [p.name for p in only] == ["monomer_stick"]


def test_h2_above_places_the_molecule_over_the_pair():
    from sitewindow.atomistic.structures.builders import add_hydrogens, graphene_supercell

    g = graphene_supercell(3)
    gh = add_hydrogens(g, [(0, "up"), (1, "up")])
    n = np.array([0.0, 0.0, 1.0])
    mol = h2_above(gh, [len(g), len(g) + 1], [0, 1], n, 3.0)
    assert abs(mol.get_distance(len(g), len(g) + 1) - 0.75) < 1e-9
    mid = 0.5 * (g.positions[0] + g.positions[1])
    assert (
        abs(0.5 * (mol.positions[len(g), 2] + mol.positions[len(g) + 1, 2]) - mid[2] - 3.0) < 1e-9
    )
    assert np.allclose(mol.positions[: len(g)], gh.positions[: len(g)])


def test_primary_method_prefers_scan_for_sticking_and_converged_neb_for_recombination():
    both = {
        CONSTRAINED_SCAN: {"barrier_eV": 1.0, "converged": True},
        CI_NEB: {"barrier_eV": 0.9, "converged": True},
    }
    assert primary_method("stick", both) == CONSTRAINED_SCAN
    assert primary_method("recomb", both) == CI_NEB
    unconverged = {
        CONSTRAINED_SCAN: {"barrier_eV": 1.0, "converged": True},
        CI_NEB: {"barrier_eV": 0.9, "converged": False},
    }
    assert primary_method("recomb", unconverged) == CONSTRAINED_SCAN
    assert primary_method("recomb", {CI_NEB: {"barrier_eV": 0.9, "converged": False}}) == CI_NEB
    failed = {CI_NEB: {"barrier_eV": None, "converged": False}}
    assert primary_method("recomb", failed) is None
    assert primary_method("stick", {}) is None
    neither = {
        CONSTRAINED_SCAN: {"barrier_eV": 2.5, "converged": True, "saddle_located": False},
        CI_NEB: {"barrier_eV": 31.0, "converged": False, "saddle_located": False},
    }
    assert primary_method("recomb", neither) == CONSTRAINED_SCAN
    wrong_product = {
        CONSTRAINED_SCAN: {"barrier_eV": 1.0, "converged": True, "saddle_located": True},
        CI_NEB: {"barrier_eV": 0.9, "converged": True, "saddle_located": False},
    }
    assert primary_method("recomb", wrong_product) == CONSTRAINED_SCAN
    assert primary_method("stick", wrong_product) == CONSTRAINED_SCAN


def test_product_distance_and_saddle_caveat_name_the_wrong_product():
    h2 = hydrogen_molecule()
    d, ok = product_distance("recomb", h2, (0, 1), [0, 1], 3.6, 0.75)
    assert ok and d == pytest.approx(h2.get_distance(0, 1))
    far = set_coordinate(h2, (0, 1), 2.9, "recomb")
    d, ok = product_distance("recomb", far, (0, 1), [0, 1], 3.6, 0.75)
    assert not ok and d == pytest.approx(2.9)
    assert product_distance("stick", far, (0, 1), [1], 3.6, 0.75)[1]
    assert not product_distance("stick", h2, (0, 1), [1], 3.6, 0.75)[1]
    assert saddle_caveat("recomb", True, True, False, 0.75) is None
    text = saddle_caveat("recomb", False, False, True, 2.9)
    assert text is not None
    assert "not every point converged" in text
    assert "H-H distance 2.90 A" in text and "end of the grid" in text
    assert (
        saddle_caveat("stick", True, True, True, 3.6, CI_NEB) == "the highest image is an end point"
    )
    assert saddle_caveat("stick", False, True, False, 3.6, CI_NEB) == "the band did not converge"
    smooth = saddle_caveat("recomb", True, True, False, 0.75, discontinuity=(0.3, 1.9, 2.0))
    assert smooth is None
    jump = saddle_caveat("recomb", True, True, False, 0.75, discontinuity=(2.69, 2.16, 2.18))
    assert jump == "the constrained path is discontinuous (2.69 eV between 2.16 and 2.18 A)"
    assert saddle_caveat("stick", True, True, False, 3.6, far_offset_eV=0.012) is None
    off = saddle_caveat("stick", True, True, False, 3.6, far_offset_eV=0.47)
    assert off == "the far end of the path lies +0.47 eV from the separated reactants"
    pts = [ScanPoint(1.1, -1.0, True, None), ScanPoint(1.2, -0.5, True, None)]
    pts += [ScanPoint(1.3, float("nan"), False, None), ScanPoint(1.4, 1.5, True, None)]
    assert largest_step(pts) == (2.0, 1.2, 1.4)
    assert largest_step(pts[:1]) is None


def test_build_barriers_writes_validated_tables_with_both_methods(tmp_path):
    spec = BarrierSpec(
        tier="verify",
        pathways=["monomer_stick", "ortho_recomb"],
        recipes={"basal": {"supercell": 3}},
        methods=(CONSTRAINED_SCAN, CI_NEB),
        fmax=0.3,
        steps=60,
        grid_step_A=0.3,
        refine_step_A=0.1,
        refine_halfwidth_A=0.3,
        far_distance_A=2.6,
        recomb_distance_A=2.4,
        neb_images=3,
        neb_steps=30,
        stick_calculator_settings={"extra": {"variant": "stick"}},
    )
    out = build_barriers(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    frame = out["records"]
    recs = frame_to_records(frame)
    assert len(recs) == len(frame) == 2
    assert set(frame["quantity"]) == {"E_barrier_stick", "E_barrier_recomb"}
    assert set(frame["reference_state"]) == {"H_atom", "adsorbed_state"}
    methods = dict(zip(frame["quantity"], frame["identifier"].str.split("@").str[1]))
    assert methods["E_barrier_stick"] == CONSTRAINED_SCAN
    table, meta = read_table(tmp_path / "out" / "barriers_verify.parquet")
    assert len(table) == 2 and meta["tier"] == "verify"
    profiles = out["profiles"]
    assert len(profiles) == 2
    for prof in profiles:
        assert set(prof["methods"]) == {CONSTRAINED_SCAN, CI_NEB}
        assert prof["primary_method"] in prof["methods"]
        scan = prof["methods"][CONSTRAINED_SCAN]
        pts = scan["points"]
        assert all(p["coordinate_A"] < q["coordinate_A"] for p, q in zip(pts, pts[1:]))
        assert pts[0]["energy_eV"] == pytest.approx(prof["E_initial_eV"])
        assert scan["transition_structure"].splitlines()[0].strip().isdigit()
        assert len(prof["methods"][CI_NEB]["image_energies_eV"]) == 3
    for prof in profiles:
        assert prof["methods"][CONSTRAINED_SCAN]["points"][0]["coordinate_A"] == pytest.approx(
            prof["initial_coordinate_A"]
        )
    recomb = next(p for p in profiles if p["kind"] == "recomb")
    assert recomb["initial_H_H_distance_A"] > 1.0 and len(recomb["hydrogen_indices"]) == 2
    for prof in profiles:
        for r in prof["methods"].values():
            assert isinstance(r["product_reached"], bool)
            assert isinstance(r["saddle_located"], bool)
            assert (r["saddle_caveat"] is None) == r["saddle_located"]
            assert r["final_distance_A"] is None or r["final_distance_A"] > 0.0
        primary = prof["methods"][prof["primary_method"]]
        row = frame[frame["identifier"] == f"{prof['pathway']}@{prof['primary_method']}"].iloc[0]
        assert bool(row["converged"]) == primary["saddle_located"]
    stick = next(p for p in profiles if p["kind"] == "stick")
    stick_row = frame[frame["quantity"] == "E_barrier_stick"].iloc[0]
    recomb_row = frame[frame["quantity"] == "E_barrier_recomb"].iloc[0]
    assert stick_row["settings_hash"] != recomb_row["settings_hash"]
    assert out["manifest"]["stick_settings_hash"] == stick_row["settings_hash"]
    assert out["manifest"]["calculator"]["settings_hash"] == recomb_row["settings_hash"]
    assert stick["methods"][CONSTRAINED_SCAN]["n_excluded_points"] == 0
    assert stick["methods"][CONSTRAINED_SCAN]["far_point_offset_eV"] is not None
    assert stick["methods"][CI_NEB]["far_point_offset_eV"] is not None
    assert stick["methods"][CONSTRAINED_SCAN]["largest_step_eV"] >= 0.0
    assert len(stick["methods"][CONSTRAINED_SCAN]["largest_step_between_A"]) == 2
    assert recomb["methods"][CONSTRAINED_SCAN]["far_point_offset_eV"] is None
    assert (tmp_path / "out" / "barrier_profiles_verify.json").exists()
    assert (tmp_path / "out" / "barriers_manifest_verify.json").exists()
    rerun = build_barriers(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    scan_vals = [p["methods"][CONSTRAINED_SCAN]["barrier_eV"] for p in rerun["profiles"]]
    assert np.allclose(scan_vals, [p["methods"][CONSTRAINED_SCAN]["barrier_eV"] for p in profiles])


def test_resample_path_keeps_the_end_points_and_spaces_images_evenly():
    h2 = hydrogen_molecule()
    geoms = []
    for d in (0.75, 0.9, 1.3, 2.0, 2.1):
        geoms.append(set_coordinate(h2, (0, 1), d, "recomb"))
    path = resample_path(geoms, 5)
    assert len(path) == 5
    assert np.allclose(path[0].positions, geoms[0].positions)
    assert np.allclose(path[-1].positions, geoms[-1].positions)
    dists = [p.get_distance(0, 1) for p in path]
    gaps = np.diff(dists)
    assert np.allclose(gaps, gaps[0], atol=1e-9)
    with pytest.raises(ValueError):
        resample_path([h2], 3)


def test_relax_retries_a_failed_scf_with_the_robust_calculator(tmp_path):
    from ase.calculators.calculator import CalculationFailed

    calls = {"n": 0}

    @dataclass
    class FlakyCalculator(EMTCalculator):
        name: str = "flaky"

        def attach(self, atoms):
            from ase.calculators.emt import EMT

            calls["n"] += 1
            if calls["n"] == 1:

                class Failing(EMT):
                    def calculate(self, *args, **kwargs):
                        raise CalculationFailed("SCF not converged")

                atoms.calc = Failing()
            else:
                atoms.calc = EMT()

        def robust(self):
            return self

    h2 = hydrogen_molecule()
    relaxed, res = relax(h2, FlakyCalculator(), tmp_path, fmax=0.05)
    assert res.scf_retried and res.converged and calls["n"] == 2
    again, res2 = relax(h2, FlakyCalculator(), tmp_path, fmax=0.05)
    assert res2.cached and res2.scf_retried


def test_profile_seed_returns_the_located_transition_structure(tmp_path):
    import io
    import json

    from ase.io import write as ase_write

    from sitewindow.atomistic.workflows.barriers import profile_seed

    def xyz(a):
        buf = io.StringIO()
        ase_write(buf, a, format="extxyz")
        return buf.getvalue()

    scan_ts = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.0]], cell=[8, 8, 8], pbc=True)
    neb_ts = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.5]], cell=[8, 8, 8], pbc=True)
    profiles = [
        {
            "pathway": "basal/x/ortho_recomb",
            "methods": {
                "constrained_scan": {"saddle_located": False, "transition_structure": xyz(scan_ts)},
                "ci_neb": {"saddle_located": True, "transition_structure": xyz(neb_ts)},
            },
        }
    ]
    p = tmp_path / "profiles.json"
    p.write_text(json.dumps(profiles))
    cache = {}
    seed = profile_seed(p, "basal/x/ortho_recomb", 2, cache)
    assert len(seed) == 1 and abs(seed[0].get_distance(0, 1) - 1.5) < 1e-6
    assert profile_seed(p, "basal/x/ortho_recomb", 3, cache) == []
    assert profile_seed(p, "basal/x/other", 2, cache) == []
    assert profile_seed(None, "basal/x/ortho_recomb", 2) == []


def test_refinement_single_points_add_a_refined_scan_record(tmp_path):
    from sitewindow.atomistic.workflows.barriers import REFINED_SCAN

    spec = BarrierSpec(
        tier="verify",
        pathways=["monomer_stick"],
        recipes={"basal": {"supercell": 3}},
        methods=(CONSTRAINED_SCAN,),
        fmax=0.3,
        steps=60,
        grid_step_A=0.3,
        refine_step_A=0.1,
        far_distance_A=2.2,
        refine_calculator=EMTCalculator(extra={"role": "refine"}),
    )
    out = build_barriers(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    frame = out["records"]
    methods = frame["identifier"].str.split("@").str[1].tolist()
    assert CONSTRAINED_SCAN in methods and REFINED_SCAN in methods
    scan = frame[frame["identifier"].str.endswith(CONSTRAINED_SCAN)].iloc[0]
    refined = frame[frame["identifier"].str.endswith(REFINED_SCAN)].iloc[0]
    # the same calculator evaluates both, so the refined barrier reproduces the scan barrier
    assert refined["value_eV"] == pytest.approx(scan["value_eV"], abs=1e-6)
    assert refined["settings_hash"] != scan["settings_hash"]
    profile = out["profiles"][0]["methods"][REFINED_SCAN]
    assert profile["structures_from"] == CONSTRAINED_SCAN and "E_transition_eV" in profile
    assert out["manifest"]["refinement"]["method"] == REFINED_SCAN
