"""Barrier workflow (C2, B2): coordinate scan, climbing-image NEB, pathway catalogue, tables."""

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
    neb_barrier,
    pathways,
    primary_method,
    resample_path,
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
        assert again.cached and again.n_steps == 0
        assert again.energies_eV == res.energies_eV
    else:
        assert not again.cached


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
