"""End-to-end landscape construction on a reduced class set with the EMT test calculator."""

from dataclasses import dataclass, field

import pytest
from ase import Atoms

from sitewindow.atomistic.calculators.gpaw_adapter import GPAWCalculator, scaled_kmesh
from sitewindow.atomistic.structures.builders import graphene_supercell, hydrogen_atom, swcnt
from sitewindow.atomistic.workflows.landscape import LandscapeSpec, build_landscape
from sitewindow.atomistic.workflows.tier1 import Tier1Settings, check_h2_bond, convergence_variants
from sitewindow.io import read_table
from sitewindow.schema import Quantity, frame_to_records


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


def test_build_landscape_writes_validated_tables(tmp_path):
    spec = LandscapeSpec(
        tier="verify",
        classes=["basal", "monovacancy", "cnt_armchair"],
        recipes={
            "basal": {"supercell": 3},
            "monovacancy": {"supercell": 3, "max_hydrogens": 2},
            "cnt_armchair": {"diameters_nm": [0.55], "length": 2},
        },
        fmax=0.3,
        steps=60,
    )
    out = build_landscape(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    radial = out["manifest"]["relaxation"]["tube_radius_optimization"]
    assert set(radial) == {"swcnt_4_4_x2"}
    assert radial["swcnt_4_4_x2"]["dE_eV"] <= 0.0 and radial["swcnt_4_4_x2"]["n_single_points"] >= 3
    frame = out["records"]
    recs = frame_to_records(frame)
    assert len(recs) == len(frame) > 0
    quantities = set(frame["quantity"])
    assert {"E_bind_H_half_H2", "E_bind_H_atom", "E_sequential", "E_formation_defect"} <= quantities
    assert (frame["tier"] == "T3").all() and frame["settings_hash"].notna().all()
    seq = frame[(frame.quantity == Quantity.E_sequential.value) & (frame.site_class == "basal")]
    assert set(seq["n_H"]) == {2}
    table, meta = read_table(tmp_path / "out" / "landscape_raw_verify.parquet")
    assert meta["sitewindow_schema_version"] and meta["tier"] == "verify"
    assert len(table) == len(frame)
    cfgs, _ = read_table(tmp_path / "out" / "landscape_configurations_verify.parquet")
    assert set(cfgs["config_hash"]) <= set(frame["config_hash"])
    assert (tmp_path / "out" / "landscape_manifest_verify.json").exists()
    rerun = build_landscape(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    assert rerun["configurations"]["cached"].all()


def test_lower_energy_prefers_converged_then_lower_energy():
    from sitewindow.atomistic.workflows.adsorption import RelaxResult
    from sitewindow.atomistic.workflows.landscape import _lower_energy

    def result(energy: float, converged: bool) -> RelaxResult:
        return RelaxResult(energy, 0.01, 1, converged, "s", "t", 0.0, False)

    a = (Atoms("H"), result(-1.0, True), 0)
    b = (Atoms("H"), result(-2.0, True), 0)
    c = (Atoms("H"), result(-5.0, False), 0)
    assert _lower_energy(a, b) is b and _lower_energy(b, a) is b
    assert _lower_energy(a, c) is a and _lower_energy(c, a) is a
    assert _lower_energy(a, (Atoms("H"), result(-1.0, True), 0)) is a


def test_restart_table_gives_a_second_start(tmp_path):
    spec = LandscapeSpec(
        tier="verify", classes=["basal"], recipes={"basal": {"supercell": 3}}, fmax=0.3, steps=60
    )
    build_landscape(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    table = tmp_path / "out" / "landscape_configurations_verify.parquet"
    key = "basal/graphene_3x3/monomer"
    spec2 = LandscapeSpec(
        tier="verify",
        classes=["basal"],
        recipes={"basal": {"supercell": 3}},
        fmax=0.3,
        steps=60,
        restarts=table,
        restart_keys=(key,),
    )
    second = build_landscape(spec2, EMTCalculator(), tmp_path / "cache", tmp_path / "out2")
    cfgs = second["configurations"]
    assert set(cfgs["start"]) <= {"builder", "restart"}
    assert (cfgs[cfgs.configuration != "monomer"]["start"] == "builder").all()
    restarts = second["manifest"]["relaxation"]["restarts"]
    assert restarts["source"] == str(table) and restarts["keys"] == [key]
    assert list(restarts["chosen"]) == [key]
    assert abs(restarts["chosen"][key]["dE_restart_minus_builder_eV"]) < 1e-3
    assert restarts["chosen"][key]["restart_converged"]


def test_scaled_kmesh_reproduces_manual_protocol():
    assert scaled_kmesh(graphene_supercell(5)) == (4, 4, 1)
    assert scaled_kmesh(graphene_supercell(6)) == (3, 3, 1)
    assert scaled_kmesh(swcnt(5, 5, 3))[:2] == (1, 1) and scaled_kmesh(swcnt(5, 5, 3))[2] > 1
    assert GPAWCalculator().kpoints(hydrogen_atom()) == (1, 1, 1)


def test_tier1_settings_and_variants():
    base = Tier1Settings(supercell=5)
    v = convergence_variants(base)
    assert set(v) == {"base", "cell_plus_one", "kmesh_dense", "cutoff_600", "vacuum_20"}
    assert v["cell_plus_one"].supercell == 6 and v["cutoff_600"].ecut_eV == 600.0
    assert v["kmesh_dense"].kmesh == (6, 6, 1) and v["vacuum_20"].vacuum_A == 20.0
    pw = base.pw()
    assert pw.settings()["ecut_eV"] == 500.0 and pw.settings()["dispersion"] == "d3bj"
    assert base.molecule_pw().kmesh == (1, 1, 1)


def test_gpaw_adapter_reports_missing_dependency():
    pytest.importorskip("ase")
    try:
        import gpaw  # noqa: F401

        pytest.skip("gpaw installed")
    except ImportError:
        pass
    from sitewindow.errors import CalculatorUnavailable

    with pytest.raises(CalculatorUnavailable):
        GPAWCalculator().attach(graphene_supercell(2))


def test_h2_reference_guard_and_spin_pairing():
    from sitewindow.atomistic.structures.builders import hydrogen_molecule, set_unpaired_electrons
    from sitewindow.errors import ConvergenceError

    check_h2_bond(0.75)
    with pytest.raises(ConvergenceError):
        check_h2_bond(3.61)
    h2 = hydrogen_molecule(pbc=True)
    set_unpaired_electrons(h2, 0)
    assert h2.get_initial_magnetic_moments().sum() == 0.0


def test_relax_tube_radius_lowers_the_energy_and_leaves_sheets_alone(tmp_path):
    import numpy as np

    from sitewindow.atomistic.workflows.adsorption import relax_tube_radius, single_point

    calc = EMTCalculator()
    t = swcnt(4, 4, length=2)
    scaled, rec = relax_tube_radius(t, calc, tmp_path / "cache")
    e0 = single_point(t, calc, tmp_path / "cache").energy_eV
    e1 = single_point(scaled, calc, tmp_path / "cache").energy_eV
    assert rec["dE_eV"] <= 0.0 and abs(e1 - e0 - rec["dE_eV"]) < 1e-9
    assert abs(rec["scale"] - 1.0) < 0.1 and rec["n_single_points"] >= 3
    centre = np.array([t.cell[0, 0], t.cell[1, 1]]) / 2.0
    r0 = np.linalg.norm(t.positions[:, :2] - centre, axis=1)
    r1 = np.linalg.norm(scaled.positions[:, :2] - centre, axis=1)
    assert np.allclose(r1 / r0, rec["scale"])
