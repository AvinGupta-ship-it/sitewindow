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
        classes=["basal", "monovacancy"],
        recipes={"basal": {"supercell": 3}, "monovacancy": {"supercell": 3, "max_hydrogens": 2}},
        fmax=0.3,
        steps=60,
    )
    out = build_landscape(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
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
