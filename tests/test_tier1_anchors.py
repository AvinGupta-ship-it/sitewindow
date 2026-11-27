"""Tier-1 anchor-set workflow on a test engine: records, checkpointing, table and sequential energies."""

from dataclasses import dataclass, field

import numpy as np
import pytest
import yaml
from ase import Atoms

from sitewindow.atomistic.workflows import tier1
from sitewindow.atomistic.workflows.tier1 import (
    ANCHOR_REGISTRY,
    Tier1Settings,
    anchor_configuration,
    collect_anchor_table,
    load_anchor_registry,
    resolve_anchors,
    run_anchor_set,
)
from sitewindow.io import read_table
from sitewindow.landscape.classes import catalogue
from sitewindow.schema import frame_to_records


@dataclass
class EMTCalc:
    name: str = "emt"
    tier: str = "T1"
    label: str = "lcao"
    extra: dict = field(default_factory=dict)

    def attach(self, atoms: Atoms) -> None:
        from ase.calculators.emt import EMT

        atoms.calc = EMT()

    def settings(self) -> dict:
        return {"model": "emt", "mode": self.label, **self.extra}

    def kpoints(self, atoms: Atoms) -> tuple[int, int, int]:
        return (1, 1, 1)


@dataclass(frozen=True)
class EMTEngine:
    """Effective-medium stand-in for the GPAW engine: no pre-relaxation, no dispersion."""

    def prerelax(self):
        return None

    def lcao(self, txt=None):
        return EMTCalc(label="lcao")

    def pw(self, txt=None):
        return EMTCalc(label="pw")

    def molecule_pw(self, txt=None):
        return EMTCalc(label="pw")

    def dispersion(self, atoms: Atoms) -> float:
        return 0.0

    def name(self) -> str:
        return "emt"


def _config(cls: str, name: str, recipe: dict):
    return next(c for c in catalogue()[cls].configurations("full", recipe) if c.name == name)


def test_anchor_registry_resolves_every_class_once_per_round_one():
    reg = load_anchor_registry(ANCHOR_REGISTRY)
    resolved = resolve_anchors(reg)
    assert 40 <= len(resolved) <= 60
    round_one = [cfg.site_class.value for rnd, cfg in resolved if rnd == 1]
    chem = [k for k, d in catalogue().items() if "T2" in d.tiers]
    assert sorted(round_one) == sorted(chem)
    sizes = [len(cfg.substrate) + cfg.n_H for rnd, cfg in resolved if rnd == 1]
    assert sizes == sorted(sizes)
    assert reg["protocol"]["ecut_eV"] == 500 and reg["protocol"]["fmax_eV_A"] == 0.02


def test_anchor_configuration_writes_a_record_with_binding_and_formation(tmp_path):
    settings = Tier1Settings(supercell=3, fmax_eV_A=0.3, relax_steps=40, prerelax_tier2=False)
    cfg = _config("monovacancy", "seq1", {"supercell": 3})
    res = anchor_configuration(
        cfg, settings, tmp_path / "cache", tmp_path / "out", engine=EMTEngine()
    )
    assert res["n_H"] == 1 and res["site_class"] == "monovacancy"
    b = res["binding_eV"]
    assert b["E_bind_H_atom_PBE_D3"] == pytest.approx(b["E_bind_H_atom_PBE"])
    assert res["formation_eV"] is not None and "E_formation_PBE_D3" in res["formation_eV"]
    assert res["descriptors"]["coordination"] >= 2
    assert len(res["spin_guess"]["substrate_sites"]) == 3 and res["spin_guess"]["moment_muB"] == 1.0
    assert set(res["magnetic_moments_pw"]) == {"substrate", "configuration"}
    assert (tmp_path / "out" / "anchors" / f"{res['tag']}.json").exists()
    assert res["structure"].splitlines()[0].strip() == str(res["n_atoms"])


def test_run_anchor_set_checkpoints_and_builds_sequential_records(tmp_path, monkeypatch):
    registry = {
        "protocol": {"ecut_eV": 500},
        "anchors": [
            {
                "round": 1,
                "site_class": "basal",
                "configuration": "monomer",
                "recipe": {"supercell": 3},
            },
            {
                "round": 2,
                "site_class": "basal",
                "configuration": "ortho",
                "recipe": {"supercell": 3},
            },
        ],
    }
    reg_path = tmp_path / "anchors.yaml"
    reg_path.write_text(yaml.safe_dump(registry))
    settings = Tier1Settings(fmax_eV_A=0.3, relax_steps=40, prerelax_tier2=False)
    plan = run_anchor_set(settings, tmp_path / "cache", tmp_path / "out", reg_path, dry_run=True)
    assert [p["status"] for p in plan["plan"]] == ["pending", "pending"]
    first = run_anchor_set(
        settings, tmp_path / "cache", tmp_path / "out", reg_path, rounds=[1], engine=EMTEngine()
    )
    assert [p["status"] for p in first["plan"]] == ["done", "pending"]
    assert first["manifest"]["n_anchors"] == 1
    second = run_anchor_set(
        settings, tmp_path / "cache", tmp_path / "out", reg_path, engine=EMTEngine()
    )
    assert all(p["status"] == "done" for p in second["plan"])
    assert "wall_s" not in second["plan"][0]
    frame = second["records"]
    recs = frame_to_records(frame)
    assert len(recs) == len(frame) > 0
    assert (frame["tier"] == "T1").all()
    assert set(frame["engine"]) == {"gpaw:PBE-D3(BJ)", "gpaw:PBE"}
    seq = frame[frame["quantity"] == "E_sequential"]
    assert set(seq["n_H"]) == {2} and len(seq) == 4
    table, meta = read_table(tmp_path / "out" / "anchors_tier1.parquet")
    assert meta["tier"] == "T1" and len(table) == len(frame)
    again = collect_anchor_table(tmp_path / "out")
    assert again["manifest"]["n_records"] == len(frame)


def test_gpaw_engine_requires_gpaw_or_reports_it():
    settings = Tier1Settings()
    engine = tier1.GPAWEngine(settings)
    assert engine.name() == "gpaw"
    assert engine.pw().settings()["ecut_eV"] == 500.0
    assert engine.lcao().settings()["basis"] == "dzp"


class _FailingEngine(EMTEngine):
    """Raises a convergence error for every structure with more than ``limit`` atoms."""

    def lcao(self, txt=None):
        return _FailingCalc(label="lcao")

    def pw(self, txt=None):
        return _FailingCalc(label="pw")


@dataclass
class _FailingCalc(EMTCalc):
    limit: int = 20

    def attach(self, atoms: Atoms) -> None:
        from ase.calculators.emt import EMT

        from sitewindow.errors import ConvergenceError

        if len(atoms) > self.limit:
            raise ConvergenceError("self-consistent field did not converge")
        atoms.calc = EMT()


def test_run_anchor_set_records_a_failed_anchor_and_continues(tmp_path):
    registry = {
        "anchors": [
            {
                "round": 1,
                "site_class": "basal",
                "configuration": "monomer",
                "recipe": {"supercell": 3},
            },
            {
                "round": 1,
                "site_class": "basal",
                "configuration": "monomer",
                "recipe": {"supercell": 4},
            },
        ],
    }
    reg_path = tmp_path / "anchors.yaml"
    reg_path.write_text(yaml.safe_dump(registry))
    settings = Tier1Settings(fmax_eV_A=0.3, relax_steps=40, prerelax_tier2=False)
    res = run_anchor_set(
        settings, tmp_path / "cache", tmp_path / "out", reg_path, engine=_FailingEngine()
    )
    statuses = {p["key"]: p["status"] for p in res["plan"]}
    assert statuses["basal/graphene_3x3/monomer"] == "done"
    assert statuses["basal/graphene_4x4/monomer"] == "failed"
    failed = list((tmp_path / "out" / "anchors_failed").glob("*.json"))
    assert len(failed) == 1 and res["manifest"]["failed"][0]["key"] == "basal/graphene_4x4/monomer"
    again = run_anchor_set(
        settings, tmp_path / "cache", tmp_path / "out", reg_path, engine=_FailingEngine()
    )
    assert {p["key"]: p["status"] for p in again["plan"]} == statuses


def test_gpaw_robust_variant_keeps_the_settings_hash():
    from sitewindow.atomistic.calculators.base import CalculatorSpec
    from sitewindow.atomistic.calculators.gpaw_adapter import ROBUST_MIXER, GPAWCalculator

    calc = GPAWCalculator(mode="pw")
    robust = calc.robust()
    assert robust.extra["mixer"] == ROBUST_MIXER and robust.maxiter >= 400
    assert (
        CalculatorSpec("gpaw", "T1", calc.settings()).settings_hash
        == CalculatorSpec("gpaw", "T1", robust.settings()).settings_hash
    )


def test_scf_failure_detection_matches_gpaw_errors_by_name():
    from ase.calculators.calculator import CalculationFailed

    from sitewindow.atomistic.workflows.adsorption import is_scf_failure

    class ConvergenceError(Exception):
        pass

    class KohnShamConvergenceError(ConvergenceError):
        pass

    assert is_scf_failure(KohnShamConvergenceError())
    assert is_scf_failure(CalculationFailed("x"))
    assert not is_scf_failure(ValueError("x"))


def test_final_magnetic_moments_reads_the_dft_part_of_a_sum_calculator():
    from ase.calculators.calculator import PropertyNotImplementedError

    from sitewindow.atomistic.workflows.adsorption import final_magnetic_moments

    class Dft:
        def get_property(self, name, atoms=None, allow_calculation=True):
            assert name == "magmoms"
            return np.array([0.7, -0.1])

    class Dispersion:
        def get_property(self, name, atoms=None, allow_calculation=True):
            raise PropertyNotImplementedError(name)

    class Mixer:
        calcs = [Dispersion(), Dft()]

    class Sum:
        mixer = Mixer()

        def get_magnetic_moments(self, atoms):
            raise PropertyNotImplementedError("magmoms")

    atoms = Atoms("CH", positions=[[0, 0, 0], [0, 0, 1.1]])
    atoms.calc = Sum()
    assert final_magnetic_moments(atoms) == pytest.approx([0.7, -0.1])
    atoms.calc = None
    assert final_magnetic_moments(atoms) is None


def test_broken_symmetry_guess_alternates_with_the_sublattice():
    from sitewindow.atomistic.structures.builders import (
        broken_symmetry_moments,
        set_unpaired_electrons,
        sublattices,
        undercoordinated_carbons,
    )
    from sitewindow.atomistic.workflows.adsorption import moment_tag, unpaired_electrons

    ribbon = _config("edge_zigzag", "edge_CH", {"width": 6, "length": 4}).substrate.copy()
    parity = int(np.sum(ribbon.numbers)) % 2
    sites = broken_symmetry_moments(ribbon, parity)
    edge = undercoordinated_carbons(ribbon)
    assert sites == edge and len(edge) == 8
    canonical = ribbon.copy()
    set_unpaired_electrons(canonical, parity)
    guess = ribbon.get_initial_magnetic_moments() - canonical.get_initial_magnetic_moments()
    colour = sublattices(ribbon)
    assert all(guess[k] == (1.0 if colour[k] == 0 else -1.0) for k in edge)
    assert all(guess[k] == 0.0 for k in range(len(ribbon)) if k not in edge)
    assert sum(guess[k] > 0 for k in edge) == sum(guess[k] < 0 for k in edge)
    assert unpaired_electrons(ribbon) == parity and moment_tag(ribbon).startswith("_m")

    vacancy = _config("monovacancy", "seq1", {"supercell": 3}).substrate.copy()
    sites = broken_symmetry_moments(vacancy, 0)
    moments = vacancy.get_initial_magnetic_moments()
    assert len(sites) == 3 and len({np.sign(moments[k]) for k in sites}) == 1
    assert unpaired_electrons(vacancy) == 0

    basal = _config("basal", "monomer", {"supercell": 3}).substrate.copy()
    assert broken_symmetry_moments(basal, 0) == [] and moment_tag(basal) == ""


def test_cache_keys_separate_the_symmetry_broken_start_from_the_canonical_one():
    from sitewindow.atomistic.calculators.base import CalculatorSpec
    from sitewindow.atomistic.structures.builders import (
        broken_symmetry_moments,
        set_unpaired_electrons,
    )
    from sitewindow.atomistic.workflows.adsorption import (
        _cache_key,
        moment_tag,
        unpaired_electrons,
    )

    vacancy = _config("monovacancy", "seq1", {"supercell": 3}).substrate.copy()
    set_unpaired_electrons(vacancy, 0)
    spec = CalculatorSpec("emt", "T1", {"model": "emt"})
    key_canonical = _cache_key(vacancy, spec, 0.05, True)
    broken_symmetry_moments(vacancy, 0)
    key_broken = _cache_key(vacancy, spec, 0.05, True)
    assert key_canonical.endswith("_u0")
    assert key_broken == key_canonical + moment_tag(vacancy) and unpaired_electrons(vacancy) == 0
    assert np.sum(vacancy.get_initial_magnetic_moments()) == pytest.approx(3.0)
