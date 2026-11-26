"""Tier-1 anchor-set workflow on a test engine: records, checkpointing, table and sequential energies."""

from dataclasses import dataclass, field

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
