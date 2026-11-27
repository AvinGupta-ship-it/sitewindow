"""Electronic-entropy decomposition of tblite energies and the extrapolated-energy option."""

import numpy as np
import pytest
from ase.units import kB

from sitewindow.atomistic.calculators.tblite_adapter import fermi_entropy, smearing_terms
from sitewindow.atomistic.workflows.landscape import ENERGY_KINDS, LandscapeSpec, build_landscape


def test_fermi_entropy_counts_spin_orbitals_and_handles_restricted_occupations():
    assert fermi_entropy(np.array([[1.0, 0.0], [1.0, 0.0]])) == 0.0
    half = np.array([[0.5], [0.5]])
    assert abs(fermi_entropy(half) - 2.0 * np.log(2.0)) < 1e-12
    restricted = np.array([[1.0]])  # spin-summed occupation 1 -> 0.5 per spin
    assert abs(fermi_entropy(restricted) - 2.0 * np.log(2.0)) < 1e-12
    third = np.full((2, 3), 1.0 / 3.0)
    s = -6.0 * (np.log(1 / 3) / 3 + 2 / 3 * np.log(2 / 3))
    assert abs(fermi_entropy(third) - s) < 1e-12


class _Result:
    def __init__(self, occ):
        self._occ = occ

    def get(self, key):
        assert key == "orbital-occupations"
        return self._occ


class _Calc:
    def __init__(self, energy, occ):
        self.results = {"energy": energy}
        self._res = _Result(occ)


def test_smearing_terms_reconstruct_internal_energy_and_extrapolation():
    occ = np.array([[1.0, 0.5, 0.0], [1.0, 0.5, 0.0]])
    T = 3000.0
    terms = smearing_terms(_Calc(-10.0, occ), T)
    ts = kB * T * 2.0 * np.log(2.0)
    assert abs(terms["TS_eV"] - ts) < 1e-12
    assert abs(terms["E_internal_eV"] - (-10.0 + ts)) < 1e-12
    assert abs(terms["E0_eV"] - (-10.0 + 0.5 * ts)) < 1e-12
    assert terms["n_fractional_orbitals"] == 2
    assert terms["F_eV"] == -10.0


def test_landscape_rejects_the_extrapolation_for_engines_without_smearing(tmp_path):
    from test_landscape_pipeline import EMTCalculator

    spec = LandscapeSpec(
        tier="verify", classes=["basal"], recipes={"basal": {"supercell": 3}}, energy="E0"
    )
    with pytest.raises(ValueError, match="tblite"):
        build_landscape(spec, EMTCalculator(), tmp_path / "cache", tmp_path / "out")
    with pytest.raises(ValueError, match="energy must be"):
        build_landscape(
            LandscapeSpec(tier="verify", classes=["basal"], energy="X"),
            EMTCalculator(),
            tmp_path / "cache",
            tmp_path / "out",
        )
    assert set(ENERGY_KINDS) == {"F", "E0"}
