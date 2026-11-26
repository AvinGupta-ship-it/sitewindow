"""GFN-xTB (tblite) calculator adapter: fidelity tier T2."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from ase import Atoms
from ase.units import kB
from scipy.special import xlogy

from sitewindow.atomistic.calculators.base import CalculatorUnavailable


@dataclass
class TBLiteCalculator:
    method: str = "GFN2-xTB"
    accuracy: float = 1.0
    electronic_temperature: float = 300.0
    max_iterations: int = 300
    spin_polarized: bool = True
    name: str = "tblite"
    tier: str = "T2"
    extra: dict = field(default_factory=dict)
    _report_settings: dict | None = field(default=None, repr=False, compare=False)

    def attach(self, atoms: Atoms) -> None:
        try:
            from tblite.ase import TBLite
        except ImportError as e:  # pragma: no cover
            raise CalculatorUnavailable("tblite is not installed; install sitewindow[xtb]") from e
        kwargs = dict(
            method=self.method,
            accuracy=self.accuracy,
            electronic_temperature=self.electronic_temperature,
            max_iterations=self.max_iterations,
            verbosity=0,
        )
        if self.spin_polarized:
            kwargs["spin_polarization"] = 1.0
        kwargs.update(self.extra)
        atoms.calc = TBLite(**kwargs)

    def robust(self) -> TBLiteCalculator:
        """The same Hamiltonian with a damped mixer and more self-consistency iterations, used to
        retry a point whose self-consistent field did not converge; ``settings`` is unchanged so
        that the result is cached under the original key."""
        return TBLiteCalculator(
            method=self.method,
            accuracy=self.accuracy,
            electronic_temperature=self.electronic_temperature,
            max_iterations=max(self.max_iterations, 1000),
            spin_polarized=self.spin_polarized,
            name=self.name,
            tier=self.tier,
            extra={**self.extra, "mixer_damping": 0.2},
            _report_settings=self.settings(),
        )

    def settings(self) -> dict:
        if self._report_settings is not None:
            return dict(self._report_settings)
        return {
            "method": self.method,
            "accuracy": self.accuracy,
            "electronic_temperature": self.electronic_temperature,
            "max_iterations": self.max_iterations,
            "spin_polarized": self.spin_polarized,
            **self.extra,
        }


def fermi_entropy(occupations: np.ndarray) -> float:
    """Electronic entropy S/k_B of Fermi-Dirac occupations.

    ``occupations`` has one row per spin channel with values in [0, 1] (the per-spin orbital
    occupations that tblite returns for a spin-polarized calculation); a single row with values in
    [0, 2] is the spin-summed occupation of a restricted calculation and is halved and counted
    twice. S/k_B = -sum_i [f_i ln f_i + (1 - f_i) ln(1 - f_i)] over spin orbitals.
    """
    occ = np.atleast_2d(np.asarray(occupations, dtype=float))
    if occ.shape[0] == 1:
        occ = np.vstack([occ / 2.0, occ / 2.0])
    f = np.clip(occ, 0.0, 1.0)
    return float(-np.sum(xlogy(f, f) + xlogy(1.0 - f, 1.0 - f)))


def smearing_terms(ase_calculator: Any, electronic_temperature_K: float) -> dict[str, float]:
    """Decompose the last tblite single point into its finite-temperature terms.

    tblite adds the electronic entropy term -T S to the band energy inside its self-consistent
    loop, so the energy it reports is the Mermin free energy F = E - T S, where E is the internal
    energy at the smeared occupations. The occupations of the converged solution give S, hence E
    and the zero-temperature extrapolation E_0 = F + T S / 2 = E - T S / 2, which is exact to
    second order in T for Fermi-Dirac smearing of a smooth density of states (Gillan 1989).
    """
    res = ase_calculator._res
    occ = np.atleast_2d(np.asarray(res.get("orbital-occupations"), dtype=float))
    s_k = fermi_entropy(occ)
    f = np.clip(occ / 2.0 if occ.shape[0] == 1 else occ, 0.0, 1.0)
    ts = kB * electronic_temperature_K * s_k
    free = float(ase_calculator.results["energy"])
    return {
        "F_eV": free,
        "TS_eV": float(ts),
        "E_internal_eV": free + float(ts),
        "E0_eV": free + 0.5 * float(ts),
        "S_over_kB": s_k,
        "n_fractional_orbitals": int(np.sum((f > 0.01) & (f < 0.99))),
    }


__all__ = ["TBLiteCalculator", "fermi_entropy", "smearing_terms"]
