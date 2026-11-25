"""GFN-xTB (tblite) calculator adapter: fidelity tier T2."""

from __future__ import annotations

from dataclasses import dataclass, field

from ase import Atoms

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

    def settings(self) -> dict:
        return {
            "method": self.method,
            "accuracy": self.accuracy,
            "electronic_temperature": self.electronic_temperature,
            "max_iterations": self.max_iterations,
            "spin_polarized": self.spin_polarized,
            **self.extra,
        }


__all__ = ["TBLiteCalculator"]
