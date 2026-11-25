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


__all__ = ["TBLiteCalculator"]
