"""MACE-MP-0 universal interatomic potential adapter: fidelity tier T3 (manual C2).

Results of this tier are used only after per-class validation against Tier 1 in C3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ase import Atoms

from sitewindow.errors import CalculatorUnavailable


@dataclass
class MACECalculator:
    model: str = "medium"
    device: str = "cpu"
    default_dtype: str = "float64"
    dispersion: bool = False
    name: str = "mace"
    tier: str = "T3"
    extra: dict[str, Any] = field(default_factory=dict)

    def attach(self, atoms: Atoms) -> None:
        try:
            from mace.calculators import mace_mp
        except ImportError as e:  # pragma: no cover - optional dependency
            raise CalculatorUnavailable(
                "mace-torch is not installed; install sitewindow[mace]"
            ) from e
        atoms.calc = mace_mp(
            model=self.model,
            device=self.device,
            default_dtype=self.default_dtype,
            dispersion=self.dispersion,
            **self.extra,
        )

    def settings(self) -> dict[str, Any]:
        return {
            "model": self.model,
            "device": self.device,
            "default_dtype": self.default_dtype,
            "dispersion": self.dispersion,
            **self.extra,
        }


__all__ = ["MACECalculator"]
