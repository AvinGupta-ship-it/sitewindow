"""Calculator protocol and registry for the fidelity tiers."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from importlib.metadata import entry_points
from typing import Protocol, runtime_checkable

from ase import Atoms


class CalculatorUnavailable(RuntimeError):
    """Raised when a calculator's dependency is not installed."""


@runtime_checkable
class Calculator(Protocol):
    name: str
    tier: str

    def attach(self, atoms: Atoms) -> None:
        """Attach an ASE calculator to ``atoms`` in place."""
        ...

    def settings(self) -> dict:
        """Serializable settings that determine the result."""
        ...


@dataclass(frozen=True)
class CalculatorSpec:
    name: str
    tier: str
    settings: dict

    @property
    def settings_hash(self) -> str:
        return hashlib.sha256(
            json.dumps({"name": self.name, **self.settings}, sort_keys=True).encode()
        ).hexdigest()[:16]


def available_calculators() -> dict[str, type]:
    """Discover calculators registered under the ``sitewindow.calculators`` entry point."""
    found: dict[str, type] = {}
    for ep in entry_points(group="sitewindow.calculators"):
        try:
            found[ep.name] = ep.load()
        except Exception:  # pragma: no cover - missing optional dependency
            continue
    return found


__all__ = ["Calculator", "CalculatorSpec", "CalculatorUnavailable", "available_calculators"]
