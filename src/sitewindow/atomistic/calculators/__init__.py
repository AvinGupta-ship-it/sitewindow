"""Calculator adapters for the fidelity tiers and a factory keyed by calculator name."""

from __future__ import annotations

from typing import Any

from sitewindow.atomistic.calculators.base import (
    Calculator,
    CalculatorSpec,
    CalculatorUnavailable,
    available_calculators,
)


def make_calculator(name: str, **settings: Any) -> Calculator:
    """Instantiate a registered calculator (``tblite``, ``gpaw``, ``mace`` or a plugin) by name."""
    registry = available_calculators()
    if name not in registry:
        raise CalculatorUnavailable(f"unknown calculator {name!r}; available: {sorted(registry)}")
    calc: Calculator = registry[name](**settings)
    return calc


__all__ = [
    "Calculator",
    "CalculatorSpec",
    "CalculatorUnavailable",
    "available_calculators",
    "make_calculator",
]
