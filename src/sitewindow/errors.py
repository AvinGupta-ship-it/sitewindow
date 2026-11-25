"""Typed exception hierarchy (manual Section 10.5).

Every error raised by ``sitewindow`` derives from :class:`SiteWindowError`, so callers can
distinguish program errors from library errors. Subclasses also derive from the closest
built-in type (``ValueError``, ``RuntimeError``) for compatibility with generic handlers.
"""

from __future__ import annotations


class SiteWindowError(Exception):
    """Base class of all exceptions raised by sitewindow."""


class UnitError(SiteWindowError, ValueError):
    """A unit string is unknown or a conversion needs information that was not supplied."""


class SchemaError(SiteWindowError, ValueError):
    """A record violates the canonical data schema or a physical range check."""


class ConvergenceError(SiteWindowError, RuntimeError):
    """A numerical procedure (SCF, relaxation, quadrature, sampler) did not converge."""


class CalculatorUnavailable(SiteWindowError, RuntimeError):
    """A calculator's optional dependency is not installed."""


class BenchmarkFailure(SiteWindowError):
    """A benchmark comparison failed its pre-registered threshold."""


__all__ = [
    "SiteWindowError",
    "UnitError",
    "SchemaError",
    "ConvergenceError",
    "CalculatorUnavailable",
    "BenchmarkFailure",
]
