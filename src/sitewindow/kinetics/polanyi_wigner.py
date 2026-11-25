"""Polanyi-Wigner release kinetics: isothermal solutions, linear ramps, Redhead relations.

Implements DERIVATIONS 4.1-4.6. Barriers in eV, prefactors in s^-1, temperatures in K.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from sitewindow.constants import K_B_EV


def rate_constant(E_d_eV: ArrayLike, nu_s: ArrayLike, T: float) -> NDArray[np.float64]:
    """k = nu exp(-E_d / k_B T), evaluated in the log domain."""
    E = np.asarray(E_d_eV, dtype=float)
    nu = np.asarray(nu_s, dtype=float)
    return np.exp(np.log(nu) - E / (K_B_EV * T))


def release_fraction_isothermal(
    E_d_eV: ArrayLike,
    nu_s: ArrayLike,
    order: ArrayLike,
    T: float,
    tau_s: float,
    theta0: ArrayLike = 1.0,
) -> NDArray[np.float64]:
    """Fraction released after time tau at temperature T (DERIVATIONS 4.2)."""
    k = rate_constant(E_d_eV, nu_s, T)
    n = np.asarray(order)
    th0 = np.asarray(theta0, dtype=float)
    kt = k * tau_s
    f1 = -np.expm1(-kt)
    f2 = kt * th0 / (1.0 + kt * th0)
    return np.where(n == 2, f2, f1)


def t90_first_order(E_d_eV: ArrayLike, nu_s: ArrayLike, tau_s: float) -> NDArray[np.float64]:
    """Temperature at which 90 % is released within tau for first-order kinetics (4.6)."""
    E = np.asarray(E_d_eV, dtype=float)
    nu = np.asarray(nu_s, dtype=float)
    return E / (K_B_EV * np.log(nu * tau_s / np.log(10.0)))


def t_release(
    E_d_eV: float,
    nu_s: float,
    order: int,
    tau_s: float,
    fraction: float = 0.9,
    theta0: float = 1.0,
    T_bounds: tuple[float, float] = (10.0, 3000.0),
) -> float:
    """Temperature at which ``fraction`` is released within tau (any order), by root finding."""
    g = lambda T: (
        float(release_fraction_isothermal(E_d_eV, nu_s, order, T, tau_s, theta0)) - fraction
    )  # noqa: E731
    lo, hi = T_bounds
    if g(lo) > 0:
        return lo
    if g(hi) < 0:
        return float("inf")
    return float(brentq(g, lo, hi, xtol=1e-6))


def tpd_linear_ramp(
    E_d_eV: float,
    nu_s: float,
    order: int,
    beta_K_s: float,
    theta0: float = 1.0,
    T0: float = 100.0,
    T1: float = 1200.0,
    n_points: int = 1101,
    rtol: float = 1e-8,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Coverage and desorption rate (per K) along a linear ramp (DERIVATIONS 4.3).

    Returns ``(T, theta, rate_per_K)``.
    """
    T_grid = np.linspace(T0, T1, n_points)

    def rhs(T: float, y: NDArray[np.float64]) -> list[float]:
        th = max(y[0], 0.0)
        k = float(rate_constant(E_d_eV, nu_s, T))
        return [-(k / beta_K_s) * th**order]

    sol = solve_ivp(rhs, (T0, T1), [theta0], t_eval=T_grid, method="BDF", rtol=rtol, atol=atol)
    theta = np.clip(sol.y[0], 0.0, None)
    k = rate_constant(E_d_eV, nu_s, T_grid)
    rate = (k / beta_K_s) * theta**order
    return T_grid, theta, rate


def redhead_peak_temperature(
    E_d_eV: float, nu_s: float, beta_K_s: float, order: int = 1, theta0: float = 1.0
) -> float:
    """Peak temperature from the Redhead relations (DERIVATIONS 4.4).

    First order: E/(k T_p^2) = (nu/beta) exp(-E/k T_p).
    Second order: E/(k T_p^2) = (2 nu theta_p/beta) exp(-E/k T_p) with theta_p ~ theta0/2,
    the standard approximation for the coverage at the second-order peak.
    """
    E = E_d_eV

    def g(T: float) -> float:
        pre = nu_s / beta_K_s if order == 1 else 2.0 * nu_s * (theta0 / 2.0) / beta_K_s
        return E / (K_B_EV * T * T) - pre * np.exp(-E / (K_B_EV * T))

    return float(brentq(g, 20.0, 5000.0, xtol=1e-8))


__all__ = [
    "rate_constant",
    "release_fraction_isothermal",
    "t90_first_order",
    "t_release",
    "tpd_linear_ramp",
    "redhead_peak_temperature",
]
