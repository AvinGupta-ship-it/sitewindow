"""Kinetic admissibility of sites under an operating window (DERIVATIONS 4.5)."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from sitewindow.constants import K_B_EV
from sitewindow.kinetics.polanyi_wigner import release_fraction_isothermal


def release_admissible(
    E_d_eV: ArrayLike,
    nu_s: ArrayLike,
    order: ArrayLike,
    T_r: float,
    tau_r: float,
    theta0: ArrayLike = 1.0,
    eta_r: float = 0.9,
) -> tuple[NDArray[np.bool_], NDArray[np.float64]]:
    """Release-side admissibility indicator and release fraction."""
    frac = release_fraction_isothermal(E_d_eV, nu_s, order, T_r, tau_r, theta0)
    return frac >= eta_r, frac


def charge_admissible_atomic(
    E_a_eV: ArrayLike,
    s0: float,
    flux_per_site_s: float,
    T_c: float,
    tau_c: float,
    eta_c: float = 0.9,
) -> tuple[NDArray[np.bool_], NDArray[np.float64]]:
    """Charging admissibility for sticking of atomic hydrogen (DERIVATIONS 3.3)."""
    E = np.asarray(E_a_eV, dtype=float)
    k = s0 * flux_per_site_s * np.exp(-E / (K_B_EV * T_c))
    theta = -np.expm1(-k * tau_c)
    return theta >= eta_c, theta


def smooth_admissibility(
    k_tau: ArrayLike, k_tau_ref: float = np.log(10.0), kappa: float = 4.0
) -> NDArray[np.float64]:
    """Logistic alternative to the indicator, centered where k tau = ln(1/(1-eta)) for eta = 0.9."""
    x = np.log(np.asarray(k_tau, dtype=float) / k_tau_ref)
    return 1.0 / (1.0 + np.exp(-kappa * x))


__all__ = ["release_admissible", "charge_admissible_atomic", "smooth_admissibility"]
