"""Hydrogen gas properties for the corpus inference (DERIVATIONS 9.2 and 12.1).

The excess-to-absolute conversion needs the molar density of the bulk gas at the isotherm
conditions, the admissibility rules need the saturation pressure below the critical
temperature, and the Langmuir kernel needs the fugacity above 50 bar. CoolProp (the Leachman
reference equation of state for normal hydrogen) is used when installed; otherwise the density
follows the second-virial equation of state solved for the density, the saturation pressure a
Clausius-Clapeyron interpolation between the normal boiling point and the critical point, and
the fugacity the correlation of the theory module. The fallbacks are labelled in the returned
provenance string so that tables record which route produced them.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from sitewindow.constants import R_GAS
from sitewindow.theory.window import _fugacity_coefficient, second_virial_lj

# Normal hydrogen (Leachman et al. 2009): critical point and normal boiling point.
T_CRIT_H2_K = 33.145
P_CRIT_H2_BAR = 12.964
T_NBP_H2_K = 20.369
P_NBP_H2_BAR = 1.01325


def _coolprop():
    try:  # pragma: no cover - optional dependency
        import CoolProp.CoolProp as cp

        return cp
    except Exception:
        return None


def eos_provenance() -> str:
    """Name of the equation-of-state route in use."""
    return "CoolProp Leachman normal hydrogen" if _coolprop() is not None else "second virial"


def h2_density_mol_m3(T: float, P_bar: ArrayLike) -> NDArray[np.float64]:
    """Molar density of bulk H2 (mol m^-3) at temperature ``T`` and pressures ``P_bar``.

    With CoolProp the Leachman equation of state is evaluated point by point; otherwise the
    second-virial equation P = rho R T (1 + B rho) is solved for rho with the Lennard-Jones
    estimate of B(T) from the theory module. Non-positive pressures return zero density.
    """
    P = np.atleast_1d(np.asarray(P_bar, dtype=float))
    out = np.zeros_like(P)
    pos = P > 0.0
    cp = _coolprop()
    if cp is not None:
        for i in np.flatnonzero(pos):
            out[i] = cp.PropsSI("Dmolar", "T", float(T), "P", float(P[i]) * 1.0e5, "Hydrogen")
        return out
    B = second_virial_lj(float(T)) * 1.0e-6  # m^3 mol^-1
    rt = R_GAS * float(T)
    x = 4.0 * B * P[pos] * 1.0e5 / rt
    if abs(B) < 1e-12:
        out[pos] = P[pos] * 1.0e5 / rt
    else:
        out[pos] = (-1.0 + np.sqrt(np.maximum(1.0 + x, 0.0))) / (2.0 * B)
    return out


def gas_in_pores_mol_kg(T: float, P_bar: ArrayLike, pore_volume_cm3_g: float) -> NDArray[np.float64]:
    """Bulk gas held in a pore volume ``pore_volume_cm3_g`` (cm^3 g^-1), in mol kg^-1 (9.2)."""
    rho = h2_density_mol_m3(T, P_bar)  # mol m^-3
    return rho * pore_volume_cm3_g * 1.0e-3


def h2_saturation_pressure_bar(T: float) -> float:
    """Vapor pressure of normal hydrogen (bar) below the critical temperature; inf above it."""
    if T >= T_CRIT_H2_K:
        return float("inf")
    cp = _coolprop()
    if cp is not None and T > 13.957:
        try:  # pragma: no cover - optional dependency
            return float(cp.PropsSI("P", "T", float(T), "Q", 0.0, "Hydrogen")) / 1.0e5
        except Exception:
            pass
    slope = np.log(P_CRIT_H2_BAR / P_NBP_H2_BAR) / (1.0 / T_CRIT_H2_K - 1.0 / T_NBP_H2_K)
    return float(P_NBP_H2_BAR * np.exp(slope * (1.0 / T - 1.0 / T_NBP_H2_K)))


def effective_pressures_bar(T: float, P_bar: ArrayLike) -> NDArray[np.float64]:
    """Fugacity (bar) above 50 bar, pressure below (DERIVATIONS 12.1), element-wise."""
    P = np.atleast_1d(np.asarray(P_bar, dtype=float))
    f = P.copy()
    for i in np.flatnonzero(P > 50.0):
        f[i] = P[i] * _fugacity_coefficient(float(T), float(P[i]))
    return f


__all__ = [
    "T_CRIT_H2_K",
    "eos_provenance",
    "h2_density_mol_m3",
    "gas_in_pores_mol_kg",
    "h2_saturation_pressure_bar",
    "effective_pressures_bar",
]
