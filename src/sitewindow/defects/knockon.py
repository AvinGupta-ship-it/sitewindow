"""Electron knock-on displacement cross sections (DERIVATIONS 6.1-6.2).

The McKinley-Feshbach expression gives the total displacement cross section for a static
target of mass M and charge Z with displacement threshold E_d. Thermal smearing integrates
over the Gaussian out-of-plane velocity distribution of the target nucleus so that the
effective cross section is finite below the static threshold.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad

from sitewindow.constants import ALPHA_FS, AMU_KG, K_B_J, M_C12_AMU, M_E_C2_EV, R_E_M

BARN_M2 = 1.0e-28
_E_D_GRID = np.linspace(18.0, 25.0, 29)
_T_EFF_GRID = np.linspace(300.0, 1500.0, 25)


def beta_squared(E_keV: ArrayLike) -> NDArray[np.float64]:
    """(v/c)^2 for an electron of kinetic energy E (keV)."""
    E = np.asarray(E_keV, dtype=float) * 1.0e3
    gamma = 1.0 + E / M_E_C2_EV
    return 1.0 - 1.0 / gamma**2


def t_max_eV(E_keV: ArrayLike, mass_amu: float = M_C12_AMU) -> NDArray[np.float64]:
    """Maximum energy transfer to the nucleus in a head-on collision (eV)."""
    E = np.asarray(E_keV, dtype=float) * 1.0e3
    M_c2 = mass_amu * AMU_KG * (2.99792458e8) ** 2 / 1.602176634e-19  # eV
    return 2.0 * E * (E + 2.0 * M_E_C2_EV) / M_c2


def mckinley_feshbach_barn(
    E_keV: ArrayLike, E_d_eV: float, Z: int = 6, mass_amu: float = M_C12_AMU
) -> NDArray[np.float64]:
    """Total displacement cross section (barn) for a static target (DERIVATIONS 6.1)."""
    E = np.asarray(E_keV, dtype=float)
    b2 = beta_squared(E)
    b = np.sqrt(b2)
    tmax = t_max_eV(E, mass_amu)
    alpha = Z * ALPHA_FS
    x = tmax / E_d_eV
    pref = np.pi * Z**2 * R_E_M**2 * (1.0 - b2) / b2**2
    with np.errstate(invalid="ignore", divide="ignore"):
        bracket = (
            (x - 1.0) - b2 * np.log(x) + np.pi * alpha * b * (2.0 * (np.sqrt(x) - 1.0) - np.log(x))
        )
    sigma = np.where(x >= 1.0, pref * bracket, 0.0)
    return sigma / BARN_M2


def mean_square_velocity_from_temperature(T_eff_K: float, mass_amu: float = M_C12_AMU) -> float:
    """<v_z^2> (m^2 s^-2) of a classical one-dimensional oscillator at effective temperature T_eff."""
    return K_B_J * T_eff_K / (mass_amu * AMU_KG)


def debye_mean_square_velocity(
    theta_D_K: float = 1287.0, T_K: float = 300.0, mass_amu: float = M_C12_AMU
) -> float:
    """One-component mean-square velocity <v_z^2> (m^2 s^-2) in the Debye model.

    <v_z^2> = (1/M) int_0^{omega_D} g(omega) (hbar omega / 2) coth(hbar omega / 2 k_B T) d omega with
    g(omega) = 3 omega^2 / omega_D^3 normalized to one mode per component. At T = 0 this gives
    3 k_B theta_D / (8 M). The out-of-plane Debye temperature of graphene, 1287 K for 12C and
    sqrt(12/13) times that for 13C, is the value used by Meyer et al. (2012).
    """
    from sitewindow.constants import HBAR

    M = mass_amu * AMU_KG
    omega_D = K_B_J * theta_D_K / HBAR

    def integrand(w: float) -> float:
        x = HBAR * w / (2.0 * K_B_J * T_K) if T_K > 0 else np.inf
        coth = 1.0 if not np.isfinite(x) else 1.0 / np.tanh(x)
        return 3.0 * w * w / omega_D**3 * (HBAR * w / 2.0) * coth

    val, _ = quad(integrand, 0.0, omega_D, limit=200)
    return float(val / M)


def thermally_smeared_cross_section_barn(
    E_keV: float,
    E_d_eV: float,
    v2_mean: float,
    Z: int = 6,
    mass_amu: float = M_C12_AMU,
    n_sigma: float = 6.0,
) -> float:
    """Effective cross section with Gaussian out-of-plane velocity smearing (DERIVATIONS 6.2).

    A nucleus moving toward the beam with velocity v_z gains an additional energy transfer in
    the laboratory frame; the maximum transfer becomes T' = (sqrt(T_max) + sqrt(E_v))^2 with
    E_v = M v_z^2 / 2 for motion toward the incoming electron (and (sqrt(T_max) - sqrt(E_v))^2
    away from it). The static cross section is evaluated with the effective threshold
    E_d * T_max / T' and averaged over the Gaussian distribution of v_z.
    """
    tmax = float(t_max_eV(E_keV, mass_amu))
    M = mass_amu * AMU_KG
    s = np.sqrt(v2_mean)

    def integrand(v: float) -> float:
        E_v = 0.5 * M * v * v / 1.602176634e-19
        sign = 1.0 if v >= 0 else -1.0
        t_eff = (np.sqrt(tmax) + sign * np.sqrt(E_v)) ** 2
        if t_eff <= 0:
            return 0.0
        E_d_eff = E_d_eV * tmax / t_eff
        sig = float(mckinley_feshbach_barn(E_keV, E_d_eff, Z, mass_amu))
        pdf = np.exp(-0.5 * (v / s) ** 2) / (s * np.sqrt(2.0 * np.pi))
        return sig * pdf

    val, _ = quad(integrand, -n_sigma * s, n_sigma * s, limit=200)
    return float(val)


def fit_threshold_and_velocity(
    E_keV: ArrayLike,
    sigma_barn: ArrayLike,
    sigma_err_barn: ArrayLike,
    E_d_grid: ArrayLike = _E_D_GRID,
    T_eff_grid: ArrayLike = _T_EFF_GRID,
    mass_amu: float = M_C12_AMU,
) -> dict[str, float]:
    """Weighted least-squares grid fit of (E_d, T_eff) to measured cross sections (benchmark B9)."""
    E = np.asarray(E_keV, dtype=float)
    y = np.asarray(sigma_barn, dtype=float)
    err = np.asarray(sigma_err_barn, dtype=float)
    best = {"chi2": np.inf, "E_d_eV": np.nan, "T_eff_K": np.nan}
    for Ed in np.asarray(E_d_grid):
        for Te in np.asarray(T_eff_grid):
            v2 = mean_square_velocity_from_temperature(float(Te), mass_amu)
            pred = np.array(
                [
                    thermally_smeared_cross_section_barn(float(e), float(Ed), v2, mass_amu=mass_amu)
                    for e in E
                ]
            )
            chi2 = float(np.sum(((pred - y) / err) ** 2))
            if chi2 < best["chi2"]:
                best = {"chi2": chi2, "E_d_eV": float(Ed), "T_eff_K": float(Te)}
    best["chi2_per_dof"] = best["chi2"] / max(len(E) - 2, 1)
    return best


__all__ = [
    "beta_squared",
    "t_max_eV",
    "mckinley_feshbach_barn",
    "thermally_smeared_cross_section_barn",
    "mean_square_velocity_from_temperature",
    "debye_mean_square_velocity",
    "fit_threshold_and_velocity",
]
