"""Mean-field rate equations for irradiation-induced defect populations (DERIVATIONS 7).

State vector (per carbon atom): x_V (monovacancies), x_I (adatoms/interstitials),
x_DV (divacancies), x_SW (Stone-Wales). Cross sections in barn are converted to per-carbon
production rates with the fluence rate in electrons per nm^2 per second.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from sitewindow.constants import K_B_EV

BARN_TO_NM2 = 1.0e-10  # 1 barn = 1e-28 m^2 = 1e-10 nm^2


@dataclass(frozen=True)
class RateParameters:
    """Thermally activated process parameters (eV) and prefactors (s^-1)."""

    E_mig_V: float = 1.3
    nu_mig_V: float = 1.0e13
    E_rec_VI: float = 0.0
    nu_rec_VI: float = 1.0e13
    E_loss_I: float = 0.4
    nu_loss_I: float = 1.0e13
    E_anneal_SW: float = 5.0
    nu_anneal_SW: float = 1.0e13
    dv_coalescence_capture_sites: float = (
        6.0  # number of neighbouring sites capturing a migrating vacancy
    )


def _k(E: float, nu: float, T: float) -> float:
    return nu * np.exp(-E / (K_B_EV * T))


def integrate_population(
    fluence_rate_e_nm2_s: float,
    t_end_s: float,
    T_K: float,
    sigma_V_barn: float,
    sigma_DV_barn: float = 0.0,
    sigma_SW_barn: float = 0.0,
    params: RateParameters | None = None,
    n_points: int = 200,
    rtol: float = 1e-9,
    atol: float = 1e-14,
) -> dict[str, NDArray[np.float64]]:
    """Integrate DERIVATIONS 7.1 and return trajectories and the carbon-conservation residual."""
    params = params or RateParameters()
    phi = fluence_rate_e_nm2_s
    sV = sigma_V_barn * BARN_TO_NM2
    sDV = sigma_DV_barn * BARN_TO_NM2
    sSW = sigma_SW_barn * BARN_TO_NM2
    k_mig = _k(params.E_mig_V, params.nu_mig_V, T_K)
    k_rec = _k(params.E_rec_VI, params.nu_rec_VI, T_K)
    k_loss = _k(params.E_loss_I, params.nu_loss_I, T_K)
    k_sw = _k(params.E_anneal_SW, params.nu_anneal_SW, T_K)
    z = params.dv_coalescence_capture_sites

    def rhs(t: float, y: NDArray[np.float64]) -> list[float]:
        xV, xI, xDV, xSW, removed = y
        x_tot = xV + 2.0 * xDV + xSW
        prod = sV * phi * max(1.0 - x_tot, 0.0)
        coalesce = z * k_mig * xV * xV
        recomb = k_rec * xV * xI
        dV = prod - recomb - 2.0 * coalesce - sDV * phi * xV
        dI = prod - recomb - k_loss * xI
        dDV = coalesce + sDV * phi * xV
        dSW = sSW * phi * max(1.0 - x_tot, 0.0) - k_sw * xSW
        # carbons removed from the lattice: production minus recombination; divacancy by beam
        # adds one more removal per event
        dRemoved = prod - recomb + sDV * phi * xV
        return [dV, dI, dDV, dSW, dRemoved]

    t_eval = np.linspace(0.0, t_end_s, n_points)
    sol = solve_ivp(
        rhs,
        (0.0, t_end_s),
        [0.0, 0.0, 0.0, 0.0, 0.0],
        t_eval=t_eval,
        method="BDF",
        rtol=rtol,
        atol=atol,
    )
    xV, xI, xDV, xSW, removed = sol.y
    # Conservation: removed carbons must equal xV + 2 xDV (coalescence conserves removals)
    residual = removed - (xV + 2.0 * xDV)
    return {
        "t_s": sol.t,
        "fluence_e_nm2": sol.t * phi,
        "x_V": xV,
        "x_I": xI,
        "x_DV": xDV,
        "x_SW": xSW,
        "x_total": xV + 2.0 * xDV + xSW,
        "conservation_residual": residual,
    }


__all__ = ["RateParameters", "integrate_population", "BARN_TO_NM2"]
