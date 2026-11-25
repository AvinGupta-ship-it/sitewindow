"""Raman defect-density relations (DERIVATIONS 8).

The activation/structural-disorder model relates the D/G intensity ratio of graphene to the
mean inter-defect distance L_D with two length scales r_A and r_S and an excitation-energy
dependent prefactor. Parameter values are configuration data read from
``data/registry/raman_parameters.yaml`` with their source; the defaults below are the values
reported for Ar+-bombarded graphene (Lucchese et al. 2010; Cancado et al. 2011) and are
overridden by the registry at run time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from sitewindow.constants import GRAPHENE_AREAL_DENSITY_NM2


@dataclass(frozen=True)
class RamanParameters:
    """Parameters of the two-length-scale model.

    Defaults follow Cancado et al. (2011): r_A = 3.1 nm, r_S = 1 nm, C_A = (160 +/- 48) E_L^-4
    (E_L in eV) and C_S = 0 (the structural term is defect-type dependent and was set to zero
    for ion-induced point defects). ``C_S = 0.87`` reproduces the earlier Lucchese et al. (2010)
    parameterization at 2.41 eV and is available as an alternative.
    """

    r_A_nm: float = 3.1
    r_S_nm: float = 1.0
    C_A_prefactor: float = 160.0  # C_A = C_A_prefactor * E_L^-4 with E_L in eV
    C_A_prefactor_sigma: float = 48.0
    C_S: float = 0.0


def LD_low_defect_regime_nm(id_ig: ArrayLike, E_L_eV: float = 2.41) -> NDArray[np.float64]:
    """Low-defect-regime relation L_D^2 (nm^2) = (4.3 +/- 1.3) x 10^3 E_L^-4 (I_D/I_G)^-1, E_L in eV.

    This is the large-L_D limit of the full model, I_D/I_G -> C_A pi (r_A^2 - r_S^2)/L_D^2, and is
    stated by the source for L_D >= 10 nm; the exponential corrections of the full model are
    about 16 % at 10 nm and below 5 % beyond 20 nm.
    """
    r = np.asarray(id_ig, dtype=float)
    return np.sqrt(4.3e3 * E_L_eV ** (-4) / r)


DEFAULT_RAMAN_PARAMETERS = RamanParameters()


def id_ig_ratio(
    L_D_nm: ArrayLike, E_L_eV: float = 2.41, params: RamanParameters = DEFAULT_RAMAN_PARAMETERS
) -> NDArray[np.float64]:
    """I_D/I_G as a function of L_D (nm) at excitation energy E_L (eV), DERIVATIONS 8.1."""
    L = np.asarray(L_D_nm, dtype=float)
    rA, rS = params.r_A_nm, params.r_S_nm
    C_A = params.C_A_prefactor * E_L_eV ** (-4)
    term_a = (
        C_A
        * (rA**2 - rS**2)
        / (rA**2 - 2.0 * rS**2)
        * (np.exp(-np.pi * rS**2 / L**2) - np.exp(-np.pi * (rA**2 - rS**2) / L**2))
    )
    term_s = params.C_S * (1.0 - np.exp(-np.pi * rS**2 / L**2))
    return term_a + term_s


def defect_density_from_LD(L_D_nm: ArrayLike) -> NDArray[np.float64]:
    """Areal defect density (nm^-2) from mean inter-defect distance: n_D = 1/(pi L_D^2)."""
    L = np.asarray(L_D_nm, dtype=float)
    return 1.0 / (np.pi * L**2)


def LD_from_defects_per_carbon(x: ArrayLike) -> NDArray[np.float64]:
    """Mean inter-defect distance (nm) for a per-carbon defect density x (DERIVATIONS 7.2)."""
    n_D = np.asarray(x, dtype=float) * GRAPHENE_AREAL_DENSITY_NM2
    return 1.0 / np.sqrt(np.pi * n_D)


def defects_per_carbon_from_LD(L_D_nm: ArrayLike) -> NDArray[np.float64]:
    return defect_density_from_LD(L_D_nm) / GRAPHENE_AREAL_DENSITY_NM2


_ID_IDP_BY_TYPE = {"sp3": 13.0, "vacancy": 7.0, "boundary": 3.5}


def id_idprime_mixture(fractions: dict[str, float]) -> float:
    """First-order linear mixture of I_D/I_D' by defect type (Eckmann et al. 2012 class values)."""
    tot = sum(fractions.values())
    if tot <= 0:
        raise ValueError("fractions must sum to a positive value")
    return sum(_ID_IDP_BY_TYPE[k] * v for k, v in fractions.items()) / tot


__all__ = [
    "RamanParameters",
    "id_ig_ratio",
    "LD_low_defect_regime_nm",
    "defect_density_from_LD",
    "LD_from_defects_per_carbon",
    "defects_per_carbon_from_LD",
    "id_idprime_mixture",
]
