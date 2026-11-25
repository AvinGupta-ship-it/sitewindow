"""Entropy-of-adsorption models (DERIVATIONS 2.3).

Each model returns the standard molar entropy change on adsorption, ``delta_S`` in
J mol^-1 K^-1, at the standard-state pressure P0 = 1 bar, as a function of temperature. The
pressure dependence of the gas chemical potential enters the Langmuir isotherm through the
product K P and must not be included here. The constant model reproduces the Bhatia-Myers
baseline (-8 R). The harmonic-oscillator and two-dimensional-gas models are derived from the
partition functions of the adsorbed state and the Sackur-Tetrode plus rigid-rotor entropy of
gaseous H2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from sitewindow.constants import AMU_KG, H_PLANCK, K_B_J, M_H2_AMU, N_A, R_GAS

# Rotational constant of H2 (cm^-1) and its energy scale in K.
_B_ROT_CM1 = 60.853
_THETA_ROT_K = 1.438776877 * _B_ROT_CM1  # hc B / k_B


class EntropyModel(Protocol):
    """Protocol for entropy-of-adsorption models."""

    name: str

    def delta_s(self, T: float, P_bar: float = 1.0) -> float:
        """Standard molar entropy change S_ads - S_gas(T, P0) in J mol^-1 K^-1; P_bar is ignored."""
        ...


def gas_translational_entropy(T: float, P_bar: float, mass_amu: float = M_H2_AMU) -> float:
    """Sackur-Tetrode translational entropy of an ideal gas, J mol^-1 K^-1."""
    m = mass_amu * AMU_KG
    P = P_bar * 1.0e5
    lam = H_PLANCK / np.sqrt(2.0 * np.pi * m * K_B_J * T)
    n_over_v = P / (K_B_J * T)
    return R_GAS * (np.log(1.0 / (n_over_v * lam**3)) + 2.5)


def gas_rotational_entropy_h2(T: float) -> float:
    """Rotational entropy of H2 in J mol^-1 K^-1 from the explicit level sum.

    Ortho/para populations follow the high-temperature equilibrium (nuclear-spin weights
    1 for even J, 3 for odd J), which is the appropriate reference for a gas equilibrated at
    the window temperature; at 77 K the sum is still evaluated explicitly.
    """
    J = np.arange(0, 40)
    g_ns = np.where(J % 2 == 0, 1.0, 3.0)
    e = J * (J + 1) * _THETA_ROT_K / T
    w = g_ns * (2 * J + 1) * np.exp(-e)
    q = w.sum()
    u_over_kT = (w * e).sum() / q
    return R_GAS * (np.log(q) + u_over_kT)


@dataclass(frozen=True)
class ConstantEntropy:
    """Constant entropy change, baseline -8 R (Bhatia and Myers)."""

    multiple_of_r: float = -8.0
    name: str = "constant"

    def delta_s(self, T: float, P_bar: float = 1.0) -> float:
        return self.multiple_of_r * R_GAS


@dataclass(frozen=True)
class HarmonicOscillatorEntropy:
    """Localized adsorbate with three frustrated-translation frequencies (cm^-1)."""

    nu_cm1: tuple[float, float, float] = (60.0, 60.0, 100.0)
    keep_rotation: bool = True
    name: str = "harmonic_oscillator"

    def delta_s(self, T: float, P_bar: float = 1.0) -> float:
        x = np.array(self.nu_cm1) * 1.438776877 / T
        s_ads = R_GAS * np.sum(x / np.expm1(x) - np.log1p(-np.exp(-x)))
        s_gas = gas_translational_entropy(T, 1.0) + gas_rotational_entropy_h2(T)
        s_rot_ads = gas_rotational_entropy_h2(T) if self.keep_rotation else 0.0
        return float(s_ads + s_rot_ads - s_gas)


@dataclass(frozen=True)
class TwoDGasEntropy:
    """Mobile two-dimensional ideal gas with area per molecule ``area_nm2``."""

    area_nm2: float = 0.10
    name: str = "two_d_gas"

    def delta_s(self, T: float, P_bar: float = 1.0) -> float:
        m = M_H2_AMU * AMU_KG
        lam = H_PLANCK / np.sqrt(2.0 * np.pi * m * K_B_J * T)
        a = self.area_nm2 * 1.0e-18
        s_ads = R_GAS * (np.log(a / lam**2) + 2.0) + gas_rotational_entropy_h2(T)
        s_gas = gas_translational_entropy(T, 1.0) + gas_rotational_entropy_h2(T)
        return float(s_ads - s_gas)


DEFAULT_MODEL_SET: tuple[EntropyModel, ...] = (
    ConstantEntropy(),
    HarmonicOscillatorEntropy(),
    TwoDGasEntropy(),
)

__all__ = [
    "EntropyModel",
    "ConstantEntropy",
    "HarmonicOscillatorEntropy",
    "TwoDGasEntropy",
    "DEFAULT_MODEL_SET",
    "gas_translational_entropy",
    "gas_rotational_entropy_h2",
    "N_A",
]
