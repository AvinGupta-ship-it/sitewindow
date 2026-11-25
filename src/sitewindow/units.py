"""Unit parsing and conversion for isotherm and energy data.

Canonical units: pressure in bar (absolute); loading in mol kg^-1 of sorbent; energy in eV.
Conversions cover the unit strings encountered in the NIST/ARPA-E ISODB carbon subset.
Loading units that require a material density (``g/l``, ``g/cm3``) are converted only when a
bulk density is supplied; otherwise :class:`UnitError` is raised so that the caller records
the record as not convertible rather than silently guessing.
"""

from __future__ import annotations

import re

import numpy as np
from numpy.typing import ArrayLike, NDArray

from sitewindow.constants import (
    ATM_TO_BAR,
    EV_TO_KJ_PER_MOL,
    KPA_TO_BAR,
    M_H2_AMU,
    MPA_TO_BAR,
    PSI_TO_BAR,
    TORR_TO_BAR,
    V_MOLAR_STP_CM3,
)


class UnitError(ValueError):
    """Raised when a unit string is unknown or a conversion needs missing information."""


_PRESSURE_FACTORS: dict[str, float] = {
    "bar": 1.0,
    "bars": 1.0,
    "mbar": 1.0e-3,
    "pa": 1.0e-5,
    "kpa": KPA_TO_BAR,
    "mpa": MPA_TO_BAR,
    "atm": ATM_TO_BAR,
    "torr": TORR_TO_BAR,
    "mmhg": TORR_TO_BAR,
    "psi": PSI_TO_BAR,
    "psia": PSI_TO_BAR,
}


def _norm(s: str) -> str:
    return re.sub(r"\s+", "", s.strip().lower())


def pressure_to_bar(values: ArrayLike, unit: str) -> NDArray[np.float64]:
    """Convert pressure values to bar.

    Parameters
    ----------
    values
        Pressure values in ``unit``.
    unit
        One of ``bar``, ``mbar``, ``Pa``, ``kPa``, ``MPa``, ``atm``, ``torr``, ``mmHg``, ``psi``.
    """
    key = _norm(unit)
    if key not in _PRESSURE_FACTORS:
        raise UnitError(f"unknown pressure unit {unit!r}")
    return np.asarray(values, dtype=float) * _PRESSURE_FACTORS[key]


def loading_to_mol_per_kg(
    values: ArrayLike,
    unit: str,
    *,
    molar_mass_g: float = M_H2_AMU,
    bulk_density_g_cm3: float | None = None,
) -> NDArray[np.float64]:
    """Convert adsorbed-amount values to mol per kg of sorbent.

    Supported units: ``mol/kg``, ``mmol/g``, ``mol/g``, ``mmol/kg``, ``wt%``, ``mg/g``,
    ``g/g``, ``g/kg``, ``cm3(STP)/g``, ``cm3/g`` (STP assumed, flagged by caller),
    ``ml(STP)/g``, ``l(STP)/kg``, ``g/l`` and ``g/cm3`` (require ``bulk_density_g_cm3``).

    ``wt%`` is interpreted as 100 * m_ads / (m_sorbent + m_ads), the convention of the
    NIST database for excess gravimetric capacity; the inverse map is used.
    """
    v = np.asarray(values, dtype=float)
    key = _norm(unit).replace("(stp)", "stp").replace("^3", "3")
    if key in {"mol/kg", "mmol/g"}:
        return v
    if key == "mol/g":
        return v * 1000.0
    if key in {"mmol/kg", "umol/g", "μmol/g"}:
        return v / 1000.0
    if key in {"wt%", "wt.%", "wt%(excess)", "%", "mass%"}:
        frac = v / 100.0
        if np.any(frac >= 1.0):
            raise UnitError("wt% value >= 100 is not physical")
        return frac / (1.0 - frac) / molar_mass_g * 1000.0
    if key in {"mg/g", "g/kg"}:
        return v / molar_mass_g
    if key == "g/g":
        return v / molar_mass_g * 1000.0
    if key in {"cm3stp/g", "cm3/g", "mlstp/g", "ml/g", "cc/g", "ccstp/g"}:
        return v / V_MOLAR_STP_CM3 * 1000.0
    if key in {"lstp/kg", "l/kg"}:
        return v / (V_MOLAR_STP_CM3 / 1000.0)
    if key in {"g/l", "g/cm3", "kg/m3"}:
        if bulk_density_g_cm3 is None:
            raise UnitError(f"loading unit {unit!r} requires a bulk density")
        g_per_cm3 = v / 1000.0 if key in {"g/l", "kg/m3"} else v
        return g_per_cm3 / bulk_density_g_cm3 / molar_mass_g * 1000.0
    raise UnitError(f"unknown loading unit {unit!r}")


def mol_per_kg_to_wt_percent(
    values: ArrayLike, molar_mass_g: float = M_H2_AMU
) -> NDArray[np.float64]:
    """Convert mol kg^-1 to wt% defined as 100 * m_ads / (m_sorbent + m_ads)."""
    v = np.asarray(values, dtype=float)
    m_ads_per_kg = v * molar_mass_g / 1000.0  # kg per kg sorbent
    return 100.0 * m_ads_per_kg / (1.0 + m_ads_per_kg)


def ev_to_kj_per_mol(x: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(x, dtype=float) * EV_TO_KJ_PER_MOL


def kj_per_mol_to_ev(x: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(x, dtype=float) / EV_TO_KJ_PER_MOL


def energy_to_ev(values: ArrayLike, unit: str) -> NDArray[np.float64]:
    """Convert energies in ``eV``, ``meV``, ``kJ/mol``, ``kcal/mol``, ``hartree``, ``ry`` to eV."""
    v = np.asarray(values, dtype=float)
    key = _norm(unit)
    factors = {
        "ev": 1.0,
        "mev": 1.0e-3,
        "kj/mol": 1.0 / EV_TO_KJ_PER_MOL,
        "kjmol-1": 1.0 / EV_TO_KJ_PER_MOL,
        "kcal/mol": 4.184 / EV_TO_KJ_PER_MOL,
        "hartree": 27.211386245988,
        "ha": 27.211386245988,
        "ry": 13.605693122994,
    }
    if key not in factors:
        raise UnitError(f"unknown energy unit {unit!r}")
    return v * factors[key]
