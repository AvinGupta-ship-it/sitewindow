"""Isosteric heats of adsorption from multi-temperature isotherm sets (manual C6).

q_st(n) = R T^2 (d ln P / dT)_n is evaluated as -R d ln P / d(1/T) at fixed loading: every
isotherm of a material is interpolated as ln P(n) by a shape-preserving cubic (PCHIP) after the
loading has been made strictly increasing, and at each loading of a common grid the slope of
ln P against 1/T is fitted by least squares over the isotherms covering that loading (two
temperatures give the Clausius-Clapeyron pair formula exactly). Uncertainty is nonparametric:
the points of every isotherm are resampled with replacement and the whole evaluation repeated.
The loadings are the reported (excess) amounts; at the low loadings where B4 is evaluated the
excess-to-absolute correction is negligible, and the AED inversion carries the full conversion.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from sitewindow.constants import R_GAS
from sitewindow.experimental.corpus import IsothermData, MaterialSet

MIN_POINTS_FOR_INTERPOLATION = 4
MIN_TEMPERATURE_SPREAD_K = 5.0
LOW_LOADING_QUANTILE = 0.25


@dataclass
class IsostericResult:
    material_id: str
    table: pd.DataFrame
    temperatures_used: list[float]
    reason: str | None = None

    @property
    def usable(self) -> bool:
        return self.reason is None and len(self.table) > 0


def _monotone(P: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """Strictly increasing loading with ln P, averaging duplicates and dropping reversals."""
    order = np.argsort(n, kind="stable")
    n, P = n[order], P[order]
    keep_n, keep_lnp = [], []
    for ni, pi in zip(n, P, strict=True):
        if pi <= 0.0:
            continue
        if keep_n and ni <= keep_n[-1] + 1e-12:
            keep_lnp[-1] = 0.5 * (keep_lnp[-1] + np.log(pi))
            continue
        keep_n.append(float(ni))
        keep_lnp.append(float(np.log(pi)))
    if len(keep_n) < MIN_POINTS_FOR_INTERPOLATION:
        return None
    return np.asarray(keep_n), np.asarray(keep_lnp)


def _interpolant(iso: IsothermData, rng: np.random.Generator | None) -> tuple | None:
    P, n = iso.pressure_bar, iso.loading_mol_kg
    if rng is not None:
        idx = rng.integers(0, len(P), len(P))
        P, n = P[idx], n[idx]
    mono = _monotone(P, n)
    if mono is None:
        return None
    x, y = mono
    return PchipInterpolator(x, y, extrapolate=False), float(x.min()), float(x.max())


def _evaluate(
    isos: list[IsothermData], grid: np.ndarray, rng: np.random.Generator | None
) -> np.ndarray:
    """q_st (kJ mol^-1) on ``grid``; nan where fewer than two temperatures cover the loading."""
    interps = []
    for iso in isos:
        f = _interpolant(iso, rng)
        if f is not None:
            interps.append((iso.temperature_K, *f))
    out = np.full(len(grid), np.nan)
    for k, n0 in enumerate(grid):
        xs, ys = [], []
        for T, f, lo, hi in interps:
            if lo <= n0 <= hi:
                xs.append(1.0 / T)
                ys.append(float(f(n0)))
        if len(set(xs)) < 2:
            continue
        xs_a, ys_a = np.asarray(xs), np.asarray(ys)
        slope = np.polyfit(xs_a, ys_a, 1)[0]
        out[k] = -R_GAS * slope / 1000.0
    return out


def isosteric_heat(
    ms: MaterialSet, n_grid: int = 25, n_boot: int = 200, seed: int = 0
) -> IsostericResult:
    """Isosteric heat versus loading with bootstrap uncertainty for one material set."""
    isos = [i for i in ms.isotherms if i.n_points >= MIN_POINTS_FOR_INTERPOLATION]
    temps = sorted({i.temperature_K for i in isos})
    if len(temps) < 2:
        return IsostericResult(ms.material_id, pd.DataFrame(), temps, "fewer than two temperatures")
    if max(temps) - min(temps) < MIN_TEMPERATURE_SPREAD_K:
        return IsostericResult(ms.material_id, pd.DataFrame(), temps, "temperature spread below 5 K")
    ranges = []
    for iso in isos:
        f = _interpolant(iso, None)
        if f is not None:
            ranges.append((iso.temperature_K, f[1], f[2]))
    if len({r[0] for r in ranges}) < 2:
        return IsostericResult(
            ms.material_id, pd.DataFrame(), temps, "fewer than two interpolable isotherms"
        )
    # The grid spans the union of the loading ranges that two isotherms of different
    # temperature cover jointly; if no such pair exists nothing is estimable.
    lo, hi = np.inf, -np.inf
    for i in range(len(ranges)):
        for j in range(i + 1, len(ranges)):
            if ranges[i][0] == ranges[j][0]:
                continue
            a = max(ranges[i][1], ranges[j][1])
            b = min(ranges[i][2], ranges[j][2])
            if b > a:
                lo, hi = min(lo, a), max(hi, b)
    if not np.isfinite(lo) or hi <= lo:
        return IsostericResult(ms.material_id, pd.DataFrame(), temps, "no overlapping loading range")
    grid = np.linspace(lo, hi, n_grid)
    central = _evaluate(isos, grid, None)
    rng = np.random.default_rng(seed)
    boots = np.array([_evaluate(isos, grid, rng) for _ in range(n_boot)])
    import warnings

    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sd = np.nanstd(boots, axis=0, ddof=1)
        q05 = np.nanquantile(boots, 0.05, axis=0)
        q95 = np.nanquantile(boots, 0.95, axis=0)
        n_eff = np.sum(np.isfinite(boots), axis=0)
    table = pd.DataFrame(
        {
            "material_id": ms.material_id,
            "material_class": ms.material_class,
            "loading_mol_kg": grid,
            "q_st_kJ_mol": central,
            "q_st_sd_kJ_mol": sd,
            "q_st_q05_kJ_mol": q05,
            "q_st_q95_kJ_mol": q95,
            "n_bootstrap_finite": n_eff,
            "n_temperatures": len(temps),
        }
    )
    table = table[np.isfinite(table["q_st_kJ_mol"])].reset_index(drop=True)
    if table.empty:
        return IsostericResult(ms.material_id, table, temps, "no loading covered by two temperatures")
    return IsostericResult(ms.material_id, table, temps)


def low_loading_summary(res: IsostericResult, ms: MaterialSet) -> dict:
    """q_st at low loading: the median over the lowest quarter of the common loading range."""
    t = res.table
    n_lo = t["loading_mol_kg"].min()
    n_hi = t["loading_mol_kg"].max()
    cut = n_lo + LOW_LOADING_QUANTILE * (n_hi - n_lo)
    low = t[t["loading_mol_kg"] <= cut + 1e-12]
    return {
        "material_id": ms.material_id,
        "doi": ms.doi,
        "adsorbent_name": ms.adsorbent_name,
        "material_class": ms.material_class,
        "graphitic": ms.graphitic,
        "temperatures_K": res.temperatures_used,
        "n_temperatures": len(res.temperatures_used),
        "loading_range_mol_kg": [float(n_lo), float(n_hi)],
        "q_st_low_kJ_mol": float(low["q_st_kJ_mol"].median()),
        "q_st_low_sd_kJ_mol": float(np.sqrt(np.nanmean(low["q_st_sd_kJ_mol"] ** 2))),
        "q_st_low_q05_kJ_mol": float(low["q_st_q05_kJ_mol"].min()),
        "q_st_low_q95_kJ_mol": float(low["q_st_q95_kJ_mol"].max()),
        "q_st_mean_kJ_mol": float(t["q_st_kJ_mol"].mean()),
        "q_st_at_max_loading_kJ_mol": float(t["q_st_kJ_mol"].iloc[-1]),
        "digitized": all(i.digitized for i in ms.isotherms),
        "generic_name": ms.generic_name,
    }


def corpus_isosteric_heats(
    sets: list[MaterialSet], n_grid: int = 25, n_boot: int = 200, seed: int = 0
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Isosteric heats for every multi-temperature material: (curves, low-loading summary, log)."""
    curves, summary, log = [], [], []
    for ms in sets:
        if ms.n_temperatures < 2:
            continue
        res = isosteric_heat(ms, n_grid=n_grid, n_boot=n_boot, seed=seed)
        if not res.usable:
            log.append({"material_id": ms.material_id, "adsorbent_name": ms.adsorbent_name,
                        "doi": ms.doi, "temperatures_K": ms.temperatures, "reason": res.reason})
            continue
        curves.append(res.table)
        summary.append(low_loading_summary(res, ms))
    curve_table = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    return curve_table, pd.DataFrame(summary), pd.DataFrame(log)


__all__ = [
    "IsostericResult",
    "corpus_isosteric_heats",
    "isosteric_heat",
    "low_loading_summary",
]
