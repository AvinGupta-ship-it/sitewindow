"""Bayesian adsorption-energy-distribution (AED) inversion of H2 isotherms (manual C6).

Model (DERIVATIONS 9.1-9.5). For a material with isotherm points k at (T_k, P_k) and reported
loading y_k:

    n_abs,k = n_m sum_j w_j theta(eps_j, T_k, P_k),   theta = sigmoid(eps_j/(R T_k) + dS/R + ln f_k),
    n_exc,k = n_abs,k - rho_gas(T_k, P_k) V_pore,
    y_k ~ N(n_exc,k, s_k^2),   s_k^2 = (a n_exc,k)^2 + b^2 + (dn/dlnP_k * dlnP_k)^2,

with eps_j > 0 the binding energy on the grid [2, 40] kJ mol^-1 (0.5 kJ mol^-1 spacing),
dS = -8 R + delta_S the entropy of adsorption (the C1 baseline plus a per-material offset that
is free when the material has isotherms at two or more temperatures and fixed at zero
otherwise, label "single_T_entropy_fixed"), f_k the fugacity above 50 bar, rho_gas the bulk gas
density from the hydrogen equation of state and V_pore the pore volume of the
excess-to-absolute conversion. The log-weights u_j carry a second-difference Gaussian
random-walk prior with scale lambda ~ HalfNormal(lambda_0), u_{j+1} - 2 u_j + u_{j-1} ~ N(0,
lambda^2), together with a weak N(0, 4^2) prior on each u_j that fixes the two free linear
degrees of freedom of the random walk; the weights are the softmax of u (the Dirichlet-like
normalization of DERIVATIONS 9.4). The centered form samples without divergences where the
non-centered cumulative-sum form (kept as an option) produces them. The
measurement-error scales (a, b) distinguish tabulated from digitized data; dlnP_k is the
pressure-axis digitization uncertainty, the larger of one percent of the isotherm's pressure
span relative to P_k and two percent of P_k, capped at 25 percent (a point far below the axis
span was read from a logarithmic axis or an inset), for digitized data and 0.5 percent for
tabulated data, propagated through the local slope of the model isotherm. The error scales of digitized
data are estimated in two stages: stage 1 fits every material with the scales free (log-normal
priors around the manual's values), stage 2 pools the posterior medians across the corpus by a
precision-weighted mean of their logarithms, leaving out the material's own source (DOI), and
refits each material with the pooled scales fixed (empirical Bayes, leave-one-source-out). The
tabulated scales keep the manual's values, since the corpus holds a single tabulated isotherm.

Outputs per material: posterior mean and 90 % band of the weights, the mean site energy, the
fraction of site capacity above 10 kJ mol^-1, the window-occupancy fraction for the standard
windows (the fraction of weight at energies whose Langmuir swing is at least half the optimum
swing of the window under the material's entropy offset), posterior predictive checks
(chi-square p-value, PIT coverage), MCMC diagnostics and an outlier classification by metadata.
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from sitewindow.constants import EV_TO_KJ_PER_MOL, R_GAS
from sitewindow.errors import SiteWindowError
from sitewindow.experimental.corpus import IsothermData, MaterialSet
from sitewindow.experimental.eos import (
    effective_pressures_bar,
    eos_provenance,
    h2_density_mol_m3,
)
from sitewindow.theory.entropy import ConstantEntropy
from sitewindow.theory.window import W1, W2, W3, OperatingWindow, swing

log = logging.getLogger(__name__)

GRID_MIN_KJ = 2.0
GRID_MAX_KJ = 40.0
GRID_STEP_KJ = 0.5
DS_BASELINE_MULTIPLE_OF_R = -8.0
ENTROPY_OFFSET_PRIOR_SD_J = 2.0 * R_GAS
LAMBDA0_DEFAULT = 1.0
LOGWEIGHT_PRIOR_SD = 4.0
CENTERED_DEFAULT = True
PORE_VOLUME_DEFAULT_CM3_G = 1.0
PORE_VOLUME_BROAD_SD = 1.0
PORE_VOLUME_REPORTED_SD = 0.25
STRONG_SITE_THRESHOLD_KJ = 10.0
STRONG_FRACTION_THRESHOLD = 0.10
STANDARD_WINDOWS: tuple[OperatingWindow, ...] = (W1, W2, W3)
N_DRAWS_SAVED = 200


class InversionUnavailable(SiteWindowError):
    """PyMC (optional dependency group ``bayes``) is not installed."""


@dataclass(frozen=True)
class ErrorScales:
    """Measurement-error scales of DERIVATIONS 9.3 by source type."""

    a_tab: float = 0.02
    b_tab: float = 0.02
    a_dig: float = 0.05
    b_dig: float = 0.05
    dig_pressure_rel: float = 0.02
    dig_pressure_axis: float = 0.01
    dig_pressure_cap: float = 0.25
    tab_pressure_rel: float = 0.005
    tau_a: float = 0.5
    tau_b: float = 0.5
    source: str = "manual defaults"


def energy_grid(step: float = GRID_STEP_KJ) -> np.ndarray:
    """Binding-energy grid in kJ mol^-1 (positive for binding)."""
    n = int(round((GRID_MAX_KJ - GRID_MIN_KJ) / step)) + 1
    return GRID_MIN_KJ + step * np.arange(n)


def log_kernel(eps_kJ: np.ndarray, T: np.ndarray, f_bar: np.ndarray, delta_s_J: float = 0.0) -> np.ndarray:
    """ln(K f) for every (point, grid energy): eps/(R T) + dS/R + ln f (DERIVATIONS 2.1, 9.1)."""
    dS = DS_BASELINE_MULTIPLE_OF_R * R_GAS + delta_s_J
    return (
        eps_kJ[None, :] * 1000.0 / (R_GAS * T[:, None]) + dS / R_GAS + np.log(f_bar)[:, None]
    )


def occupancy_matrix(eps_kJ, T, f_bar, delta_s_J=0.0) -> np.ndarray:
    x = log_kernel(np.asarray(eps_kJ, float), np.asarray(T, float), np.asarray(f_bar, float), delta_s_J)
    return 1.0 / (1.0 + np.exp(-x))


@dataclass
class MaterialData:
    """Point-wise arrays of one material, concatenated over its isotherms."""

    material_id: str
    T: np.ndarray
    P: np.ndarray
    f: np.ndarray
    rho: np.ndarray
    y: np.ndarray
    is_digitized: np.ndarray
    dlnP: np.ndarray
    isotherm_index: np.ndarray
    temperatures: list[float]
    pore_volume_prior: tuple[float, float]
    n_m_prior_median: float

    @property
    def multi_temperature(self) -> bool:
        return len(self.temperatures) >= 2


def material_data(
    ms: MaterialSet,
    scales: ErrorScales | None = None,
    pore_volume_cm3_g: float | None = None,
) -> MaterialData:
    """Assemble the point arrays, gas densities, fugacities and digitization terms."""
    scales = scales or ErrorScales()
    T, P, f, rho, y, dig, dlnP, idx = [], [], [], [], [], [], [], []
    for i, iso in enumerate(ms.isotherms):
        n = iso.n_points
        T.append(np.full(n, iso.temperature_K))
        P.append(iso.pressure_bar)
        f.append(effective_pressures_bar(iso.temperature_K, iso.pressure_bar))
        rho.append(h2_density_mol_m3(iso.temperature_K, iso.pressure_bar))
        y.append(iso.loading_mol_kg)
        dig.append(np.full(n, not iso.tabular))
        if iso.tabular:
            dlnP.append(np.full(n, scales.tab_pressure_rel))
        else:
            span = float(iso.pressure_bar.max())
            dlnP.append(
                np.minimum(
                    np.maximum(scales.dig_pressure_axis * span / iso.pressure_bar, scales.dig_pressure_rel),
                    scales.dig_pressure_cap,
                )
            )
        idx.append(np.full(n, i))
    yv = np.concatenate(y)
    if pore_volume_cm3_g is not None and pore_volume_cm3_g > 0:
        vprior = (float(pore_volume_cm3_g), PORE_VOLUME_REPORTED_SD)
    else:
        vprior = (PORE_VOLUME_DEFAULT_CM3_G, PORE_VOLUME_BROAD_SD)
    return MaterialData(
        material_id=ms.material_id,
        T=np.concatenate(T),
        P=np.concatenate(P),
        f=np.concatenate(f),
        rho=np.concatenate(rho),
        y=yv,
        is_digitized=np.concatenate(dig).astype(float),
        dlnP=np.concatenate(dlnP),
        isotherm_index=np.concatenate(idx),
        temperatures=ms.temperatures,
        pore_volume_prior=vprior,
        n_m_prior_median=max(2.0 * float(np.nanmax(yv)), 0.5),
    )


def window_mask(
    grid_kJ: np.ndarray, window: OperatingWindow, delta_s_J: float = 0.0, fraction: float = 0.5
) -> np.ndarray:
    """Grid energies inside the thermodynamic window: swing at least ``fraction`` of its maximum."""
    model = ConstantEntropy(DS_BASELINE_MULTIPLE_OF_R + delta_s_J / R_GAS)
    ds = swing(-np.asarray(grid_kJ, float) / EV_TO_KJ_PER_MOL, window, model)
    return ds >= fraction * float(np.max(ds))


def _pymc():
    try:
        import pymc as pm
        import pytensor.tensor as pt

        return pm, pt
    except Exception as e:  # pragma: no cover - optional dependency
        raise InversionUnavailable("pymc is required for the AED inversion") from e


def build_model(
    data: MaterialData,
    grid_kJ: np.ndarray,
    scales: ErrorScales,
    *,
    free_scales: bool,
    lambda0: float = LAMBDA0_DEFAULT,
    entropy_free: bool | None = None,
    centered: bool = CENTERED_DEFAULT,
):
    """The PyMC model of one material (see the module docstring)."""
    pm, pt = _pymc()
    J = len(grid_kJ)
    entropy_free = data.multi_temperature if entropy_free is None else entropy_free
    inv_rt = grid_kJ[None, :] * 1000.0 / (R_GAS * data.T[:, None])
    log_f = np.log(data.f)
    v0, v_sd = data.pore_volume_prior
    with pm.Model() as model:
        lam = pm.HalfNormal("lam", sigma=lambda0)
        if centered:
            u = pm.Normal("z", 0.0, LOGWEIGHT_PRIOR_SD, shape=J)
            pm.Potential("smoothness", pm.logp(pm.Normal.dist(0.0, lam), u[2:] - 2 * u[1:-1] + u[:-2]).sum())
        else:
            z = pm.Normal("z", 0.0, 1.0, shape=J - 2)
            c1 = pm.Normal("c1", 0.0, 1.0)
            d1 = c1 + pt.concatenate([pt.zeros(1), pt.cumsum(lam * z)])
            u = pt.concatenate([pt.zeros(1), pt.cumsum(d1)])
        w = pm.Deterministic("w", pt.special.softmax(u))
        log_nm = pm.Normal("log_nm", mu=math.log(data.n_m_prior_median), sigma=1.0)
        n_m = pm.Deterministic("n_m", pt.exp(log_nm))
        log_v = pm.Normal("log_V", mu=math.log(v0), sigma=v_sd)
        v_pore = pm.Deterministic("V_pore", pt.exp(log_v))
        if entropy_free:
            delta_s = pm.Normal("delta_s", 0.0, ENTROPY_OFFSET_PRIOR_SD_J)
        else:
            delta_s = pm.Deterministic("delta_s", pt.constant(0.0))
        x = inv_rt + (DS_BASELINE_MULTIPLE_OF_R * R_GAS + delta_s) / R_GAS + log_f[:, None]
        theta = pm.math.sigmoid(x)
        n_abs = n_m * pt.sum(w[None, :] * theta, axis=1)
        mu = pm.Deterministic("mu", n_abs - data.rho * v_pore * 1.0e-3)
        slope = n_m * pt.sum(w[None, :] * theta * (1.0 - theta), axis=1)
        if free_scales:
            # Stage 1: the manual's values with a log-sd of 0.5; stage 2: the corpus estimate
            # with its between-material spread (partial pooling, DERIVATIONS 9.3).
            a_dig = pm.LogNormal("a_dig", mu=math.log(scales.a_dig), sigma=scales.tau_a)
            b_dig = pm.LogNormal("b_dig", mu=math.log(scales.b_dig), sigma=scales.tau_b)
        else:
            a_dig, b_dig = scales.a_dig, scales.b_dig
        a = data.is_digitized * a_dig + (1.0 - data.is_digitized) * scales.a_tab
        b = data.is_digitized * b_dig + (1.0 - data.is_digitized) * scales.b_tab
        s2 = (a * mu) ** 2 + b**2 + (slope * data.dlnP) ** 2
        sigma = pm.Deterministic("sigma", pt.sqrt(s2))
        pm.Normal("y", mu=mu, sigma=sigma, observed=data.y)
    return model


@dataclass
class AEDResult:
    material_id: str
    label: str
    grid_kJ: list[float]
    w_mean: list[float]
    w_q05: list[float]
    w_q95: list[float]
    mean_energy: dict[str, float]
    strong_fraction: dict[str, float]
    window_occupancy: dict[str, dict[str, float]]
    n_m: dict[str, float]
    V_pore: dict[str, float]
    delta_s: dict[str, float]
    error_scales: dict[str, Any]
    diagnostics: dict[str, Any]
    ppc: dict[str, float]
    n_points: int
    temperatures_K: list[float]
    metadata: dict[str, Any] = field(default_factory=dict)
    draws: np.ndarray | None = None

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("draws", None)
        return d


def _summ(x: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "q05": float(np.quantile(x, 0.05)),
        "q50": float(np.quantile(x, 0.50)),
        "q95": float(np.quantile(x, 0.95)),
    }


def _occupancies(w: np.ndarray, grid: np.ndarray, delta_s: np.ndarray) -> dict[str, dict[str, float]]:
    """Window-occupancy fractions per draw; the swing band is evaluated per unique entropy offset
    rounded to 0.1 J mol^-1 K^-1 so that the cost stays bounded."""
    out = {}
    keys = np.round(delta_s, 1)
    for win in STANDARD_WINDOWS:
        masks = {k: window_mask(grid, win, float(k)) for k in np.unique(keys)}
        occ = np.array([float(np.sum(w[i][masks[keys[i]]])) for i in range(len(w))])
        out[win.name] = _summ(occ)
    return out


def summarize_posterior(
    idata,
    data: MaterialData,
    grid: np.ndarray,
    label: str,
    scales: ErrorScales,
    free_scales: bool,
    entropy_free: bool,
    centered: bool = CENTERED_DEFAULT,
) -> AEDResult:
    import arviz as az

    post = idata.posterior
    w = post["w"].stack(s=("chain", "draw")).transpose("s", ...).values
    nm = post["n_m"].stack(s=("chain", "draw")).values
    vp = post["V_pore"].stack(s=("chain", "draw")).values
    ds = post["delta_s"].stack(s=("chain", "draw")).values
    mu = post["mu"].stack(s=("chain", "draw")).transpose("s", ...).values
    sig = post["sigma"].stack(s=("chain", "draw")).transpose("s", ...).values
    mean_e = w @ grid
    strong = w[:, grid > STRONG_SITE_THRESHOLD_KJ].sum(axis=1)
    var_names = ["lam", "log_nm", "log_V"] + ([] if centered else ["c1"])
    var_names += ["delta_s"] if entropy_free else []
    if free_scales:
        var_names += ["a_dig", "b_dig"]
    rhat = az.rhat(idata, var_names=var_names + ["z"])
    ess = az.ess(idata, var_names=var_names + ["z"], method="bulk")
    rhat_max = float(max(float(rhat[v].max()) for v in rhat.data_vars))
    ess_min = float(min(float(ess[v].min()) for v in ess.data_vars))
    div = int(idata.sample_stats["diverging"].sum()) if "diverging" in idata.sample_stats else 0
    # Diagnostics of the reported quantities themselves: the log-weights of saturated or
    # unoccupied grid energies follow their prior and mix slowly without affecting them.
    import xarray as xr

    w_da = post["w"]
    grid_da = xr.DataArray(grid, dims=w_da.dims[-1:])
    derived = xr.Dataset(
        {
            "mean_energy": (w_da * grid_da).sum(w_da.dims[-1]),
            "strong_fraction": w_da.where(grid_da > STRONG_SITE_THRESHOLD_KJ, 0.0).sum(w_da.dims[-1]),
            "n_m": post["n_m"],
            "V_pore": post["V_pore"],
            "lam": post["lam"],
        }
    )
    if entropy_free:
        derived["delta_s"] = post["delta_s"]
    rhat_d = az.rhat(derived)
    ess_d = az.ess(derived, method="bulk")
    rhat_derived = float(max(float(rhat_d[v].max()) for v in rhat_d.data_vars))
    ess_derived = float(min(float(ess_d[v].min()) for v in ess_d.data_vars))
    # Posterior predictive checks.
    y = data.y
    y_rep = idata.posterior_predictive["y"].stack(s=("chain", "draw")).transpose("s", ...).values
    t_obs = np.sum(((y[None, :] - mu) / sig) ** 2, axis=1)
    t_rep = np.sum(((y_rep - mu) / sig) ** 2, axis=1)
    p_value = float(np.mean(t_rep >= t_obs))
    from scipy.stats import norm

    pit = norm.cdf((y[None, :] - mu) / sig).mean(axis=0)
    cover90 = float(np.mean(np.abs(pit - 0.5) <= 0.45))
    err = {"a_tab": scales.a_tab, "b_tab": scales.b_tab, "source": scales.source}
    if free_scales:
        a = post["a_dig"].stack(s=("chain", "draw")).values
        b = post["b_dig"].stack(s=("chain", "draw")).values
        err.update({"a_dig": _summ(a), "b_dig": _summ(b), "estimated": True})
    else:
        err.update({"a_dig": scales.a_dig, "b_dig": scales.b_dig, "estimated": False})
    rng = np.random.default_rng(0)
    keep = rng.choice(len(w), size=min(N_DRAWS_SAVED, len(w)), replace=False)
    return AEDResult(
        material_id=data.material_id,
        label=label,
        grid_kJ=grid.tolist(),
        w_mean=w.mean(axis=0).tolist(),
        w_q05=np.quantile(w, 0.05, axis=0).tolist(),
        w_q95=np.quantile(w, 0.95, axis=0).tolist(),
        mean_energy=_summ(mean_e),
        strong_fraction=_summ(strong),
        window_occupancy=_occupancies(w, grid, ds),
        n_m=_summ(nm),
        V_pore=_summ(vp),
        delta_s=_summ(ds),
        error_scales=err,
        diagnostics={
            "rhat_max": rhat_max,
            "ess_min": ess_min,
            "rhat_derived_max": rhat_derived,
            "ess_derived_min": ess_derived,
            "divergences": div,
            "n_draws": int(len(w)),
        },
        ppc={"chi2_p_value": p_value, "pit_coverage_90": cover90},
        n_points=int(len(y)),
        temperatures_K=data.temperatures,
        draws=w[keep].astype(np.float32),
    )


def fit_material(
    ms: MaterialSet,
    *,
    scales: ErrorScales | None = None,
    free_scales: bool = False,
    pore_volume_cm3_g: float | None = None,
    grid_step: float = GRID_STEP_KJ,
    lambda0: float = LAMBDA0_DEFAULT,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    cores: int = 2,
    seed: int = 0,
    target_accept: float = 0.9,
    entropy_free: bool | None = None,
    centered: bool = CENTERED_DEFAULT,
) -> AEDResult:
    """Fit the AED model of one material and summarize the posterior."""
    pm, _ = _pymc()
    scales = scales or ErrorScales()
    data = material_data(ms, scales, pore_volume_cm3_g)
    grid = energy_grid(grid_step)
    entropy_free = data.multi_temperature if entropy_free is None else entropy_free
    label = "multi_T" if entropy_free else "single_T_entropy_fixed"
    model = build_model(
        data, grid, scales, free_scales=free_scales, lambda0=lambda0, entropy_free=entropy_free,
        centered=centered,
    )
    import warnings

    with model, warnings.catch_warnings():
        # The smoothness Potential acts on the prior of the log-weights only; the posterior
        # predictive draws of y use the sampled mu and sigma, so the PyMC warning does not apply.
        warnings.filterwarnings("ignore", message="The effect of Potentials")
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=cores,
            random_seed=seed,
            target_accept=target_accept,
            progressbar=False,
            compute_convergence_checks=False,
        )
        pm.sample_posterior_predictive(
            idata, var_names=["y"], progressbar=False, random_seed=seed, extend_inferencedata=True
        )
    res = summarize_posterior(idata, data, grid, label, scales, free_scales, entropy_free, centered)
    res.metadata = {
        "doi": ms.doi,
        "adsorbent_name": ms.adsorbent_name,
        "material_class": ms.material_class,
        "n_isotherms": len(ms.isotherms),
        "digitized_only": all(not i.tabular for i in ms.isotherms),
        "generic_name": ms.generic_name,
        "quality_flags": sorted({f for i in ms.isotherms for f in i.quality_flags.split(",") if f}),
        "max_pressure_bar": float(max(i.pressure_bar.max() for i in ms.isotherms)),
        "max_loading_mol_kg": float(max(i.loading_mol_kg.max() for i in ms.isotherms)),
        "pore_volume_reported_cm3_g": pore_volume_cm3_g,
        "lambda0": lambda0,
        "grid_step_kJ": grid_step,
        "eos": eos_provenance(),
        "sampler": {"draws": draws, "tune": tune, "chains": chains, "target_accept": target_accept},
    }
    return res


# ----------------------------------------------------------------------------------------------
# Synthetic data and recovery tests (acceptance criterion F2)


def synthetic_material(
    weights: np.ndarray,
    grid_kJ: np.ndarray,
    *,
    n_m: float = 10.0,
    pore_volume_cm3_g: float = 0.8,
    delta_s_J: float = 0.0,
    temperatures: tuple[float, ...] = (77.0, 87.0),
    pressures_bar: np.ndarray | None = None,
    scales: ErrorScales | None = None,
    digitized: bool = True,
    seed: int = 0,
    name: str = "synthetic",
) -> MaterialSet:
    """Isotherms generated from a planted distribution under the measurement-error model."""
    scales = scales or ErrorScales()
    rng = np.random.default_rng(seed)
    P = np.logspace(-2, np.log10(50.0), 20) if pressures_bar is None else np.asarray(pressures_bar, float)
    w = np.asarray(weights, float) / np.sum(weights)
    isos = []
    for k, T in enumerate(temperatures):
        f = effective_pressures_bar(T, P)
        theta = occupancy_matrix(grid_kJ, np.full(len(P), T), f, delta_s_J)
        n_abs = n_m * theta @ w
        n_exc = n_abs - h2_density_mol_m3(T, P) * pore_volume_cm3_g * 1.0e-3
        slope = n_m * (theta * (1 - theta)) @ w
        if digitized:
            dlnP = np.minimum(
                np.maximum(scales.dig_pressure_axis * P.max() / P, scales.dig_pressure_rel),
                scales.dig_pressure_cap,
            )
            s = np.sqrt((scales.a_dig * n_exc) ** 2 + scales.b_dig**2 + (slope * dlnP) ** 2)
        else:
            s = np.sqrt((scales.a_tab * n_exc) ** 2 + scales.b_tab**2 + (slope * scales.tab_pressure_rel) ** 2)
        y = n_exc + s * rng.standard_normal(len(P))
        isos.append(
            IsothermData(
                filename=f"{name}.Isotherm{k}",
                temperature_K=float(T),
                pressure_bar=P.copy(),
                loading_mol_kg=np.maximum(y, 0.0),
                digitized=digitized,
                tabular=not digitized,
            )
        )
    return MaterialSet(
        material_id=f"synthetic__{name}",
        doi="synthetic",
        adsorbent_hashkey="synthetic",
        adsorbent_name=name,
        material_class="synthetic",
        isotherms=isos,
    )


def planted_weights(grid_kJ: np.ndarray, modes: list[tuple[float, float]]) -> np.ndarray:
    """Weights of planted modes (energy, fraction) placed on the nearest grid point each."""
    w = np.zeros(len(grid_kJ))
    for e, frac in modes:
        w[int(np.argmin(np.abs(grid_kJ - e)))] += frac
    return w / w.sum()


def mode_energies(draws: np.ndarray, grid_kJ: np.ndarray, splits: list[float]) -> np.ndarray:
    """Weight-mean energy of each mode region separated by ``splits`` for every draw."""
    edges = [-np.inf, *splits, np.inf]
    out = np.zeros((len(draws), len(edges) - 1))
    for m in range(len(edges) - 1):
        sel = (grid_kJ > edges[m]) & (grid_kJ <= edges[m + 1])
        mass = draws[:, sel].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[:, m] = (draws[:, sel] @ grid_kJ[sel]) / mass
    return out


def resolvable(eps_kJ: float, temperatures, pressures_bar, delta_s_J: float = 0.0) -> bool:
    """Whether a measurement design resolves a site energy: some measured point has an
    occupancy between 0.05 and 0.95 at that energy, so that the data see the Langmuir step of
    the site rather than a constant (saturated) or vanishing contribution."""
    P = np.asarray(pressures_bar, float)
    for T in temperatures:
        th = occupancy_matrix(np.array([eps_kJ]), np.full(len(P), float(T)), effective_pressures_bar(T, P), delta_s_J)[:, 0]
        if np.any((th > 0.05) & (th < 0.95)):
            return True
    return False


def recovery_test(
    kind: str = "single",
    *,
    replicates: int = 5,
    tolerance_kJ: float = 1.0,
    level: float = 0.90,
    draws: int = 800,
    tune: int = 800,
    chains: int = 2,
    cores: int = 2,
    temperatures: tuple[float, ...] = (77.0, 87.0),
    entropy_free: bool | None = None,
    seed: int = 0,
    lambda0: float = LAMBDA0_DEFAULT,
    grid_step: float = GRID_STEP_KJ,
) -> dict[str, Any]:
    """Synthetic recovery of planted single-site or bimodal distributions.

    A replicate passes when, for every planted mode, the posterior probability that the mode's
    weight-mean energy lies within ``tolerance_kJ`` of the planted energy is at least ``level``
    (the planted energy within 1 kJ mol^-1 at 90 % posterior, acceptance criterion F2).
    """
    grid = energy_grid(grid_step)
    if kind == "single":
        modes = [(6.0, 1.0)]
        splits: list[float] = []
    elif kind == "bimodal":
        modes = [(5.0, 0.7), (12.0, 0.3)]
        splits = [8.5]
    else:
        raise ValueError(kind)
    w0 = planted_weights(grid, modes)
    pressures = np.logspace(-2, np.log10(50.0), 20)
    resolvable_modes = [resolvable(e, temperatures, pressures) for e, _ in modes]
    results = []
    for r in range(replicates):
        ms = synthetic_material(w0, grid, temperatures=temperatures, seed=seed + r, name=f"{kind}_{r}")
        res = fit_material(
            ms, scales=ErrorScales(), free_scales=False, draws=draws, tune=tune, chains=chains,
            cores=cores, seed=seed + r, entropy_free=entropy_free, lambda0=lambda0,
            grid_step=grid_step,
        )
        assert res.draws is not None
        energies = mode_energies(res.draws, grid, splits)
        probs, passes = [], []
        for m, (e, _) in enumerate(modes):
            p = float(np.mean(np.abs(energies[:, m] - e) <= tolerance_kJ))
            probs.append(p)
            passes.append(p >= level)
        fracs = [float(np.mean(res.draws[:, grid > splits[0]].sum(axis=1))) if splits else 1.0]
        results.append(
            {
                "replicate": r,
                "planted": modes,
                "mode_energy_posterior_mean": energies.mean(axis=0).tolist(),
                "p_within_tolerance": probs,
                "strong_mode_weight_posterior_mean": fracs,
                "passed": bool(all(passes)),
                "passed_resolvable": bool(all(p for p, ok in zip(passes, resolvable_modes, strict=True) if ok)),
                "rhat_max": res.diagnostics["rhat_max"],
                "divergences": res.diagnostics["divergences"],
                "ppc_p_value": res.ppc["chi2_p_value"],
            }
        )
    return {
        "kind": kind,
        "temperatures_K": list(temperatures),
        "entropy": "free" if (len(temperatures) >= 2 if entropy_free is None else entropy_free) else "fixed",
        "lambda0": lambda0,
        "grid_step_kJ": grid_step,
        "tolerance_kJ_mol": tolerance_kJ,
        "level": level,
        "resolvable_modes": resolvable_modes,
        "replicates": results,
        "pass_fraction": float(np.mean([r["passed"] for r in results])),
        "passed": bool(all(r["passed"] for r in results)),
        "pass_fraction_resolvable": float(np.mean([r["passed_resolvable"] for r in results])),
    }


# ----------------------------------------------------------------------------------------------
# Corpus run with checkpoints


def pooled_scales(stage1: list[dict[str, Any]], exclude_doi: str | None = None) -> ErrorScales:
    """Precision-weighted mean of the log posterior medians of (a_dig, b_dig) across materials."""
    rows = [r for r in stage1 if r["metadata"].get("doi") != exclude_doi and r["error_scales"].get("estimated")]
    if not rows:
        return ErrorScales(source="manual defaults (no stage-1 estimates)")
    vals, taus = {}, {}
    for key in ("a_dig", "b_dig"):
        med = np.array([r["error_scales"][key]["q50"] for r in rows])
        q05 = np.array([r["error_scales"][key]["q05"] for r in rows])
        q95 = np.array([r["error_scales"][key]["q95"] for r in rows])
        sd_log = np.maximum((np.log(q95) - np.log(q05)) / (2 * 1.645), 0.05)
        wgt = 1.0 / sd_log**2
        vals[key] = float(np.exp(np.sum(wgt * np.log(med)) / np.sum(wgt)))
        # Between-material spread of the log scale: the variance of the medians minus the
        # mean within-material posterior variance (method of moments), floored at 0.1.
        between = np.var(np.log(med), ddof=1) - np.mean(sd_log**2) if len(med) > 1 else 0.0
        taus[key] = float(np.sqrt(max(between, 0.01)))
    return ErrorScales(
        a_dig=vals["a_dig"],
        b_dig=vals["b_dig"],
        tau_a=taus["a_dig"],
        tau_b=taus["b_dig"],
        source=f"pooled over {len(rows)} stage-1 fits"
        + (f" excluding {exclude_doi}" if exclude_doi else ""),
    )


PORE_VOLUME_PHYSICAL_MAX_CM3_G = 3.0


def degenerate_fit(res: dict[str, Any]) -> bool:
    """A fit that explains the excess branch with an unphysical pore volume or places the mean
    site energy at the lower grid edge: the independent-Langmuir mixture does not describe the
    isotherm (typically cryogenic high-pressure data past the excess maximum)."""
    grid = res["grid_kJ"]
    step = grid[1] - grid[0] if len(grid) > 1 else 0.5
    return bool(
        res["V_pore"]["q50"] > PORE_VOLUME_PHYSICAL_MAX_CM3_G
        or res["mean_energy"]["mean"] <= grid[0] + step
    )


def classify_outlier(res: dict[str, Any]) -> list[str]:
    """Metadata features of a material whose strong-site fraction exceeds the H5 statement."""
    flags = []
    if degenerate_fit(res):
        flags.append("excess_branch_misfit")
    if res["strong_fraction"]["q50"] <= STRONG_FRACTION_THRESHOLD:
        return flags
    md = res["metadata"]
    if res["label"] == "single_T_entropy_fixed":
        flags.append("single_T_entropy_fixed")
    if md.get("digitized_only"):
        flags.append("digitized_only")
    if md.get("generic_name"):
        flags.append("generic_adsorbent_name")
    if md.get("max_pressure_bar", 0.0) < 1.5 and min(res["temperatures_K"]) > 100.0:
        flags.append("low_coverage_only")
    if res["ppc"]["chi2_p_value"] < 0.05 or res["ppc"]["chi2_p_value"] > 0.95:
        flags.append("posterior_predictive_misfit")
    if md.get("quality_flags"):
        flags.append("quality_flags:" + "+".join(md["quality_flags"]))
    if md.get("material_class") == "carbon_other":
        flags.append("unclassified_carbon")
    return flags or ["unexplained"]


def _save(res: AEDResult, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res.as_dict(), indent=1, default=float))
    if res.draws is not None:
        np.savez_compressed(path.with_suffix(".npz"), w=res.draws, grid_kJ=np.asarray(res.grid_kJ))


def run_corpus(
    sets: list[MaterialSet],
    out_dir: Path,
    *,
    pore_volumes: dict[str, float] | None = None,
    stage: str = "both",
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    cores: int = 2,
    stage1_draws: int = 500,
    lambda0: float = LAMBDA0_DEFAULT,
    grid_step: float = GRID_STEP_KJ,
    seed: int = 0,
    only: set[str] | None = None,
    partial_pooling: bool = True,
) -> dict[str, Any]:
    """Two-stage inversion of every material set with per-material checkpoints.

    Stage 1 (``out_dir/aed_stage1``) fits each material with free digitized error scales under
    the manual's prior; stage 2 (``out_dir/aed``) refits each material with its scales drawn
    from the hierarchical prior estimated over the other sources (centre and between-material
    spread by the method of moments; partial pooling), or, with ``partial_pooling=False``, with
    the pooled scales fixed. Finished materials are skipped on rerun. Returns the assembled
    report (also written as ``aed_report.json``).
    """
    out_dir = Path(out_dir)
    pore_volumes = pore_volumes or {}
    todo = [s for s in sets if only is None or s.material_id in only]
    stage1_dir, stage2_dir = out_dir / "aed_stage1", out_dir / "aed"
    stage1: list[dict[str, Any]] = []
    if stage in ("1", "both"):
        for ms in todo:
            p = stage1_dir / f"{ms.material_id}.json"
            if p.exists():
                stage1.append(json.loads(p.read_text()))
                continue
            log.info("stage 1 %s (%d isotherms, %d temperatures)", ms.material_id, len(ms.isotherms), ms.n_temperatures)
            res = fit_material(
                ms, scales=ErrorScales(), free_scales=True, pore_volume_cm3_g=pore_volumes.get(ms.material_id),
                grid_step=grid_step, lambda0=lambda0, draws=stage1_draws, tune=stage1_draws, chains=chains,
                cores=cores, seed=seed,
            )
            res.draws = None
            _save(res, p)
            stage1.append(res.as_dict())
    else:
        stage1 = [json.loads(p.read_text()) for p in sorted(stage1_dir.glob("*.json"))]
    results: list[dict[str, Any]] = []
    if stage in ("2", "both"):
        for ms in todo:
            p = stage2_dir / f"{ms.material_id}.json"
            if p.exists():
                results.append(json.loads(p.read_text()))
                continue
            scales = pooled_scales(stage1, exclude_doi=ms.doi)
            log.info(
                "stage 2 %s hierarchical scales a=%.3f (tau %.2f) b=%.3f (tau %.2f)",
                ms.material_id, scales.a_dig, scales.tau_a, scales.b_dig, scales.tau_b,
            )
            res = fit_material(
                ms, scales=scales, free_scales=partial_pooling, pore_volume_cm3_g=pore_volumes.get(ms.material_id),
                grid_step=grid_step, lambda0=lambda0, draws=draws, tune=tune, chains=chains, cores=cores,
                seed=seed,
            )
            _save(res, p)
            results.append(res.as_dict())
    report = assemble_report(results, stage1, out_dir)
    return report


def assemble_report(results: list[dict[str, Any]], stage1: list[dict[str, Any]], out_dir: Path) -> dict[str, Any]:
    """Material table, distribution table and the corpus report from the checkpoint records."""
    from sitewindow.io.tables import write_table

    out_dir = Path(out_dir)
    rows, dist = [], []
    for r in results:
        flags = classify_outlier(r)
        rows.append(
            {
                "material_id": r["material_id"],
                "doi": r["metadata"].get("doi"),
                "adsorbent_name": r["metadata"].get("adsorbent_name"),
                "material_class": r["metadata"].get("material_class"),
                "label": r["label"],
                "n_points": r["n_points"],
                "n_temperatures": len(r["temperatures_K"]),
                "temperatures_K": json.dumps(r["temperatures_K"]),
                "mean_energy_kJ_mol": r["mean_energy"]["mean"],
                "mean_energy_q05": r["mean_energy"]["q05"],
                "mean_energy_q95": r["mean_energy"]["q95"],
                "strong_fraction_q50": r["strong_fraction"]["q50"],
                "strong_fraction_q05": r["strong_fraction"]["q05"],
                "strong_fraction_q95": r["strong_fraction"]["q95"],
                "p_strong_fraction_above_0.10": None,
                **{
                    f"occupancy_{w}_mean": r["window_occupancy"][w]["mean"] for w in r["window_occupancy"]
                },
                **{f"occupancy_{w}_q95": r["window_occupancy"][w]["q95"] for w in r["window_occupancy"]},
                "n_m_q50": r["n_m"]["q50"],
                "V_pore_q50": r["V_pore"]["q50"],
                "V_pore_q05": r["V_pore"]["q05"],
                "V_pore_q95": r["V_pore"]["q95"],
                "delta_s_q50": r["delta_s"]["q50"],
                "delta_s_sd": r["delta_s"]["sd"],
                "a_dig": r["error_scales"]["a_dig"] if not isinstance(r["error_scales"]["a_dig"], dict) else r["error_scales"]["a_dig"]["q50"],
                "b_dig": r["error_scales"]["b_dig"] if not isinstance(r["error_scales"]["b_dig"], dict) else r["error_scales"]["b_dig"]["q50"],
                "rhat_max": r["diagnostics"]["rhat_max"],
                "ess_min": r["diagnostics"]["ess_min"],
                "rhat_derived_max": r["diagnostics"].get("rhat_derived_max"),
                "ess_derived_min": r["diagnostics"].get("ess_derived_min"),
                "divergences": r["diagnostics"]["divergences"],
                "degenerate_fit": degenerate_fit(r),
                "ppc_p_value": r["ppc"]["chi2_p_value"],
                "pit_coverage_90": r["ppc"]["pit_coverage_90"],
                "outlier": bool(r["strong_fraction"]["q50"] > STRONG_FRACTION_THRESHOLD),
                "outlier_features": ";".join(flags),
                "digitized_only": bool(r["metadata"].get("digitized_only")),
                "generic_name": bool(r["metadata"].get("generic_name")),
            }
        )
        for e, m, lo, hi in zip(r["grid_kJ"], r["w_mean"], r["w_q05"], r["w_q95"], strict=True):
            dist.append({"material_id": r["material_id"], "eps_kJ_mol": e, "w_mean": m, "w_q05": lo, "w_q95": hi})
    table = pd.DataFrame(rows)
    dist_table = pd.DataFrame(dist)
    if not table.empty:
        write_table(table, out_dir / "aed_corpus.parquet")
        write_table(dist_table, out_dir / "aed_distributions.parquet")
    n = len(table)
    ppc_ok = int(((table["ppc_p_value"] >= 0.05) & (table["ppc_p_value"] <= 0.95)).sum()) if n else 0
    pooled = pooled_scales(stage1)
    report = {
        "n_materials": n,
        "n_multi_temperature": int((table["label"] == "multi_T").sum()) if n else 0,
        "n_single_temperature_entropy_fixed": int((table["label"] == "single_T_entropy_fixed").sum()) if n else 0,
        "pooled_error_scales": asdict(pooled),
        "stage1_scale_spread": {
            "a_dig_q50_range": [float(min(r["error_scales"]["a_dig"]["q50"] for r in stage1)), float(max(r["error_scales"]["a_dig"]["q50"] for r in stage1))] if stage1 else None,
            "b_dig_q50_range": [float(min(r["error_scales"]["b_dig"]["q50"] for r in stage1)), float(max(r["error_scales"]["b_dig"]["q50"] for r in stage1))] if stage1 else None,
        },
        "ppc_p_in_[0.05,0.95]": {"count": ppc_ok, "fraction": float(ppc_ok / n) if n else None, "criterion_0.80_met": bool(n and ppc_ok / n >= 0.80)},
        "diagnostics": {
            "rhat_max": float(table["rhat_max"].max()) if n else None,
            "n_rhat_above_1.01": int((table["rhat_max"] > 1.01).sum()) if n else 0,
            "ess_min": float(table["ess_min"].min()) if n else None,
            "n_with_divergences": int((table["divergences"] > 0).sum()) if n else 0,
            "rhat_derived_max": float(table["rhat_derived_max"].max()) if n and table["rhat_derived_max"].notna().any() else None,
            "n_rhat_derived_above_1.01": int((table["rhat_derived_max"] > 1.01).sum()) if n else 0,
            "ess_derived_min": float(table["ess_derived_min"].min()) if n and table["ess_derived_min"].notna().any() else None,
            "n_ess_derived_below_400": int((table["ess_derived_min"] < 400).sum()) if n else 0,
        },
        "degenerate_fits": {
            "n": int(table["degenerate_fit"].sum()) if n else 0,
            "materials": table.loc[table["degenerate_fit"], ["material_id", "adsorbent_name", "V_pore_q50", "mean_energy_kJ_mol", "ppc_p_value"]].to_dict(orient="records") if n else [],
        },
        "strong_sites": {
            "threshold_kJ_mol": STRONG_SITE_THRESHOLD_KJ,
            "fraction_threshold": STRONG_FRACTION_THRESHOLD,
            "n_above": int(table["outlier"].sum()) if n else 0,
            "fraction_above": float(table["outlier"].mean()) if n else None,
            "n_above_excluding_degenerate": int((table["outlier"] & ~table["degenerate_fit"]).sum()) if n else 0,
            "n_above_q05": int((table["strong_fraction_q05"] > STRONG_FRACTION_THRESHOLD).sum()) if n else 0,
            "features": table.loc[table["outlier"], ["material_id", "adsorbent_name", "material_class", "label", "strong_fraction_q50", "outlier_features"]].to_dict(orient="records") if n else [],
        },
        "window_occupancy_corpus": {
            w.name: {
                "median_of_means": float(table[f"occupancy_{w.name}_mean"].median()) if n else None,
                "max_q95": float(table[f"occupancy_{w.name}_q95"].max()) if n else None,
            }
            for w in STANDARD_WINDOWS
        },
        "mean_energy_corpus_kJ_mol": {
            "median": float(table["mean_energy_kJ_mol"].median()) if n else None,
            "q25": float(table["mean_energy_kJ_mol"].quantile(0.25)) if n else None,
            "q75": float(table["mean_energy_kJ_mol"].quantile(0.75)) if n else None,
        },
        "eos": eos_provenance(),
    }
    (out_dir / "aed_report.json").write_text(json.dumps(report, indent=2, default=float))
    return report


__all__ = [
    "AEDResult",
    "ErrorScales",
    "InversionUnavailable",
    "MaterialData",
    "STANDARD_WINDOWS",
    "assemble_report",
    "build_model",
    "classify_outlier",
    "degenerate_fit",
    "energy_grid",
    "fit_material",
    "log_kernel",
    "material_data",
    "mode_energies",
    "occupancy_matrix",
    "planted_weights",
    "pooled_scales",
    "recovery_test",
    "resolvable",
    "run_corpus",
    "synthetic_material",
    "window_mask",
]
