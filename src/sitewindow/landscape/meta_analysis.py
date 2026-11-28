"""Hierarchical meta-analysis of published site energies (manual C3).

For a published value y_i of class c(i) with method descriptors m(i):

    y_i = mu_c(i) + beta_func[m_i] + beta_disp[m_i] + gamma / N_cell,i + u_study(i) + e_i,

with mu_c ~ N(mu_0, tau_c^2) (non-centered), the method effects beta ~ N(0, sigma_beta^2) partially
pooled, study effects u ~ N(0, tau_study^2) and e_i ~ N(0, s_i^2) with s_i the stated uncertainty or
the registry default. Inference by NUTS (PyMC), four chains; diagnostics R-hat, bulk ESS,
divergences, posterior predictive coverage and leave-one-study-out cross-validation. The model is
fitted per quantity and reference state, since binding energies are never averaged across
reference states (manual 6.9). Without PyMC the fit raises ``MetaAnalysisUnavailable``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from sitewindow.errors import SiteWindowError


class MetaAnalysisUnavailable(SiteWindowError):
    """PyMC (optional dependency group ``bayes``) is not installed."""


@dataclass
class MetaAnalysisResult:
    quantity: str
    reference_state: str
    n_values: int
    classes: list[str]
    class_means: pd.DataFrame
    method_effects: pd.DataFrame
    heterogeneity: dict[str, float]
    diagnostics: dict[str, Any]
    posterior_predictive: dict[str, Any]
    loso: list[dict[str, Any]] = field(default_factory=list)
    idata: Any = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "quantity": self.quantity,
            "reference_state": self.reference_state,
            "n_values": self.n_values,
            "classes": self.classes,
            "class_means": self.class_means.to_dict(orient="records"),
            "method_effects": self.method_effects.to_dict(orient="records"),
            "heterogeneity": self.heterogeneity,
            "diagnostics": self.diagnostics,
            "posterior_predictive": self.posterior_predictive,
            "leave_one_study_out": self.loso,
        }


def _codes(values: pd.Series) -> tuple[np.ndarray, list[str]]:
    cats = sorted(str(v) for v in values.unique())
    index = {c: i for i, c in enumerate(cats)}
    return np.array([index[str(v)] for v in values]), cats


def _require_pymc() -> tuple[Any, Any]:
    try:
        import arviz as az
        import pymc as pm
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise MetaAnalysisUnavailable(
            "pymc and arviz are required for the meta-analysis (pip install 'sitewindow[bayes]')"
        ) from exc
    return pm, az


def build_model(frame: pd.DataFrame, prior_scale: float = 1.0) -> tuple[Any, dict[str, Any]]:
    """The PyMC model of the registry frame; ``prior_scale`` multiplies every prior width for
    the sensitivity study of the manual (three prior scales)."""
    pm, _ = _require_pymc()
    c_idx, classes = _codes(frame["site_class"])
    f_idx, functionals = _codes(frame["functional"])
    d_idx, dispersions = _codes(frame["dispersion"])
    s_idx, studies = _codes(frame["study"])
    y = frame["value_eV"].to_numpy(dtype=float)
    s = frame["sigma_eV"].to_numpy(dtype=float)
    inv_n = frame["inv_n_cell"].to_numpy(dtype=float)
    coords = {
        "cls": classes,
        "func": functionals,
        "disp": dispersions,
        "study": studies,
        "obs": np.arange(len(y)),
    }
    with pm.Model(coords=coords) as model:
        mu0 = pm.Normal("mu0", mu=float(np.mean(y)), sigma=2.0 * prior_scale)
        tau_c = pm.HalfNormal("tau_c", sigma=1.5 * prior_scale)
        z_c = pm.Normal("z_c", 0.0, 1.0, dims="cls")
        mu_c = pm.Deterministic("mu_c", mu0 + tau_c * z_c, dims="cls")
        sigma_beta = pm.HalfNormal("sigma_beta", sigma=0.3 * prior_scale)
        z_f = pm.Normal("z_f", 0.0, 1.0, dims="func")
        z_d = pm.Normal("z_d", 0.0, 1.0, dims="disp")
        beta_func = pm.Deterministic("beta_func", sigma_beta * z_f, dims="func")
        beta_disp = pm.Deterministic("beta_disp", sigma_beta * z_d, dims="disp")
        gamma = pm.Normal("gamma", 0.0, 10.0 * prior_scale)
        tau_study = pm.HalfNormal("tau_study", sigma=0.2 * prior_scale)
        z_u = pm.Normal("z_u", 0.0, 1.0, dims="study")
        u = pm.Deterministic("u", tau_study * z_u, dims="study")
        mean = mu_c[c_idx] + beta_func[f_idx] + beta_disp[d_idx] + gamma * inv_n + u[s_idx]
        pm.Normal("y", mu=mean, sigma=s, observed=y, dims="obs")
    return model, coords


def fit_meta_analysis(
    frame: pd.DataFrame,
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    seed: int = 0,
    target_accept: float = 0.95,
    prior_scale: float = 1.0,
    loso: bool = False,
    loso_draws: int = 400,
) -> MetaAnalysisResult:
    """Fit the hierarchical model to one quantity and reference state of the registry frame."""
    pm, az = _require_pymc()
    frame = frame.reset_index(drop=True)
    if frame.empty:
        raise ValueError("the meta-analysis frame is empty")
    quantities = frame["quantity"].unique()
    references = frame["reference_state"].unique()
    if len(quantities) != 1 or len(references) != 1:
        raise ValueError("fit one quantity and one reference state at a time (manual 6.9)")
    model, coords = build_model(frame, prior_scale)
    with model:
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            random_seed=seed,
            target_accept=target_accept,
            progressbar=False,
            compute_convergence_checks=False,
        )
        pm.compute_log_likelihood(idata, progressbar=False)
        ppc = pm.sample_posterior_predictive(idata, progressbar=False, random_seed=seed)
    var_names = [
        "mu0",
        "tau_c",
        "mu_c",
        "sigma_beta",
        "beta_func",
        "beta_disp",
        "gamma",
        "tau_study",
        "u",
    ]
    summary = az.summary(idata, var_names=var_names, hdi_prob=0.90)
    rhat_max = float(np.nanmax(summary["r_hat"].to_numpy(dtype=float)))
    ess_min = float(np.nanmin(summary["ess_bulk"].to_numpy(dtype=float)))
    n_div = int(np.asarray(idata.sample_stats["diverging"]).sum())
    post = idata.posterior
    mu_c = post["mu_c"].stack(sample=("chain", "draw")).to_numpy()
    class_means = pd.DataFrame(
        {
            "site_class": coords["cls"],
            "mean_eV": mu_c.mean(axis=1),
            "sd_eV": mu_c.std(axis=1, ddof=1),
            "hdi5_eV": np.quantile(mu_c, 0.05, axis=1),
            "hdi95_eV": np.quantile(mu_c, 0.95, axis=1),
            "n_values": [int((frame["site_class"] == c).sum()) for c in coords["cls"]],
        }
    )
    effects: list[dict[str, Any]] = []
    for name, dim, labels in (
        ("beta_func", "func", coords["func"]),
        ("beta_disp", "disp", coords["disp"]),
    ):
        arr = post[name].stack(sample=("chain", "draw")).to_numpy()
        for i, label in enumerate(labels):
            effects.append(
                {
                    "effect": name,
                    "level": label,
                    "mean_eV": float(arr[i].mean()),
                    "sd_eV": float(arr[i].std(ddof=1)),
                    "hdi5_eV": float(np.quantile(arr[i], 0.05)),
                    "hdi95_eV": float(np.quantile(arr[i], 0.95)),
                    "n_values": int(
                        (frame["functional" if dim == "func" else "dispersion"] == label).sum()
                    ),
                }
            )
    func = post["beta_func"].stack(sample=("chain", "draw")).to_numpy()
    fl = list(coords["func"])
    pairwise: list[dict[str, Any]] = []
    for i, a in enumerate(fl):
        for j, b in enumerate(fl):
            if j <= i:
                continue
            d = func[i] - func[j]
            pairwise.append(
                {
                    "effect": "beta_func_difference",
                    "level": f"{a} - {b}",
                    "mean_eV": float(d.mean()),
                    "sd_eV": float(d.std(ddof=1)),
                    "hdi5_eV": float(np.quantile(d, 0.05)),
                    "hdi95_eV": float(np.quantile(d, 0.95)),
                    "n_values": int(
                        ((frame["functional"] == a) | (frame["functional"] == b)).sum()
                    ),
                }
            )
    gamma = post["gamma"].stack(sample=("chain", "draw")).to_numpy()
    effects.append(
        {
            "effect": "gamma_cell_size",
            "level": "eV per (1/N_cell)",
            "mean_eV": float(gamma.mean()),
            "sd_eV": float(gamma.std(ddof=1)),
            "hdi5_eV": float(np.quantile(gamma, 0.05)),
            "hdi95_eV": float(np.quantile(gamma, 0.95)),
            "n_values": int((frame["inv_n_cell"] > 0).sum()),
        }
    )
    method_effects = pd.DataFrame(effects + pairwise)
    tau_c = post["tau_c"].to_numpy().ravel()
    tau_s = post["tau_study"].to_numpy().ravel()
    sig_b = post["sigma_beta"].to_numpy().ravel()
    heterogeneity = {
        "tau_class_eV": float(tau_c.mean()),
        "tau_class_sd_eV": float(tau_c.std(ddof=1)),
        "tau_study_eV": float(tau_s.mean()),
        "tau_study_sd_eV": float(tau_s.std(ddof=1)),
        "sigma_beta_eV": float(sig_b.mean()),
        "sigma_beta_sd_eV": float(sig_b.std(ddof=1)),
    }
    yrep = ppc.posterior_predictive["y"].stack(sample=("chain", "draw")).to_numpy()
    y = frame["value_eV"].to_numpy(dtype=float)
    lo, hi = np.quantile(yrep, [0.05, 0.95], axis=1)
    posterior_predictive = {
        "coverage_90": float(np.mean((y >= lo) & (y <= hi))),
        "mean_abs_residual_eV": float(np.mean(np.abs(y - yrep.mean(axis=1)))),
    }
    diagnostics = {
        "rhat_max": rhat_max,
        "ess_bulk_min": ess_min,
        "divergences": n_div,
        "draws": draws,
        "tune": tune,
        "chains": chains,
        "target_accept": target_accept,
        "prior_scale": prior_scale,
        "passes_rhat_1.01": bool(rhat_max < 1.01),
        "passes_ess_400": bool(ess_min > 400),
    }
    try:
        loo = az.loo(idata, pointwise=True)
        diagnostics["loo_elpd"] = float(loo.elpd_loo)
        diagnostics["loo_pareto_k_max"] = float(np.max(np.asarray(loo.pareto_k)))
    except Exception as exc:  # pragma: no cover - arviz version differences
        diagnostics["loo_error"] = str(exc).splitlines()[0]
    result = MetaAnalysisResult(
        str(quantities[0]),
        str(references[0]),
        len(frame),
        list(coords["cls"]),
        class_means,
        method_effects,
        heterogeneity,
        diagnostics,
        posterior_predictive,
        idata=idata,
    )
    if loso:
        result.loso = leave_one_study_out(
            frame, draws=loso_draws, chains=min(chains, 2), seed=seed, prior_scale=prior_scale
        )
    return result


def leave_one_study_out(
    frame: pd.DataFrame, draws: int = 400, chains: int = 2, seed: int = 0, prior_scale: float = 1.0
) -> list[dict[str, Any]]:
    """Refit without each study and predict its values from the class means and method effects
    of the remaining studies (the study effect integrates to its prior)."""
    pm, _ = _require_pymc()
    out: list[dict[str, Any]] = []
    for study in sorted(frame["study"].unique()):
        held = frame[frame["study"] == study]
        rest = frame[frame["study"] != study]
        if rest.empty or not set(held["site_class"]) <= set(rest["site_class"]):
            out.append(
                {
                    "study": study,
                    "n_values": int(len(held)),
                    "note": "class absent from the other studies; not predictable",
                }
            )
            continue
        model, coords = build_model(rest, prior_scale)
        with model:
            idata = pm.sample(
                draws=draws,
                tune=draws,
                chains=chains,
                random_seed=seed,
                progressbar=False,
                compute_convergence_checks=False,
                target_accept=0.95,
            )
        post = idata.posterior
        mu = post["mu_c"].stack(sample=("chain", "draw")).to_numpy()
        bf = post["beta_func"].stack(sample=("chain", "draw")).to_numpy()
        bd = post["beta_disp"].stack(sample=("chain", "draw")).to_numpy()
        g = post["gamma"].stack(sample=("chain", "draw")).to_numpy()
        ts = post["tau_study"].stack(sample=("chain", "draw")).to_numpy()
        rng = np.random.default_rng(seed)
        rows = []
        for _, r in held.iterrows():
            ci = coords["cls"].index(str(r["site_class"]))
            pred = mu[ci] + g * float(r["inv_n_cell"]) + rng.normal(0.0, 1.0, len(ts)) * ts
            if str(r["functional"]) in coords["func"]:
                pred = pred + bf[coords["func"].index(str(r["functional"]))]
            if str(r["dispersion"]) in coords["disp"]:
                pred = pred + bd[coords["disp"].index(str(r["dispersion"]))]
            pred = pred + rng.normal(0.0, float(r["sigma_eV"]), len(ts))
            lo, hi = np.quantile(pred, [0.05, 0.95])
            rows.append(
                {
                    "id": r["id"],
                    "value_eV": float(r["value_eV"]),
                    "predicted_eV": float(pred.mean()),
                    "sd_eV": float(pred.std(ddof=1)),
                    "covered_90": bool(lo <= r["value_eV"] <= hi),
                }
            )
        out.append(
            {
                "study": study,
                "n_values": int(len(held)),
                "coverage_90": float(np.mean([x["covered_90"] for x in rows])),
                "values": rows,
            }
        )
    return out


def synthetic_recovery(
    frame: pd.DataFrame,
    n_replicates: int = 100,
    draws: int = 500,
    chains: int = 2,
    seed: int = 0,
    level: float = 0.90,
) -> dict[str, Any]:
    """Acceptance criterion (iii) for the meta-analysis: class means are planted on the design of
    the registry frame (its classes, methods, studies, cell sizes and standard errors), values
    are drawn from the model and the fitted posterior intervals are checked against the planted
    means."""
    pm, _ = _require_pymc()
    rng = np.random.default_rng(seed)
    frame = frame.reset_index(drop=True)
    c_idx, classes = _codes(frame["site_class"])
    f_idx, functionals = _codes(frame["functional"])
    d_idx, dispersions = _codes(frame["dispersion"])
    s_idx, studies = _codes(frame["study"])
    hits: list[bool] = []
    per_class: dict[str, list[bool]] = {c: [] for c in classes}
    for rep in range(n_replicates):
        mu = rng.normal(float(frame["value_eV"].mean()), 0.8, len(classes))
        bf = rng.normal(0.0, 0.15, len(functionals))
        bd = rng.normal(0.0, 0.05, len(dispersions))
        g = rng.normal(0.0, 3.0)
        u = rng.normal(0.0, 0.05, len(studies))
        mean = mu[c_idx] + bf[f_idx] + bd[d_idx] + g * frame["inv_n_cell"].to_numpy() + u[s_idx]
        y = mean + rng.normal(0.0, frame["sigma_eV"].to_numpy())
        synthetic = frame.assign(value_eV=y)
        model, coords = build_model(synthetic)
        with model:
            idata = pm.sample(
                draws=draws,
                tune=draws,
                chains=chains,
                random_seed=seed + rep,
                progressbar=False,
                compute_convergence_checks=False,
                target_accept=0.95,
            )
        post = idata.posterior["mu_c"].stack(sample=("chain", "draw")).to_numpy()
        # the planted class mean is identified up to the pooled method offset; compare the
        # class contrasts to the grand mean, which the model identifies
        planted = mu - mu.mean()
        fitted = post - post.mean(axis=0, keepdims=True)
        lo = np.quantile(fitted, 0.5 - level / 2, axis=1)
        hi = np.quantile(fitted, 0.5 + level / 2, axis=1)
        for i, c in enumerate(coords["cls"]):
            hit = bool(lo[i] <= planted[i] <= hi[i])
            hits.append(hit)
            per_class[c].append(hit)
    frac = float(np.mean(hits))
    return {
        "n_replicates": n_replicates,
        "level": level,
        "coverage_of_planted_class_contrasts": frac,
        "coverage_by_class": {c: float(np.mean(v)) for c, v in per_class.items()},
        "passes_0.90": bool(frac >= 0.90),
        "draws": draws,
        "chains": chains,
        "seed": seed,
        "note": "class means are identified relative to the grand mean because every method effect is partially pooled around zero",
    }


__all__ = [
    "MetaAnalysisResult",
    "MetaAnalysisUnavailable",
    "build_model",
    "fit_meta_analysis",
    "leave_one_study_out",
    "synthetic_recovery",
]
