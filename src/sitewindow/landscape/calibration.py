"""Multi-fidelity calibration of the Tier-2 and Tier-3 landscapes against the Tier-1 anchors
(manual C3).

For each class c and tier t the transfer E_T1 = a_c,t + b_c,t E_t + δ_t(x) + e maps a tier-t energy
to the Tier-1 scale, where δ is a Gaussian process over the site descriptors x (pyramidalization
angle, coordination, neighbour-hydrogen count, distance to the defect, curvature) with a
Matérn-5/2 kernel. The process is the pooling device across classes: it is fitted to the residuals
of the pooled transfer of the tier over every anchor, so that a class with one or two anchors is
calibrated by the process interpolating its anchors (status ``offset_calibrated``; its slope is the
pooled slope, the manual's fallback of fixing b and estimating an offset), a class without anchors
is predicted by the pooled transfer and the process extrapolating in descriptor space with its
prior variance restored (``awaiting_anchors``), and a class with at least three anchors receives
its own (a, b) with the process fitted to its own residuals (``calibrated``). Predictive standard
deviations are inflated by the Student-t factor of the degrees of freedom behind the plug-in
variance. The Tier-0 offset, the Tier-1 bias against the DMC and CCSD(T) anchors, is estimated
from the available anchor pairs and carried as a class-independent shift with its uncertainty.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy import optimize, stats

DESCRIPTORS: tuple[str, ...] = (
    "pyramidalization_deg",
    "coordination",
    "n_neighbor_H",
    "distance_to_defect_nm",
    "curvature_inv_nm",
)
T1_ENGINE = "gpaw:PBE-D3(BJ)"
T1_SIGMA_EV = 0.02
"""Numerical uncertainty of Tier-1 chemisorption energies (manual C2, Tier-1 convergence)."""
TIER_SIGMA_EV: dict[str, float] = {"T2": 0.03, "T3": 0.05, "T1": T1_SIGMA_EV}
"""Per-configuration numerical uncertainty of each tier's energies (manual C2): the Tier-2
extrapolated protocol converges across cell families within 0.03 eV; Tier 3 carries the
tolerance default."""
T2_EXTRAPOLATION_SPREAD_EV: dict[str, float] = {"basal": 0.08, "monovacancy": 0.13}
"""Spread of the Tier-2 zero-temperature extrapolation between 4500 and 6000 K in the 12x12 cell
(manual C2, Convergence and sensitivity): a class-level systematic of the Tier-2 protocol, which
the class offsets of the calibration absorb for anchored classes and the between-class variance
carries for the others; it is not a per-configuration noise."""
POOLED = "__pooled__"
MIN_ANCHORS_FOR_SLOPE = 3
PRIOR_OFFSET_SD_EV = 0.5
"""Scale of the class-offset prior used before any anchor is seen."""
TIER0_SIGMA_FLOOR_EV = 0.05
"""Uncertainty floor carried with the Tier-0 offset while it rests on a single anchor pair, the
scale the manual states for the transferability assumption."""
HELD_OUT_GROUPS: dict[str, tuple[str, ...]] = {
    "curvature": ("cnt_armchair", "cnt_zigzag", "mixed_vacancy_cnt", "mixed_sw_cnt"),
    "oxygen": ("epoxide", "hydroxyl", "carbonyl_pair", "carboxyl"),
}


def configuration_name(identifier: Any) -> str:
    """The configuration part of a record identifier ``class/substrate/configuration``."""
    return str(identifier).split("/")[-1]


def t_inflation(dof: int, level: float = 0.90) -> float:
    """Factor turning a plug-in standard deviation estimated with ``dof`` degrees of freedom into
    the normal-equivalent width whose ``level`` interval has the Student-t coverage."""
    q = 0.5 + level / 2
    return float(stats.t.ppf(q, max(int(dof), 1)) / stats.norm.ppf(q))


# --------------------------------------------------------------------------- anchor pairs

_SIZE_TOKENS = re.compile(r"_(\d+x\d+|x\d+)(?=_|$)")


def substrate_family(name: Any) -> str:
    """A substrate name without its cell, width or length tokens, so that the same site in
    different cells compares (``graphene_5x5`` and ``graphene_12x12`` give ``graphene``,
    ``swcnt_5_5_x3`` and ``swcnt_5_5_x6`` give ``swcnt_5_5``, the ribbons keep their edge)."""
    return _SIZE_TOKENS.sub("", str(name))


PAIR_COLUMNS = (
    "site_class",
    "substrate_family",
    "configuration",
    "n_H",
    "E_T1",
    "E_t",
    "supercell_T1",
    "supercell_t",
    "tier",
    *DESCRIPTORS,
)


def anchor_pairs(
    anchors: pd.DataFrame,
    landscapes: pd.DataFrame | Sequence[pd.DataFrame],
    quantity: str = "E_bind_H_atom",
    reference_state: str = "H_atom",
    engine: str = T1_ENGINE,
) -> pd.DataFrame:
    """Anchor configurations that a landscape table also holds, matched by class, substrate family
    (the same sheet defect, ribbon edge or tube chirality in any cell), configuration name and
    hydrogen count, with the landscape descriptors of the site. Several tables are searched in
    order, so that a tube of the verify tier can pair with an anchor the full tier lacks."""
    tables = [landscapes] if isinstance(landscapes, pd.DataFrame) else list(landscapes)
    a = anchors[
        (anchors["quantity"] == quantity)
        & (anchors["reference_state"] == reference_state)
        & (anchors["engine"] == engine)
    ].copy()
    if a.empty:
        return pd.DataFrame(columns=list(PAIR_COLUMNS))
    a["configuration"] = a["identifier"].map(configuration_name)
    a["substrate_family"] = a["supercell"].map(substrate_family)
    keys = ["site_class", "substrate_family", "configuration", "n_H"]
    rows: list[dict[str, Any]] = []
    for _, anchor in a.iterrows():
        for t in tables:
            sel = t[(t["quantity"] == quantity) & (t["reference_state"] == reference_state)]
            if sel.empty:
                continue
            sel = sel[
                (sel["site_class"] == anchor["site_class"])
                & (sel["n_H"] == anchor["n_H"])
                & (sel["identifier"].map(configuration_name) == anchor["configuration"])
                & (sel["supercell"].map(substrate_family) == anchor["substrate_family"])
            ]
            if sel.empty:
                continue
            row = sel.iloc[0]
            rows.append(
                {
                    **{k: anchor[k] for k in keys},
                    "E_T1": float(anchor["value_eV"]),
                    "E_t": float(row["value_eV"]),
                    "supercell_T1": anchor["supercell"],
                    "supercell_t": row["supercell"],
                    "tier": row["tier"],
                    **{d: float(row[d]) if pd.notna(row[d]) else float("nan") for d in DESCRIPTORS},
                }
            )
            break
    return pd.DataFrame(rows, columns=list(PAIR_COLUMNS))


# --------------------------------------------------------------------------- linear transfer


@dataclass
class Transfer:
    """Linear transfer E_T1 = a + b E_t of one class (or of the pooled anchors) and tier."""

    site_class: str
    tier: str
    n: int
    a: float
    b: float
    cov: np.ndarray
    resid_sd: float
    b_fixed: bool
    note: str | None = None
    dof: int = 1
    """Degrees of freedom behind the variance estimate; predictive intervals use Student-t
    quantiles with this many degrees of freedom, as the plug-in variance of few anchors demands."""

    def mean(self, e_t: np.ndarray) -> np.ndarray:
        return self.a + self.b * np.asarray(e_t, dtype=float)

    def var(self, e_t: np.ndarray) -> np.ndarray:
        e = np.asarray(e_t, dtype=float)
        x = np.stack([np.ones_like(e), e], axis=-1)
        return np.einsum("ij,jk,ik->i", x, self.cov, x)

    def as_dict(self) -> dict[str, Any]:
        return {
            "site_class": self.site_class,
            "tier": self.tier,
            "n_anchors": self.n,
            "a_eV": self.a,
            "b": self.b,
            "sigma_a_eV": math.sqrt(max(self.cov[0, 0], 0.0)),
            "sigma_b": math.sqrt(max(self.cov[1, 1], 0.0)),
            "residual_sd_eV": self.resid_sd,
            "b_fixed": self.b_fixed,
            "dof": self.dof,
            "note": self.note,
        }


def fit_transfer(
    e_t: np.ndarray,
    e_t1: np.ndarray,
    sigma: np.ndarray,
    site_class: str,
    tier: str,
    pooled_sd: float | None = None,
    pooled_dof: int = 1,
) -> Transfer:
    """Weighted least squares of E_T1 on E_t; with fewer than three anchors the slope is fixed to
    one and only the offset is estimated, its uncertainty taken from the larger of the class
    scatter and the pooled residual scatter of the tier (or the offset prior when no scatter is
    known), with the degrees of freedom of that scatter estimate recorded."""
    e_t = np.asarray(e_t, dtype=float)
    e_t1 = np.asarray(e_t1, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    n = len(e_t)
    if n == 0:
        raise ValueError("a transfer needs at least one anchor")
    w = 1.0 / sigma**2
    if n >= MIN_ANCHORS_FOR_SLOPE and np.ptp(e_t) > 1e-6:
        x = np.stack([np.ones(n), e_t], axis=1)
        xtw = x.T * w
        cov = np.linalg.inv(xtw @ x)
        beta = cov @ (xtw @ e_t1)
        r = e_t1 - x @ beta
        chi2 = float(np.sum(w * r**2))
        scale = max(1.0, chi2 / (n - 2))
        resid_sd = math.sqrt(float(np.sum(r**2) / (n - 2)))
        return Transfer(
            site_class,
            tier,
            n,
            float(beta[0]),
            float(beta[1]),
            cov * scale,
            resid_sd,
            False,
            dof=n - 2,
        )
    offsets = e_t1 - e_t
    a = float(np.sum(w * offsets) / np.sum(w))
    pooled = pooled_sd if pooled_sd is not None and pooled_sd > 0 else None
    if n >= 2:
        resid_sd = float(np.std(offsets, ddof=1))
        scatter = max(resid_sd, float(np.sqrt(1.0 / np.sum(w))), pooled or 0.0)
        var_a = scatter**2 / n
        dof = pooled_dof if pooled is not None and pooled >= resid_sd else n - 1
        note = "slope fixed to 1: two anchors"
    else:
        resid_sd = float("nan")
        prior = pooled if pooled is not None else PRIOR_OFFSET_SD_EV
        var_a = float(1.0 / np.sum(w)) + prior**2
        dof = pooled_dof if pooled is not None else 1
        note = "slope fixed to 1: one anchor; offset uncertainty from the pooled scatter"
    cov = np.array([[var_a, 0.0], [0.0, 0.0]])
    return Transfer(site_class, tier, n, a, 1.0, cov, resid_sd, True, note, dof=max(dof, 1))


# --------------------------------------------------------------------------- Gaussian process


def matern52(r: np.ndarray) -> np.ndarray:
    s = math.sqrt(5.0) * r
    return (1.0 + s + s**2 / 3.0) * np.exp(-s)


def matern32(r: np.ndarray) -> np.ndarray:
    s = math.sqrt(3.0) * r
    return (1.0 + s) * np.exp(-s)


KERNELS = {"matern52": matern52, "matern32": matern32}


@dataclass
class Standardizer:
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, x: np.ndarray) -> Standardizer:
        x = np.asarray(x, dtype=float)
        mean = np.nanmean(x, axis=0) if len(x) else np.zeros(x.shape[1])
        scale = np.nanstd(x, axis=0) if len(x) else np.ones(x.shape[1])
        mean = np.where(np.isfinite(mean), mean, 0.0)
        scale = np.where(np.isfinite(scale) & (scale > 1e-9), scale, 1.0)
        return cls(mean, scale)

    def __call__(self, x: np.ndarray) -> np.ndarray:
        z = (np.asarray(x, dtype=float) - self.mean) / self.scale
        return np.where(np.isfinite(z), z, 0.0)


@dataclass
class GPHyper:
    amplitude: float
    lengthscales: np.ndarray
    noise: float
    kernel: str = "matern52"
    log_marginal_likelihood: float = float("nan")
    optimized: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "amplitude_eV": self.amplitude,
            "lengthscales_standardized": [float(v) for v in self.lengthscales],
            "noise_eV": self.noise,
            "kernel": self.kernel,
            "log_marginal_likelihood": self.log_marginal_likelihood,
            "optimized": self.optimized,
        }


def _distance(za: np.ndarray, zb: np.ndarray, lengthscales: np.ndarray) -> np.ndarray:
    d = (za[:, None, :] - zb[None, :, :]) / lengthscales
    return np.sqrt(np.sum(d**2, axis=-1))


def kernel_matrix(za: np.ndarray, zb: np.ndarray, hyper: GPHyper) -> np.ndarray:
    return hyper.amplitude**2 * KERNELS[hyper.kernel](_distance(za, zb, hyper.lengthscales))


def _block_lml(z: np.ndarray, r: np.ndarray, sigma: np.ndarray, hyper: GPHyper) -> float:
    k = kernel_matrix(z, z, hyper) + np.diag(sigma**2 + hyper.noise**2)
    try:
        chol = np.linalg.cholesky(k)
    except np.linalg.LinAlgError:
        return -1e30
    alpha = np.linalg.solve(chol.T, np.linalg.solve(chol, r))
    return float(
        -0.5 * r @ alpha - np.sum(np.log(np.diag(chol))) - 0.5 * len(r) * math.log(2 * math.pi)
    )


def fit_gp_hyper(
    blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    n_dim: int,
    kernel: str = "matern52",
    amplitude0: float = 0.1,
    noise0: float = T1_SIGMA_EV,
) -> GPHyper:
    """Shared hyperparameters (amplitude, per-descriptor lengthscales, noise) maximizing the sum of
    the marginal likelihoods of independent blocks of residuals, with weak log-normal penalties
    that keep the optimum finite where the residuals are few. Fewer than three residuals in total
    leave the starting values (not optimized)."""
    n_total = sum(len(b[1]) for b in blocks)
    start = GPHyper(amplitude0, np.ones(n_dim), noise0, kernel)
    if n_total < 3:
        return start

    def unpack(theta: np.ndarray) -> GPHyper:
        return GPHyper(
            math.exp(theta[0]), np.exp(theta[1 : 1 + n_dim]), math.exp(theta[1 + n_dim]), kernel
        )

    def objective(theta: np.ndarray) -> float:
        h = unpack(theta)
        lml = sum(_block_lml(z, r, s, h) for z, r, s in blocks if len(r))
        penalty = (
            0.5 * ((theta[0] - math.log(amplitude0)) / 1.5) ** 2
            + 0.5 * float(np.sum((theta[1 : 1 + n_dim] / 1.5) ** 2))
            + 0.5 * ((theta[1 + n_dim] - math.log(noise0)) / 1.0) ** 2
        )
        return -lml + penalty

    theta0 = np.concatenate([[math.log(amplitude0)], np.zeros(n_dim), [math.log(noise0)]])
    bounds = [(math.log(1e-3), math.log(5.0))] + [(math.log(0.1), math.log(50.0))] * n_dim
    bounds += [(math.log(1e-3), math.log(1.0))]
    res = optimize.minimize(objective, theta0, method="L-BFGS-B", bounds=bounds)
    hyper = unpack(res.x)
    hyper.log_marginal_likelihood = float(
        sum(_block_lml(z, r, s, hyper) for z, r, s in blocks if len(r))
    )
    hyper.optimized = bool(res.success)
    return hyper


def gp_predict(
    z_train: np.ndarray,
    r_train: np.ndarray,
    sigma_train: np.ndarray,
    z_test: np.ndarray,
    hyper: GPHyper,
) -> tuple[np.ndarray, np.ndarray]:
    """Posterior mean and covariance of the discrepancy at ``z_test``; with no training residuals
    the prior (zero mean, amplitude² variance) is returned."""
    n_test = len(z_test)
    if len(r_train) == 0:
        return np.zeros(n_test), kernel_matrix(z_test, z_test, hyper)
    k = kernel_matrix(z_train, z_train, hyper) + np.diag(sigma_train**2 + hyper.noise**2)
    ks = kernel_matrix(z_test, z_train, hyper)
    chol = np.linalg.cholesky(k)
    alpha = np.linalg.solve(chol.T, np.linalg.solve(chol, r_train))
    v = np.linalg.solve(chol, ks.T)
    mean = ks @ alpha
    cov = kernel_matrix(z_test, z_test, hyper) - v.T @ v
    return mean, cov


# --------------------------------------------------------------------------- calibration


@dataclass
class ClassState:
    """Transfer and the residuals (with their descriptor coordinates) the discrepancy process is
    conditioned on when this state predicts."""

    transfer: Transfer
    z: np.ndarray
    residuals: np.ndarray
    sigma: np.ndarray
    n_anchors: int
    class_offset: float = 0.0
    """Posterior mean of the class random effect u_c (zero for a class with its own transfer)."""
    offset_var: float = 0.0
    """Posterior variance of u_c: the between-class variance for a class without anchors, shrunk
    by the anchors of the class otherwise."""


@dataclass
class Calibration:
    """Per-class transfers and the pooled discrepancy process of one tier."""

    tier: str
    kernel: str = "matern52"
    transfers: dict[str, ClassState] = field(default_factory=dict)
    pooled: ClassState | None = None
    hyper: GPHyper | None = None
    standardizer: Standardizer | None = None
    pairs: pd.DataFrame = field(default_factory=pd.DataFrame)
    sigma_t: float = TIER_SIGMA_EV["T2"]
    tau: float = 0.0
    """Between-class standard deviation of the class random effect, from the fitted noise of the
    discrepancy process after the per-configuration numerical noise is removed."""

    @property
    def anchored_classes(self) -> list[str]:
        return sorted(self.transfers)

    def n_anchors(self, site_class: str) -> int:
        st = self.transfers.get(site_class)
        return 0 if st is None else st.n_anchors

    def status(self, site_class: str) -> str:
        n = self.n_anchors(site_class)
        if n >= MIN_ANCHORS_FOR_SLOPE:
            return "calibrated"
        if n >= 1:
            return "offset_calibrated"
        return "awaiting_anchors"

    def fit(self, pairs: pd.DataFrame) -> Calibration:
        pairs = pairs.reset_index(drop=True)
        self.pairs = pairs
        x = (
            pairs[list(DESCRIPTORS)].to_numpy(dtype=float)
            if len(pairs)
            else np.zeros((0, len(DESCRIPTORS)))
        )
        self.standardizer = Standardizer.fit(x)
        z_all = self.standardizer(x)
        sig = np.full(len(pairs), math.hypot(T1_SIGMA_EV, self.sigma_t))
        if len(pairs) == 0:
            self.hyper = GPHyper(
                PRIOR_OFFSET_SD_EV,
                np.ones(len(DESCRIPTORS)),
                math.hypot(T1_SIGMA_EV, self.sigma_t),
                self.kernel,
            )
            self.pooled = None
            return self
        e_t = pairs["E_t"].to_numpy(dtype=float)
        e_t1 = pairs["E_T1"].to_numpy(dtype=float)
        sigma_e2 = T1_SIGMA_EV**2 + self.sigma_t**2
        pooled_tr = fit_transfer(e_t, e_t1, sig, POOLED, self.tier)
        pooled_res = e_t1 - pooled_tr.mean(e_t)
        pooled_sd = float(np.std(pooled_res, ddof=1)) if len(pairs) > 1 else PRIOR_OFFSET_SD_EV
        pooled_dof = max(len(pairs) - 2, 1)
        pooled_tr.dof = pooled_dof
        groups = {
            str(c): np.asarray(list(idx)) for c, idx in pairs.groupby("site_class").groups.items()
        }
        # between-class variance tau^2 of the class random effect, by the method of moments on
        # the class mean residuals of the pooled transfer (the prior scale with a single class)
        means = np.array([float(np.mean(pooled_res[ii])) for ii in groups.values()])
        counts = np.array([len(ii) for ii in groups.values()], dtype=float)
        if len(means) > 1:
            tau2 = max(float(np.var(means, ddof=1)) - float(np.mean(sigma_e2 / counts)), 1e-6)
        else:
            tau2 = PRIOR_OFFSET_SD_EV**2
        self.tau = math.sqrt(tau2)
        self.pooled = ClassState(pooled_tr, z_all[:0], pooled_res[:0], sig[:0], 0, 0.0, tau2)
        # the discrepancy process is fitted to the within-class residuals (around the class
        # transfer, or around the pooled transfer with the class offset), the structure the
        # descriptors can explain; while no class has two anchors its amplitude is the prior
        # scale tau, since within-class variation is then unobserved
        blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        states: dict[str, tuple[Transfer, np.ndarray, int, float, float]] = {}
        for cls, ii in groups.items():
            n = len(ii)
            if n >= MIN_ANCHORS_FOR_SLOPE:
                tr = fit_transfer(e_t[ii], e_t1[ii], sig[ii], cls, self.tier, pooled_sd, pooled_dof)
                res = e_t1[ii] - tr.mean(e_t[ii])
                states[cls] = (tr, res, n, 0.0, 0.0)
            else:
                res = pooled_res[ii]
                states[cls] = (pooled_tr, res, n, float(np.mean(res)), 0.0)
            if n >= 2:
                blocks.append((z_all[ii], res - np.mean(res), sig[ii]))
        within_dof = sum(len(b[1]) - 1 for b in blocks)
        if within_dof >= 3:
            self.hyper = fit_gp_hyper(
                blocks,
                len(DESCRIPTORS),
                self.kernel,
                amplitude0=max(self.tau, 0.05),
                noise0=math.sqrt(sigma_e2),
            )
        else:
            self.hyper = GPHyper(
                self.tau, np.ones(len(DESCRIPTORS)), math.sqrt(sigma_e2), self.kernel
            )
        amp2 = self.hyper.amplitude**2
        for cls, (tr, res, n, mean_res, _) in states.items():
            ii = groups[cls]
            if n >= MIN_ANCHORS_FOR_SLOPE:
                self.transfers[cls] = ClassState(tr, z_all[ii], res, sig[ii], n)
                continue
            # posterior of u_c ~ N(0, tau^2) from n anchors whose residuals scatter by the
            # numerical noise and the process amplitude around the class offset
            per_conf = sigma_e2 + amp2
            u_var = 1.0 / (n / per_conf + 1.0 / tau2)
            u_mean = u_var * n * mean_res / per_conf
            self.transfers[cls] = ClassState(
                pooled_tr, z_all[ii], res - u_mean, sig[ii], n, u_mean, u_var
            )
        return self

    def _state(self, site_class: str) -> ClassState | None:
        return self.transfers.get(site_class, self.pooled)

    def _predict_arrays(
        self, site_class: str, e_t: np.ndarray, z: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """Mean, covariance (before the Student-t inflation) of the Tier-1 energies of one class,
        the linear-transfer variance and the degrees of freedom."""
        assert self.hyper is not None
        st = self._state(site_class)
        if st is None:
            cov = kernel_matrix(z, z, self.hyper) + PRIOR_OFFSET_SD_EV**2 * np.ones(
                (len(z), len(z))
            )
            return e_t.copy(), cov, np.zeros(len(e_t)), 1
        lin_mean = st.transfer.mean(e_t)
        lin_var = st.transfer.var(e_t)
        gp_mean, gp_cov = gp_predict(st.z, st.residuals, st.sigma, z, self.hyper)
        sigma_e2 = T1_SIGMA_EV**2 + self.sigma_t**2
        n = len(e_t)
        # the class random effect is common to every configuration of the class, so its variance
        # enters every element of the covariance and does not average out over configurations
        cov = gp_cov + np.diag(lin_var + sigma_e2) + st.offset_var * np.ones((n, n))
        return lin_mean + gp_mean + st.class_offset, cov, lin_var, st.transfer.dof

    def predict(self, frame: pd.DataFrame) -> pd.DataFrame:
        """Calibrated mean and standard deviation on the Tier-1 scale for rows holding
        ``site_class``, ``value_eV`` and the descriptors, with the per-row status and transfer."""
        assert self.hyper is not None and self.standardizer is not None
        out = frame.copy()
        mean = np.full(len(out), np.nan)
        sd = np.full(len(out), np.nan)
        status = np.array(["awaiting_anchors"] * len(out), dtype=object)
        a_col = np.full(len(out), np.nan)
        b_col = np.full(len(out), np.nan)
        n_col = np.zeros(len(out), dtype=int)
        if len(out) == 0:
            return out.assign(
                calibrated_mean_eV=mean,
                calibrated_sd_eV=sd,
                calibration_status=status,
                transfer_a_eV=a_col,
                transfer_b=b_col,
                n_anchors_class=n_col,
            )
        x = out[list(DESCRIPTORS)].to_numpy(dtype=float)
        z = self.standardizer(x)
        e_t = out["value_eV"].to_numpy(dtype=float)
        for cls in out["site_class"].unique():
            ii = np.flatnonzero(out["site_class"].to_numpy() == cls)
            m, cov, _, dof = self._predict_arrays(str(cls), e_t[ii], z[ii])
            mean[ii] = m
            sd[ii] = np.sqrt(np.clip(np.diag(cov), 0.0, None)) * t_inflation(dof)
            status[ii] = self.status(str(cls))
            st = self._state(str(cls))
            if st is not None:
                a_col[ii] = st.transfer.a + st.class_offset
                b_col[ii] = st.transfer.b
            n_col[ii] = self.n_anchors(str(cls))
        return out.assign(
            calibrated_mean_eV=mean,
            calibrated_sd_eV=sd,
            calibration_status=status,
            transfer_a_eV=a_col,
            transfer_b=b_col,
            n_anchors_class=n_col,
        )

    def predict_class_mean(self, frame: pd.DataFrame) -> tuple[float, float]:
        """Predictive mean and standard deviation of the average Tier-1 energy over the rows of
        one class, with the covariance of the discrepancy process and of the transfer included."""
        assert self.hyper is not None and self.standardizer is not None
        cls = str(frame["site_class"].iloc[0])
        e_t = frame["value_eV"].to_numpy(dtype=float)
        z = self.standardizer(frame[list(DESCRIPTORS)].to_numpy(dtype=float))
        m, cov, lin_var, dof = self._predict_arrays(cls, e_t, z)
        n = len(e_t)
        ones = np.ones(n) / n
        var = float(ones @ cov @ ones)
        st = self._state(cls)
        if st is not None:
            # the transfer parameters are common to every row: their variance does not average out
            x_bar = np.array([1.0, float(np.mean(e_t))])
            var += float(x_bar @ st.transfer.cov @ x_bar) - float(np.mean(lin_var)) / n
        return float(np.mean(m)), math.sqrt(max(var, 0.0)) * t_inflation(dof)

    def report(self) -> dict[str, Any]:
        assert self.hyper is not None
        classes: dict[str, Any] = {}
        for cls, st in sorted(self.transfers.items()):
            d = st.transfer.as_dict()
            d.update(
                site_class=cls,
                n_anchors=st.n_anchors,
                class_offset_eV=st.class_offset,
                status=self.status(cls),
            )
            if st.n_anchors < MIN_ANCHORS_FOR_SLOPE:
                d["note"] = (
                    "pooled transfer; the discrepancy process interpolates the class anchors"
                    if st.n_anchors
                    else d["note"]
                )
            classes[cls] = d
        return {
            "tier": self.tier,
            "n_anchor_pairs": int(len(self.pairs)),
            "classes": classes,
            "pooled_transfer": None if self.pooled is None else self.pooled.transfer.as_dict(),
            "gp": self.hyper.as_dict(),
        }


def loo_coverage(
    pairs: pd.DataFrame,
    tier: str,
    level: float = 0.90,
    kernel: str = "matern52",
    sigma_t: float | None = None,
) -> dict[str, Any]:
    """Leave-one-anchor-out predictive coverage: each anchor is predicted from a calibration
    fitted to the others; the fraction of anchors inside the ``level`` interval is reported with
    the standardized residuals."""
    zq = stats.norm.ppf(0.5 + level / 2)
    sigma_t = TIER_SIGMA_EV.get(tier, 0.05) if sigma_t is None else sigma_t
    hits: list[bool] = []
    zs: list[float] = []
    rows: list[dict[str, Any]] = []
    for i in range(len(pairs)):
        train = pairs.drop(index=pairs.index[i])
        test = pairs.iloc[[i]].rename(columns={"E_t": "value_eV"})
        cal = Calibration(tier, kernel, sigma_t=sigma_t).fit(train)
        pred = cal.predict(test)
        m = float(pred["calibrated_mean_eV"].iloc[0])
        s = float(pred["calibrated_sd_eV"].iloc[0])
        truth = float(pairs["E_T1"].iloc[i])
        hit = abs(truth - m) <= zq * s
        hits.append(bool(hit))
        zs.append((truth - m) / s if s > 0 else float("nan"))
        rows.append(
            {
                "site_class": pairs["site_class"].iloc[i],
                "configuration": pairs["configuration"].iloc[i],
                "E_T1": truth,
                "predicted": m,
                "sd": s,
                "covered": bool(hit),
                "status_without": cal.status(str(pairs["site_class"].iloc[i])),
            }
        )
    n = len(hits)
    return {
        "level": level,
        "n": n,
        "coverage": float(np.mean(hits)) if n else float("nan"),
        "rmse_eV": float(np.sqrt(np.mean([(r["E_T1"] - r["predicted"]) ** 2 for r in rows])))
        if n
        else float("nan"),
        "standardized_residual_sd": float(np.nanstd(zs)) if n > 1 else float("nan"),
        "passes_0.85_0.95": bool(n and 0.85 <= float(np.mean(hits)) <= 0.95),
        "anchors": rows,
    }


def held_out_class_errors(
    pairs: pd.DataFrame,
    tier: str,
    groups: dict[str, tuple[str, ...]] | None = None,
    kernel: str = "matern52",
) -> dict[str, Any]:
    """Out-of-distribution check of the manual (8.3): the anchors of a class group are withheld,
    the calibration is fitted to the rest and the withheld anchors are predicted."""
    groups = groups or HELD_OUT_GROUPS
    out: dict[str, Any] = {}
    for name, classes in groups.items():
        mask = pairs["site_class"].isin(classes)
        if not mask.any():
            out[name] = {"n_held_out": 0, "note": "no anchors in these classes yet"}
            continue
        cal = Calibration(tier, kernel).fit(pairs[~mask])
        pred = cal.predict(pairs[mask].rename(columns={"E_t": "value_eV"}))
        err = pred["calibrated_mean_eV"].to_numpy() - pairs.loc[mask, "E_T1"].to_numpy()
        sd = pred["calibrated_sd_eV"].to_numpy()
        out[name] = {
            "n_held_out": int(mask.sum()),
            "classes": sorted(pairs.loc[mask, "site_class"].unique()),
            "rmse_eV": float(np.sqrt(np.mean(err**2))),
            "mean_error_eV": float(np.mean(err)),
            "coverage_90": float(np.mean(np.abs(err) <= 1.645 * sd)),
            "errors_eV": [float(e) for e in err],
        }
    return out


def synthetic_recovery(
    n_replicates: int = 100,
    seed: int = 0,
    anchors_per_class: tuple[int, ...] = (1, 2, 3, 4, 6, 8),
    tier: str = "T2",
    level: float = 0.90,
) -> dict[str, Any]:
    """Acceptance criterion (iii) of the manual for the calibration: in each replicate class
    transfers (a_c, b_c) and a Matérn discrepancy are planted, anchors are drawn with Tier-1 and
    tier noise, the calibration is fitted and the planted class mean of E_T1 over the class's
    configurations is checked against the ``level`` predictive interval of the class mean."""
    rng = np.random.default_rng(seed)
    zq = stats.norm.ppf(0.5 + level / 2)
    covered: list[bool] = []
    per_class: dict[int, list[bool]] = {n: [] for n in anchors_per_class}
    for _ in range(n_replicates):
        hyper = GPHyper(0.08, rng.uniform(1.0, 3.0, len(DESCRIPTORS)), 0.0)
        frames = []
        truth: dict[str, float] = {}
        for k, n in enumerate(anchors_per_class):
            cls = f"class_{k}"
            a, b = rng.normal(0.0, 0.3), rng.normal(1.0, 0.1)
            n_pred = 6
            e_t = rng.normal(-1.5 - 0.5 * k, 0.5, n + n_pred)
            x = rng.normal(0.0, 1.0, (n + n_pred, len(DESCRIPTORS)))
            cov = kernel_matrix(x, x, hyper) + 1e-10 * np.eye(n + n_pred)
            delta = rng.multivariate_normal(np.zeros(n + n_pred), cov)
            e_t1 = a + b * e_t + delta
            df = pd.DataFrame(x, columns=list(DESCRIPTORS))
            df["site_class"] = cls
            df["configuration"] = [f"c{i}" for i in range(n + n_pred)]
            df["n_H"] = 1
            df["E_t"] = e_t + rng.normal(0.0, TIER_SIGMA_EV["T2"], n + n_pred)
            df["E_T1"] = e_t1 + rng.normal(0.0, T1_SIGMA_EV, n + n_pred)
            df["truth"] = e_t1
            df["anchor"] = [True] * n + [False] * n_pred
            frames.append(df)
            truth[cls] = float(np.mean(e_t1[n:]))
        data = pd.concat(frames, ignore_index=True)
        cal = Calibration(tier).fit(data[data["anchor"]].drop(columns=["truth", "anchor"]))
        for k, n in enumerate(anchors_per_class):
            cls = f"class_{k}"
            test = data[(data["site_class"] == cls) & (~data["anchor"])].rename(
                columns={"E_t": "value_eV"}
            )
            m, s = cal.predict_class_mean(test)
            hit = abs(truth[cls] - m) <= zq * s
            covered.append(bool(hit))
            per_class[n].append(bool(hit))
    frac = float(np.mean(covered))
    return {
        "n_replicates": n_replicates,
        "level": level,
        "coverage_of_planted_class_means": frac,
        "coverage_by_anchor_count": {str(n): float(np.mean(v)) for n, v in per_class.items()},
        "passes_0.90": bool(frac >= 0.90),
        "seed": seed,
    }


# --------------------------------------------------------------------------- Tier-0 offset

TIER0_SITE_MAP: dict[tuple[str, str], tuple[str, int]] = {
    ("basal", "top_monomer"): ("monomer", 1),
}
"""Published site identifiers of Tier-0 anchors mapped to the catalogue configuration and
hydrogen count they correspond to."""


def tier0_offset(
    anchors: pd.DataFrame, registry: pd.DataFrame, engine: str = T1_ENGINE
) -> dict[str, Any]:
    """Tier-1 bias against the Tier-0 (DMC, CCSD(T)) anchors: the mean of E_T1 minus E_T0 over the
    anchor pairs that exist, its statistical uncertainty (Tier-1 numerical and stated Tier-0
    uncertainties, plus the spread when more than one pair exists) and the uncertainty carried,
    which does not fall below the manual's transferability scale while one pair decides it."""
    pairs: list[dict[str, Any]] = []
    t0 = registry[registry["tier0"]] if "tier0" in registry else registry.iloc[0:0]
    for _, row in t0.iterrows():
        key = (row["site_class"], row["site_id"])
        if key not in TIER0_SITE_MAP:
            continue
        configuration, n_h = TIER0_SITE_MAP[key]
        sel = anchors[
            (anchors["site_class"] == row["site_class"])
            & (anchors["quantity"] == row["quantity"])
            & (anchors["reference_state"] == row["reference_state"])
            & (anchors["engine"] == engine)
            & (anchors["n_H"] == n_h)
            & (anchors["identifier"].map(configuration_name) == configuration)
        ]
        if sel.empty:
            continue
        e_t1 = float(sel["value_eV"].iloc[0])
        pairs.append(
            {
                "published_id": row["id"],
                "site_class": row["site_class"],
                "configuration": configuration,
                "level": row["level"],
                "E_T0_eV": float(row["value_eV"]),
                "sigma_T0_eV": float(row["sigma_eV"]),
                "E_T1_eV": e_t1,
                "offset_eV": e_t1 - float(row["value_eV"]),
            }
        )
    if not pairs:
        return {
            "n_pairs": 0,
            "offset_eV": 0.0,
            "sigma_statistical_eV": float("nan"),
            "sigma_carried_eV": TIER0_SIGMA_FLOOR_EV,
            "note": "no Tier-0 anchor pair available; offset set to zero with the "
            "transferability floor",
            "pairs": [],
        }
    offsets = np.array([p["offset_eV"] for p in pairs])
    stat = math.sqrt(float(np.mean([p["sigma_T0_eV"] ** 2 for p in pairs])) + T1_SIGMA_EV**2)
    if len(pairs) > 1:
        stat = math.hypot(stat, float(np.std(offsets, ddof=1)))
    carried = stat if len(pairs) > 1 else max(stat, TIER0_SIGMA_FLOOR_EV)
    return {
        "n_pairs": len(pairs),
        "offset_eV": float(np.mean(offsets)),
        "sigma_statistical_eV": stat,
        "sigma_carried_eV": carried,
        "note": (
            "estimated from the anchor pairs; carried uncertainty held at the transferability "
            "floor while one pair decides it"
            if len(pairs) == 1
            else "estimated from the anchor pairs"
        ),
        "pairs": pairs,
    }


__all__ = [
    "DESCRIPTORS",
    "HELD_OUT_GROUPS",
    "Calibration",
    "GPHyper",
    "T1_ENGINE",
    "T1_SIGMA_EV",
    "TIER_SIGMA_EV",
    "Transfer",
    "anchor_pairs",
    "configuration_name",
    "fit_gp_hyper",
    "fit_transfer",
    "gp_predict",
    "held_out_class_errors",
    "kernel_matrix",
    "loo_coverage",
    "matern32",
    "matern52",
    "substrate_family",
    "synthetic_recovery",
    "t_inflation",
    "tier0_offset",
]
