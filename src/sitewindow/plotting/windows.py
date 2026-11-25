"""Figure 1 of Manuscript 1: swing functions, optimum and threshold energies, admissibility boundary."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from sitewindow.constants import EV_TO_KJ_PER_MOL  # noqa: E402
from sitewindow.kinetics.polanyi_wigner import t90_first_order  # noqa: E402
from sitewindow.theory.entropy import DEFAULT_MODEL_SET, ConstantEntropy  # noqa: E402
from sitewindow.theory.window import (  # noqa: E402
    W1,
    W2,
    W3,
    OperatingWindow,
    optimum_energy,
    swing,
)

OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def figure_window_theory(
    out_dir: Path, windows: tuple[OperatingWindow, ...] = (W1, W2, W3)
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 8,
            "font.family": "sans-serif",
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "figure.dpi": 150,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4), constrained_layout=True)
    eps_kJ = np.linspace(-40.0, -0.5, 800)
    eps_eV = eps_kJ / EV_TO_KJ_PER_MOL
    ax = axes[0]
    for w, c in zip(windows, OKABE_ITO):
        s = swing(eps_eV, w, ConstantEntropy())
        ax.plot(-eps_kJ, s, color=c, lw=1.2, label=w.name)
        e_star, eta = optimum_energy(w, ConstantEntropy())
        ax.plot(-e_star * EV_TO_KJ_PER_MOL, eta, "o", color=c, ms=3.5)
    ax.axvspan(4.0, 8.0, color="0.85", lw=0, label="carbon physisorption")
    ax.set_xlabel("$-\\varepsilon$ (kJ mol$^{-1}$ per H$_2$)")
    ax.set_ylabel("swing $\\Delta\\theta$")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=6.5)
    ax.text(0.02, 0.95, "(a)", transform=ax.transAxes, va="top")
    # (b) entropy-model dependence of the optimum for W1 and W2
    ax = axes[1]
    for w, c in zip(windows[:2], OKABE_ITO):
        for m, ls in zip(DEFAULT_MODEL_SET, ("-", "--", ":")):
            s = swing(eps_eV, w, m)
            ax.plot(
                -eps_kJ, s, color=c, lw=1.0, ls=ls, label=f"{w.name}, {m.name.replace('_', ' ')}"
            )
    ax.set_xlabel("$-\\varepsilon$ (kJ mol$^{-1}$ per H$_2$)")
    ax.set_ylabel("swing $\\Delta\\theta$")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=5.5, ncol=1)
    ax.text(0.02, 0.95, "(b)", transform=ax.transAxes, va="top")
    # (c) release admissibility boundary: T90 vs barrier for tau_r = 1e3 s, prefactor band
    ax = axes[2]
    E = np.linspace(0.3, 3.0, 300)
    for nu, c, lab in (
        (1e12, OKABE_ITO[5], "$\\nu=10^{12}$ s$^{-1}$"),
        (1e13, OKABE_ITO[0], "$\\nu=10^{13}$ s$^{-1}$"),
        (1e14, OKABE_ITO[6], "$\\nu=10^{14}$ s$^{-1}$"),
    ):
        ax.plot(E, t90_first_order(E, nu, 1.0e3), color=c, lw=1.0, label=lab)
    for T, name in ((358.15, "358 K"), (473.15, "473 K"), (573.15, "573 K")):
        ax.axhline(T, color="0.6", lw=0.6, ls="--")
        ax.text(0.32, T + 8, name, fontsize=6, color="0.3")
    ax.set_xlabel("release barrier $E_d^\\ddagger$ (eV)")
    ax.set_ylabel("$T_{90}$ for $\\tau_r = 10^3$ s (K)")
    ax.set_ylim(200, 1200)
    ax.legend(frameon=False, fontsize=6)
    ax.text(0.02, 0.95, "(c)", transform=ax.transAxes, va="top")
    paths = []
    for ext in ("pdf", "png"):
        p = out_dir / f"fig1_window_theory_verify.{ext}"
        fig.savefig(p, dpi=600 if ext == "png" else None)
        paths.append(p)
    plt.close(fig)
    return paths


__all__ = ["figure_window_theory"]
