"""Deliverable-capacity functional, mass accounting, and the variational bound.

Implements DERIVATIONS 5.1-5.5. A landscape is a table of discrete sites with weights
(sites per carbon atom), channel, binding energy, release barrier, prefactor, multiplicity,
maximum occupancy and the lattice mass change of the defect that creates the site.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from sitewindow.constants import GRAPHENE_SSA_M2_G, M_C_AMU, M_H2_AMU, M_H_AMU, N_A
from sitewindow.kinetics.admissibility import release_admissible
from sitewindow.theory.entropy import ConstantEntropy, EntropyModel
from sitewindow.theory.window import OperatingWindow, occupancy, optimum_energy


@dataclass
class SiteTable:
    """Discrete landscape.

    Attributes
    ----------
    eps_eV
        Binding energy per entity (H2 for molecular, H relative to 1/2 H2 for atomic).
    E_d_eV
        Release barrier (eV); for molecular sites the desorption barrier equals ``-eps``
        when ``E_d_eV`` is NaN.
    nu_s
        Prefactor (s^-1).
    order
        Kinetic order of release (1 or 2).
    sites_per_C
        Number of such sites per carbon atom (x_c * g_s).
    n_max
        Maximum entities per site.
    dmass_per_C_amu
        Lattice mass change per carbon attributable to the defect that creates the site,
        distributed over the sites it creates.
    channel
        0 for molecular, 1 for atomic.
    """

    eps_eV: NDArray[np.float64]
    E_d_eV: NDArray[np.float64]
    nu_s: NDArray[np.float64]
    order: NDArray[np.int64]
    sites_per_C: NDArray[np.float64]
    n_max: NDArray[np.float64]
    dmass_per_C_amu: NDArray[np.float64]
    channel: NDArray[np.int64]
    labels: list[str] | None = None

    def __post_init__(self) -> None:
        n = len(self.eps_eV)
        for name in (
            "E_d_eV",
            "nu_s",
            "order",
            "sites_per_C",
            "n_max",
            "dmass_per_C_amu",
            "channel",
        ):
            arr = np.asarray(getattr(self, name))
            if arr.shape != (n,):
                raise ValueError(f"{name} must have shape ({n},), got {arr.shape}")
            setattr(self, name, arr)
        self.eps_eV = np.asarray(self.eps_eV, dtype=float)
        if np.any(self.sites_per_C < 0) or np.any(self.n_max < 0):
            raise ValueError("site densities and occupancies must be non-negative")


@dataclass(frozen=True)
class CapacityResult:
    D_wt_percent: float
    charged_wt_percent: float
    per_site_delivered_H_per_C: NDArray[np.float64]
    swing: NDArray[np.float64]
    admissible: NDArray[np.bool_]
    occupancy_deficit: float
    admissibility_deficit: float
    mass_deficit: float


def mass_fraction(h_per_C: float, dmass_per_C_amu: float, entity_mass_amu: float) -> float:
    """Mass fraction of hydrogen entities per carbon (DERIVATIONS 5.2-5.3)."""
    m_h = h_per_C * entity_mass_amu
    return m_h / (M_C_AMU + dmass_per_C_amu + m_h)


def deliverable(
    table: SiteTable,
    window: OperatingWindow,
    entropy: EntropyModel | None = None,
    *,
    delta_s_atomic_eV_K: float = 0.0,
    saturate_atomic_charge: bool | None = None,
) -> CapacityResult:
    """Deliverable capacity (wt%) for a discrete landscape under a window (DERIVATIONS 5.3).

    For the atomic channel with ``window.charging_mode == 'atomic_H_source'`` the charged
    occupancy is taken as 1 (saturation by the atomic source) and the released occupancy as the
    kinetic residual; otherwise the dissociative Langmuir equilibrium is used on both sides.
    """
    entropy = entropy or ConstantEntropy()
    mol = table.channel == 0
    atm = ~mol
    th_c = np.zeros_like(table.eps_eV)
    th_r = np.zeros_like(table.eps_eV)
    if mol.any():
        th_c[mol] = occupancy(
            table.eps_eV[mol], window.T_c, window.P_c, entropy, use_fugacity=window.fugacity
        )
        th_r[mol] = occupancy(
            table.eps_eV[mol], window.T_r, window.P_r, entropy, use_fugacity=window.fugacity
        )
    if atm.any():
        if saturate_atomic_charge is None:
            saturate_atomic_charge = window.charging_mode == "atomic_H_source"
        if saturate_atomic_charge:
            th_c[atm] = 1.0
            th_r[atm] = 0.0  # released fraction enters through admissibility/release fraction
        else:
            th_c[atm] = occupancy(
                table.eps_eV[atm],
                window.T_c,
                window.P_c,
                channel="atomic",
                delta_s_atomic_eV_K=delta_s_atomic_eV_K,
                use_fugacity=window.fugacity,
            )
            th_r[atm] = occupancy(
                table.eps_eV[atm],
                window.T_r,
                window.P_r,
                channel="atomic",
                delta_s_atomic_eV_K=delta_s_atomic_eV_K,
                use_fugacity=window.fugacity,
            )
    dtheta = th_c - th_r
    E_d = np.where(np.isnan(table.E_d_eV), -table.eps_eV, table.E_d_eV)
    adm, frac = release_admissible(
        E_d, table.nu_s, table.order, window.T_r, window.tau_r, th_c, window.eta_r
    )
    # Released amount: thermodynamic swing where admissible; for atomic saturated charging use the
    # kinetic release fraction directly (thermodynamic residual is ~0 at these energies).
    delivered = np.where(adm, dtheta, 0.0)
    if atm.any() and saturate_atomic_charge:
        delivered[atm] = np.where(adm[atm], frac[atm] * th_c[atm], 0.0)
    entity_mass = np.where(mol, M_H2_AMU, M_H_AMU)
    h_del = table.sites_per_C * table.n_max * delivered
    h_charged = table.sites_per_C * table.n_max * th_c
    m_del = float(np.sum(h_del * entity_mass))
    m_charged = float(np.sum(h_charged * entity_mass))
    dmass = float(np.sum(table.dmass_per_C_amu))
    D = 100.0 * m_del / (M_C_AMU + dmass + m_charged)
    charged = 100.0 * m_charged / (M_C_AMU + dmass + m_charged)
    # Deficit decomposition relative to the bound for the dominant channel.
    cap_sites = float(np.sum(table.sites_per_C * table.n_max * entity_mass))
    if cap_sites > 0:
        eta = _window_efficiency(window, entropy, mol.any(), atm.any(), table)
        ideal = 100.0 * cap_sites * eta / (M_C_AMU + cap_sites * eta)
        thermo = (
            100.0
            * float(np.sum(table.sites_per_C * table.n_max * dtheta * entity_mass))
            / (M_C_AMU + dmass + m_charged)
        )
        occ_def = max(ideal - thermo, 0.0)
        adm_def = max(thermo - 100.0 * m_del / (M_C_AMU + m_charged), 0.0)
        mass_def = max(100.0 * m_del / (M_C_AMU + m_charged) - D, 0.0)
    else:
        occ_def = adm_def = mass_def = 0.0
    return CapacityResult(D, charged, h_del, dtheta, adm, occ_def, adm_def, mass_def)


def _window_efficiency(
    window: OperatingWindow, entropy: EntropyModel, has_mol: bool, has_atm: bool, table: SiteTable
) -> float:
    if has_mol:
        _, eta = optimum_energy(window, entropy)
        return float(eta)
    return 1.0


def window_efficiency(
    window: OperatingWindow, entropy: EntropyModel | None = None
) -> tuple[float, float]:
    """eta(W) = sup over eps of the molecular swing, and the optimizing eps* (DERIVATIONS 5.4)."""
    entropy = entropy or ConstantEntropy()
    eps_star, eta = optimum_energy(window, entropy)
    return eta, eps_star


def monolayer_capacity_mol_per_kg(
    areal_density_nm2: float = 10.0, ssa_m2_g: float = GRAPHENE_SSA_M2_G
) -> float:
    """Monolayer H2 capacity (mol kg^-1) for a given areal density and specific surface area."""
    molecules_per_g = areal_density_nm2 * 1.0e18 * ssa_m2_g
    return molecules_per_g / N_A * 1000.0


def bound_wt_percent(
    window: OperatingWindow,
    channel: str,
    entropy: EntropyModel | None = None,
    *,
    areal_density_nm2: float = 10.0,
    ssa_m2_g: float = GRAPHENE_SSA_M2_G,
    release_fraction: float = 1.0,
) -> dict[str, float]:
    """Variational bound D_max(W) (DERIVATIONS 5.4) for one channel.

    Molecular: N_max from the monolayer capacity on ``ssa_m2_g``; the bound is the supremum over
    the site energy of the mass-accounted deliverable fraction of a homogeneous saturated
    landscape, attained at ``eps_D`` close to the swing optimum ``eps_star`` (the two differ by
    the mass of the charged hydrogen in the denominator, about 1e-3 eV).
    Atomic: N_max = 1 H per C (graphane) and eta given by the kinetic release fraction
    ``release_fraction`` (thermodynamic swing is ~1 for saturated charging).
    """
    entropy = entropy or ConstantEntropy()
    if channel == "molecular":
        from scipy.optimize import minimize_scalar

        eta, eps_star = window_efficiency(window, entropy)
        n_mol_per_kg = monolayer_capacity_mol_per_kg(areal_density_nm2, ssa_m2_g)
        m_h_per_kg = n_mol_per_kg * M_H2_AMU / 1000.0  # kg H2 per kg C

        def deliverable_fraction(eps: float) -> float:
            # Mass-accounted deliverable fraction of a homogeneous saturated landscape (5.4):
            # the bound is its supremum over eps, which the mass of the charged hydrogen shifts
            # slightly off the swing optimum eps*.
            th_c = float(
                occupancy(eps, window.T_c, window.P_c, entropy, use_fugacity=window.fugacity)
            )
            th_r = float(
                occupancy(eps, window.T_r, window.P_r, entropy, use_fugacity=window.fugacity)
            )
            return m_h_per_kg * (th_c - th_r) / (1.0 + m_h_per_kg * th_c)

        res = minimize_scalar(
            lambda e: -deliverable_fraction(e),
            bracket=(eps_star - 0.05, eps_star, eps_star + 0.05),
            tol=1e-12,
        )
        eps_d = float(res.x) if res.success else float(eps_star)
        D = 100.0 * max(deliverable_fraction(eps_d), deliverable_fraction(float(eps_star)))
        return {
            "D_max_wt_percent": float(D),
            "eta": float(eta),
            "eps_star_eV": float(eps_star),
            "eps_D_eV": eps_d,
            "N_max_mol_per_kg": float(n_mol_per_kg),
        }
    # atomic
    m_h_per_C = M_H_AMU / M_C_AMU
    D = 100.0 * m_h_per_C * release_fraction / (1.0 + m_h_per_C)
    return {
        "D_max_wt_percent": float(D),
        "eta": float(release_fraction),
        "eps_star_eV": float("nan"),
        "N_max_H_per_C": 1.0,
    }


def random_landscape_bound_check(
    window: OperatingWindow,
    rng: np.random.Generator,
    n_trials: int = 1000,
    entropy: EntropyModel | None = None,
) -> bool:
    """Randomized verification that no molecular landscape exceeds the bound (acceptance A4)."""
    entropy = entropy or ConstantEntropy()
    bound = bound_wt_percent(window, "molecular", entropy)
    n_sites_max = monolayer_capacity_mol_per_kg() / 1000.0 * M_C_AMU  # sites per C equivalent
    for _ in range(n_trials):
        k = int(rng.integers(1, 12))
        eps = rng.uniform(-1.0, -0.01, size=k)
        w = rng.dirichlet(np.ones(k)) * n_sites_max
        table = SiteTable(
            eps,
            np.full(k, np.nan),
            np.full(k, 1e13),
            np.ones(k, dtype=int),
            w,
            np.ones(k),
            np.zeros(k),
            np.zeros(k, dtype=int),
        )
        res = deliverable(table, window, entropy)
        if res.D_wt_percent > bound["D_max_wt_percent"] * (1 + 1e-9) + 1e-12:
            return False
    return True


__all__ = [
    "SiteTable",
    "CapacityResult",
    "deliverable",
    "window_efficiency",
    "bound_wt_percent",
    "monolayer_capacity_mol_per_kg",
    "mass_fraction",
    "random_landscape_bound_check",
]
