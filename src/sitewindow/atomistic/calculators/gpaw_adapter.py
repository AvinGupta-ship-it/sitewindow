"""GPAW density-functional calculator adapter: fidelity tier T1 (manual C2).

Settings follow the Tier-1 protocol: PBE with D3(BJ) dispersion, spin polarization, LCAO dzp
mode for relaxations and plane-wave single points for energies, a k-point mesh scaled from
4x4x1 on the 5x5 graphene cell, and the dipole-layer correction for slabs. The dispersion
energy is added through the ``dftd3-python`` ASE calculator so that the PBE and D3 parts can be
reported separately (benchmark B1).

Parallel execution uses ``mpirun -np N gpaw python --module sitewindow.workflows.cli ...``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

from ase import Atoms

from sitewindow.constants import GRAPHENE_LATTICE_A
from sitewindow.errors import CalculatorUnavailable

K_REF = 4
L_REF_A = 5 * GRAPHENE_LATTICE_A
ROBUST_MIXER: dict[str, Any] = {"name": "pulay", "beta": 0.02, "nmaxold": 8, "weight": 100.0}


def scaled_kmesh(
    atoms: Atoms, k_ref: int = K_REF, length_ref_A: float = L_REF_A
) -> tuple[int, int, int]:
    """k-mesh with the same reciprocal-space density as ``k_ref`` points per ``length_ref_A``."""
    lengths = atoms.cell.lengths()
    mesh = []
    for axis in range(3):
        if atoms.pbc[axis] and lengths[axis] > 0:
            mesh.append(max(1, int(round(k_ref * length_ref_A / lengths[axis]))))
        else:
            mesh.append(1)
    return mesh[0], mesh[1], mesh[2]


@dataclass
class GPAWCalculator:
    """Tier-1 GPAW adapter implementing the :class:`Calculator` protocol."""

    mode: str = "lcao"
    basis: str = "dzp"
    ecut_eV: float = 500.0
    xc: str = "PBE"
    dispersion: str | None = "d3bj"
    k_ref: int = K_REF
    kmesh: tuple[int, int, int] | None = None
    spin_polarized: bool = True
    dipole_correction: bool = True
    smearing_eV: float = 0.05
    grid_spacing_A: float = 0.18
    convergence_energy_eV: float = 5.0e-4
    maxiter: int = 400
    txt: str | None = None
    name: str = "gpaw"
    tier: str = "T1"
    extra: dict[str, Any] = field(default_factory=dict)
    _report_settings: dict[str, Any] | None = field(default=None, repr=False, compare=False)

    def robust(self) -> GPAWCalculator:
        """The same Hamiltonian with a damped Pulay mixer, a heavier long-wavelength metric and
        the full iteration budget, used to retry a self-consistent field that sloshes (ribbons and
        tubes with a large vacuum cross-section); ``settings`` is unchanged so that the result is
        cached under the original key."""
        return replace(
            self,
            maxiter=max(self.maxiter, 400),
            extra={**self.extra, "mixer": dict(ROBUST_MIXER)},
            _report_settings=self.settings(),
        )

    def kpoints(self, atoms: Atoms) -> tuple[int, int, int]:
        if self.kmesh is not None:
            return self.kmesh
        if atoms.info.get("geometry") == "molecule" or not any(atoms.pbc):
            return (1, 1, 1)
        return scaled_kmesh(atoms, self.k_ref)

    def attach(self, atoms: Atoms) -> None:
        try:
            from gpaw import GPAW, PW, FermiDirac
        except ImportError as e:  # pragma: no cover - optional dependency
            raise CalculatorUnavailable(
                "gpaw is not installed; create the conda-forge environment in environment.yml"
            ) from e
        kwargs: dict[str, Any] = dict(
            xc=self.xc,
            kpts={"size": self.kpoints(atoms), "gamma": True},
            occupations=FermiDirac(self.smearing_eV),
            spinpol=self.spin_polarized,
            convergence={"energy": self.convergence_energy_eV},
            maxiter=self.maxiter,
            txt=self.txt,
        )
        if self.mode == "pw":
            kwargs["mode"] = PW(self.ecut_eV)
        elif self.mode == "lcao":
            kwargs["mode"] = "lcao"
            kwargs["basis"] = self.basis
            kwargs["h"] = self.grid_spacing_A
        else:
            raise ValueError("mode must be 'lcao' or 'pw'")
        slab = bool(atoms.pbc[0] and atoms.pbc[1] and not atoms.pbc[2])
        if self.dipole_correction and slab:
            kwargs["poissonsolver"] = {"dipolelayer": "xy"}
        kwargs.update(self.extra)
        # Initial magnetic moments are the caller's responsibility (set_unpaired_electrons):
        # overriding them here would start even-electron systems such as H2 in a high-spin state.
        gpaw = GPAW(**kwargs)
        if self.dispersion is None:
            atoms.calc = gpaw
            return
        atoms.calc = _with_dispersion(gpaw, self.xc, self.dispersion)

    def settings(self) -> dict[str, Any]:
        if self._report_settings is not None:
            return dict(self._report_settings)
        s: dict[str, Any] = {
            "mode": self.mode,
            "xc": self.xc,
            "dispersion": self.dispersion,
            "k_ref": self.k_ref,
            "kmesh": None if self.kmesh is None else list(self.kmesh),
            "spin_polarized": self.spin_polarized,
            "dipole_correction": self.dipole_correction,
            "smearing_eV": self.smearing_eV,
            "convergence_energy_eV": self.convergence_energy_eV,
            **self.extra,
        }
        if self.mode == "pw":
            s["ecut_eV"] = self.ecut_eV
        else:
            s["basis"] = self.basis
            s["grid_spacing_A"] = self.grid_spacing_A
        return s


def _with_dispersion(base: Any, xc: str, dispersion: str) -> Any:
    try:
        from ase.calculators.mixing import SumCalculator
        from dftd3.ase import DFTD3
    except ImportError as e:  # pragma: no cover - optional dependency
        raise CalculatorUnavailable("dftd3-python is not installed") from e
    if dispersion != "d3bj":
        raise ValueError("only D3(BJ) dispersion is supported")
    return SumCalculator([base, DFTD3(method=xc, damping="d3bj")])


def dispersion_energy_eV(atoms: Atoms, xc: str = "PBE", damping: str = "d3bj") -> float:
    """Stand-alone D3 dispersion energy (eV) of a structure, for separate reporting."""
    try:
        from dftd3.ase import DFTD3
    except ImportError as e:  # pragma: no cover - optional dependency
        raise CalculatorUnavailable("dftd3-python is not installed") from e
    a = atoms.copy()
    a.calc = DFTD3(method=xc, damping=damping)
    return float(a.get_potential_energy())


__all__ = [
    "GPAWCalculator",
    "ROBUST_MIXER",
    "scaled_kmesh",
    "dispersion_energy_eV",
    "K_REF",
    "L_REF_A",
]
