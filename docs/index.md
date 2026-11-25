# SiteWindow

SiteWindow evaluates the landscape of hydrogen binding sites in defective, curved and
oxygen-functionalized sp² carbon against the thermodynamic, kinetic and stoichiometric
constraints of a storage cycle. The scientific program, components, data sources, benchmarks and
acceptance criteria are specified in the research manual and the linked specifications.

## Components

| Subpackage | Content |
|---|---|
| `sitewindow.theory` | Langmuir and dissociative occupancy, swing functions, entropy models, optimum and threshold energies, deliverable-capacity functional, variational bound |
| `sitewindow.kinetics` | Polanyi–Wigner kinetics, Redhead relations, linear-ramp desorption, kinetic admissibility |
| `sitewindow.defects` | McKinley–Feshbach knock-on cross sections with Debye-model smearing, defect rate equations, Raman defect relations |
| `sitewindow.experimental` | NIST ISODB mirror loader, carbon selection, unit normalization |
| `sitewindow.atomistic` | structure builders, calculator adapters, relaxation and binding-energy workflows with caching |
| `sitewindow.validation` | benchmark runners reading thresholds from the registry |
| `sitewindow.workflows` | command-line interface |

## Interpretation of outputs

- `theory_window_tables`: for each window and entropy model, the optimal site energy ε* (kJ mol⁻¹ per H₂), the window efficiency η (maximum fraction of site capacity deliverable), the molecular-channel bound D_max (wt%) for a monolayer on 2630 m² g⁻¹, and the heterogeneity threshold ε_c for a −5.8 kJ mol⁻¹ reference site.
- `validation_report.json`: one record per benchmark with predicted and reference values, the metric, and pass/fail against the registry threshold; `passed: null` marks runners that report without a pass criterion at the current tier.
- Verification-tier outputs (`data/processed/verify/`) use reduced configurations and are not evidence for publication claims.

## Reference hardware

The reference workstation for this program is one Apple Silicon machine:

- Processor: Apple M4, 10 cores (4 performance and 6 efficiency cores), arm64.
- Memory: 24 GB unified memory.
- Operating system: macOS 15.6.1 (Darwin 24.6.0).
- Package environment: Python 3.12.16 in the project virtual environment (`.venv`) with the
  `xtb`, `bayes`, `uq` and `dev` extras; tblite 0.8.0 provides the Tier-2 GFN2-xTB engine.
- Tier-1 density-functional engine: GPAW 26.7.0 from conda-forge in the `sitewindow` conda
  environment defined by `environment.yml` (Open MPI 5.0.11, OpenMP, libxc 7.1.2, libvdwxc, FFTW,
  ScaLAPACK, ELPA, PAW datasets `gpaw_data` 1.0.1), with D3(BJ) dispersion from `dftd3-python`
  1.7.0. GPAW is the active Tier-1 engine; the PySCF fallback is not installed. The installation
  was verified by LCAO (dzp) calculations of H and H₂ and a plane-wave calculation of a periodic
  graphene cell.

## Known limitations at this version

- Tier-2 (GFN2-xTB) chemisorption energies on the basal plane are offset by about +0.5 eV from the diffusion Monte Carlo reference; the calibration component that models this offset is not yet implemented.
- The Debye-smeared knock-on model reproduces the static onset and the 100 keV cross section within a factor of two but underestimates the published 80 keV tail by about two orders of magnitude; alternative smearing kinematics remain to be compared.
- Corpus inference is limited to data harmonization; isosteric heats and adsorption-energy-distribution inversion are not yet implemented.
