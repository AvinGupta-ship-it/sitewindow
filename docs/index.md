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
| `sitewindow.atomistic` | structure builders for every site class, site descriptors, calculator adapters (tblite T2, GPAW T1, MACE T3), relaxation and single-point workflows with content-addressed caching, landscape construction, Tier-1 anchors and convergence, pairwise physisorption potentials |
| `sitewindow.landscape` | site-class catalogue: recipes, enumerated sites, multiplicities, mass changes and hydrogen configurations per tier |
| `sitewindow.schema` | pydantic models of every record type with provenance, uncertainty and physical range checks; flat record tables |
| `sitewindow.io` | Parquet tables carrying the schema version; SHA-256 manifests |
| `sitewindow.validation` | benchmark runners reading thresholds from the registry |
| `sitewindow.workflows` | command-line interface |

## Interpretation of outputs

- `theory_window_tables`: for each window and entropy model, the optimal site energy ε* (kJ mol⁻¹ per H₂), the window efficiency η (maximum fraction of site capacity deliverable), the molecular-channel bound D_max (wt%) for a monolayer on 2630 m² g⁻¹, and the heterogeneity threshold ε_c for a −5.8 kJ mol⁻¹ reference site.
- `validation_report.json`: one record per benchmark with predicted and reference values, the metric, and pass/fail against the registry threshold; `passed: null` marks runners that report without a pass criterion at the current tier.
- `landscape_raw_<tier>.parquet`: one row per energy record of the atomistic landscape (binding energies per H with the ½H₂ and the free-atom references, sequential binding energies, defect formation energies) with tier, engine, settings hash, descriptors of the site carbon and provenance; `landscape_configurations_<tier>.parquet` holds the relaxed structures (extended XYZ) and spin multiplicities, `landscape_raw_<tier>_excluded.csv` the non-converged configurations with reasons, and `landscape_manifest_<tier>.json` the calculator settings, reference energies and file hashes.
- `physisorption_<tier>.parquet`: strongest H₂ site energy (eV per H₂, negative for binding) of each physisorption geometry for each pair-potential parameter set, classically and with the Feynman–Hibbs correction at 77, 233 and 298 K.
- `tier1/tier1_anchor_basal_<n>x<n>.json`: Tier-1 basal-monomer anchor with the PBE and D3 parts of the binding energy reported separately, the relaxation record and the comparison with benchmark B1; `tier1_convergence_basal_<n>x<n>.json` the convergence variants of acceptance criterion B1.
- Verification-tier outputs (`data/processed/verify/`) use reduced configurations and are not evidence for publication claims.

## Running the atomistic workflows

```bash
sitewindow landscape build --tier verify            # GFN2-xTB landscape of the verify-tier class set
sitewindow landscape physisorption --tier verify    # pair-potential physisorption set
mpirun -np 4 gpaw python --module sitewindow.workflows.cli tier1 anchor --supercell 5   # Tier-1 anchor (conda environment)
mpirun -np 4 gpaw python --module sitewindow.workflows.cli tier1 converge --supercell 5 # convergence series
```

Tier-1 runs use the `gpaw python` wrapper under `mpirun`; every rank computes and rank 0 writes. Results are cached per structure, calculator settings and spin state, so interrupted runs resume from the cache.

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

- Tier-2 (GFN2-xTB) chemisorption energies on the basal plane are offset by about +0.5 eV from the diffusion Monte Carlo reference and by about +0.6 eV from the Tier-1 PBE value; the calibration component that models these offsets is not yet implemented. The Tier-1 PBE value itself lies 0.09 eV below the diffusion Monte Carlo reference, consistent with the Tier-0 offset the manual anticipates.
- The Debye-smeared knock-on model reproduces the static onset and the 100 keV cross section within a factor of two but underestimates the published 80 keV tail by about two orders of magnitude; alternative smearing kinematics remain to be compared.
- Corpus inference is limited to data harmonization; isosteric heats and adsorption-energy-distribution inversion are not yet implemented.
- The single-site Lennard-Jones H₂–C potential cannot satisfy benchmarks B3 and B4 jointly: the literature parameter sets underestimate the CCSD(T) coronene binding by about 1.4 kJ mol⁻¹, and the set calibrated to B3 places the 0.7 nm slit-pore site above the isosteric-heat interval of B4. Both parameterizations are carried as model-form alternatives; their spread (about 45 % of the site energy) is the physisorption model-form uncertainty until an anisotropic or many-body potential is adopted.
- The Tier-2 landscape at the verify tier covers 17 configurations (4×4 cells, three nanotube diameters); the full-tier enumeration (125 configurations in 6×6 and 7×7 cells, eight diameters per chirality, edges, oxygen groups, mixed classes), the Tier-1 convergence series and the barrier calculations have not been executed.
