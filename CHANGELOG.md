# Changelog

All notable changes to this project are documented in this file. The format follows
Keep a Changelog, and the project adheres to Semantic Versioning.

## [Unreleased]

### Added
- Window theory: Langmuir and dissociative occupancy, swing functions, isothermal optimum and
  efficiency, heterogeneity threshold, entropy model set, fugacity fallback (`theory/`).
- Deliverable-capacity functional with mass accounting and the variational bound (`theory/deliverable.py`).
- Polanyi-Wigner kinetics, Redhead relations, linear-ramp desorption, admissibility (`kinetics/`).
- McKinley-Feshbach knock-on cross sections with Debye-model thermal smearing, defect rate
  equations, Raman defect relations (`defects/`).
- NIST ISODB mirror loader, carbon selection rules, unit normalization, harmonized carbon H2 table (`experimental/`).
- Structure builders for graphene supercells, vacancies, Stone-Wales, nanotubes; GFN2-xTB
  calculator adapter; relaxation and binding-energy workflows with content-addressed caching (`atomistic/`).
- Benchmark registry (`data/registry/benchmarks.yaml`) and runners for B8, B9 (endpoints), B10, B1 (tier-2 report).
- Command-line interface: `theory-tables`, `fetch-isodb`, `knockon`, `verify-atomistic`, `validate`.
- Research manual and linked specifications.
- Reference hardware description and active Tier-1 engine in `docs/index.md`; `dftd3-python` in
  `environment.yml` for PBE-D3(BJ) dispersion with GPAW.
- Canonical data schema as pydantic models with provenance and physical range checks, flat
  Parquet tables carrying the schema version, JSON Schema export (`schema/`, `io/`,
  `sitewindow schema export`); typed exception hierarchy (`errors.py`).
- Structure builders for divacancies (5-8-5 and 555-777 starting geometries), nanoribbons with
  zigzag and armchair edges, nanotubes with vacancies and Stone-Wales defects, oxygen groups
  (epoxide, hydroxyl, carbonyl pair, carboxyl), coronene and slit pores; hydrogen placement
  towards vacancy holes and ribbon edges; sublattice colouring and graphane placements.
- Site descriptors: pyramidalization angle (pi-orbital axis vector), coordination, neighbour
  hydrogen count, distance to the defect centre, curvature (`atomistic/descriptors.py`).
- GPAW Tier-1 adapter (PBE-D3(BJ), LCAO and plane-wave modes, scaled k-mesh, dipole-layer
  correction, spin polarization) and MACE Tier-3 adapter; calculator factory by name.
- Site-class catalogue with recipes, sites, multiplicities, mass changes and hydrogen
  configurations for the fifteen chemisorption and four physisorption classes
  (`landscape/classes.py`), extensible through the `sitewindow.site_classes` entry point.
- Landscape workflow writing `landscape_raw_<tier>.parquet`, configuration tables, exclusion
  logs and manifests (`sitewindow landscape build`); verify-tier GFN2-xTB landscape executed.
- Tier-1 basal-monomer anchor and convergence-series workflow with Tier-2 pre-relaxation, LCAO
  relaxation and plane-wave single points, D3 reported separately (`sitewindow tier1`); the
  full-protocol 5x5 anchor is executed (`data/processed/tier1/`) and meets benchmark B1; the
  convergence series (6x6 cell, 6x6x1 k-mesh, 600 eV, 20 A vacuum) is executed and meets the B1
  convergence criterion, with the cell-size change at the 0.02 eV bound.
- Pairwise H2-C physisorption potentials (Steele-Buch and Steele-Darkrim-Levesque Lorentz-Berthelot
  sets, B3-calibrated set) with Feynman-Hibbs correction; sheet, pore, tube and bundle site
  energies (`sitewindow landscape physisorption`); parameter registry
  `data/registry/physisorption_potentials.yaml`.
- Validation runners for B1 (Tier-1 anchor and Tier-1 convergence series), B3 (pair potential on
  coronene) and B4 (computed basal and pore site energies).
- `sitewindow validate --markdown` renders the human-readable validation summary from the
  machine-readable report (`validation/report.py`).

### Changed
- The benchmark registry carries the B1 convergence bound (0.02 eV) transcribed from the
  acceptance criteria, so the convergence runner reads its threshold from the registry.
- Verification-tier configuration lists the manual's reduced class set (basal, monovacancy,
  Stone-Wales, three armchair nanotubes) and the reduced Tier-1 anchor settings.
- `environment.yml` installs tblite from conda-forge (`tblite-python`) so that GPAW and tblite
  share one OpenMP runtime in the same process.
- `fetch-isodb` checks out the ISODB commit recorded in `data/sources.lock` instead of the mirror's
  current HEAD, so the corpus is pinned as specified in manual Section 6.2; `--update` moves the pin
  to the current HEAD. The lock file records the mirror URL together with the commit.

### Fixed
- Type annotations in `theory` and `kinetics` satisfy mypy strict with current numpy stubs: the
  `EntropyModel` protocol exposes `name` as a read-only property so that frozen dataclass models
  conform, array-returning functions declare their dtype, and `rate_constant` accepts array
  temperatures.
