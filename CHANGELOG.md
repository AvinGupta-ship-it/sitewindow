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
  logs and manifests (`sitewindow landscape build`); verify-tier and full-tier GFN2-xTB landscapes
  executed (`data/processed/verify/`, 17 configurations; `data/processed/full/`, 125 configurations,
  no exclusions; both computed from breathing-optimized tube substrates).
- Tier-1 basal-monomer anchor and convergence-series workflow with Tier-2 pre-relaxation, LCAO
  relaxation and plane-wave single points, D3 reported separately (`sitewindow tier1`); the
  full-protocol 5x5 anchor is executed (`data/processed/tier1/`) and meets benchmark B1; the
  convergence series (6x6 cell, 6x6x1 k-mesh, 600 eV, 20 A vacuum) is executed and meets the B1
  convergence criterion, with the cell-size change at the 0.02 eV bound.
- Pairwise H2-C physisorption potentials (Steele-Buch and Steele-Darkrim-Levesque Lorentz-Berthelot
  sets, B3-calibrated set) with Feynman-Hibbs correction; sheet, pore, tube and bundle site
  energies (`sitewindow landscape physisorption`); parameter registry
  `data/registry/physisorption_potentials.yaml`.
- Barrier workflow (`sitewindow landscape barriers`): constrained relaxation along the reaction
  coordinate on a 0.1 A grid refined to 0.02 A, and climbing-image NEB with seven images and
  resumable bands, for monomer sticking, ortho- and para-dimer recombination, dangling-bond pair
  recombination at the monovacancy and ortho-dimer recombination on armchair nanotubes at two
  diameters; tables `barriers_<tier>.parquet` with reaction profiles and transition structures.
  `relax` accepts ASE constraints. Recombination bands start from the scan geometries, a failed
  self-consistent field is retried with a damped mixer, and calculations that still fail are
  excluded with their reason instead of aborting the run. At Tier 2 the GFN2-xTB self-consistent
  field does not converge at the hydrogen-abstraction saddle of dimer recombination, so Tier-2
  recombination barriers are constrained-scan maxima flagged as unconverged. Every method result
  is checked against its product (the separated H atom for sticking, H2 for recombination): a
  value is a saddle-point energy only when every point or image converged, the path ends in the
  product, the maximum is interior, the scan is continuous and, for sticking, the far end of the
  path reproduces the separated reactants, and the table's `converged` flag is that test; bands
  whose scan did not reach H2 run to an H2 molecule relaxed above the pair (its centre of mass
  held when it would fall back, and the desorbed molecule 6 A above the pair when it dissociates
  at 3 A) from an IDPP start; a band whose final state is not the product is not run, a band with
  an image whose energy did not converge is recorded as a failure without a barrier value, and an
  unconverged band resumes once with a fresh step budget. Runs of scan points that converged to a self-consistent-field solution
  above the ground state of their neighbours are excluded from the path. Sticking pathways take
  their own cell and calculator settings (`barriers.stick` in the tier configuration): 300 K in
  the 10x10 cell at Tier 2, where the far hydrogen reproduces the separated fragments, since
  Fermi smearing at 3000 K displaces it by 0.5 eV in every cell. The schema gains the reference
  state `adsorbed_state` for release barriers (schema version 0.2.0).
- Radial breathing-mode minimization of clean nanotube substrates before their relaxation
  (`relax_tube_radius`), recorded in the landscape and barrier manifests: a wide tube rolled at
  the graphene bond length meets the per-atom force criterion unrelaxed while its hydrogenated
  configurations release the collective strain, which shifted tube binding energies by up to
  0.34 eV in the first full-tier series; at Tier 3 the unrelaxed tubes fail the force criterion
  and the step changes the verify-tier values by at most 0.004 eV.
- Tier-3 (MACE-MP-0 medium) verify-tier landscape and barrier tables under
  `data/processed/verify/tier3/`, computed on the Tier-2 configuration set; the B1 runner reports
  the Tier-3 basal monomer and the B2 runner reports Tier-3 barriers and marks values that are
  path maxima rather than saddle points. The Tier-3 tube recombination bands do not converge
  within twice the step budget and are recorded as flagged path maxima.
- Validation runners for B1 (Tier-1 anchor, Tier-1 convergence series and Tier-3 report), B2 (sticking barrier
  per tier and method), B3 (pair potential on coronene) and B4 (computed basal and pore site
  energies).
- `sitewindow validate --markdown` renders the human-readable validation summary from the
  machine-readable report (`validation/report.py`).

### Changed
- Tier-2 protocol: Fermi smearing at an electronic temperature of 3000 K and 12x12 (full) or
  9x9 (verify) graphene cells for every sheet class, after the convergence series showed that
  the Gamma-only sampling of tblite makes chemisorption energies oscillate by up to 1 eV between
  cell sizes that fold the Dirac point onto Gamma and those that do not. Calculator settings are
  read from the tier configuration (`calculator_settings_<engine>`); armchair nanotubes use six unit
  cells and zigzag tubes four. The Tier-2 B1 report reads the landscape table and names the cell
  and the smearing; the verify-tier landscape is recomputed and the full-tier landscape computed under
  the revised protocol.
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
