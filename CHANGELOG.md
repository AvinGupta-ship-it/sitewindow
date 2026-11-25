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
