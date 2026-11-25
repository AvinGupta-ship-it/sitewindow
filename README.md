# SiteWindow

Site-energy landscapes and the deliverable-capacity limit of defect-engineered sp² carbon for hydrogen storage.

SiteWindow is a scientific Python package for constructing, calibrating, and evaluating the
landscape of hydrogen binding sites in defective, curved, and oxygen-functionalized sp² carbon
against the thermodynamic, kinetic, and stoichiometric constraints of a storage cycle. It
implements a window theory of deliverable capacity for energetically heterogeneous adsorbents
(molecular and dissociative channels, kinetic admissibility, mass accounting, and a variational
bound), multi-fidelity atomistic workflows for site energies and barriers, electron knock-on
and defect-population models with Raman observables, forward thermal-desorption models, and
Bayesian inference of adsorption-energy distributions from the public NIST/ARPA-E isotherm
corpus.

The complete scientific and engineering specification is `RESEARCH_MANUAL.md`, with linked
specifications in `docs/specs/` (source registry, benchmark registry, data schema, derivations,
acceptance criteria).

## Installation

```bash
pip install .                 # core: theory, kinetics, defects, corpus tools
pip install ".[xtb]"          # GFN2-xTB tier (tblite)
pip install ".[mace]"         # universal interatomic potential tier (CPU)
pip install ".[bayes,uq]"     # PyMC/ArviZ inference and Sobol sensitivity
pip install ".[all,dev]"      # everything plus the test and documentation toolchain
```

Tier-1 density-functional anchors use GPAW (conda-forge: `conda env create -f environment.yml`)
or, where GPAW cannot be installed, PySCF (`pip install ".[pyscf]"`).

## Quick start

```bash
sitewindow theory-tables       # window efficiencies, optimum and threshold energies, bounds
sitewindow fetch-isodb         # clone the NIST ISODB mirror and harmonize carbon H2 isotherms
sitewindow knockon             # electron knock-on cross sections for graphene
sitewindow verify-atomistic    # GFN2-xTB verification of H binding on graphene and at a vacancy
sitewindow validate            # run the benchmark runners and write validation_report.json
```

```python
from sitewindow.theory.window import W1, optimum_energy, swing
from sitewindow.theory.deliverable import bound_wt_percent

eps_star, eta = optimum_energy(W1)       # optimal site energy (eV per H2) and window efficiency
bound_wt_percent(W1, "molecular")        # variational bound D_max(W1) for the molecular channel
```

## Implementation status

| Component (manual section) | Status |
|---|---|
| C1 Window theory, deliverable functional, bound (`theory/`) | implemented; analytic and randomized tests pass; benchmark B8 passes |
| C2 Atomistic landscape (`atomistic/`) | structure builders, tier-2 calculator adapter, relaxation and binding-energy workflows with caching implemented; verification-tier runs executed; full class enumeration and tier-1/tier-3 engines pending |
| C3 Registry, meta-analysis, calibration (`landscape/`) | pending |
| C4 Kinetics and forward desorption (`kinetics/`) | Polanyi–Wigner, Redhead, ramps, admissibility implemented and tested; multi-species pathways pending |
| C5 Defect populations and Raman (`defects/`) | knock-on with Debye smearing, rate equations, Raman relations implemented and tested; benchmark B9 endpoint check fails at 80 keV (model-form difference under investigation), B10 passes |
| C6 Corpus inference (`experimental/`) | ISODB loader, carbon selection, unit normalization, harmonized table implemented and executed (306 isotherms); isosteric heats and AED inversion pending |
| C7 Evaluation, sensitivity, hypothesis tests (`evaluate/`) | pending |
| C8 Validation runners (`validation/`) | B8, B9 (endpoints), B10, B1 (tier-2 report) implemented; remaining runners pending |

The verification-tier outputs in `data/processed/verify/` are labeled as such and do not
support publication claims.

## Documentation

`docs/` contains the specifications; the API reference uses MkDocs
(`mkdocs serve`).

## Citation

See `CITATION.cff`.

## License

MIT License; derived datasets under CC BY 4.0 (see `data/LICENSE-DATA`).
