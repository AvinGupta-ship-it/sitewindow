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
sitewindow landscape build     # relax every catalogue configuration of a tier (GFN2-xTB by default)
sitewindow landscape physisorption  # pair-potential H2 site energies for the physisorption set
sitewindow landscape barriers       # sticking and recombination barriers along the designated pathways
sitewindow tier2 smearing      # free-energy, internal-energy and extrapolated Tier-2 binding energies by electronic temperature
sitewindow validate            # run the benchmark runners and write validation_report.json
```

Tier-1 anchors run in the conda-forge environment under MPI; the anchor set over the registry
`data/registry/tier1_anchors.yaml` is checkpointed per structure and per anchor and resumes
where it stopped:

```bash
mpirun -np 4 gpaw python --module sitewindow.workflows.cli tier1 anchor --supercell 5
mpirun -np 5 gpaw python --module sitewindow.workflows.cli tier1 anchors
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
| C2 Atomistic landscape (`atomistic/`, `landscape/classes.py`) | builders for every catalogue class, descriptors, calculator adapters for tblite (T2), GPAW (T1, PBE-D3(BJ), MPI) and MACE (T3), site-class catalogue (15 chemisorption and 4 physisorption classes; 17 verify-tier and 125 full-tier configurations), landscape workflow writing schema-validated tables; verify-tier Tier-2 landscape executed under the revised Tier-2 protocol (Fermi smearing at 3000 K, 9×9 sheets, six-cell armchair tubes, adopted after the convergence series exposed a Γ-only finite-size error of up to 1 eV; clean tubes minimized along their radial breathing coordinate before relaxation; 17 configurations, no exclusions); full-tier Tier-2 landscape (12×12 sheets, eight diameters per chirality with four-cell zigzag tubes, edges, oxygen groups, mixed classes; 125 configurations, no exclusions) executed from breathing-optimized tubes, with the tube exterior/interior mean within 0.06 eV of the sheet monomer across 0.5–2.5 nm; Tier-3 (MACE-MP-0) landscape and barriers executed for the verify set, binding about 1 eV more strongly than Tier 2 on sp² sites; pairwise physisorption potentials with Feynman–Hibbs correction executed for the verify set; Tier-1 anchor and convergence workflows implemented; the full-protocol basal-monomer anchor (5×5 cell, 500 eV, 4×4×1 k-mesh, 15 Å vacuum) is executed and meets the B1 criterion (PBE −0.781 eV against −0.821 eV); the Tier-1 convergence series (6×6 cell, 6×6×1 k-mesh, 600 eV cutoff, 20 Å vacuum) is executed and meets the B1 convergence criterion, the largest change being −0.020 eV for the 6×6 cell against the 0.02 eV bound; barrier workflow (constrained scans along a C–H distance, climbing-image NEB seeded from the scan or from an IDPP path to H₂ above the pair, saddle validation of every result, SCF retry and exclusion logging) executed for the verify-tier pathway set at Tier 2 and Tier 3: basal sticking 0.25 eV at Tier 2 (300 K, 10×10 cell; scan 0.252, NEB 0.251 eV) and 0.24 eV at Tier 3 (0.238/0.235 eV), both within the B2 band that is decided at Tier 1 only; para-dimer recombination 1.34 eV (Tier 2) and 2.18 eV (Tier 3); ortho-dimer recombination 1.12 eV (Tier 3; the Tier-2 band fails at the abstraction saddle and its scan maximum of 1.64 eV is flagged); ortho recombination on the (5,5) and (10,10) tubes 2.46 and 2.62 eV (Tier 2, located saddles; the Tier-3 bands did not converge within twice their step budget: the (5,5) band ends at a flagged path maximum of 1.72 eV after 1200 steps and the (10,10) climbing image leaves the recombination channel for the two separated atoms at 4.50 eV, so the Tier-3 tube barriers are pending a band seeded from the Tier-2 transition structure); H₂ release from the vacancy dangling-bond pair without a located saddle at either tier (flagged scan maxima 2.22 and 2.50 eV as lower bounds); the full-tier physisorption set executed (312 records: seven slit widths, seven tube diameters with bundle grooves and interstitials at three, the basal and vacancy-adjacent sheets and the sites adjacent to the epoxide, hydroxyl, carbonyl-pair and carboxyl groups through OPLS-AA oxygen parameters, which strengthen the Steele–Buch site energy from −3.9 kJ mol⁻¹ on the basal plane to −4.2 to −4.9 kJ mol⁻¹ beside the groups); the Tier-1 anchor set of 47 configurations spanning every chemisorption class (`data/registry/tier1_anchors.yaml`, PBE-D3(BJ) at the converged protocol, checkpointed per structure and per anchor, `sitewindow tier1 anchors`) implemented and running; the electronic-temperature study of the Tier-2 protocol (`sitewindow tier2 smearing`) executed: tblite reports the Mermin free energy, the far-hydrogen offset of the sticking reference is not a smearing artifact (internal-energy offset +0.8 to +1.0 eV in 3n cells at every temperature), and the zero-temperature extrapolation at 6000 K converges the basal-monomer, vacancy and epoxide binding energies across the cell families within 0.03 eV where the 3000 K free energies spread by 0.1 eV, so the Tier-2 protocol is revised to E₀ = F + TS/2 at 6000 K (numerical uncertainty 0.03 eV, meeting the 0.05 eV target) and both tiers of the Tier-2 landscape are being recomputed under it (the 3000 K tables remain as first-protocol calibration input with 0.1 eV uncertainty); NEB bands can be seeded from the transition structure another tier located; the Tier-2 sticking-cell convergence check executed (11×11 cell: scan 0.244 eV, band 0.242 eV against 0.252 and 0.251 eV in the 10×10 cell, converged within 0.01 eV); Tier-1 barrier anchors, the Tier-3 tube recombination bands seeded from the Tier-2 saddles, the NEB image-count convergence, full-tier barriers and the full-tier Tier-3 landscape pending |
| C3 Registry, meta-analysis, calibration (`landscape/`) | pending |
| C4 Kinetics and forward desorption (`kinetics/`) | Polanyi–Wigner, Redhead, ramps, admissibility implemented and tested; multi-species pathways pending |
| C5 Defect populations and Raman (`defects/`) | knock-on with Debye smearing, rate equations, Raman relations implemented and tested; benchmark B9 endpoint check fails at 80 keV (model-form difference under investigation), B10 passes |
| C6 Corpus inference (`experimental/`) | ISODB loader, carbon selection, unit normalization, harmonized table implemented and executed (306 isotherms); isosteric heats and AED inversion pending |
| C7 Evaluation, sensitivity, hypothesis tests (`evaluate/`) | pending |
| C8 Validation runners (`validation/`) | B8, B9 (endpoints), B10, B1 (tier-2 and tier-3 reports from the landscape tables, tier-1 anchor and tier-1 convergence), B2 (sticking barrier per tier and method, path maxima flagged), B3 (pair potential), B4 (computed site energies), B11 (monovacancy and Stone–Wales formation energies per tier from the landscape and anchor tables; the Tier-2 values pass) implemented; remaining runners pending |
| Data schema and I/O (`schema/`, `io/`) | pydantic models of every record type with provenance and physical range checks, flat Parquet tables with schema-version metadata, JSON Schema export implemented and tested |

The verification-tier outputs in `data/processed/verify/` are labeled as such and do not
support publication claims.

## Documentation

`docs/` contains the specifications; the API reference uses MkDocs
(`mkdocs serve`).

## Citation

See `CITATION.cff`.

## License

MIT License; derived datasets under CC BY 4.0 (see `data/LICENSE-DATA`).
