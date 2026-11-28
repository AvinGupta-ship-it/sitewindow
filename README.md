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
| C2 Atomistic landscape (`atomistic/`, `landscape/classes.py`) | builders for every catalogue class, descriptors, calculator adapters for tblite (T2), GPAW (T1, PBE-D3(BJ), MPI) and MACE (T3), site-class catalogue (15 chemisorption and 4 physisorption classes; 17 verify-tier and 125 full-tier configurations), landscape workflow writing schema-validated tables, with an optional second relaxation start from the minimum an earlier table found for the same configuration; verify-tier Tier-2 landscape executed under the revised Tier-2 protocol (Fermi smearing at 3000 K, 9×9 sheets, six-cell armchair tubes, adopted after the convergence series exposed a Γ-only finite-size error of up to 1 eV; clean tubes minimized along their radial breathing coordinate before relaxation; 17 configurations, no exclusions); full-tier Tier-2 landscape (12×12 sheets, eight diameters per chirality with four-cell zigzag tubes, edges, oxygen groups, mixed classes; 125 configurations, no exclusions) executed from breathing-optimized tubes, with the tube exterior/interior mean within 0.06 eV of the sheet monomer across 0.5–2.5 nm; Tier-3 (MACE-MP-0) landscape and barriers executed for the verify set, binding about 1 eV more strongly than Tier 2 on sp² sites; full-tier Tier-3 landscape executed (123 of 125 configurations, the fifth and sixth vacancy hydrogens excluded for non-convergence; class-dependent offsets against Tier 2 from −0.3 eV at the vacancy dangling bond to −1.0 eV at the Stone–Wales site and −2.1 eV at the carbonyl-pair hole carbon, defect formation energies 0.5 to 1.9 eV below Tier 2); pairwise physisorption potentials with Feynman–Hibbs correction executed for the verify set; Tier-1 anchor and convergence workflows implemented; the full-protocol basal-monomer anchor (5×5 cell, 500 eV, 4×4×1 k-mesh, 15 Å vacuum) is executed and meets the B1 criterion (PBE −0.781 eV against −0.821 eV); the Tier-1 convergence series (6×6 cell, 6×6×1 k-mesh, 600 eV cutoff, 20 Å vacuum) is executed and meets the B1 convergence criterion, the largest change being −0.020 eV for the 6×6 cell against the 0.02 eV bound; barrier workflow (constrained scans along a C–H distance, climbing-image NEB seeded from the scan or from an IDPP path to H₂ above the pair, saddle validation of every result, SCF retry and exclusion logging) executed for the verify-tier pathway set at Tier 2 and Tier 3: basal sticking 0.25 eV at Tier 2 (300 K, 10×10 cell; scan 0.252, NEB 0.251 eV) and 0.24 eV at Tier 3 (0.238/0.235 eV), both within the B2 band that is decided at Tier 1 only; para-dimer recombination 1.34 eV (Tier 2) and 2.18 eV (Tier 3); ortho-dimer recombination 1.12 eV (Tier 3; the Tier-2 band fails at the abstraction saddle and its scan maximum of 1.64 eV is flagged); ortho recombination on the (5,5) and (10,10) tubes 2.46 and 2.62 eV (Tier 2, located saddles; the Tier-3 bands do not converge: from an IDPP start the (5,5) band ends at a flagged path maximum of 1.72 eV after 1200 steps and the (10,10) climbing image leaves the recombination channel for the two separated atoms at 4.50 eV, while bands seeded from the Tier-2 transition structures stay in the channel but end unconverged after 600 steps at 2.06 and 1.29 eV, so the Tier-3 tube barriers remain flagged path maxima); H₂ release from the vacancy dangling-bond pair without a located saddle at either tier (flagged scan maxima 2.22 and 2.50 eV as lower bounds); the full-tier physisorption set executed (312 records: seven slit widths, seven tube diameters with bundle grooves and interstitials at three, the basal and vacancy-adjacent sheets and the sites adjacent to the epoxide, hydroxyl, carbonyl-pair and carboxyl groups through OPLS-AA oxygen parameters, which strengthen the Steele–Buch site energy from −3.9 kJ mol⁻¹ on the basal plane to −4.2 to −4.9 kJ mol⁻¹ beside the groups); the Tier-1 anchor set of 47 configurations spanning every chemisorption class (`data/registry/tier1_anchors.yaml`, PBE-D3(BJ) at the converged protocol, checkpointed per structure and per anchor, `sitewindow tier1 anchors`) implemented and running, with the symmetry-broken spin start adopted after the first monovacancy and zigzag-edge records showed non-magnetic plane-wave states (both anchors recomputed: antiferromagnetic ribbon edges, magnetic monovacancy with a formation energy inside the B11 interval); the Stone–Wales rotated-bond anchor executed (non-magnetic monomer, −2.34 eV PBE-D3; clean-defect formation energy 5.41 eV PBE-D3 in the 6×6 cell, 0.41 eV above the B11 reference, where the Tier-2 cell series of the same quantity falls by 0.29 eV from the 6×6 to the 9×9 cell, so the Tier-1 value in a 9×9 cell decides the benchmark and is queued as a round-2 anchor); the electronic-temperature study of the Tier-2 protocol (`sitewindow tier2 smearing`) executed: tblite reports the Mermin free energy, the far-hydrogen offset of the sticking reference is not a smearing artifact (internal-energy offset +0.8 to +1.0 eV in 3n cells at every temperature), and the zero-temperature extrapolation at 6000 K converges the basal-monomer, vacancy and epoxide binding energies across the cell families within 0.03 eV where the 3000 K free energies spread by 0.1 eV, so the Tier-2 protocol is revised to E₀ = F + TS/2 at 6000 K (cell-family spread 0.03 eV; the extrapolation itself changes by 0.06 to 0.13 eV between 4500 and 6000 K in the 12×12 cell, which is recorded as the class-level Tier-2 protocol uncertainty that the C3 calibration absorbs) and both tiers of the Tier-2 landscape are recomputed under it: the verify tier is complete (17 configurations, no exclusions; sp² sites move by 0.00 to 0.02 eV and the vacancy dangling bonds by up to 0.1 eV against the first-protocol table) and the full tier is complete (125 configurations, 368 records, no exclusions; basal monomer −0.819 eV against −0.821 eV PBE; sheet sp² and tube records within −0.08 to +0.03 eV of the first-protocol table, vacancy dangling bonds within −0.11 eV, formation energies +0.02 to +0.11 eV higher; three records had relaxed into different minima than the first-protocol run and were resolved by a second relaxation start from the first-protocol minimum, which lowered the (8,8) tube vacancy's first and second hydrogen by 1.87 and 0.25 eV per atom and confirmed the divacancy records); NEB bands can be seeded from the transition structure another tier located; the Tier-2 sticking-cell convergence check executed (11×11 cell: scan 0.244 eV, band 0.242 eV against 0.252 and 0.251 eV in the 10×10 cell, converged within 0.01 eV); the NEB image-count convergence executed at Tier 3 (basal sticking 0.235 eV with 5, 7 and 9 images, spread 0.001 eV); the Tier-3 tube recombination bands seeded from the Tier-2 saddles executed (unconverged, flagged); Tier-1 barrier anchors (monomer sticking, ortho and para recombination on the 5×5 sheet, ortho recombination on the (5,5) tube; LCAO scans with plane-wave single points on the stationary structures, `sitewindow tier1 barriers`) queued after round 1 of the anchor set; full-tier barriers pending |
| C3 Registry, meta-analysis, calibration (`landscape/`) | published-value registry of 63 records (54 site energies from twelve open-access sources, each with DOI or arXiv identifier, method metadata, reference state, extraction and the place where the number was read; the Boukhvalov–Katsnelson defect values registered but excluded from pooling for sign ambiguity); hierarchical meta-analysis (PyMC NUTS, site means nested in class means, method and study effects, cell-size term) executed per quantity and reference state with R̂ ≤ 1.005 and ESS ≥ 1078 in every group (DMC − PBE +0.134 ± 0.037 eV for the basal monomer); multi-fidelity calibration (pooled transfer, class random effect, Matérn-5/2 discrepancy over the site descriptors, Student-t widths, Tier-0 offset −0.129 eV from the DMC anchor) executed on the seven round-1 anchors with leave-one-anchor-out coverage 7 of 7, held-out curvature and oxygen classes covered and synthetic recovery 0.96 over 100 replicates; `landscape_calibrated.parquet` (1040 records, both tiers on the Tier-1 and Tier-0 scales) and the calibration report assembled by `sitewindow landscape calibrate` (seven classes offset-calibrated, eight awaiting anchors, none with its own slope yet; the Tier-3 route is less precise and does not transfer to the held-out oxygen classes); meta-analysis synthetic recovery executed (planted site means recovered in 96 % of 100 replicates); prior-scale study executed (site means move by at most 0.18 eV, 0.43 posterior standard deviations, across prior widths of 0.5 to 2 times the reference; method offsets stable to 0.003 eV) and kernel study executed (Matérn-3/2 leaves the calibration unchanged while no class holds two anchors), `sitewindow landscape sensitivity` |
| C4 Kinetics and forward desorption (`kinetics/`) | Polanyi–Wigner, Redhead, ramps, admissibility implemented and tested; multi-species pathways pending |
| C5 Defect populations and Raman (`defects/`) | knock-on with Debye smearing, rate equations, Raman relations implemented and tested; benchmark B9 endpoint check fails at 80 keV (model-form difference under investigation), B10 passes |
| C6 Corpus inference (`experimental/`) | ISODB loader with corrected carbon selection (metal–organic framework composites and alkali-intercalated carbons removed), unit normalization and harmonized table executed (280 isotherms, below the 300 of criterion (i) at the pinned mirror commit); material sets by adsorbent and publication with logged admissibility rules and a material metadata registry (simulated records, ambiguous names, the two verified BET areas); hydrogen equation of state (CoolProp); Clausius–Clapeyron isosteric heats with bootstrap uncertainty executed for 21 multi-temperature materials (cryogenic sets 3.3 to 7.5 kJ mol⁻¹ at low loading; room-temperature-only sets ill-conditioned); Bayesian adsorption-energy-distribution inversion with the measurement-error model, two-stage error-scale estimation, window-occupancy fractions and per-material checkpoints implemented, synthetic recovery and the corpus run in progress |
| C7 Evaluation, sensitivity, hypothesis tests (`evaluate/`) | pending |
| C8 Validation runners (`validation/`) | B8, B9 (endpoints), B10, B1 (tier-2 and tier-3 reports from the landscape tables, tier-1 anchor and tier-1 convergence), B2 (sticking barrier per tier and method, path maxima flagged), B3 (pair potential), B4 (computed site energies, and the corpus median of the low-loading isosteric heat, which fails as registered at 10.8 kJ mol⁻¹ over six activated carbons, four of them room-temperature-only sets, while the graphitic, cryogenic and all-carbon medians of 4.4 to 6.2 kJ mol⁻¹ lie inside the interval), B7 (surface-area scaling of the 77 K excess uptake; reported, not decided, with two materials of verified BET area where three are needed), B11 (monovacancy and Stone–Wales formation energies per tier from the landscape and anchor tables; the Tier-2 values and the Tier-1 monovacancy value pass; the Tier-1 Stone–Wales value in the 6×6 anchor cell, 5.41 eV PBE-D3, exceeds the 0.3 eV threshold by 0.11 eV and the benchmark fails at Tier 1 until the 9×9 Tier-1 value is computed) implemented; remaining runners pending |
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
