# Research Manual

## Site-Energy Landscapes and the Deliverable-Capacity Limit of Defect-Engineered sp² Carbon for Hydrogen Storage

**Software framework:** SiteWindow (Python package `sitewindow`)

**Document role.** This manual is the authoritative scientific and engineering specification of the research program. Every component, equation, data source, benchmark, acceptance criterion, and deliverable required to execute the program is defined here or in the linked specifications under `docs/specs/`. The program is organized by scientific component and logical dependency, never by schedule.

**Linked specifications**

| Document | Content |
|---|---|
| `docs/specs/SOURCE_REGISTRY.md` | Verified public data and literature sources, identifiers, access mechanisms, licensing, and verification status |
| `docs/specs/BENCHMARK_REGISTRY.md` | Validation benchmarks with quantities, conditions, comparison metrics, thresholds, and independence from calibration |
| `docs/specs/DATA_SCHEMA.md` | Canonical data models, units, provenance, uncertainty representation, and serialization formats |
| `docs/specs/DERIVATIONS.md` | Mathematical derivations of the window theory, deliverable-capacity functional, bounds, kinetic admissibility, and inference models |
| `docs/specs/ACCEPTANCE_CRITERIA.md` | Consolidated, measurable acceptance criteria for every component and deliverable |

---

## 1. Project identity and scientific abstract

### 1.1 Title and naming

- **Project title:** Site-energy landscapes and the deliverable-capacity limit of defect-engineered sp² carbon for hydrogen storage.
- **Framework name:** SiteWindow. The name denotes the central object of the program: the *window* in site-energy space within which a hydrogen binding site contributes deliverable, reversible capacity under a specified operating cycle, and the *site-energy landscape* of defective carbon that is tested against that window.
- **Package name:** `sitewindow` (PyPI-style distribution name `sitewindow`). The name was checked against PyPI (no existing distribution) and GitHub (no scientific software of that name; unrelated web-development repositories named "SiteWindows" exist and do not conflict with a scientific Python package).
- The title names the object of investigation (a limit and the landscapes tested against it), not a result. The program is designed to determine the limit; it does not presuppose its value.

### 1.2 Scientific abstract

Hydrogen storage in sp² carbon is governed by a binding-energy landscape that is, to first approximation, bimodal: molecular physisorption at 4–8 kJ mol⁻¹, which is too weak for near-ambient storage, and atomic chemisorption through C–H bonds at 0.7–4 eV per atom, which is too strong and too slow to release. Defect engineering (vacancies, divacancies, Stone–Wales rearrangements, edges, curvature, and oxygen functionalities, produced for example by electron irradiation or by reduction of graphene oxide) is widely proposed as the route to intermediate, usable binding. Whether any physically realizable defect population actually places a substantial density of sites in the usable range has never been quantified with uncertainty, because the required ingredients exist only in isolation: scattered first-principles site energies with unquantified method dependence; a deliverable-capacity theory (Bhatia–Myers) that treats homogeneous molecular physisorption without kinetics; experimental isotherms and thermal-desorption data that have not been inverted into site-energy distributions at corpus scale; and no framework connecting defect populations, site energies, kinetic admissibility, and mass stoichiometry into a single deliverable-capacity functional.

This program constructs that framework. It (i) derives a window theory of deliverable capacity for energetically heterogeneous adsorbents with molecular and dissociative channels, kinetic admissibility constraints, and stoichiometric mass accounting, including a variational upper bound; (ii) builds a calibrated, uncertainty-quantified site-energy landscape for hydrogen on defective, curved, and oxygen-functionalized sp² carbon by multi-fidelity atomistic computation (tight-binding, universal interatomic potentials, local density-functional anchors) and hierarchical Bayesian meta-analysis of published high-level results; (iii) maps irradiation and reduction conditions to defect populations through displacement physics and rate equations, with Raman-observable defect density as the experimental coordinate; (iv) infers adsorption-energy distributions from the public hydrogen-isotherm corpus for carbons (NIST/ARPA-E ISODB) and compares predicted thermal-desorption spectra with published benchmarks; and (v) evaluates the deliverable-capacity functional across defect classes, densities, and operating windows to determine the attainable limit, the defect populations that approach it, and the circumstances under which the limit is or is not compatible with practical targets. The result is either a quantitative, falsifiable statement of a fundamental limit for carbon-only defect engineering, or the identification of specific site classes that break it, with uncertainty and with spectroscopic predictions by which the conclusion can be tested experimentally.

### 1.3 Central scientific question

> Over the complete landscape of hydrogen binding sites that defective, curved, and oxygen-functionalized sp² carbon can present, from molecular physisorption to atomic chemisorption, what is the maximum deliverable hydrogen capacity attainable under practical charge–discharge windows once thermodynamic, kinetic, and stoichiometric constraints are imposed jointly, and do any physically realizable defect populations, including those produced by electron irradiation and chemical reduction, approach that maximum?

### 1.4 Principal scientific contribution being investigated

A quantitative, uncertainty-bounded determination of the deliverable-capacity limit for defect-engineered sp² carbon, resolved by site class, defect density, and operating window, together with the identification of which site classes populate the usable window and which do not. The investigation is designed so that both outcomes are informative: a demonstrated limit below practical targets is a field-shaping negative result that redirects defect-engineering effort; a demonstrated window-populating site class is a field-shaping positive result with a concrete synthesis target and spectroscopic signature.

### 1.5 Expected categories of scientific outputs

1. A theory of deliverable capacity on heterogeneous adsorbents with kinetic constraints (closed-form results, a variational bound, and numerical maps).
2. A calibrated site-energy landscape dataset for hydrogen on sp² carbon defects with per-site posterior distributions and provenance.
3. Deliverable-capacity maps and bounds across defect classes, densities, and operating windows, with global sensitivity analysis.
4. Corpus-level adsorption-energy-distribution inference for carbons from public isotherms.
5. Forward-predicted thermal-desorption and Raman observables for defect populations, benchmarked against published measurements.
6. The `sitewindow` software package, DOI-ready releases, and preprint-ready manuscripts.

### 1.6 Relationship to hydrogen energy and carbon nanomaterials

The question is the central question of carbon-based solid-state hydrogen storage: can the binding landscape of carbon be engineered into the usable range without heteroatoms or metals? The materials (graphene, reduced graphene oxide, single-walled carbon nanotubes), the defect-generation routes (electron irradiation, chemical reduction, atomic-hydrogen exposure), and the observables (Raman D/G and D/D′ ratios, temperature-programmed desorption, hydrogen isotherms) are exactly those of experimental carbon hydrogen-storage research.

### 1.7 Distinction from existing work

- Bhatia and Myers (Langmuir 2006, 22, 1688) established the optimal homogeneous adsorption enthalpy for Langmuir physisorption delivery and examined heterogeneity for a two-site model without kinetic constraints, dissociative channels, defect chemistry, or stoichiometric mass accounting. SiteWindow generalizes this to arbitrary site-energy distributions over both channels, adds kinetic admissibility and mass accounting, derives a variational bound, and connects the distribution to defect populations.
- First-principles studies of hydrogen at carbon defects (for example Casolo et al. 2009; Boukhvalov and Katsnelson 2008; Dumi et al. 2022) report site energies for individual configurations with method-dependent spreads of tens of percent. SiteWindow treats these as data for a hierarchical meta-analysis and as anchors for a multi-fidelity landscape; it does not merely compute additional values.
- Experimental work on diameter-dependent hydrogenation of nanotubes (Nikitin et al. 2005, 2008) and on hydrogen on graphite (Zecho et al. 2002; Hornekær et al. 2006) provides benchmarks. SiteWindow predicts these observables from the landscape; it does not fit to them.
- Adsorption-energy-distribution inversion of isotherms is a classical technique; SiteWindow applies it at corpus scale, with a Bayesian measurement-error model that distinguishes tabulated from digitized data, to a question (window occupancy in carbons) it has not been used to answer.
- The previously released independent projects RamanUQ (Raman uncertainty quantification), tpd-sbi (Bayesian inversion of thermal-desorption spectra), H2STAR (thermodynamic hydrogen storage modeling), and HyCAN-DB (literature corpus) address different questions. SiteWindow uses only *forward* spectroscopic models (it predicts desorption spectra from the landscape; it does not invert measured spectra), does not model storage-system thermodynamics, and builds its data layer from public databases rather than a literature corpus. Interfaces by which those projects could consume SiteWindow outputs are defined in Section 10.6 without any dependency on them.

### 1.8 Consequential assumptions stated at the outset

- Operating windows are defined by the U.S. DOE onboard light-duty targets as published (minimum delivery pressure 5 bar, maximum delivery pressure 12 bar, delivery temperature −40 to 85 °C, fill time 3–5 min) and by common adsorbent test conditions (100 bar maximum charging pressure, 77 K cryogenic charging). Windows are configuration data, not hard-coded assumptions, and sensitivity to them is part of the program.
- The atomic (chemisorption) channel is charged by an external atomic-hydrogen source (irradiation-driven fragmentation of adsorbates, thermal cracking, or plasma), consistent with the experimental practice that motivates the program; charging by dissociation of H₂ on pure carbon is kinetically negligible at the conditions considered and is treated as such, with the barrier quantified rather than assumed.
- Pure carbon and oxygen-functionalized carbon are in scope; heteroatom dopants and metal decoration are out of scope for the core program and are the first extension direction.
## 2. Scientific case

### 2.1 The open problem

Solid-state hydrogen storage requires a material whose hydrogen binding free energy places a large fraction of capacity in the usable range between the charging and discharging states. For adsorptive storage of molecular hydrogen delivered between 30 bar and 1.5 bar at 298 K, Bhatia and Myers derived an optimal adsorption enthalpy of approximately 15 kJ mol⁻¹ under an entropy-change assumption of about −8R, and noted that carbons, with typical enthalpies near 6 kJ mol⁻¹, reach their optimum only near 115 K (Bhatia and Myers, Langmuir 2006, 22, 1688–1700). Estimates of the useful range are usually quoted as 15–25 kJ mol⁻¹ per H₂ (0.15–0.26 eV) for ambient operation and depend on the enthalpy–entropy correlation assumed.

Atomic hydrogen chemisorbs on the graphene basal plane with a binding energy of roughly 0.7–1.0 eV per atom relative to a free H atom (Sha and Jackson, Surf. Sci. 2002, 496, 318; Casolo et al., J. Chem. Phys. 2009, 130, 054704; diffusion Monte Carlo value −0.691 ± 0.019 eV versus PBE −0.821 eV, Dumi et al., J. Chem. Phys. 2022, 156, 144702), after a sticking barrier of about 0.2 eV and with substantial lattice puckering. Relative to ½ H₂, chemisorption of isolated atoms on the pristine basal plane is endothermic or nearly thermoneutral; clustering into ortho and para dimers and larger clusters stabilizes the chemisorbed state. Thermal desorption of H₂ from hydrogenated graphite shows two peaks near 445 K and 560 K (Zecho et al., J. Chem. Phys. 2002, 117, 8486), assigned to dimer species with different recombination pathways (Hornekær et al., Phys. Rev. Lett. 2006, 96, 156104). Hydrogenated single-walled carbon nanotubes release hydrogen at temperatures that depend on tube diameter, and C–H bonds are completely broken by 600 °C (Nikitin et al., Phys. Rev. Lett. 2005, 95, 225507; Nano Lett. 2008, 8, 162). Defects change the chemisorption energetics strongly: dangling bonds at monovacancies bind hydrogen by several electronvolts, Stone–Wales defects lower the chemisorption energy relative to perfect graphene over an extended neighborhood (Boukhvalov and Katsnelson, Nano Lett. 2008, 8, 4373), and curvature lowers both binding and desorption barriers (Tozzini and Pellegrini, Phys. Rev. B 2010, 81, 113404; Phys. Chem. Chem. Phys. 2013, 15, 80).

These facts define the problem. Between the physisorption mode (too weak) and the chemisorption mode (too strong, too slow to release, and not chargeable from H₂ without a catalyst or an atomic-hydrogen source) lies the usable window. Defect engineering is proposed to populate it, and electron irradiation is one practical means to create defects and, in the presence of adsorbates, to hydrogenate graphene (Jones et al., Carbon 2010, 48, 2335–2340). What has never been established is a quantitative answer to three coupled questions:

1. **Occupancy.** For each defect class and for realistic mixtures, what fraction of binding sites lies inside the usable window, with uncertainty that reflects the known spread of first-principles methods?
2. **Admissibility.** Of the sites in the thermodynamic window, which are kinetically admissible (charge within the fill time, release within the discharge time at the discharge temperature) given the barriers, the recombinative (second-order) nature of H₂ formation, and dimer- or cluster-mediated pathways?
3. **Accounting.** When the mass of the defect-bearing lattice, the stoichiometric limit on hydrogen per site, the stability limit on defect density, and the trade-off between site strength and site count are included, what deliverable capacity results, and how does it compare with practical targets?

### 2.2 Why the problem is significant to hydrogen energy

Carbon sorbents remain the most scalable, lowest-cost, lightest class of storage media, and the conclusion that defect engineering cannot reach the window would redirect a substantial research effort toward heteroatom doping, metal decoration, or system-level cryogenic solutions. Conversely, a demonstrated window-populating site class in pure carbon would be a synthesis target of unusual value, because it would be compatible with irradiation- or reduction-based processing of bulk powders. Either result changes how defect engineering is pursued. The question also underlies the long-standing reproducibility problems in carbon hydrogen storage: claims of enhanced ambient uptake in carbons have repeatedly failed to replicate, and a site-resolved theory with uncertainty provides the quantitative expectation against which such claims can be judged.

### 2.3 Why existing approaches have not resolved it

- First-principles studies compute individual configurations. Reported values for the same site differ by 20–40% between functionals, supercell sizes, and dispersion treatments (for the basal-plane monomer, published PBE/PW91 values range from about −0.79 to −0.98 eV, against −0.69 eV from diffusion Monte Carlo). No study has treated this spread as a quantity to be modeled, nor has any mapped configuration-level energies onto a density of sites per defect for a realistic population.
- The deliverable-capacity theory exists for a single Langmuir site energy and a two-site extension; the heterogeneous, two-channel, kinetically constrained generalization with stoichiometric accounting has not been derived.
- Thermal-desorption benchmarks have been modeled individually (for example kinetic Monte Carlo of H₂ desorption from graphite), but not predicted from a calibrated landscape and not used to constrain deliverable capacity.
- The public isotherm corpus for carbons has not been inverted into adsorption-energy distributions with a measurement-error model that distinguishes tabulated from digitized data, nor queried for window occupancy.
- Experimental irradiation work demonstrates hydrogenation and defect creation but cannot by itself separate site classes or determine what fraction of the resulting sites are usable.

### 2.4 Underlying physics and materials chemistry

The program rests on the following physics, each of which is modeled explicitly in Section 4.

- **Molecular physisorption** of H₂ on sp² carbon is dominated by dispersion; the H₂–coronene CCSD(T)/CBS reference binding is −4.7 kJ mol⁻¹ (Kocman et al., Phys. Chem. Chem. Phys. 2015, 17, 6423), and quantum delocalization of H₂ is significant below about 100 K (Feynman–Hibbs corrections). Curvature, slit-pore overlap of wall potentials, and polar oxygen functionalities raise the well depth by amounts of order 1–5 kJ mol⁻¹.
- **Atomic chemisorption** forms an sp³ C–H bond with lattice puckering, a sticking barrier, strong site dependence on local coordination (dangling bonds at vacancies), sublattice-dependent pairing (ortho/para dimers), and curvature dependence through pyramidalization. Desorption of H₂ is recombinative and generally second order in coverage with pathway-specific barriers.
- **Defect populations** created by electron irradiation follow knock-on displacement physics above a threshold near 80 keV for carbon in graphene, with cross sections measured directly for single-layer graphene (Meyer et al., Phys. Rev. Lett. 2012, 108, 196102) and reviewed for nanostructures (Krasheninnikov and Nordlund, J. Appl. Phys. 2010, 107, 071301); below the knock-on threshold, hydrogenation proceeds through fragmentation of adsorbates by backscattered and secondary electrons (Jones et al. 2010). Defect density is observable through the Raman D/G intensity ratio (Lucchese et al., Carbon 2010, 48, 1592; Cançado et al., Nano Lett. 2011, 11, 3190), and defect type through the D/D′ ratio (Eckmann et al., Nano Lett. 2012, 12, 3925).
- **Deliverable capacity** for a cycle between a charging state and a discharging state is the difference in equilibrium loading, restricted to sites whose kinetics permit charge and release within the specified times.

### 2.5 Limitations of previous investigations that this program addresses

| Limitation | Treatment in this program |
|---|---|
| Method-dependent site energies reported as point values | Hierarchical Bayesian meta-analysis with method random effects; multi-fidelity calibration with per-class discrepancy |
| Homogeneous or two-site deliverable-capacity theory | Functional over arbitrary site distributions with two channels, kinetic admissibility, and mass accounting; variational bound |
| Thermodynamic treatments of chemisorption without release kinetics | Explicit recombinative desorption kinetics with pathway barriers; admissibility defined by release fraction within a specified time |
| Isolated configuration energies without site counting | Site-density accounting per defect, including maximum hydrogen per defect and lattice mass |
| Defect density used without defect type | Defect-class-resolved populations from irradiation physics and from reduction chemistry; Raman D/D′ as a type-sensitive observable |
| Claims of enhanced binding from single isotherms | Corpus-scale Bayesian inversion with measurement-error modeling |

### 2.6 Scientific consequences of resolving the problem

- A quantitative limit on deliverable capacity for carbon-only defect engineering under each operating window, with posterior uncertainty.
- Identification of the defect classes that contribute usable sites and of the classes that only add mass or trap hydrogen irreversibly.
- A curvature–kinetics relation for release temperature of chemisorbed hydrogen as a function of nanotube diameter, with experimental validation.
- A corpus-level statement on site heterogeneity in real carbons, which bears directly on reproducibility of enhanced-uptake claims.
- Testable spectroscopic predictions (desorption peak temperatures and shapes; Raman-coordinate dependence of deliverable capacity) for irradiated and reduced carbons.

### 2.7 What would constitute a monumental result

Either of the following, established with the full validation and uncertainty program of this manual:

- **A limit.** A demonstration that, under every near-ambient window considered and for every physically admissible defect population of pure or oxygen-functionalized sp² carbon, the deliverable capacity attributable to defect engineering is bounded by a value well below practical targets, with the bound derived variationally and confirmed numerically, and with the binding-energy gap between the two channels quantified as a property of carbon rather than of any particular defect.
- **A window-populating class.** A demonstration that a specific, synthetically accessible defect class (for example a particular curvature range combined with a particular vacancy-derived configuration) places a large fraction of its sites in the kinetically admissible window, yielding predicted deliverable capacities that approach practical targets at a defect density below the amorphization limit, with a specified Raman signature and desorption spectrum.

### 2.8 Why the result is plausible to investigate computationally

Every ingredient is computable at workstation scale or available publicly: site energies by tight-binding and universal interatomic potentials with local density-functional anchors on supercells of 50–200 atoms; high-level anchors from the literature; kinetics by transition-state estimates and published barriers; thermodynamics by closed-form statistical mechanics; populations by closed-form cross sections and rate equations; and experimental constraints from a public isotherm database (thousands of H₂ isotherms, several hundred on carbons) and from published desorption spectra. The central conclusions depend on distributions and bounds, which are robust to the residual uncertainty of any single configuration energy; that robustness is itself quantified.

### 2.9 Why the question can be investigated using public evidence and computationally generated data

Inputs are limited to the following classes: public databases (NIST/ARPA-E ISODB, mirrored on GitHub), published first-principles values (extracted from open-access articles and preprints with full provenance), published measurements of displacement cross sections, Raman defect relations, desorption spectra and release temperatures (used as benchmarks, not as calibration data), and data generated by the program's own atomistic and statistical workflows. No unpublished or laboratory-internal data are used.

### 2.10 Boundaries of the scientific claims

- Claims concern deliverable capacity per unit mass of sorbent material under specified windows; system-level gravimetric and volumetric capacities, tank engineering, and heat management are out of scope.
- Claims concern pure and oxygen-functionalized sp² carbon; heteroatom-doped and metal-decorated carbons are out of scope for the core program.
- Claims about chemisorbed hydrogen distinguish *dischargeable* capacity (release-side only, with charging by an atomic-hydrogen source) from *cyclable* capacity (charge and release under the window conditions).
- Computed site energies are reported with the fidelity tier and calibration uncertainty; surrogate and low-fidelity values are never described as first-principles results.
- Correlations observed in the corpus analysis are reported as correlations; mechanistic statements are reserved for results supported by the landscape and kinetics components.

### 2.11 Relationship to the researcher's laboratory and prior computational work

The program addresses the question that motivates irradiation-based defect engineering of carbon nanomaterials for hydrogen storage, and its observables (Raman D/G and D/D′, desorption spectra, hydrogen content) are those used in that research. It uses only published irradiation results (for example Jones et al. 2010; Meyer et al. 2012) and requires no laboratory data. Relative to the prior independent projects, SiteWindow's forward desorption model predicts spectra from the landscape rather than inverting measured spectra (tpd-sbi), does not quantify Raman measurement uncertainty (RamanUQ), does not model storage-system thermodynamics (H2STAR), and does not build a literature corpus (HyCAN-DB). The data contracts in Section 10.6 would allow those projects to consume SiteWindow outputs through files, with no code dependency in either direction.
## 3. Scientific hypotheses and falsifiability

Notation used throughout: ε denotes a site binding energy (per H₂ for the molecular channel, per H relative to ½ H₂ for the atomic channel, both negative for exothermic binding and reported in kJ mol⁻¹ or eV as stated); E_a† a charging (sticking or dissociation) barrier; E_d† a release (desorption or recombination) barrier; W an operating window specified by (T_c, P_c) for charging, (T_r, P_r) for release, fill time τ_c, and release time τ_r; ρ(ε, E_a†, E_d†) the density of sites per carbon atom; D[ρ; W] the deliverable capacity in mass fraction. Formal definitions are in `docs/specs/DERIVATIONS.md`.

Operating windows used for hypothesis tests (all are configuration data; Section 4.4):

| Window | Charge | Release | τ_c | τ_r | Purpose |
|---|---|---|---|---|---|
| W1 | 298 K, 100 bar | 298 K, 5 bar | 300 s | 10³ s | Ambient pressure swing |
| W2 | 233 K, 100 bar | 358 K, 5 bar | 300 s | 10³ s | DOE temperature envelope, pressure and temperature swing |
| W3 | 77 K, 100 bar | 77 K, 5 bar | 300 s | 10³ s | Cryogenic pressure swing (reference for Chahine-type behavior) |
| W4(T_r) | atomic-H source (coverage set by site saturation) | T_r ∈ {358, 423, 473, 573} K, 5 bar, release side only | — | 10³ s | Dischargeable capacity of chemisorbed hydrogen |

### H1. Window-occupancy hypothesis (confirmatory)

- **Statement.** For every defect class of pure sp² carbon in the landscape (basal plane, monovacancy, divacancy 5-8-5 and 555-777, Stone–Wales, zigzag and armchair edges, nanotube curvature series with diameters 0.6–2.5 nm, and their hydrogenated derivatives), the fraction of sites with binding energy inside the thermodynamic window of W1 or W2 and simultaneously kinetically admissible is small, and the resulting deliverable capacity at the maximum stable defect density is below 1 wt% for W1 and W2.
- **Physical rationale.** The physisorption well (4–8 kJ mol⁻¹) sits below the window and the chemisorption well (0.7–4 eV per atom, i.e. 70–400 kJ mol⁻¹ per H) sits far above it; defect-induced shifts are of order 1–5 kJ mol⁻¹ for physisorption and change chemisorption energies by electronvolts without bringing release barriers below the admissibility threshold at 298–358 K.
- **Quantities.** Posterior site-energy distributions per class; admissibility indicator per site; D[ρ_class; W1], D[ρ_class; W2]; the posterior probability P(D < 1 wt%).
- **Null / competing explanation.** Some class (for example curvature-plus-vacancy configurations, or hydrogenated vacancy sites with low recombination barriers) populates the window with a fraction large enough that D exceeds 1 wt% at admissible defect densities.
- **Computational experiments.** E1 (landscape construction), E3 (deliverable maps), E4 (class-resolved capacity versus defect density).
- **Data and benchmarks.** Landscape anchors B1–B4, B11; stability limit from B10.
- **Criterion.** H1 is supported for a class if P(D_class(W1) < 1 wt%) > 0.95 and P(D_class(W2) < 1 wt%) > 0.95 under the calibrated landscape; H1 is falsified for a class if P(D_class(W) > 1 wt%) > 0.95 for either window. Intermediate cases are reported as undetermined with the full posterior.
- **Uncertainty.** Dominated by multi-fidelity calibration uncertainty in chemisorption energies (expected ±0.1–0.2 eV) and by barrier uncertainty (±0.2 eV); propagated by Monte Carlo over the posterior.
- **Limitations of the test.** The landscape enumerates finite configuration sets per class; coverage-dependent lateral interactions are included through dimer and cluster configurations but not through a full lattice-gas simulation at high coverage. The result is therefore a statement about the enumerated classes with explicit coverage of the configuration space reported.

### H2. Curvature–kinetics hypothesis (confirmatory)

- **Statement.** For chemisorbed hydrogen at saturating coverage on single-walled nanotubes, the characteristic release temperature T₉₀ (temperature at which 90% of chemisorbed hydrogen is released within τ_r = 10³ s under a linear ramp) decreases monotonically with curvature 1/d over 0.6 ≤ d ≤ 2.5 nm, and the calibrated landscape reproduces, without fitting, the published diameter dependence of hydrogenation stability and the complete C–H bond breaking by about 873 K.
- **Rationale.** Pyramidalization lowers the sp²→sp³ rehybridization cost, lowering both binding energy and the recombination barrier.
- **Quantities.** ε(d), E_d†(d) per site type; forward desorption spectra; T₉₀(d).
- **Null.** No monotonic trend, or a trend of the wrong sign, or predicted T₉₀ differing from benchmark release temperatures by more than the acceptance tolerance.
- **Experiments.** E2 (curvature series), E6 (forward desorption benchmark).
- **Benchmarks.** B6 (Nikitin et al. 2005, 2008), B5 (Zecho 2002 and Hornekær 2006 for the flat limit).
- **Criterion.** Supported if the posterior predictive T₉₀ for the flat limit lies within ±50 K of the benchmark high-temperature desorption peak region (560 K for H on graphite, with the ramp-rate correction defined in the benchmark registry) and the predicted ordering of release temperatures across 1 nm and 2 nm tubes matches the published ordering; falsified if the flat-limit prediction misses by more than 100 K or the diameter ordering is reversed with posterior probability > 0.9.
- **Uncertainty.** Barrier uncertainty ±0.2 eV translates to about ±60 K at these temperatures; prefactor uncertainty (10¹² – 10¹⁴ s⁻¹) adds about ±30 K.
- **Limitations.** Desorption from real nanotube films involves bundles, defects, and diameter distributions; the benchmark comparison uses the published mean diameters and treats the diameter distribution as a nuisance parameter.

### H3. Heterogeneity-threshold hypothesis (confirmatory, theoretical with numerical test)

- **Statement.** For the molecular channel, when the mean site enthalpy lies below the homogeneous optimum, adding a minority of stronger sites increases deliverable capacity only if their enthalpy lies below a computable threshold ε_c(W); for pure-carbon physisorption sites (basal plane, slit pores, nanotube interstitials and grooves, vacancy-adjacent and oxygen-adjacent sites) the strong-site enthalpies lie below ε_c for W1 and W2, so that the deliverable-capacity gain from defect-enhanced physisorption is bounded by a quantity that is computed and reported.
- **Rationale.** The Langmuir swing Δθ(ε; W) is unimodal in ε; sites stronger than the enthalpy at which Δθ falls below the mean swing reduce rather than increase delivery per site.
- **Quantities.** Δθ(ε; W), ε_c(W), the physisorption landscape per class.
- **Null.** The physisorption landscape contains site populations above ε_c with site density sufficient to raise D by more than 0.5 wt% in W1 or W2.
- **Experiments.** E3, E5 (physisorption landscape), E8 (corpus inference as an independent empirical check).
- **Criterion.** The derivation is verified by analytic limits and numerical tests (Section 7); the empirical part is supported if the posterior for the window-occupying site density across classes gives P(ΔD < 0.5 wt%) > 0.95 in W1 and W2.
- **Uncertainty.** Physisorption energies from empirical potentials validated against B3 (CCSD(T) coronene) and B4 (isosteric heats): ±1 kJ mol⁻¹.
- **Limitations.** Entropy of adsorption is modeled (Section 4.3), with the assumed form varied in sensitivity analysis.

### H4. Dischargeable-capacity hypothesis for atomic-hydrogen charging (confirmatory)

- **Statement.** Under W4 with release at T_r ≤ 473 K, the dischargeable capacity of hydrogen chemisorbed on flat defective graphene is limited by dimer- and cluster-mediated recombination barriers to a small fraction of the saturation hydrogen content, whereas curved sites (d ≲ 1.5 nm) and specific vacancy-derived configurations release a majority of their hydrogen; the site classes that satisfy T₉₀ ≤ 473 K are enumerated with posterior probabilities.
- **Rationale.** Published desorption peaks on graphite (445 K and 560 K) and the 450 °C annealing required to restore graphene from graphane indicate that flat-sheet release at 473 K is incomplete within 10³ s; curvature lowers barriers.
- **Quantities.** Release fraction f(T_r, τ_r) per class; dischargeable capacity.
- **Null.** Flat defective graphene releases > 50% of saturation hydrogen at 473 K within 10³ s; or no curved class releases > 50%.
- **Experiments.** E4, E6.
- **Benchmarks.** B5, B6, B12.
- **Criterion.** Supported if posterior P(f_flat(473 K) < 0.5) > 0.9 and at least one curved or vacancy-derived class has P(f(473 K) > 0.5) > 0.9; falsified otherwise.
- **Limitations.** Charging by an atomic-H source is assumed to saturate the sites; the actual coverage achieved by irradiation-driven hydrogenation depends on dose and adsorbate supply and is treated parametrically.

### H5. Corpus-heterogeneity hypothesis (exploratory)

- **Statement.** Across H₂ isotherms on carbon materials in the NIST/ARPA-E ISODB, inferred adsorption-energy distributions place less than 10% of site capacity above 10 kJ mol⁻¹ for the large majority of materials; materials that appear to violate this are associated with identifiable features (low-coverage data only, digitized sources, non-carbon components, or inconsistent multi-temperature data).
- **Rationale.** Carbon physisorption energetics are dominated by dispersion with modest enhancements from pore overlap and functional groups.
- **Quantities.** Posterior adsorption-energy distributions per material; isosteric heat versus loading; window-occupancy fraction.
- **Exploratory method.** Bayesian inversion with a measurement-error model; posterior predictive checks; classification of outliers by metadata; no pre-registered threshold beyond the descriptive 10%/10 kJ mol⁻¹ statement, which is reported with its posterior rather than tested as a hypothesis.
- **Benchmarks.** B4, B7 (Chahine-type surface-area scaling), B8 (AX-21 reference).
- **Limitations.** Excess-to-absolute conversion requires pore volume or skeletal density, which the database does not always provide; the inference reports sensitivity to this conversion explicitly.

### H6. Irradiation-population hypothesis (confirmatory, conditional on H1–H4)

- **Statement.** For electron irradiation of graphene and nanotubes, the deliverable capacity per unit Raman-observable defect density, dD/d(I_D/I_G), is bounded by the class composition that knock-on and electron-stimulated processes produce, and no irradiation condition yields a population whose deliverable capacity exceeds that of the best single class at the same defect density.
- **Rationale.** Irradiation produces vacancy-dominated populations with divacancy and Stone–Wales fractions set by rate equations; the population is a mixture, and mixtures cannot exceed the best component per site.
- **Quantities.** Class fractions versus dose and energy; I_D/I_G and I_D/I_D′ predicted; D versus I_D/I_G.
- **Null.** A synergistic combination (for example vacancy adjacent to curvature) exceeds the best single class per site, which would require the landscape to contain non-additive site classes; this is tested by including mixed configurations explicitly.
- **Experiments.** E7 (populations), E4.
- **Benchmarks.** B9 (Meyer 2012 cross sections), B10 (Raman relations), Jones et al. 2010 as a qualitative benchmark for hydrogenation under low-energy irradiation.
- **Criterion.** Supported if the mixture bound holds within posterior uncertainty for all simulated conditions; falsified if a mixed class with P > 0.9 exceeds the best pure class by more than 20% per site.
- **Limitations.** Low-energy electron-stimulated hydrogenation lacks absolute cross sections; this branch is modeled phenomenologically with explicit parameter uncertainty and is labeled as such.

### 3.1 Statement on hypothesis independence from outcome

Each hypothesis is formulated with its contradicting result defined before any landscape computation is executed. A falsified hypothesis is reported as such. The variational bound (Section 4.4) is a theorem whose correctness is verified by analytic limits and numerical tests independent of the landscape; the landscape determines where real materials sit relative to the bound, which is the empirical content of the program.
## 4. Comprehensive research program

The program consists of eight scientific components (C1–C8) and the software framework (Section 10). Each component is specified below with purpose, inputs, outputs, governing equations, assumptions, algorithms, statistics, implementation, dependencies, data sources, benchmarks, resources, uncertainty, convergence and sensitivity requirements, limiting cases, failure modes, fallbacks, relationships, independently interpretable outputs, and acceptance criteria. Benchmarks B1–B12 are defined in `docs/specs/BENCHMARK_REGISTRY.md`; data sources S1–S9 in `docs/specs/SOURCE_REGISTRY.md`; computational experiments E1–E9 in Section 9.

Unit conventions (Section 5 and `docs/specs/DATA_SCHEMA.md`): energies stored in eV per adsorbed entity, with kJ mol⁻¹ as a derived display unit (1 eV = 96.485 kJ mol⁻¹); temperature in K; pressure in bar (absolute); loading in mol kg⁻¹ of sorbent, with wt% = 100 × m_H / (m_sorbent + m_H) as a derived quantity; site densities in sites per carbon atom; areal densities in nm⁻².

### C1. Window theory and deliverable-capacity functional

**Scientific purpose.** Provide the thermodynamic–kinetic theory that maps a site-energy distribution to deliverable capacity under an operating window, including a variational upper bound. This component answers the "accounting" question and defines the quantity every other component feeds.

**Relationship to the central question.** D[ρ; W] is the quantity whose maximum is sought; the window is the object against which landscapes are tested.

**Inputs.** Site distribution ρ over (channel, ε, E_a†, E_d†, g, n_max, Δm) where g is the number of equivalent sites per defect, n_max the maximum hydrogen entities per site, Δm the lattice mass change per defect; window W; entropy model parameters; kinetic prefactors.

**Outputs.** Per-site swing Δθ, admissibility A, per-site deliverable mass; D[ρ; W]; dischargeable capacity D_dis[ρ; W4]; window efficiency η(W) = sup_ε Δθ·A; threshold ε_c(W); the variational bound D_max(W) for each channel; closed-form optimum ε*(W).

**Governing equations.**

Molecular channel (non-dissociative, Langmuir sites):

- Equilibrium coverage: θ(ε, T, P) = K P / (1 + K P), K = (1/P°) exp[−ΔG_ads(ε, T)/(R T)], ΔG_ads = ΔH_ads − T ΔS_ads, with ΔH_ads(T) = ε − (c_p corrections, Section C1 assumptions) and ΔS_ads modeled as ΔS° + ΔS_conf(T) (baseline ΔS° = −8R following Bhatia and Myers; alternatives: localized harmonic oscillator and mobile two-dimensional gas models, both derived in `DERIVATIONS.md` §2).
- Quantum correction for light adsorbate at low T: effective potential by the Feynman–Hibbs second-order expansion, U_FH(r) = U(r) + (ħ²/(24 μ k_B T)) ∇²U(r), applied in C2 when generating physisorption site energies; C1 consumes temperature-dependent ε(T) where supplied.
- Swing: Δθ(ε; W) = θ(ε, T_c, P_c) − θ(ε, T_r, P_r).
- For an isothermal pressure swing with ratio r = P_c/P_r, the maximum swing over ε is Δθ_max = (√r − 1)/(√r + 1), attained at K = 1/√(P_c P_r). For W1 (r = 20), Δθ_max = 0.634. The optimum enthalpy with ΔS° = −8R for T = 298 K, P_c = 30 bar, P_r = 1.5 bar is 15.1 kJ mol⁻¹ (benchmark B8, analytic).

Atomic channel (dissociative; 2 H per H₂; one H per site):

- Equilibrium with gas-phase H₂ for a site energy ε per H relative to ½ H₂: θ = √(K P) / (1 + √(K P)), K = (1/P°) exp[−2(ε − T Δs)/(R T)] where Δs is the per-H entropy change (derivation in `DERIVATIONS.md` §3). For ε of order −1 eV per H, θ ≈ 1 at all window conditions; the swing is kinetically, not thermodynamically, limited.
- Release-side kinetics (Polanyi–Wigner): dθ/dt = −ν_n θⁿ exp[−E_d†/(k_B T)], n ∈ {1, 2}; isothermal release fraction f(T_r, τ_r) = 1 − θ(τ_r)/θ₀ with closed forms for n = 1 (f = 1 − exp(−k τ)) and n = 2 (f = k θ₀ τ / (1 + k θ₀ τ)).
- Charging-side kinetics for sticking of atomic H: k_c = s₀ Φ_H exp(−E_a†/(k_B T_c)) with incident flux Φ_H set by the atomic-H source (parametric); for H₂ dissociative charging on pure carbon the barrier (≥ 2 eV on the basal plane in the landscape; verified in C2) renders k_c τ_c ≪ 1, and the channel is labeled non-cyclable under W1–W3.

Admissibility: A_c = 1 if k_c τ_c ≥ −ln(1 − η_c), A_r = 1 if f(T_r, τ_r) ≥ η_r, with η_c = η_r = 0.9 baseline and smooth (logistic in log k τ) alternative used in sensitivity analysis.

Mass accounting and the functional:

- For defect class c with density x_c (defects per C), sites s ∈ c with multiplicity g_s and maximum occupancy n_s, the deliverable hydrogen per carbon is h_c = x_c Σ_s g_s n_s Δθ_s A_s and the lattice mass per carbon is m_C + x_c Δm_c.
- D[ρ; W] = (Σ_c h_c m_H₂-equivalent) / (m_C + Σ_c x_c Δm_c + Σ_c x_c Σ_s g_s n_s θ_s,c m_H) expressed as a mass fraction, with the hydrogen mass in the denominator evaluated at the charged state.
- Variational bound (theorem, `DERIVATIONS.md` §4): for fixed window and channel, D ≤ D_max(W) = w(n_tot,max) · η(W), where n_tot,max is the maximum admissible site density for the channel (molecular: monolayer capacity on the maximum accessible specific surface area of 2630 m² g⁻¹ with the H₂ areal density from the landscape; atomic: 1 H per C) and w(·) converts site density to mass fraction. Equality requires all sites at ε* and admissible; heterogeneity lowers D unless mean enthalpy is below ε* and added sites lie between the mean and ε_c (H3).

**Physical and chemical assumptions.** Sites are independent (lateral interactions enter only through the site enumeration in C2, which includes dimers and clusters as distinct sites); adsorbed hydrogen entropy follows a specified model; gas-phase H₂ is treated with the Leachman reference equation of state for fugacity above 50 bar (implemented via CoolProp if available, with ideal-gas and virial fallbacks that agree within 2% below 100 bar at 298 K; the fugacity correction is included in K P → K f).

**Numerical algorithms.** Closed forms where available; adaptive quadrature over ε for continuous distributions (scipy `quad` with absolute tolerance 10⁻⁸); root finding for ε*(W) and ε_c(W) by Brent's method with bracketing; Monte Carlo propagation over posterior samples of ρ (vectorized NumPy).

**Statistical methods.** None beyond propagation; this component is deterministic given inputs.

**Implementation.** `sitewindow.theory` (`langmuir.py`, `dissociative.py`, `entropy.py`, `window.py`, `deliverable.py`, `bounds.py`, `fugacity.py`), typed, pure NumPy/SciPy.

**Dependencies.** numpy, scipy; optional CoolProp.

**Data sources.** DOE target tables (S9) for window definitions.

**Benchmarks.** B8 (Bhatia–Myers optimum, analytic); B4 (isosteric heats, as consistency of the entropy model with measured temperature dependence on AX-21 and graphite-like carbons); internal analytic limits (Section 7).

**Resources.** Negligible (milliseconds per evaluation; posterior propagation with 10⁴ samples in seconds).

**Uncertainty.** Entropy-model structural uncertainty treated by model set {−8R, harmonic oscillator, 2D gas}; prefactor priors ν ∈ 10¹² – 10¹⁴ s⁻¹ (log-uniform); η thresholds varied.

**Convergence and sensitivity.** Quadrature tolerance study; Monte Carlo sample-size study (10³–10⁵) to confirm posterior summaries stable to < 1% of their value; Sobol indices of D with respect to (ε, E_d†, ν, ΔS°, window parameters).

**Limiting cases.** Homogeneous Langmuir recovers Bhatia–Myers; r → 1 gives Δθ → 0; T_r → ∞ with τ_r fixed gives A_r → 1; E_d† → 0 gives f → 1; ε → −∞ gives θ → 1 and Δθ → 0 (irreversible trap).

**Failure modes and fallbacks.** Numerical overflow for extreme K handled by log-domain evaluation; non-unimodal Δθ for exotic entropy models detected by multi-start root finding and reported.

**Relationships.** Consumes landscapes from C3, kinetics from C4, populations from C5; feeds C7.

**Independently interpretable outputs.** η(W), ε*(W), ε_c(W), D_max(W) for all windows; published as a table and as `theory_window_tables.parquet`.

**Acceptance criteria.** (i) B8 reproduced to 0.1 kJ mol⁻¹; (ii) all limiting cases pass automated tests at 10⁻⁶ relative tolerance; (iii) Δθ_max(r) formula verified numerically for r ∈ {2, 5, 20, 100} to 10⁻⁸; (iv) bound theorem verified by randomized tests (10⁴ random distributions, no violation); (v) documentation of every equation with symbol table.

### C2. Site-energy landscape by multi-fidelity atomistic computation

**Scientific purpose.** Generate site energies and barriers for hydrogen on an enumerated set of defect classes, curvatures, coverages, and oxygen functionalities using workstation-feasible methods, with explicit fidelity tiers for calibration in C3.

**Relationship to the central question.** This is the occupancy ingredient: the raw landscape before statistical calibration.

**Inputs.** Structure specifications (defect class, supercell, curvature, functionalization, hydrogen configuration); calculator settings per tier; anchor list.

**Outputs.** For each configuration: relaxed structure, total energy, binding energy per H or per H₂ with explicit reference (½ H₂ or H atom; H₂ molecule), sticking and recombination barriers for designated pathways, local descriptors (coordination, pyramidalization angle, nearest-defect distance, neighbor hydrogen count), fidelity tier, settings hash. Stored as `landscape_raw.parquet` with provenance.

**Site classes enumerated (minimum set).**

| Class | Structures | Hydrogen configurations |
|---|---|---|
| Basal plane | 5×5 and 6×6 graphene | monomer; ortho, para, meta dimers; 3- and 4-atom clusters; graphane limit |
| Monovacancy | 6×6 with 1 vacancy (Jahn–Teller reconstructed) | 1–6 H sequential (dangling-bond saturation then basal) |
| Divacancy 5-8-5 and 555-777 | 7×7 | 1–4 H at reconstruction sites and neighbors |
| Stone–Wales | 6×6 | H at the four inequivalent pentagon/heptagon sites |
| Zigzag and armchair edges | nanoribbon cells (H-terminated edges as reference; bare-edge variants) | edge C–H; sub-edge basal |
| Curvature series | SWCNT (n,n) and (n,0) with d = 0.6–2.5 nm (at least 8 diameters), exterior and interior sites | monomer; ortho/para dimers; high-coverage stripes |
| Oxygen functionalities (GO/rGO) | epoxide, hydroxyl, carbonyl pair, carboxyl on 6×6 graphene | H adjacent; H₂ physisorption near group |
| Mixed | vacancy on (10,0) and (8,8) tubes; SW on tubes | monomer and dimers |
| Physisorption set | basal, slit pores 0.6–1.2 nm, (n,n) tube interior/interstitial/groove in bundles, vacancy-adjacent, oxygen-adjacent | single H₂ |

**Fidelity tiers.**

- *Tier 0 (anchors from literature).* Published high-level values: DMC basal monomer (B1), CCSD(T)/CBS H₂–coronene (B3), published DFT values with method metadata (S4). Used only in C3.
- *Tier 1 (local DFT).* GGA-PBE with D3(BJ) dispersion for chemisorption energetics and the vdW-DF family (optB88-vdW) for physisorption cross-checks, spin-polarized, PAW or norm-conserving pseudopotentials. Primary engine GPAW (conda-forge binary) in LCAO dzp mode for relaxations and plane-wave single points (cutoff 500 eV) for energies, with the D3(BJ) dispersion energy and forces added through the `dftd3-python` interface to simple-dftd3 as an ASE calculator (included in `environment.yml`); k-mesh 4×4×1 for 5×5 cells (scaled by cell), vacuum ≥ 15 Å, dipole correction; forces < 0.02 eV Å⁻¹. Fallback engine when system libraries are unavailable: PySCF periodic DFT with GTH pseudopotentials and gth-dzvp basis, PBE-D3, same supercells, Γ-centered meshes, used for the verification tier and for anchors where GPAW is unavailable. Anchor set: 40–60 configurations spanning all classes.
- *Tier 2 (semi-empirical tight binding).* GFN2-xTB via `tblite` with periodic boundary conditions (and GFN1-xTB as a cross-check), for all chemisorption configurations including sequential hydrogenation and barriers by constrained relaxation. GFN2-xTB is known to describe C–H bond formation and sp³ rehybridization qualitatively; its systematic deviation per class is the object of calibration in C3.
- *Tier 3 (universal interatomic potential).* MACE-MP-0 (medium, CPU) through ASE for screening of curvature series, cluster configurations, and nudged-elastic-band (climbing image, 7 images) barrier estimates. Transferability to hydrogenated, defective carbon is not assumed: MACE results are used only after per-class validation against Tier 1 in C3, and MACE results for classes failing validation are excluded.
- *Physisorption potentials.* Pairwise H₂–C potentials with Feynman–Hibbs correction (Lennard-Jones with Crowell–Brown-type parameters, and a second parameter set from the literature as a model-form alternative), validated against B3 and B4; used for slit pores, tube geometries, and defect-adjacent sites by numerical integration over carbon positions from the relaxed structures.

**Governing principles.** Binding energy per H: E_b(H) = E(sub + H) − E(sub) − ½ E(H₂) (storage reference) and E_b′(H) = E(sub + H) − E(sub) − E(H) (sticking reference), both stored. Sequential binding: ΔE_n = E(sub + nH) − E(sub + (n−1)H) − ½ E(H₂). Physisorption: E_b(H₂) = E(sub + H₂) − E(sub) − E(H₂) with counterpoise correction for localized-basis engines. Barriers: constrained minimization along the C–H distance (Tier 2) and climbing-image NEB (Tier 3) with Tier-1 anchors for designated pathways (monomer sticking; ortho-dimer and para-dimer recombination; vacancy-site recombination; curvature series recombination at two diameters).

**Assumptions.** Periodic supercells with defect separation ≥ 1 nm (defect–defect interaction negligible for chemisorption energies at this separation, consistent with Boukhvalov and Katsnelson 2008); ground-state spin multiplicities searched (singlet/triplet for vacancy systems); zero-point energy corrections estimated from finite-difference C–H stretch frequencies for Tier-1 anchors and applied as class-level corrections; no thermal expansion.

**Numerical algorithms.** BFGS/FIRE geometry optimization (ASE); climbing-image NEB with 7 images and spring constant 0.1 eV Å⁻¹; constrained relaxation with fixed C–H distance on a 0.1 Å grid refined to 0.02 Å near the saddle; convergence tests on supercell size, k-mesh, cutoff, vacuum.

**Implementation.** `sitewindow.atomistic` (`structures/` builders using ASE; `calculators/` adapters for tblite, MACE, GPAW, PySCF; `workflows/` relax, adsorb, sequential, barrier; `cache.py` content-addressed caching keyed by structure hash and settings hash; `descriptors.py`).

**Dependencies.** ase, tblite, mace-torch (optional extra), torch (CPU), gpaw or pyscf (optional extras), numpy, scipy.

**Data sources.** None external beyond anchor literature (S4); structures are generated.

**Benchmarks.** B1 (DMC basal monomer), B2 (sticking barrier), B3 (H₂–coronene), B11 (vacancy formation and migration energies as structural sanity checks of the engines).

**Required computational resources.** Tier 2: ≈ 1–10 s per single point and ≈ 1–5 min per relaxation for 50–100 atoms on 1 core; the full enumerated set (≈ 2000 configurations including sequential hydrogenation) ≈ 50–150 core-hours. Tier 3: ≈ 0.5–2 s per force call for 100–250 atoms on 1 core; relaxations ≈ 2–5 min; NEB ≈ 20–40 min; ≈ 100–300 core-hours for the curvature and cluster sets. Tier 1 (GPAW LCAO + PW single point): ≈ 1–3 wall-hours per anchor on 8 cores; 40–60 anchors ≈ 60–180 wall-hours, embarrassingly parallel across anchors and checkpointed per structure. Memory < 8 GB for all tiers. Storage < 3 GB for trajectories and caches.

**Uncertainty.** Per-tier systematic deviations estimated against anchors (C3); numerical uncertainty from convergence studies (target: ≤ 0.02 eV for chemisorption energies, ≤ 0.5 kJ mol⁻¹ for physisorption at Tier 1; ≤ 0.05 eV at Tiers 2–3 relative to their own converged settings).

**Convergence and sensitivity.** Supercell 4×4/5×5/6×6/7×7 for monomer and vacancy; k-mesh 2/4/6; cutoff 400/500/600 eV; vacuum 10/15/20 Å; NEB image count 5/7/9; basis dzp vs. plane-wave agreement for 10 anchors.

**Limiting cases.** Large-diameter tube → graphene values (monotone convergence required); isolated H₂ far from surface → zero binding; graphane limit → published graphane formation energy range.

**Failure modes and fallbacks.** GPAW unavailable → PySCF; tblite failure for specific spin states → GFN1 or MACE with flagged tier; MACE instability on a class → class excluded from Tier 3 and computed at Tier 2 with reduced enumeration; NEB non-convergence → constrained-relaxation barrier with flagged method.

**Relationships.** Feeds C3 and C4; consumes structure definitions from C5 for defect geometries.

**Independently interpretable outputs.** The raw multi-tier landscape dataset with descriptors (publishable as a data product).

**Acceptance criteria.** (i) convergence targets met and documented; (ii) Tier-1 basal monomer within 0.1 eV of the published PBE range and Tier-1 H₂–coronene (cluster calculation) within 1.0 kJ mol⁻¹ of B3 with the vdW-DF functional; (iii) every class has ≥ 1 Tier-1 anchor and all Tier-2/3 values carry tier metadata; (iv) all configurations relaxed to the force tolerance; (v) curvature series converges to graphene values within 0.05 eV at d = 2.5 nm.

### C3. Published-value registry, hierarchical meta-analysis, and multi-fidelity calibration

**Scientific purpose.** Convert heterogeneous evidence (published high-level and DFT values; Tier 1–3 computed values) into calibrated posterior site-energy distributions per class with quantified epistemic uncertainty, so that method dependence is modeled rather than ignored.

**Relationship to the central question.** The occupancy question is answered with distributions, not point values; this component produces them.

**Inputs.** Registry of published values (S4; fields: DOI, system, site, functional, dispersion treatment, code, basis/cutoff, supercell, k-mesh, spin treatment, reference state, value, stated uncertainty, extraction method); C2 landscape across tiers.

**Outputs.** Posterior class means and dispersions; method effect estimates (functional, dispersion, cell size); calibrated predictive distributions for every Tier-2/3 configuration; `landscape_calibrated.parquet`; diagnostic report.

**Statistical model (hierarchical meta-analysis of published values).** For published value y_i with class c(i), method descriptors m(i):

y_i = μ_c(i) + β_func[m_i] + β_disp[m_i] + γ · (1/N_cell,i) + u_i + e_i,

with μ_c ~ N(μ₀, τ_c²), β ~ N(0, σ_β²) partially pooled, u_i ~ N(0, τ_study²) for study-level effects, e_i ~ N(0, s_i²) with s_i from stated uncertainty or a default reflecting convergence practice (0.05 eV); priors weakly informative; inference by NUTS (PyMC), 4 chains, R̂ < 1.01, ESS > 400 per parameter; posterior predictive checks and leave-one-study-out cross-validation.

**Multi-fidelity calibration.** For each class and tier t ∈ {2, 3}: E_T1 = a_c,t + b_c,t E_t + δ_c,t(x) + e, where δ is a Gaussian process over descriptors x (pyramidalization angle, coordination, neighbor-H count, distance to defect) with a Matérn-5/2 kernel and partially pooled hyperparameters across classes; fitted on anchors; LOO-CV predictive coverage checked; predictive distribution assigned to every non-anchor configuration. Tier-0 anchors enter as the top of the hierarchy: the Tier-1 (PBE-D3) functional bias is itself calibrated against DMC/CCSD(T) where available (basal monomer, coronene) and carried as a class-independent offset with uncertainty for chemisorption (−0.13 ± 0.05 eV scale, estimated, not assumed) and for physisorption.

**Assumptions.** Additivity of method effects; exchangeability of studies within a method; GP stationarity within class; Tier-0 offset transferable across classes (tested by comparison with the two available high-level anchors and flagged as an assumption in the manuscripts).

**Algorithms.** NUTS; GP by exact inference (anchor counts are small); identifiability diagnostics by prior-posterior overlap and by profile likelihood for (a, b).

**Implementation.** `sitewindow.landscape` (`registry.py`, `meta_analysis.py`, `calibration.py`, `assemble.py`).

**Dependencies.** pymc, arviz, numpy, pandas; scikit-learn or GPyTorch-free custom GP (NumPy) to avoid heavy dependencies.

**Data sources.** S4 (published values with provenance), C2 outputs.

**Benchmarks.** B1, B3 as hold-out anchors in cross-validation; synthetic-data recovery tests.

**Resources.** Minutes to an hour on 2–8 cores.

**Uncertainty.** Fully Bayesian; reported as posterior intervals per class; epistemic (method) and numerical components separated by the model structure.

**Convergence and sensitivity.** MCMC diagnostics; prior sensitivity (three prior scales); kernel sensitivity (Matérn-3/2 vs 5/2); exclusion of any single study.

**Limiting cases.** With one method and no anchors, posterior reverts to the prior-widened computed value (tested); with many concordant studies, τ → 0.

**Failure modes and fallbacks.** Non-identifiable (a, b) for a class with < 3 anchors → fix b = 1 and estimate a only, flagged; divergent transitions → reparameterization (non-centered) and increased target acceptance.

**Relationships.** Consumes C2; feeds C1, C4, C7.

**Independently interpretable outputs.** Method-effect table (for example, PBE versus DMC offsets) of standalone value to the first-principles community.

**Acceptance criteria.** (i) R̂ < 1.01 and ESS > 400; (ii) LOO-CV 90% predictive intervals cover 85–95% of held-out anchors; (iii) synthetic recovery of planted class means within posterior 90% intervals in ≥ 90% of 100 replicates; (iv) every class has a documented calibration status.

### C4. Kinetics and forward desorption model

**Scientific purpose.** Determine kinetic admissibility of sites and predict thermal-desorption observables from the landscape for benchmarking.

**Relationship to the central question.** The admissibility question.

**Inputs.** Calibrated barriers and binding energies with uncertainty (C3); pathway definitions (direct monomer desorption as H atom, which is negligible; ortho-dimer direct recombination; para-dimer via meta transformation; vacancy-site recombination; curvature-series recombination; molecular physisorption desorption); prefactor priors; ramp rate β; initial coverages.

**Outputs.** Release fraction f(T_r, τ_r) per site; T₉₀ per class; simulated TPD spectra (rate versus T for given β) with posterior bands; peak temperatures and widths.

**Equations.** Polanyi–Wigner per pathway; for a linear ramp T = T₀ + βt, the rate r(T) = ν θⁿ exp(−E/k_B T)/β; coupled ordinary differential equations for multi-species (monomer, ortho, para, meta populations with diffusion/transformation barriers) integrated with a stiff solver (scipy `solve_ivp`, BDF, rtol 10⁻⁸); Redhead peak relation used as an analytic check. Isotope mass enters through ν (√(m_H/m_D) scaling as a baseline) and zero-point corrections to barriers (class-level).

**Assumptions.** Mean-field kinetics within a class; no readsorption; prefactor independent of coverage; transformation between dimer species included for the basal plane following the published mechanism (ortho → meta → para), with barriers from the landscape.

**Algorithms.** Stiff ODE integration; posterior propagation by sampling barriers and prefactors; peak detection.

**Implementation.** `sitewindow.kinetics` (`polanyi_wigner.py`, `pathways.py`, `tpd_forward.py`, `admissibility.py`).

**Dependencies.** numpy, scipy.

**Benchmarks.** B5 (graphite TPD peaks at 445 and 560 K for H, 490 and 580 K for D, with the ramp rate from the source), B6 (nanotube release temperatures), B12 (graphane restoration by annealing).

**Resources.** Seconds per spectrum; 10⁴ posterior spectra in minutes.

**Uncertainty.** Barrier and prefactor posteriors; model-form (n = 1 vs 2; with/without transformation) compared by the benchmark residuals.

**Convergence and sensitivity.** ODE tolerance study; ramp-rate dependence; initial-coverage dependence; Sobol indices of peak temperature.

**Limiting cases.** Single first-order species recovers Redhead's relation to 1 K; β → 0 recovers isothermal results.

**Failure modes and fallbacks.** Pathway omission detected by systematic benchmark residual → add pathway from the landscape with documented justification (never by fitting barriers to the benchmark).

**Relationships.** Consumes C3; feeds C1 admissibility and C7; provides the forward model that would be the natural counterpart of an inversion tool (external interface, Section 10.6).

**Independently interpretable outputs.** Predicted TPD spectra per class and per irradiation population (C5).

**Acceptance criteria.** (i) Redhead check passes; (ii) B5 reproduced within the stated tolerance (peak positions within ±50 K at 90% posterior) *without* fitting; (iii) B6 ordering reproduced; (iv) documented pathway set.
### C5. Defect-population model and Raman observables

**Scientific purpose.** Map processing conditions (electron irradiation energy, flux, fluence, substrate; chemical reduction state of graphene oxide) to defect-class densities, bound those densities by the amorphization limit, and express the population in the Raman coordinate used experimentally.

**Relationship to the central question.** The "physically realizable population" clause.

**Inputs.** Beam energy E (keV), fluence Φ (e⁻ nm⁻²), temperature, material (graphene, SWCNT of diameter d, rGO with O/C ratio), displacement threshold distribution; published cross sections (B9); Raman relation parameters (B10).

**Outputs.** Class densities x_c(Φ, E) with uncertainty; I_D/I_G and I_D/I_D′ predictions; the admissible population set used in C7.

**Equations.**

- Knock-on cross section (McKinley–Feshbach) for a target of mass M and charge Z with displacement threshold E_d, incident electron kinetic energy E, β = v/c, maximum energy transfer T_max = 2E(E + 2m_e c²)/(M c²):

  σ_d = (π Z² e⁴ (1 − β²))/(m_e² c⁴ β⁴) · { (T_max/E_d − 1) − β² ln(T_max/E_d) + π α β [ 2(√(T_max/E_d) − 1) − ln(T_max/E_d) ] }, α = Z/137.

- Thermal smearing: the target atom's out-of-plane velocity distribution (Gaussian, with mean-square velocity from the Debye model using the out-of-plane Debye temperature of graphene, 1287 K for ¹²C, as in Meyer et al. 2012) is convolved into the energy transfer so that the effective cross section is nonzero below the static threshold (static onset at about 108 keV for the first-principles displacement threshold of 22 eV, which is used without adjustment, as in the source). The smearing kinematics is the object of validation against B9; a fitted threshold is reported only as a sensitivity analysis, never used as calibration. Angular dependence of the threshold on nanotubes (emission direction relative to the tube axis and curvature) follows the functional form reviewed by Krasheninnikov and Nordlund.

- Rate equations for defect classes (per carbon atom basis):

  dx_V/dt = σ_d Φ̇ (1 − x_tot) − k_mig x_V (x_V + x_I) − k_dv x_V², dx_DV/dt = k_dv x_V² + σ_dv Φ̇ x_V, dx_SW/dt = σ_SW Φ̇ (1 − x_tot) − k_SW,anneal x_SW, with vacancy migration barrier ≈ 1.3 eV and interstitial recombination from the literature (B11), SW formation treated as a low-probability electronic-excitation channel with a cross section parameter bounded by published observation frequencies.

- Low-energy electron-stimulated hydrogenation (sub-threshold, adsorbate-mediated): dθ_H/dt = σ_H,eff Φ̇ n_ads (1 − θ_H/θ_sat) − σ_ESD Φ̇ θ_H, with effective cross sections as free parameters with wide priors; this branch produces hydrogen on sites of the classes present, not new carbon defects, and is labeled phenomenological.

- Amorphization limit: the mean inter-defect distance L_D = (π n_D)⁻¹ᐟ² must exceed about 1–2 nm for the lattice to remain sp²-dominated (stage 1 of the Raman disorder trajectory); the admissible set is x_tot ≤ x_amorph with x_amorph derived from L_D,min and the areal carbon density (38.2 nm⁻²), carried as an interval (L_D,min ∈ [1, 2] nm).

- Raman relation: I_D/I_G(L_D, E_L) per Cançado et al. 2011 (two-length-scale activation/structural model with r_A ≈ 3.1 nm, r_S ≈ 1 nm, and the excitation-energy-dependent prefactor as published), and I_D/I_D′ by defect type (≈ 13 for sp³, ≈ 7 for vacancy-like, ≈ 3.5 for boundaries; Eckmann et al. 2012) combined linearly by class fraction as a first-order model with its limitation stated.

**Assumptions.** Mean-field rate equations (no spatial correlation beyond the amorphization bound); knock-on cross sections for graphene applied to tube walls with the angular threshold correction; sub-threshold hydrogenation phenomenological.

**Algorithms.** Numerical integration of the McKinley–Feshbach expression with thermal convolution (adaptive quadrature); ODE integration (BDF); posterior sampling of parameters.

**Implementation.** `sitewindow.defects` (`knockon.py`, `thermal.py`, `rates.py`, `populations.py`, `raman.py`).

**Dependencies.** numpy, scipy.

**Data sources.** B9 values (extracted from the open-access preprint arXiv:1203.2372 with digitization uncertainty where tabular values are unavailable), B10 parameters (open-access arXiv:1105.0175), B11 energies (S4).

**Benchmarks.** B9 (cross section versus beam energy: the published endpoints of about 10⁻⁴ barn at 80 keV and about 0.2 barn at 100 keV for ¹²C, the ¹²C/¹³C ordering, and the digitized points of the source's figure with digitization uncertainty), B10 (Raman relation reproduces the published I_D/I_G versus L_D curve at 2.41 eV), qualitative: Jones et al. 2010 (hydrogenation observed under sub-threshold irradiation with adsorbates).

**Resources.** Negligible.

**Uncertainty.** Parameter posteriors; interval on L_D,min; the sub-threshold branch carries wide priors and is reported separately.

**Convergence and sensitivity.** Quadrature and ODE tolerances; Sobol indices of class fractions with respect to E_d, velocity parameter, k_mig, σ_SW.

**Limiting cases.** E below the smeared threshold → σ_d → 0; Φ → ∞ → x_tot → x_amorph (saturation enforced by the (1 − x_tot) factor and by reporting beyond-limit populations as inadmissible).

**Failure modes and fallbacks.** If B9 cannot be reproduced within uncertainty with the published functional form, the discrepancy is reported and the population model uses the measured cross sections directly by interpolation.

**Relationships.** Consumes C2 structures (classes); feeds C7; links to C4 (predicted TPD of irradiated populations).

**Independently interpretable outputs.** Class-fraction maps versus (E, Φ); Raman-coordinate conversion tables.

**Acceptance criteria.** (i) B9: computed curve within the digitization error bars for ≥ 80% of digitized points and within a factor of 2 of the published 100 keV value, with the 80 keV tail within one order of magnitude; (ii) B10 curve reproduced to within 5% over L_D ∈ [1, 20] nm; (iii) rate-equation conservation (total carbon) holds to 10⁻¹⁰; (iv) documented parameter priors.

### C6. Experimental corpus inference (hydrogen isotherms on carbons)

**Scientific purpose.** Determine, from public measured isotherms, the adsorption-energy distributions of real carbons and whether any population of sites lies in the window; provide an independent empirical constraint on the physisorption landscape and on corpus-scale relations (surface-area scaling).

**Relationship to the central question.** Independent evidence on occupancy from experiment, for the molecular channel.

**Inputs.** S1: NIST/ARPA-E ISODB GitHub mirror (`NIST-ISODB/isodb-library`, JSON per isotherm, with `digitizer` and `tabular_data` flags, units, temperature, adsorbent name, DOI). Carbon selection by adsorbent-name rules and by cross-reference to the adsorbent registry; metadata (BET area, pore volume) extracted from the isotherm records where present.

**Outputs.** Harmonized carbon H₂ isotherm table (`isodb_carbon_h2.parquet`); per-material multi-temperature sets; isosteric heats versus loading with uncertainty; posterior adsorption-energy distributions; window-occupancy fractions; corpus-level regressions (uptake at 77 K versus surface area) with uncertainty.

**Equations.**

- Unit normalization to mol kg⁻¹ (wt%, mmol g⁻¹, cm³(STP) g⁻¹, g L⁻¹ with stated bulk density only) and bar.
- Excess-to-absolute: n_abs = n_exc + ρ_gas(T, P) V_pore, with V_pore a parameter per material (prior from reported pore volume, else a broad prior), ρ_gas from the H₂ equation of state.
- Adsorption integral equation (heterogeneous Langmuir): n_abs(T, P) = n_m ∫ ρ(ε) θ(ε, T, P) dε on a grid ε ∈ [2, 40] kJ mol⁻¹ (spacing 0.5 kJ mol⁻¹), with ΔS model as in C1 and a per-material entropy offset parameter.
- Measurement-error model: y = n_abs + e, e ~ N(0, s²), s² = s_tab² for tabular data (2% relative plus 0.02 mol kg⁻¹) and s_dig² for digitized data (5% relative plus 0.05 mol kg⁻¹ plus a pressure-axis digitization term propagated through the local slope), with the scale factors estimated hierarchically across the corpus.
- Prior on ρ(ε): log-weights with a second-difference Gaussian random-walk (smoothness) prior and a Dirichlet-like normalization; smoothness scale given a half-normal hyperprior.
- Isosteric heat: q_st = R T² (∂ ln P/∂T)_n from pairs of temperatures, by local spline interpolation in ln P at fixed n, with bootstrap uncertainty.

**Assumptions.** Independent Langmuir sites (Tóth-type pore-filling effects treated as heterogeneity); negligible H₂–H₂ lateral interactions at the loadings considered (checked against virial coefficients); quantum effects at 77 K absorbed into effective ε (the comparison to the C2 landscape uses Feynman–Hibbs-corrected energies at the same temperature).

**Algorithms.** NUTS over per-material models with shared hyperparameters (two-stage: per-material fits with empirical-Bayes hyperparameters from a corpus-level fit); synthetic recoverability tests (planted bimodal distributions); posterior predictive checks; outlier classification by metadata.

**Implementation.** `sitewindow.experimental` (`isodb.py` loader, `select.py` carbon rules, `units.py`, `excess.py`, `isosteric.py`, `aed.py`).

**Dependencies.** pandas, pyarrow, numpy, scipy, pymc, arviz; optional CoolProp.

**Data sources.** S1; S2 (NIST adsorbent registry for material identity); S9.

**Benchmarks.** B7 (surface-area scaling of 77 K excess uptake; slope consistent with the Chahine-type rule within corpus scatter), B8 (AX-21 isotherms as the Bhatia–Myers reference carbon), B4 (isosteric heats 4–7 kJ mol⁻¹ for graphitic and activated carbons at low loading).

**Resources.** Corpus load seconds; per-material NUTS minutes; full corpus (≈ 100 materials) a few hours on 2–8 cores; checkpointed per material.

**Uncertainty.** Fully Bayesian; measurement-error scales estimated; V_pore sensitivity reported.

**Convergence and sensitivity.** Grid spacing 1.0/0.5/0.25 kJ mol⁻¹; smoothness hyperprior scales; entropy model alternatives; MCMC diagnostics.

**Limiting cases.** Synthetic single-site data recover a delta-like distribution at the planted energy within the grid resolution; synthetic Langmuir data at two temperatures recover the planted q_st.

**Failure modes and fallbacks.** Materials with single-temperature data only: energy and entropy are not jointly identifiable; such materials are analyzed under the fixed entropy model and labeled "single-T, entropy-fixed". Name-based carbon selection errors: a manual review list is generated and cross-checked against the registry; ambiguous materials are excluded with the exclusion logged.

**Relationships.** Independent of C2–C5 (no shared inputs); compared with the C2 physisorption landscape in C7; feeds H3 and H5.

**Independently interpretable outputs.** Corpus table, isosteric heats, and AEDs for carbons (data product).

**Acceptance criteria.** (i) ≥ 300 carbon H₂ isotherms harmonized with unit audit (every record's unit conversion covered by a test); (ii) synthetic recovery tests pass (planted energies within 1 kJ mol⁻¹ at 90% posterior); (iii) B7 slope within the published scatter; (iv) posterior predictive p-values uniform in [0.05, 0.95] for ≥ 80% of materials; (v) all exclusions logged with reasons.

### C7. Deliverable-capacity evaluation, bounds, sensitivity, and hypothesis tests

**Scientific purpose.** Combine C1–C6 to compute deliverable and dischargeable capacities per class, per population, and per window with full uncertainty; evaluate the bound; identify the regime structure; test H1–H6.

**Inputs.** Calibrated landscape (C3), admissibility (C4), populations (C5), window definitions (C1), corpus posteriors (C6).

**Outputs.** Capacity maps D(ε, E_d†; W); class-resolved D versus x_c; population-resolved D versus (E, Φ) and versus I_D/I_G; Pareto frontiers (D versus stability margin x_amorph − x_tot, and D versus release temperature); counterfactual maps (D if barriers reduced by δ, if entropy model changed, if window changed); hypothesis-test summaries with posterior probabilities; regime diagram in (ε, E_d†) space labeled by limiting mechanism (thermodynamic, kinetic-release, kinetic-charge, mass-limited).

**Methods.** Monte Carlo over the joint posterior (≥ 10⁴ samples); Sobol global sensitivity (Saltelli sampling, first- and total-order indices with bootstrap confidence intervals) of D with respect to all landscape, kinetic, entropy, population, and window parameters; conditional analyses by class; mixture-bound test (H6) by explicit comparison of mixed-configuration classes with pure classes.

**Implementation.** `sitewindow.evaluate` (`maps.py`, `populations.py`, `pareto.py`, `counterfactual.py`, `hypotheses.py`, `sobol.py`).

**Dependencies.** numpy, scipy, SALib, pandas.

**Resources.** Minutes to hours on 2–8 cores depending on sample counts.

**Convergence.** Sobol index stability versus sample size (2¹⁰–2¹⁴); MC sample size.

**Acceptance criteria.** (i) Every hypothesis has a computed posterior statement in the prescribed form; (ii) Sobol indices with confidence intervals narrower than 0.1 for the top three parameters; (iii) regime diagram and all maps regenerated by one command from cached inputs; (iv) the bound is never violated by any computed population (automated check).

### C8. Validation and benchmarking framework

**Scientific purpose.** Make every benchmark comparison an executable, versioned test with pre-registered tolerances.

**Specification.** `docs/specs/BENCHMARK_REGISTRY.md` defines B1–B12 with source, identifier, system, quantity, conditions, acquisition method, prediction method, metric, uncertainty treatment, independence status, acceptance threshold, and limitations. `sitewindow.validation` implements one runner per benchmark producing a machine-readable report (`validation_report.json`) and a human-readable summary; the thresholds are read from the registry file, which is under version control, and may not be altered after a benchmark run without a recorded justification in the changelog.

**Levels.** Mathematical correctness (analytic identities), numerical correctness (convergence, tolerances), analytic limiting cases, reproduction of published computational results (B1, B3, B11 at the appropriate tier), comparison with experiment (B4–B10, B12), cross-material validation (graphene versus nanotube versus activated carbon), out-of-distribution validation (classes held out from calibration), robustness to uncertainty (posterior coverage), and integrated-prediction validation (TPD spectra and corpus AEDs from the full pipeline).

**Discrepancy policy.** A failed benchmark is reported with the residual, its posterior probability, and a mechanistic investigation (which component contributes most, by the sensitivity machinery); thresholds are never relaxed to pass. If a discrepancy reflects a limitation of the landscape or of the kinetic model, that limitation is stated in the manuscripts.

**Acceptance criteria.** All benchmark runners execute in CI (fast tier) or in the reproducibility workflow (full tier) and produce reports; the registry is complete for every benchmark used in any manuscript claim.
## 5. Integration of research components

### 5.1 Scientific dependency graph

```
C1 Window theory ──────────────────────────────┐
                                               │
C2 Atomistic landscape ──► C3 Calibration ──► C4 Kinetics ──► C7 Evaluation ──► Manuscripts
         │                      │                  │               ▲
         └──► C5 Populations ───┴──────────────────┘               │
                                                                   │
C6 Corpus inference (independent inputs) ──────────────────────────┘
C8 Validation (cross-cutting; consumes every component's outputs)
```

Logical dependencies: C3 requires C2 and the registry; C4 requires C3; C5 requires C2 structures; C7 requires C1, C3, C4, C5, and C6; C8 requires all. C1 and C6 have no upstream scientific dependencies and are executed first in the reproducibility workflow because their outputs validate the theory and establish the empirical baseline independently of the atomistic landscape.

### 5.2 Computational data flow

```
registry/published_values.yaml ─┐
atomistic cache (parquet) ──────┼─► landscape_raw.parquet ─► landscape_calibrated.parquet ─┐
                                 │                                                           │
isodb-library (git mirror) ─► isodb_carbon_h2.parquet ─► aed_posteriors.nc ──────────────────┤
                                                                                             ▼
windows.yaml ─► theory_window_tables.parquet ───────────────────────────────► capacity_maps.nc
populations (knockon + rates) ─► populations.parquet ───────────────────────► hypotheses.json
                                                                              figures/*.pdf, tables/*.csv
```

All intermediates are content-addressed (SHA-256 of inputs and settings) and regenerated by `sitewindow reproduce` (Section 12).

### 5.3 Interfaces and exchanged quantities

| From → To | Quantity | Units | Uncertainty representation |
|---|---|---|---|
| C2 → C3 | site energies, barriers, descriptors, tier | eV, dimensionless | numerical tolerance per tier |
| C3 → C1/C4/C7 | posterior samples of (ε, E_a†, E_d†) per site | eV | NetCDF posterior arrays (ArviZ InferenceData) |
| C4 → C1/C7 | admissibility A, release fraction f, T₉₀ | dimensionless, K | posterior samples |
| C5 → C7 | class densities x_c, Raman ratios | per C, dimensionless | posterior samples; admissible-set interval |
| C6 → C7 | AED ρ(ε), q_st(n), window-occupancy | kJ mol⁻¹, mol kg⁻¹ | posterior samples |
| C1 → C7 | Δθ, η(W), ε*, ε_c, D_max | dimensionless, kJ mol⁻¹, wt% | deterministic given inputs |

### 5.4 Microscopic to macroscopic quantities

Site energies (microscopic, per configuration) → site classes with multiplicities and maximum occupancies (mesoscopic, per defect) → class densities (per carbon) → mass fractions (macroscopic, per sorbent mass). Specific surface area enters the molecular channel through the areal site density; for the atomic channel the carbon lattice itself is the site reservoir. Volumetric capacity is out of scope (Section 2.10).

### 5.5 Integration of data sources and fidelities

Published high-level values anchor local DFT; local DFT anchors tight-binding and universal potentials; the calibrated landscape is compared, for the molecular channel, with corpus-inferred AEDs (independent experimental evidence). Disagreements between the computed physisorption landscape and corpus AEDs are analyzed by material class and are reported as findings about either the landscape (missing site types) or the corpus (measurement artifacts), with the diagnosis supported by the sensitivity machinery.

### 5.6 From experiments to conclusions

Each computational experiment (Section 9) produces one or more quantities that enter a hypothesis test (Section 3) or a benchmark (C8). The final conclusions are statements about D_max(W), class-resolved window occupancy, admissible populations, and the corpus, each carrying posterior uncertainty. The combined program yields what no component yields alone: a limit (C1) evaluated on real site distributions (C2–C3) with kinetic admissibility (C4) and realistic populations (C5), checked against independent experimental evidence (C6, C8).

### 5.7 Uncertainty propagation across components

Posterior samples are passed as arrays, never summarized to point estimates between components. Correlations induced by shared calibration parameters (for example the Tier-0 offset) are preserved by using the same posterior draw index across components in each Monte Carlo replicate.

## 6. Data architecture and provenance

### 6.1 Public source registry

`docs/specs/SOURCE_REGISTRY.md` lists every external source (S1–S9) with identity, persistent identifier, availability, acquisition mechanism, authentication (none required for all sources used), licensing and redistribution terms, quantities, units, material characteristics, suitability, limitations, and overlap notes.

### 6.2 Acquisition mechanisms

- S1 (ISODB): `git clone --depth 1 https://github.com/NIST-ISODB/isodb-library` pinned to a commit hash recorded in `data/sources.lock`; JSON parsed by `sitewindow.experimental.isodb`. The clone is cached under `data/raw/isodb-library/` and is not redistributed in the release (NIST license terms permit use; the release contains the derived carbon table with per-record DOI and the commit hash).
- S4 (published values): transcribed from open-access full texts or preprints by the program (programmatic fetch of the open-access PDF where available, text extraction, and value capture with page/table references), stored in `data/registry/published_values.yaml` with one record per value and a `verification` field (`tabular_text`, `digitized_figure`, `abstract_statement`). Values from paywalled articles are included only when the value is stated in a public abstract or in an open-access citing source, with the field recording that provenance.
- S5/S6 (benchmark data from open-access preprints: Meyer 2012 cross sections; Cançado 2011 Raman parameters): extracted from arXiv full texts; where only figures are available, digitized with explicit uncertainty by the program's digitization routine (`sitewindow.experimental.digitize`) with the digitization uncertainty stored.
- DOE targets (S9): transcribed with the URL and retrieval date.

### 6.3 Raw and processed formats

Raw: JSON (ISODB), PDF/TXT (literature), YAML (registry). Processed: Apache Parquet tables with pandas-compatible schemas; posteriors as NetCDF via ArviZ; configuration as YAML; figures as PDF (vector) and PNG (600 dpi); tables as CSV and LaTeX.

### 6.4 Canonical schemas

Defined in `docs/specs/DATA_SCHEMA.md` and implemented as pydantic v2 models in `sitewindow.schema`: `Provenance`, `SiteClass`, `SiteConfiguration`, `SiteEnergyRecord`, `PublishedValue`, `Isotherm`, `Material`, `DefectPopulation`, `OperatingWindow`, `CapacityResult`, `BenchmarkRecord`. Every record carries `kind ∈ {experimental_observation, published_computation, computed_here, model_prediction, derived}` so that the five provenance types are never merged.

### 6.5 Identifiers, validation, deduplication, units

- Identifiers: ISODB filename and DOI; registry record IDs (`PV-####`); configuration hashes (SHA-256 of canonical structure and settings); benchmark IDs (`B#`); source IDs (`S#`).
- Validation: pydantic validation at load; unit fields are enumerated types; physical range checks (loading ≥ 0, temperature > 0, energies within plausible bounds) raise typed exceptions.
- Deduplication: isotherms are deduplicated by (DOI, adsorbent, temperature, adsorbate, data hash); registry values by (DOI, system, site, method); duplicates across digitized and tabular sources prefer tabular.
- Units: stored canonical units per field; conversions in `sitewindow.units` with tests for every supported unit string in ISODB.

### 6.6 Material, structure, and defect representations

Materials: registry name, class (graphene, rGO, GO, SWCNT, MWCNT, activated carbon, templated carbon, other), reported BET area, pore volume, oxygen content where available. Structures: ASE `Atoms` serialized to extended XYZ with cell and periodic flags. Defect classes: enumerated type with a canonical construction recipe (supercell, atoms removed/rotated, functional groups added), so that a class name maps deterministically to a structure. Hydrogen sites: index of carbon atom(s), sublattice, local descriptors, occupancy count. Coverage: hydrogen per carbon in the cell.

### 6.7 Uncertainty representation

Scalars with uncertainty are stored as (value, sigma, distribution_type) or as posterior sample arrays; digitized values carry (value, sigma_digitization); computed values carry (value, numerical_tolerance, tier).

### 6.8 Missing data, versioning, licensing, integrity, caching, regeneration

Missing fields are explicit nulls with reasons where known. Processed datasets are versioned with the package (semantic versioning) and carry the source lock file. Licensing: NIST data (NIST open license terms), arXiv-extracted values (facts, with citation), DOE targets (public domain); the package is released under the MIT License; the derived datasets under CC BY 4.0. Integrity: SHA-256 manifests for every processed artifact (`data/MANIFEST.sha256`). Caching: content-addressed cache directories with a `cache.lock` recording settings; `sitewindow cache verify` recomputes hashes. Regeneration: `sitewindow reproduce --tier {verify, full}` regenerates all processed data from raw sources.

### 6.9 Comparability of heterogeneous quantities

Binding energies are stored with their reference state and never averaged across references; physisorption energies at different temperatures are compared only after Feynman–Hibbs correction to a common temperature; excess and absolute loadings are distinct fields; computed benchmark values are labeled `published_computation` and never used as experimental ground truth.

### 6.10 Treatment of disagreement

Disagreement between sources is modeled (C3 method effects; C6 measurement-error scales) and reported with its attribution; it is not treated as noise to be averaged away. Where two experiments on nominally the same material differ beyond their stated uncertainty, both are retained with a `conflict_group` tag, and the inference reports results with and without each.

## 7. Validation and benchmarking framework

The complete benchmark registry is `docs/specs/BENCHMARK_REGISTRY.md`. Summary:

| ID | Benchmark | Quantity | Used for |
|---|---|---|---|
| B1 | H chemisorption on graphene, DMC and PBE (Dumi et al. 2022) | −0.691 ± 0.019 eV (DMC); −0.821 eV (PBE) | Tier-0 anchor; Tier-1 check |
| B2 | H sticking barrier on graphene (Sha and Jackson 2002; Casolo et al. 2009) | ≈ 0.2 eV | Tier-1/2/3 check |
| B3 | H₂–coronene CCSD(T)/CBS (Kocman et al. 2015) | −4.7 kJ mol⁻¹ | Physisorption anchor |
| B4 | Isosteric heats of H₂ on carbons (AX-21 average 5.8 kJ mol⁻¹ as used by Bhatia and Myers; activated carbons, Knight et al. 2020) | 4–7 kJ mol⁻¹ | Entropy/physisorption consistency |
| B5 | TPD of H/D from graphite (Zecho et al. 2002; Hornekær et al. 2006) | peaks 445/560 K (H), 490/580 K (D) | Forward kinetics |
| B6 | Hydrogenated SWCNT stability versus diameter; complete C–H breaking by 600 °C (Nikitin et al. 2005, 2008) | ordering and release temperature | Curvature kinetics |
| B7 | Surface-area scaling of 77 K excess uptake (Chahine-type rule; corpus) | slope | Corpus consistency |
| B8 | Bhatia–Myers optimum (15.1 kJ mol⁻¹ at 298 K, 30/1.5 bar, ΔS = −8R) | analytic | Theory |
| B9 | Electron knock-on cross sections in graphene (Meyer et al. 2012) | σ(E) for 80–100 keV | Population model |
| B10 | Raman I_D/I_G versus L_D (Lucchese 2010; Cançado 2011); I_D/I_D′ by type (Eckmann 2012) | curve; ratios | Raman coordinate |
| B11 | Vacancy formation/migration energies and defect energetics in graphene (Krasheninnikov and Nordlund 2010 and sources therein) | eV | Engine sanity; rate equations |
| B12 | Graphane formation by atomic H and restoration by annealing at 450 °C (Elias et al. 2009) | reversibility at ≈ 723 K | High-coverage kinetics |

Independence: B1–B3, B11 may enter calibration (C3) and are therefore also tested in leave-one-out mode; B4–B10, B12 are never used for calibration and constitute the external validation set. Thresholds are stated in the registry and fixed before evaluation.
## 8. Statistical analysis and uncertainty quantification

### 8.1 Selected methods and their roles

| Method | Purpose | Model | Assumptions | Inputs | Outputs | Computation | Diagnostics | Effect on conclusions |
|---|---|---|---|---|---|---|---|---|
| Hierarchical Bayesian meta-analysis (C3) | Pool published site energies with method effects | Random-effects regression with partial pooling | Additive method effects; exchangeable studies | Registry values, metadata | Class means, method offsets, heterogeneity τ | NUTS (PyMC) | R̂, ESS, PPC, LOO | Sets the epistemic width of the landscape |
| Multi-fidelity calibration with GP discrepancy (C3) | Map Tier-2/3 energies to Tier-1/0 scale | Linear transfer plus GP over descriptors | Smooth discrepancy in descriptor space | Anchor pairs | Predictive distributions for all configurations | Exact GP, MAP hyperparameters with posterior check | LOO coverage, residual plots | Determines per-site posterior used in every test |
| Measurement-error AED inversion (C6) | Infer site-energy distributions from isotherms | Discretized adsorption integral equation with smoothness prior and error model | Independent sites; known kernel; error model form | Isotherms, metadata | ρ(ε), V_pore, error scales | NUTS | PPC, synthetic recovery | Empirical occupancy evidence (H3, H5) |
| Identifiability analysis | Detect non-identifiable parameters | Profile likelihood; prior–posterior overlap | — | Each model | Flags and reparameterizations | Deterministic | Overlap statistics | Prevents over-interpretation (for example single-T materials) |
| Posterior predictive validation | Check models against held-out data | Replicated data under the posterior | — | Benchmarks | Coverage, p-values | Sampling | Calibration curves (PIT histograms) | Determines whether predictive intervals are trusted |
| Global sensitivity (Sobol) (C7) | Attribute output variance | Variance decomposition | Independent input distributions (correlated inputs handled by sampling the joint posterior and reporting conditional indices) | Posterior/prior samples | S_i, S_Ti with bootstrap CIs | Saltelli (SALib) | Index stability versus N | Identifies which uncertainties matter for D |
| Model comparison | Choose among entropy models, kinetic orders, pathway sets | LOO-ELPD (ArviZ) and benchmark residuals | — | Fitted models | ELPD differences with SE | PSIS-LOO | Pareto-k | Selects the reported model and quantifies model-form uncertainty |
| Bootstrap | Uncertainty of corpus regressions and isosteric heats | Nonparametric | Exchangeable materials | Corpus tables | CIs | Resampling | — | Corpus-level statements |
| Multiple-comparison control | Class-level hypothesis statements | Posterior probabilities reported jointly; decisions by pre-registered thresholds, with the number of classes stated | — | Hypothesis outputs | Posterior tables | — | — | Avoids selective reporting |

### 8.2 Uncertainty taxonomy

- Aleatoric: measurement noise in isotherms; sample-to-sample variation across materials (modeled hierarchically).
- Epistemic: method dependence of site energies (meta-analysis τ, Tier-0 offset); barrier and prefactor uncertainty; entropy model form; population parameters.
- Numerical: convergence tolerances from C2 studies, propagated as additive Gaussian terms.
- Systematic measurement: digitization uncertainty; excess-to-absolute conversion.
- Model discrepancy: GP discrepancy in calibration; residual benchmark discrepancies reported, never absorbed.

### 8.3 Validation partitions

Calibration uses Tier-0/1 anchors; out-of-distribution validation holds out entire classes (for example, calibrate without the curvature series and predict it; calibrate without oxygen classes and predict them). Corpus inference uses leave-one-source-out (by DOI) to test the error model. No benchmark in the external set (B4–B10, B12) is used to tune any parameter.

### 8.4 Inferences the evidence cannot support

If posterior widths for chemisorption barriers exceed ±0.3 eV after calibration, the admissibility classification of sites near the threshold is reported as undetermined for those sites and the hypothesis statements are made conditional. If the corpus contains too few multi-temperature datasets for a material class, AEDs for that class are reported under the fixed-entropy model only. Such delimitations are stated in the manuscripts' limitations sections.

## 9. Computational experiments and scientific discovery

Each experiment specifies purpose, controlled and independent variables, dependent quantities, comparisons, statistical treatment, and acceptance criteria. Experiments are implemented as workflow commands with configuration files under `configs/experiments/`.

**E1. Landscape construction and calibration.** Purpose: produce the calibrated landscape. Controlled: calculator settings per tier, supercell protocol. Independent: class, configuration, tier. Dependent: ε, barriers, descriptors, posterior predictive distributions. Comparisons: tier-versus-anchor residuals by class; held-out class predictions. Statistics: C3. Acceptance: C2 and C3 criteria.

**E2. Curvature series.** Purpose: quantify curvature dependence of chemisorption energy and recombination barrier; test H2. Controlled: site type (exterior top, ortho/para dimers), tier, coverage. Independent: diameter (≥ 8 values), chirality (armchair/zigzag). Dependent: ε(d), E_d†(d), T₉₀(d). Comparisons: B6 ordering; graphene limit. Statistics: posterior monotonicity probability (fraction of posterior draws with monotone decreasing T₉₀ over d). Acceptance: H2 criteria; convergence to graphene within 0.05 eV at 2.5 nm.

**E3. Window maps and bounds.** Purpose: compute Δθ, A, η(W), ε*, ε_c, D_max for W1–W4 and a two-dimensional grid of windows (T_r ∈ [233, 600] K, P_c/P_r ∈ [2, 100]). Controlled: entropy model; prefactor prior. Dependent: maps. Comparisons: B8; model-set spread. Acceptance: C1 criteria; maps generated as vector figures.

**E4. Class-resolved deliverable and dischargeable capacity versus defect density.** Purpose: evaluate D_c(x_c; W) for each class up to the amorphization bound; test H1, H4, H6. Controlled: window; admissibility threshold. Independent: class, x_c, mixed-class configurations. Dependent: D, dischargeable capacity, limiting mechanism. Statistics: posterior probability statements per hypothesis; Sobol indices. Acceptance: hypothesis statements computed; bound never violated.

**E5. Physisorption landscape and heterogeneity threshold.** Purpose: compute physisorption site energies for basal, pores, tube geometries, defect- and oxygen-adjacent sites with Feynman–Hibbs corrections at 77, 233, 298 K; compute ε_c(W); test H3. Controlled: potential parameter set (two alternatives). Dependent: ε by site; ΔD gain. Comparisons: B3, B4. Acceptance: H3 criteria; model-form spread reported.

**E6. Forward desorption benchmarks.** Purpose: predict TPD of H/D on graphite and release temperatures of hydrogenated nanotubes from the landscape; test C4 and H2/H4. Controlled: ramp rate per source; initial coverage per source. Dependent: spectra, peak temperatures, release fractions. Comparisons: B5, B6, B12. Statistics: posterior predictive residuals. Acceptance: C4 criteria.

**E7. Irradiation populations.** Purpose: compute class densities versus (E, Φ) for graphene and nanotubes, Raman coordinates, and D versus I_D/I_G; test H6. Controlled: temperature; substrate. Independent: E ∈ {5, 20, 30, 60, 80, 100, 200, 300} keV; Φ up to the amorphization limit. Dependent: x_c, I_D/I_G, I_D/I_D′, D. Comparisons: B9, B10. Acceptance: C5 criteria; mixture-bound check.

**E8. Corpus inference.** Purpose: AEDs and isosteric heats for all carbon materials with sufficient data; corpus regressions; test H5. Controlled: error-model hyperparameters (shared). Independent: material. Dependent: ρ(ε), q_st(n), occupancy, outlier class. Comparisons: B4, B7, B8; comparison with E5 landscape. Acceptance: C6 criteria.

**E9. Counterfactual and regime analysis.** Purpose: identify what would have to change for carbon to reach targets (barrier reductions, entropy changes, site densities) and map the regime diagram; negative controls (pristine graphene, which must give near-zero deliverable capacity at W1; a hypothetical homogeneous site at ε* which must recover D_max). Dependent: D under counterfactuals; regime labels. Acceptance: negative controls pass; regime diagram produced with posterior-labeled boundaries.

**Screening statement.** Where classes or populations are ranked, computationally promising configurations are labeled as such and distinguished from experimentally demonstrated materials (only graphite, graphene, SWCNT, and activated carbons have experimental benchmarks in this program). Apparent optima are tested for survival under posterior uncertainty, the amorphization constraint, mass accounting, and alternative entropy and kinetic models before being reported as optima. Derived relations (for example T₉₀(d)) are reported with their domain of validity and the diameter range where the landscape is converged.
## 10. Scientific software architecture

### 10.1 Package and boundaries

One package, `sitewindow`, with internally bounded subpackages whose public APIs are stable enough that any of them could be extracted into a separate distribution without changing callers. The scientific purpose of the package is the construction, calibration, and evaluation of site-energy landscapes for hydrogen in sp² carbon against operating windows.

```
sitewindow/
  __init__.py            version, public re-exports
  constants.py           physical constants (CODATA 2018), unit factors
  units.py               unit parsing and conversion with tests per ISODB unit string
  schema/                pydantic v2 data models (Section 6.4)
  theory/                langmuir, dissociative, entropy, fugacity, window, deliverable, bounds
  kinetics/              polanyi_wigner, pathways, tpd_forward, admissibility
  landscape/             registry, meta_analysis, calibration, assemble, classes (site-class catalogue)
  atomistic/             structures/, calculators/, workflows/, descriptors, cache
  defects/               knockon, thermal, rates, populations, raman
  experimental/          isodb, select, units, excess, isosteric, aed, digitize
  evaluate/              maps, populations, pareto, counterfactual, hypotheses, sobol
  validation/            registry loader, runners (one per benchmark), report
  plotting/              figure functions (matplotlib), style, export
  workflows/             CLI (typer): fetch, build, calibrate, evaluate, validate, reproduce, cache
  io/                    parquet/netcdf/yaml readers and writers, manifests
```

### 10.2 Public API (stable)

- `sitewindow.theory.swing(eps, window, entropy_model) -> float | ndarray`
- `sitewindow.theory.deliverable(landscape, window, admissibility) -> CapacityResult`
- `sitewindow.theory.bound(window, channel) -> BoundResult`
- `sitewindow.kinetics.release_fraction(E_d, nu, order, T, tau, theta0)`
- `sitewindow.kinetics.tpd(spectrum_spec, landscape_samples) -> TPDResult`
- `sitewindow.landscape.build(config) -> Landscape`; `Landscape.sample(n) -> LandscapeSamples`
- `sitewindow.defects.population(E_keV, fluence, material, params) -> DefectPopulation`
- `sitewindow.defects.raman(population, excitation_eV) -> RamanRatios`
- `sitewindow.experimental.load_isodb(path) -> IsothermTable`; `aed.fit(material_isotherms, config) -> AEDPosterior`
- `sitewindow.evaluate.hypotheses(results) -> HypothesisReport`
- `sitewindow.validation.run(benchmark_ids, tier) -> ValidationReport`

### 10.3 Internal interfaces and extension points

- `Calculator` protocol (`atomistic.calculators.base`): `energy_forces(atoms) -> (E, F)`, `tier`, `settings_hash`; implementations for tblite, MACE, GPAW, PySCF; new engines register through an entry point group `sitewindow.calculators`.
- `SiteClass` catalogue (`landscape.classes`): construction recipe, descriptors, maximum occupancy, mass change; new classes register through `sitewindow.site_classes`.
- `EntropyModel` protocol (`theory.entropy`) and `Pathway` protocol (`kinetics.pathways`).
- `Benchmark` protocol (`validation`): `load()`, `predict(pipeline)`, `metric()`, `threshold`.
- Data contracts: Parquet schemas and NetCDF variable names versioned (`schema_version` field); backward-compatible readers for prior minor versions.

### 10.4 Numerical kernels, precision, determinism

Float64 throughout; log-domain evaluation for Boltzmann factors; explicit tolerances in configuration; deterministic RNG by seed (NumPy `default_rng(seed)`, PyMC `random_seed`), with seeds recorded in every output's provenance block; torch set to deterministic single-threaded mode for MACE where determinism is required (tests) and multi-threaded for production with results compared at tolerance.

### 10.5 Configuration, CLI, logging, errors, resources, caching, parallelization, checkpointing

- Configuration: YAML files validated by pydantic; CLI overrides; `configs/default.yaml`, `configs/windows.yaml`, `configs/experiments/E*.yaml`, `configs/tiers/{verify,full}.yaml`.
- CLI: `sitewindow fetch isodb`, `sitewindow landscape build --tier`, `sitewindow landscape calibrate`, `sitewindow evaluate maps`, `sitewindow evaluate populations`, `sitewindow corpus fit`, `sitewindow validate --benchmarks`, `sitewindow reproduce --tier`, `sitewindow cache verify`, `sitewindow figures`.
- Logging: standard `logging` with structured records (component, configuration hash, step, wall time); no progress narration in persisted logs beyond timing and settings.
- Errors: typed exceptions (`SiteWindowError` → `UnitError`, `SchemaError`, `ConvergenceError`, `CalculatorUnavailable`, `BenchmarkFailure`); input validation at every public entry point; numerical stability checks (NaN/Inf guards, condition-number checks for GP kernels).
- Resources: `--workers` for process-level parallelism (joblib) across configurations, materials, and posterior chunks; memory-aware chunking of posterior arrays; `--max-memory-gb`.
- Caching: content-addressed directories keyed by SHA-256 of inputs and settings; cache verification and pruning commands.
- Checkpointing: per-structure results written atomically on completion; NEB and relaxations resume from trajectory files; MCMC traces saved per chain; `reproduce` is idempotent.

### 10.6 External data contracts (no code dependency)

- TPD spectra export (`tpd_forward.nc`: temperature grid, rate posterior quantiles, pathway decomposition) for use by any inversion tool.
- Landscape export (`landscape_calibrated.parquet`) for use by system-level or spectroscopic models.
- Raman predictions export (`raman_predictions.parquet`).
- Corpus tables with ISODB identifiers.

### 10.7 Dependency selection

| Dependency | Role | Rationale |
|---|---|---|
| numpy, scipy | kernels, ODEs, quadrature, optimization | standard, stable |
| pandas, pyarrow | tables | Parquet interoperability |
| pydantic ≥ 2 | schemas | validation, JSON schema export |
| ase | structures and optimizers | de facto standard |
| tblite | GFN2-xTB with PBC | pip-installable, fast |
| mace-torch, torch (CPU) | universal potential | pip-installable; optional extra `[mace]` |
| gpaw (conda-forge) / pyscf (pip) | Tier-1 DFT | optional extras `[gpaw]`, `[pyscf]` |
| dftd3-python (conda-forge) | D3(BJ) dispersion for Tier-1 energies and forces | ASE-compatible; published PBE damping parameters |
| pymc, arviz | Bayesian inference | mature NUTS, diagnostics |
| SALib | Sobol indices | standard |
| matplotlib | figures | vector output |
| typer | CLI | typed CLI |
| pyyaml | configuration | standard |
| CoolProp (optional) | H₂ equation of state | fugacity; fallback virial |
| pytest, hypothesis, ruff, mypy | testing and static analysis | standard |

### 10.8 Packaging, versioning, releases, DOI preparation

`pyproject.toml` (PEP 621), `src/` layout, semantic versioning, `CHANGELOG.md` (Keep a Changelog format), `CITATION.cff`, `codemeta.json`, `.zenodo.json`, MIT `LICENSE`, GitHub Actions for tests, lint, type-check, build, and release archives; `sitewindow --version` prints version and git describe. Release procedure in Section 15.

## 11. Software engineering standards

- Package metadata per PEP 621; dependencies pinned by lower bounds in `pyproject.toml` and exactly in `requirements-lock.txt` (generated by `pip-compile`) and `environment.yml` (conda-forge, for GPAW).
- Typed public interfaces (`py.typed`), mypy strict on `sitewindow.theory`, `kinetics`, `schema`, `units`; ruff formatting and linting; docstrings (NumPy style) for every public function with units and references to manual equations.
- Explicit input validation; typed exceptions; numerical stability checks; deterministic seeds.
- Tests (pytest): unit (analytic identities, unit conversions for every ISODB unit string, schema validation), numerical (convergence of quadrature and ODE against closed forms), physics (limiting cases in Section 4; bound theorem randomized tests; Redhead check; McKinley–Feshbach against tabulated values; Raman relation reproduction), regression (stored outputs for fixed seeds), integration (end-to-end verify-tier pipeline on a reduced configuration set), reproducibility (hash of regenerated verify-tier outputs equals stored hash). Property-based tests (hypothesis) for monotonicity and bounds of θ, Δθ, f.
- Continuous integration: GitHub Actions matrix (Ubuntu; Python 3.11, 3.12, 3.13); fast tier (< 15 min) on every push; `[mace]` extra tested on release tags and on manual dispatch; build of sdist and wheel; `twine check`; installation verification in a clean environment.
- Coverage is reported but is not an acceptance criterion; scientific test content is (Section 18).

## 12. Documentation and reproducibility

Documentation (`docs/`, MkDocs with mkdocstrings): scientific overview; installation (pip, conda-forge environment for GPAW); environment requirements; quick start (window tables and a deliverable-capacity map in under a minute); usage guides per component; API reference (auto-generated); data schema; theory (window theory, kinetics, populations, inference) with equation numbering matching `DERIVATIONS.md`; validation (benchmark registry and latest report); benchmark reproduction commands; workflow examples (E1–E9); interpretation of outputs; known limitations; extension development (adding a calculator, a site class, a benchmark, a pathway); citation; manuscript reproduction (`docs/reproduce/manuscript_1.md`).

Reproducibility workflow: `sitewindow reproduce --tier verify` regenerates all processed data, figures, and tables from raw sources on a 2-core machine in under two hours using reduced configuration sets, and labels every output "verification tier"; `sitewindow reproduce --tier full` regenerates the publication outputs with the full configuration sets, recording wall time and hardware. Outputs carry the tier label in file names and figure captions. No machine-specific paths; all paths relative to a configured `SITEWINDOW_DATA` root; the ISODB commit hash, package version, dependency lock, seeds, and settings hashes are written to `reproducibility_manifest.json` beside every output set.

## 13. Scientific figures and tables

Conventions: vector PDF plus 600-dpi PNG; colorblind-safe palettes (Okabe–Ito for categorical; viridis/cividis for sequential; a diverging map for signed quantities); units in axis labels; 90% posterior bands shaded with the median as a line; panel labels (a), (b), ...; consistent typography (Helvetica-like sans, 8–9 pt at column width); significant figures matched to posterior precision; no decorative elements. All figures generated by `sitewindow figures --manuscript 1` from cached results; captions carry the result hashes.

### 13.1 Manuscript 1 figure set

| Fig. | Purpose and question | Quantities and panels | Data/workflow | Uncertainty | Role |
|---|---|---|---|---|---|
| 1 | Define the window and the landscape gap | (a) Δθ(ε; W1, W2, W3) with ε* and ε_c marked; (b) admissibility boundary in (E_d†, T_r); (c) schematic of channels on one energy axis with class ranges | E3 | model-set spread as bands | Frames the question |
| 2 | Calibrated site-energy landscape | Per-class posterior violins of ε (chemisorption, per H vs ½H₂) and of E_d†; anchors overlaid; curvature series as a function of 1/d | E1, E2 | posterior | Occupancy evidence |
| 3 | Deliverable-capacity maps | D(ε, E_d†; W) for W1, W2, W4(473 K); regime labels; class centroids overlaid | E3, E4 | contour of posterior median with 90% boundary | Core result |
| 4 | Class-resolved D versus defect density and Raman coordinate | D_c(x_c) up to amorphization; secondary axis I_D/I_G at 2.41 eV; populations from E7 as trajectories | E4, E7 | bands | Realizability |
| 5 | Curvature–kinetics benchmark | T₉₀(d) with posterior band; benchmark release temperatures; predicted TPD for graphite versus B5 | E2, E6 | bands; benchmark error bars | Validation |
| 6 | Corpus evidence | AED posteriors for representative carbons; window-occupancy fraction histogram across the corpus; 77 K uptake versus surface area with B7 line | E8 | posterior | Independent evidence |
| 7 | Bound, Pareto, and counterfactuals | D versus stability margin; D_max(W) across windows; what-if barrier reduction | E3, E9 | bands | Synthesis |

Supplementary figures: convergence (cell, k, cutoff, vacuum, NEB images, grids, MCMC), calibration residuals and LOO coverage, Sobol indices, PIT histograms, entropy-model comparison, knock-on cross section versus B9, Raman relation versus B10, ablations (no dimer pathways; no Tier-0 offset; no measurement-error model), full class catalogue.

### 13.2 Manuscript 1 table set

T1 operating windows; T2 site-class catalogue (structure, sites, multiplicity, n_max, Δm); T3 calibrated class energies and barriers (posterior median, 90% interval, tier composition); T4 D_max(W) and η(W) per window; T5 hypothesis outcomes with posterior probabilities; T6 benchmark report (B1–B12 with residuals); ST tables: convergence, registry of published values, corpus materials with metadata and exclusions, Sobol indices.

### 13.3 Manuscript 2 figures (theory)

Analytic swing functions and optimum conditions; heterogeneity theorem illustrations; admissibility surfaces; bound tightness examples; isotope and prefactor sensitivity; correspondence with Bhatia–Myers limits.
## 14. Manuscript plan

Two manuscripts are planned on the basis of the distinct contributions: a theory of deliverable capacity for heterogeneous, kinetically constrained adsorbents (self-contained, methodological), and the application of the full program to carbon with its scientific conclusions (the major manuscript). A third, corpus-focused manuscript is contingent on E8 producing findings that stand independently of Manuscript 1 (for example a corpus-wide measurement-artifact diagnosis); it is not planned by default.

### 14.1 Manuscript 1 (major)

- **Working title.** The deliverable-capacity limit of defect-engineered sp² carbon for hydrogen storage from calibrated site-energy landscapes.
- **Central question.** Section 1.3.
- **Contribution.** A class-resolved, uncertainty-quantified determination of deliverable and dischargeable capacity for defect-engineered carbon across practical windows; identification of window-populating and window-missing classes; curvature–kinetics relation validated against experiment; corpus-level evidence; a variational bound evaluated on real landscapes.
- **Knowledge gap.** Section 2.3.
- **Relationship to literature.** Section 2.1 and the bibliography in `docs/specs/SOURCE_REGISTRY.md` (verified references only).
- **Hypotheses.** H1–H6.
- **Claim–evidence framework.**

| Anticipated claim (to be confirmed or refuted) | Experiments | Validation | Uncertainty analysis | Figures |
|---|---|---|---|---|
| The molecular and atomic channels of pure carbon leave the W1/W2 window sparsely populated; D(W1), D(W2) bounded as computed | E1, E3, E4, E5 | B1–B4, B8 | posterior P statements; Sobol | 1, 2, 3, 4 |
| Curvature lowers release temperature monotonically; flat-sheet release is dimer-limited | E2, E6 | B5, B6, B12 | posterior bands versus benchmarks | 5 |
| Irradiation populations are mixtures bounded by their best class; D per Raman unit is bounded | E7, E4 | B9, B10 | posterior; mixture test | 4 |
| Real carbons show no window-populating physisorption sites beyond the computed threshold | E8, E5 | B4, B7 | PPC, LOO | 6 |
| The bound D_max(W) is tight only for homogeneous landscapes; counterfactual requirements for carbon to reach targets are quantified | E3, E9 | analytic tests | model-set spread | 7 |

- **Structure.** Abstract (background, gap, approach, result, implication); Introduction; Results (landscape; window theory applied; class-resolved capacity; populations and Raman coordinate; kinetics benchmarks; corpus evidence; bound and counterfactuals); Discussion (what the limit means for defect engineering; what would break it; implications for irradiation-based hydrogenation; limitations); Methods (each component, concisely, with references to the software and SI); Data availability; Code availability; Author contributions; Competing interests; References; Supplementary Information (derivations summary, convergence, calibration, full tables, extended benchmarks, ablations).
- **Venue category.** Multidisciplinary or energy-materials journal of the Nature family (Nature Communications, Nature Energy) or a leading materials/physics journal; the determination is made after results exist and is not a claim of acceptance.
- **Limitations to be discussed.** Enumerated configuration space; mean-field kinetics; entropy model dependence; Tier-0 offset transferability; sub-threshold hydrogenation phenomenology; mass-only (no volumetric) capacity; pure/oxygen-functionalized carbon only.
- **Acceptance criteria.** All claims mapped to executed experiments with passing or explicitly discussed benchmarks; every figure regenerated from cached results by one command; every reference verified (Section 18); SI complete; LaTeX compiles reproducibly.

### 14.2 Manuscript 2 (theory)

- **Working title.** Deliverable capacity of energetically heterogeneous adsorbents with kinetic constraints: a window theory with application to hydrogen.
- **Contribution.** Closed-form swing functions for molecular and dissociative channels; optimum and threshold energies for general windows including temperature swings; the heterogeneity theorem (when added strong sites help); kinetic admissibility and the resulting regime diagram; variational bound; isotope and prefactor dependence; recovery of Bhatia–Myers as a special case.
- **Evidence.** Derivations (`DERIVATIONS.md`), analytic limit tests, numerical verification, illustrative evaluations on the carbon landscape of Manuscript 1.
- **Structure.** Abstract; Introduction; Theory (definitions; molecular channel; dissociative channel; kinetics; mass accounting; bound); Results (maps, regimes, theorem illustrations); Discussion; Methods (numerical verification); Appendices.
- **Venue category.** Physical chemistry or chemical engineering thermodynamics journal.
- **Limitations.** Independent-site assumption; mean-field kinetics; entropy model.

### 14.3 Manuscript files

LaTeX source tree under `manuscripts/m1/` and `manuscripts/m2/` (`main.tex`, `sections/*.tex`, `figures/` symlinked to generated outputs, `tables/` generated, `si/main_si.tex`, `references.bib` with DOI fields for every entry, `latexmkrc`), compiled by `latexmk -pdf` in CI on release tags; `references.bib` is generated from the verified reference registry so that no unverified entry can enter.

## 15. Software publication and DOI-ready releases

Repository: Zenodo (DOI-issuing, supports direct software deposit via GitHub integration or manual upload, no journal peer review required). Deposit requirements verified against Zenodo's published metadata schema: title, creators with affiliations and ORCID where available, description, license, version, related identifiers (preprint DOIs when they exist), keywords, language, upload type `software`. `.zenodo.json` carries this metadata; `codemeta.json` carries machine-readable software metadata; `CITATION.cff` carries citation metadata.

Release contents per version: tagged source archive; `sdist` and `wheel`; `README.md` (purpose, installation, quick start, documentation link, citation, license); `CHANGELOG.md`; `RELEASE_NOTES.md` per tag; `docs/` built site; `data/processed/` derived datasets with `MANIFEST.sha256`; `reproducibility_manifest.json`; `validation_report.json` from the full tier; checksums file for all archive members.

Release procedure: (1) all CI jobs pass on the tag; (2) `sitewindow reproduce --tier verify` output hashes match the stored manifest; (3) full-tier validation report attached; (4) `python -m build`, `twine check`; (5) archive assembled by `scripts/make_release.py`, which refuses to proceed if any placeholder, unresolved reference, or disallowed content pattern (Section 18) is found; (6) consolidated handoff bundle `release/handoff/` containing the archive, metadata files, manuscript PDFs and sources, figures, SI, and a checklist of the account-bound actions (Zenodo deposit confirmation, preprint submission), which are the only actions not performed by the program. No DOI is written into any file until the deposit has been made; the release metadata reserves the field.

## 16. Computational resource architecture

| Workload | Class | Complexity | Memory | Storage | Parallelism | Bottleneck | Mitigation |
|---|---|---|---|---|---|---|---|
| Theory maps (C1, E3) | essential | O(grid × samples) | < 1 GB | < 100 MB | vectorized | none | — |
| Tier-2 landscape (C2) | essential | O(N configs × SCF) | < 2 GB | < 1 GB | embarrassingly parallel per configuration | SCF for spin-polarized vacancy cells | caching; reduced enumeration in verify tier |
| Tier-3 screening (C2) | essential for curvature/cluster sets | O(N × force calls) | < 4 GB | < 1 GB | per configuration | NEB | 7 images; resume |
| Tier-1 anchors (C2) | essential (anchor subset), refinement (extended set) | O(N_anchor × DFT) | < 8 GB | < 2 GB | per anchor | DFT SCF | LCAO relaxations, PW single points; GPAW; PySCF fallback for verify tier |
| Calibration (C3) | essential | O(MCMC) | < 2 GB | < 200 MB | chains | NUTS | non-centered parameterization |
| Kinetics (C4) | essential | O(samples × ODE) | < 1 GB | < 200 MB | samples | stiff ODE | BDF; analytic forms where possible |
| Populations (C5) | essential | negligible | — | — | — | — | — |
| Corpus (C6) | essential | O(materials × MCMC) | < 2 GB | < 1 GB (raw clone 0.5 GB) | materials | NUTS | per-material checkpoints |
| Evaluation (C7) | essential | O(samples × classes × windows) | < 4 GB | < 1 GB | chunks | Sobol sample size | adaptive N with stability check |
| Validation (C8) | independent validation | varies | — | — | — | — | — |
| Extended Tier-1 set beyond anchors | optional refinement | — | — | — | — | — | not required for conclusions |

Total essential compute on an 8-core workstation: Tier 2 ≈ 50–150 core-hours; Tier 3 ≈ 100–300 core-hours; Tier 1 anchors ≈ 500–1500 core-hours (60–180 wall-hours on 8 cores); statistics and evaluation ≈ 20–50 core-hours. All essential workloads are checkpointed and resumable. Verification tier (2 cores, < 2 wall-hours) uses: a reduced class set (basal monomer and dimers, monovacancy 1–3 H, one SW, three diameters), Tier 2 only with a PySCF single anchor, reduced posterior sample counts, and a corpus subset (20 materials), and is labeled as verification in all outputs. Alternative implementations if dependencies are unavailable: PySCF for GPAW; GFN1-xTB for GFN2; exclusion of Tier 3 with wider calibration uncertainty if MACE cannot be installed (conclusions then rely on Tier 2 plus anchors, with the reduced configuration coverage reported).

## 17. Failure modes and scientific contingency design

| Failure mode | Detection | Response |
|---|---|---|
| ISODB mirror unavailable | fetch failure | use the pinned local clone from the release archive's lock file; if absent, the NIST web API (documented alternative) with the same parser; corpus component delayed but all other components unaffected |
| Insufficient multi-temperature isotherms for a material class | identifiability flags | report fixed-entropy AEDs; delimit H5 to materials with multi-T data |
| Incompatible descriptors across site classes for the GP | residual structure, poor LOO coverage | class-specific GPs; report calibration status per class |
| Non-identifiable (a, b) in calibration | profile likelihood | fix b = 1; widen uncertainty; flag |
| Poor out-of-class predictive performance | held-out class test | do not extrapolate; compute Tier-1 anchors for the class; report |
| Numerical instability (SCF non-convergence, NEB failure) | convergence errors | alternative mixing/smearing; constrained relaxation; exclude with log |
| Resource overrun | timing logs | reduce enumeration per class with documented coverage; keep anchors |
| Invalid approximation (e.g., mean-field kinetics for clustered hydrogen) | benchmark residuals (B5, B12) | add explicit multi-species pathways; if still failing, state the limitation and restrict claims to classes where the forward model validates |
| Inconsistent experimental benchmarks (e.g., spread in release temperatures) | conflict groups | report comparison against each with posterior residuals; no selection |
| Unreliable parameter extrapolation (windows outside validated range) | domain flags | restrict reported maps to validated domains; label extrapolations |
| Weak uncertainty calibration (coverage < 85% or > 95%) | PIT, coverage | inflate/deflate error model hierarchically; report pre- and post-recalibration |
| Inadequate transferability (MACE on a class) | per-class validation | exclude class from Tier 3 |
| Hypothesis contradicted | posterior statements | report falsification as a result; the bound and landscape remain findings |
| Underdetermined inverse problem (AED) | PPC fine, posterior wide | report resolution limits; use the posterior as is; do not sharpen by prior tuning |

Scientific status of failures: a negative or bounded result from any hypothesis is a finding and is reported as such. Thresholds in the benchmark registry are not changed after results are seen.

## 18. Comprehensive acceptance criteria

The consolidated criteria are in `docs/specs/ACCEPTANCE_CRITERIA.md`; summary:

- **Scientific and mathematical correctness.** All analytic identities and limiting cases pass automated tests at the stated tolerances; derivations in `DERIVATIONS.md` are verified numerically in the test suite (each numbered equation with a test).
- **Numerical stability and convergence.** Convergence studies documented for every discretization and sampler with target tolerances met (C2: ≤ 0.02 eV Tier 1; C6 grid; MCMC R̂ < 1.01).
- **Validation and benchmark accuracy.** B1–B12 executed with reports; external-set benchmarks pass pre-registered thresholds or are discussed as discrepancies with mechanistic attribution.
- **Uncertainty calibration.** LOO/held-out coverage of 90% intervals within [0.85, 0.95].
- **Physical consistency.** Carbon conservation in rate equations; θ ∈ [0, 1]; Δθ bounded by Δθ_max; bound never violated; graphene limit of curvature series.
- **Data integrity.** Manifest hashes verified; every record validated; exclusions logged.
- **Reproducibility.** `reproduce --tier verify` regenerates outputs with matching hashes on a clean machine in CI.
- **Software functionality and installation.** `pip install .` in a clean environment; CLI commands execute; extras documented.
- **Test quality.** Every public function has a test of numerical or physical content; tests of existence only are not counted.
- **Documentation.** Every CLI command and public API documented; theory pages reference equations; limitations page present.
- **Computational feasibility.** Full-tier wall times recorded and within the Section 16 envelope on the reference workstation description.
- **Manuscript completeness.** All sections present; all figures and tables generated; SI complete; LaTeX compiles.
- **Reference accuracy.** Every bibliography entry has a verified DOI or arXiv identifier recorded in the reference registry with a verification note; no entry without verification.
- **Publication metadata and release integrity.** `.zenodo.json`, `codemeta.json`, `CITATION.cff` validated; archive checksums match; release script content check passes (no placeholders, no tool attribution, no time-based plans, no delegation language).
- **Extensibility.** A new site class, calculator, pathway, and benchmark can each be added by implementing the documented protocol without modifying core modules (demonstrated by an example plugin in `examples/`).

## 19. Open-ended research and extension architecture

### 19.1 Architectural extensibility

The extension points (Section 10.3) are protocols with registered entry points: calculators, site classes, entropy models, pathways, benchmarks, and operating windows. The landscape and capacity data contracts are versioned. Each extension below connects through one or more of these protocols and through the exported data contracts (Section 10.6), without modification of the core.

### 19.2 Extension landscape

**X1. Heteroatom-doped and metal-decorated carbon landscapes.** Question: do B/N substitution or alkali/transition-metal decoration populate the window, and at what mass and stability cost? Unresolved: scattered DFT claims without kinetic admissibility or mass accounting. Methods: new site classes with Tier-1 anchors; dispersion-corrected DFT for Kubas-type binding; the same calibration and evaluation. Data: published values; ISODB isotherms for doped carbons. Validation: B3-type CCSD(T) anchors for doped coronene (Kocman et al. provide a B/Li-modified coronene reference), experimental isosteric heats for boron-substituted carbons. Contribution: extension of the limit to the full carbon-derived materials space. Deliverable: `sitewindow-dopants` plugin and a manuscript. Interfaces: `SiteClass`, `Calculator`. Independent of experiments; feasible at workstation scale with the same tiers.

**X2. Spillover kinetics on the landscape.** Question: can catalyst-generated atomic hydrogen populate carbon sites at useful rates and coverages at ambient temperature, given migration and recombination barriers? Unresolved: the long-standing reproducibility controversy. Methods: kinetic Monte Carlo on lattices with landscape-derived site energies and barriers, catalyst particles as boundary sources; comparison with published uptake kinetics. Data: published spillover uptake curves (open-access where available); landscape. Validation: diffusion barriers against published values; uptake versus time. Contribution: quantitative resolution of the achievable spillover capacity. Deliverable: `sitewindow-kmc` and a manuscript. Interfaces: landscape export; `Pathway`.

**X3. Coupled irradiation–hydrogenation kinetics with spectroscopic forward models.** Question: what dose–energy–adsorbate conditions maximize usable hydrogen in irradiated carbons, and what Raman and XPS signatures (D/G, D/D′, C 1s sp³ fraction, Auger D-parameter) result? Methods: spatially resolved kinetic Monte Carlo of displacement and adsorbate fragmentation; forward Raman and XPS models from site populations. Data: published irradiation and spectroscopy results. Validation: published dose dependence of D/G and sp³ fractions. Contribution: predictive processing maps for irradiation-based hydrogenation. Deliverable: plugin and manuscript. Interfaces: `defects` populations; data contracts to external spectroscopy tools.

**X4. Volumetric capacity in realistic porous carbon models.** Question: how do site heterogeneity and pore structure jointly determine volumetric deliverable capacity? Methods: grand canonical Monte Carlo (open-source engine or in-package) on amorphous and templated carbon models with Feynman–Hibbs potentials; coupling to the window theory through site-resolved occupancy. Data: ISODB with pore-size distributions; published models. Validation: 77 K and 298 K isotherms of AX-21, MSC-30, zeolite-templated carbons. Contribution: extends the limit to volumetric targets. Deliverable: `sitewindow-gcmc` and a manuscript. Interfaces: `EntropyModel`, landscape export.

**X5. Quantum nuclear and isotope effects.** Question: how do zero-point and tunneling corrections shift the window for H₂, D₂, and HD, and does isotope selectivity in carbons follow? Methods: path-integral corrections to physisorption; isotope-dependent barriers; validation against published D₂/H₂ selectivity in carbons. Contribution: isotope-resolved landscape. Interfaces: `EntropyModel`, kinetics.

**X6. Active-learning interatomic potential for hydrogenated defective carbon.** Question: can a dedicated potential reach Tier-1 accuracy for the full landscape at Tier-3 cost? Methods: fine-tuning of a universal potential on the Tier-1 anchor set with uncertainty-driven acquisition; validation on held-out anchors. Contribution: a publishable potential and the extension of the enumerated configuration space by orders of magnitude. Deliverable: potential files and a methods manuscript. Interfaces: `Calculator`.

**X7. Beyond-carbon two-dimensional hosts.** Question: do h-BN, MoS₂, or covalent organic frameworks present landscapes that populate the window? Methods: same pipeline with new site classes and anchors. Contribution: comparative limits across host families. Interfaces: `SiteClass`.

**X8. Operando thermal-management coupling.** Question: how do heat-of-adsorption distributions from the landscape constrain charging rates in packed beds? Methods: coupling the exported landscape to a bed-scale heat and mass transfer model. Contribution: links material limits to system limits. Interfaces: landscape export (data contract only).

Each extension is a scientific direction with its own question, methods, validation, software deliverable, and manuscript potential; none is scheduled.
