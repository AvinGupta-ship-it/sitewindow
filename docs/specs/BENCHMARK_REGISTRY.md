# Benchmark Registry

Thresholds are fixed before evaluation. The machine-readable copy is `data/registry/benchmarks.yaml`; `sitewindow validate` reads thresholds from that file. A change to any threshold requires a changelog entry with justification and is not permitted after the corresponding benchmark has been run for a reported result.

Columns: source; identifier; physical system; quantity; conditions; acquisition; prediction method; metric; uncertainty treatment; independence; acceptance; limitations.

## B1. Chemisorption energy of a single H atom on graphene (high-level reference)

- Source: Dumi et al. 2022 (S7). Identifier: DOI 10.1063/5.0085982.
- System: 5×5 graphene supercell, H on top site.
- Quantity: binding energy relative to free H: DMC −0.691 ± 0.019 eV; PBE (plane wave) −0.821 eV; literature PBE/PW91 range −0.79 to −0.98 eV (as summarized in the source).
- Conditions: 0 K electronic energies, no ZPE.
- Acquisition: tabular values from the open-access preprint.
- Prediction: Tier-1 PBE-D3 single H on 5×5 and 6×6 cells (D3 contribution reported separately); Tier-0 offset estimated as (DMC − PBE).
- Metric: absolute difference.
- Uncertainty: DMC statistical error; Tier-1 numerical tolerance.
- Independence: used as a Tier-0 anchor in calibration; also evaluated in leave-one-out mode.
- Acceptance: Tier-1 PBE value within 0.10 eV of −0.821 eV (same functional, converged settings); calibrated posterior for the basal monomer covers −0.691 eV within its 90% interval.
- Limitations: single configuration; DMC finite-size effects as stated by the source.

## B2. Sticking barrier for H on graphene

- Source: Sha and Jackson 2002; Casolo et al. 2009 (barriers 0.0–0.2 eV for sequential sticking, 0.2 eV for the monomer at PBE).
- Quantity: barrier height ≈ 0.2 eV (monomer).
- Prediction: Tier-1 constrained relaxation and Tier-2/3 NEB.
- Metric: absolute difference.
- Independence: anchor (calibration) with leave-one-out check.
- Acceptance: Tier-1 within 0.05 eV of 0.2 eV; calibrated Tier-2/3 posterior covers 0.2 eV.
- Limitations: barrier depends on substrate relaxation protocol; the comparison uses the fully relaxed protocol of the sources.

## B3. H₂ physisorption on coronene (CCSD(T)/CBS)

- Source: Kocman et al. 2015 (S8). Identifier: DOI 10.1039/C4CP04354E.
- Quantity: −4.7 kJ mol⁻¹ per H₂ (CCSD(T)/CBS, cross-checked against DMC in the source).
- Prediction: Feynman–Hibbs-uncorrected classical minimum of the pairwise H₂–C potential summed over the coronene geometry (from the ESI); and Tier-1 vdW-DF cluster calculation.
- Metric: absolute difference.
- Independence: anchor for the physisorption potential parameters; also leave-one-out.
- Acceptance: pairwise potential within 1.0 kJ mol⁻¹; Tier-1 vdW-DF within 1.5 kJ mol⁻¹.
- Limitations: cluster versus periodic; orientation averaging of H₂.

## B4. Isosteric heats of H₂ on carbons

- Sources: Bhatia and Myers 2006 (AX-21 average 5.8 kJ mol⁻¹); Knight et al. 2020 (activated carbon at room temperature, method with reduced uncertainty); corpus-derived q_st from S1 multi-temperature sets.
- Quantity: q_st at low loading 4–7 kJ mol⁻¹ for graphitic and activated carbons.
- Prediction: physisorption landscape (basal and pore sites) and corpus inference.
- Metric: interval containment and distribution overlap.
- Independence: external validation (not used for calibration).
- Acceptance: computed basal/pore q_st within [4, 8] kJ mol⁻¹; corpus median for activated carbons within [4, 8] kJ mol⁻¹.
- Limitations: q_st depends on loading and on the excess–absolute conversion.

## B5. Thermal desorption of H and D from graphite(0001)

- Sources: Zecho et al. 2002 (peaks 445 K and 560 K for H₂; 490 K and 580 K for D₂); Hornekær et al. 2006 (dimer mechanism; recombination barriers).
- Conditions: ramp rate as reported by the source (recorded in `benchmarks.yaml`); low coverage.
- Prediction: C4 forward multi-species model with landscape barriers for ortho/para dimers.
- Metric: peak temperature differences; presence of two peaks.
- Independence: external validation.
- Acceptance: both peak positions within ±50 K at 90% posterior; two-peak structure reproduced; H/D ordering reproduced.
- Limitations: coverage-dependent spectra; the comparison uses the low-coverage limit.

## B6. Hydrogenated single-walled nanotubes versus diameter

- Sources: Nikitin et al. 2005 (65 ± 15 at% hydrogenation, C–H bonds broken by 600 °C); Nikitin et al. 2008 (maximal hydrogenation depends on diameter; near-complete hydrogenation at d ≈ 2 nm stable at room temperature).
- Quantity: ordering of hydrogenation stability with diameter; complete release by ≈ 873 K.
- Prediction: curvature series (E2) and forward release (E6).
- Metric: ordering agreement; predicted T at which f > 0.99 within the stated window.
- Independence: external validation.
- Acceptance: predicted complete release (f > 0.99 under the source's thermal protocol) below 900 K for d ≈ 1–2 nm; diameter ordering of C–H stability reproduced with posterior probability > 0.9.
- Limitations: diameter distributions, bundling, and XPS quantification uncertainty in the sources.

## B7. Surface-area scaling of 77 K excess uptake in carbons

- Source: corpus S1; the empirical Chahine-type rule (about 1 wt% per 500 m² g⁻¹ at 77 K and moderate pressure) as discussed in the literature (e.g., Chahine and Bose 1994 and subsequent compilations).
- Quantity: slope of maximum excess uptake at 77 K versus BET area for materials with reported areas.
- Prediction: corpus regression.
- Metric: slope with bootstrap CI.
- Independence: external consistency check.
- Acceptance: slope CI overlaps 1.5–2.5 wt% per 1000 m² g⁻¹.
- Limitations: BET area reported for a subset only.

## B8. Bhatia–Myers optimum (analytic)

- Source: Bhatia and Myers 2006.
- Quantity: ΔH* = −15.1 kJ mol⁻¹ at 298 K between 30 and 1.5 bar with ΔS = −8R; carbons' optimum temperature ≈ 115 K for ΔH = −5.8 kJ mol⁻¹.
- Prediction: C1 closed form (DERIVATIONS 2.7) and numerical optimization.
- Acceptance: within 0.1 kJ mol⁻¹ and 5 K respectively.
- Independence: theory check.

## B9. Electron knock-on cross sections for graphene

- Source: Meyer et al. 2012 (S5).
- Quantity: measured displacement cross section versus beam energy for ¹²C (90, 95, 100 keV) and ¹³C (95, 100 keV) single-layer graphene; the source states that the cross section rises from about 10⁻⁴ barn at 80 keV to about 0.2 barn at 100 keV, that the static-lattice onset is at 108 kV, and that the thermal-vibration model uses a first-principles threshold of 22 eV without adjustment with the Debye model (out-of-plane Debye temperature 1287 K for ¹²C, scaled by √(12/13) for ¹³C), the 0 K and 300 K curves being almost identical. Individual values are available only graphically (Figs. 2 and 4).
- Prediction: McKinley–Feshbach with Debye-model velocity smearing (DERIVATIONS 6), E_d = 22 eV, no fitted parameters.
- Metric: ratio of predicted to published values at the stated endpoints; residuals against digitized points with digitization uncertainty; static onset energy.
- Independence: external validation (no parameter is fitted to this benchmark).
- Acceptance: static onset within 2 keV of 108 keV; predicted ¹²C cross section at 100 keV within a factor of 2 of 0.2 barn; predicted 80 keV value within one order of magnitude of 10⁻⁴ barn; ¹³C below ¹²C at 95 and 100 keV; ≥ 80% of digitized points within their digitization error bars. A failure of the 80 keV criterion is treated as evidence of a model-form difference in the smearing kinematics and is reported, with the alternative kinematic formulations compared.
- Limitations: digitization of figure data; the source's kinematic treatment is reproduced from its description rather than from its equations.

## B10. Raman defect relations

- Sources: Lucchese et al. 2010; Cançado et al. 2011; Eckmann et al. 2012 (S6).
- Quantity: I_D/I_G versus L_D at 2.41 eV; I_D/I_D′ ≈ 13 (sp³), ≈ 7 (vacancy-like), ≈ 3.5 (boundaries).
- Prediction: DERIVATIONS 8 with source parameters.
- Acceptance: reproduction of the published curve within 5% over L_D ∈ [1, 20] nm; ratios as tabulated.
- Independence: parameter transcription check (not a prediction).

## B11. Defect energetics in graphene

- Sources: Krasheninnikov and Nordlund 2010 and sources therein (monovacancy formation ≈ 7.4–7.7 eV; migration ≈ 1.3 eV; divacancy formation lower per vacancy; Stone–Wales formation ≈ 5 eV, barrier ≈ 9 eV); Lehtinen et al. 2004 for vacancy magnetism.
- Prediction: Tier-1 and Tier-2 formation and migration energies.
- Acceptance: Tier-1 within 0.3 eV of the published ranges; Tier-2 within 1.0 eV with the deviation recorded as a calibration input.
- Independence: engine sanity; may inform calibration.

## B12. Graphane formation and reversibility

- Source: Elias et al. 2009 (graphane formed by exposure to atomic hydrogen; original graphene properties restored by annealing at 450 °C for 24 h).
- Quantity: near-complete dehydrogenation at ≈ 723 K over ≈ 10⁵ s.
- Prediction: C4 release at high coverage with cluster/graphane-edge pathways.
- Acceptance: predicted release fraction at 723 K over 8.6 × 10⁴ s > 0.9 at 90% posterior.
- Independence: external validation.
- Limitations: substrate effects (graphene on SiO₂); double-sided versus single-sided hydrogenation.
