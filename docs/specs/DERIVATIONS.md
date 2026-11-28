# Mathematical Derivations

This document derives every equation used by `sitewindow`. Equation numbers are referenced by docstrings and by tests (`tests/test_derivations.py`). Symbols: R gas constant, k_B Boltzmann constant, T temperature (K), P pressure (bar, absolute), P° = 1 bar, ε site binding energy (negative for exothermic binding), θ fractional occupancy, Δθ swing, W = (T_c, P_c, T_r, P_r, τ_c, τ_r) an operating window.

## 1. Definitions

**(1.1) Site.** A binding site is characterized by channel (molecular or atomic), binding energy ε, charging barrier E_a†, release barrier E_d†, multiplicity g (equivalent sites per defect), maximum occupancy n_max (entities per site), and the lattice mass change Δm per defect that creates it.

**(1.2) Landscape.** ρ(ε, E_a†, E_d†) is the density of sites per carbon atom; for discrete enumerations, ρ is a weighted sum of Dirac measures with weights x_c g_s.

**(1.3) Deliverable capacity.** For a cycle between charged state (T_c, P_c) and released state (T_r, P_r), the deliverable capacity is the mass of hydrogen released per mass of sorbent in the charged state, restricted to kinetically admissible sites.

## 2. Molecular channel

**(2.1) Langmuir occupancy.** For independent sites with adsorption free energy ΔG(ε, T) = ΔH(ε, T) − TΔS(T):

θ(ε, T, P) = K P / (1 + K P), K = (1/P°) exp[−ΔG/(R T)].

With fugacity f replacing P above 50 bar (Section 2.6).

**(2.2) Enthalpy.** ΔH(ε, T) = ε + Δ_th(T), where Δ_th collects thermal corrections (difference of heat capacities between adsorbed and gas states integrated from 0 K). Baseline: Δ_th = 0 with ε understood as the isosteric heat at the stated temperature; alternative: Δ_th = −R T/2 for a localized adsorbate losing translational degrees of freedom. The choice is a model-set member in sensitivity analysis.

**(2.3) Entropy models.**
(a) Constant: ΔS = −8R (Bhatia–Myers baseline).
(b) Localized harmonic oscillator: S_ads = R Σ_{i=1}^{3} [ x_i/(e^{x_i} − 1) − ln(1 − e^{−x_i}) ], x_i = hν_i/(k_B T), with ν_i the frustrated translational frequencies of H₂ in the site (from the curvature of the site potential); S_gas = S_trans(T, P) + S_rot(T) from Sackur–Tetrode and the rigid-rotor partition function of H₂ (ortho/para equilibrium composition at T); ΔS = S_ads − S_gas.
(c) Mobile two-dimensional gas: S_ads = R [ ln(a_m/Λ²) + 2 ] + S_rot, with a_m the area per molecule and Λ = h/√(2π m k_B T); ΔS = S_ads − S_gas.

**(2.4) Swing.** Δθ(ε; W) = θ(ε, T_c, P_c) − θ(ε, T_r, P_r).

**(2.5) Isothermal optimum.** For T_c = T_r = T, maximizing Δθ over K:

d/dK [K P_c/(1 + K P_c) − K P_r/(1 + K P_r)] = P_c/(1 + K P_c)² − P_r/(1 + K P_r)² = 0 ⇒ K* = 1/√(P_c P_r).

Then K* P_c = √r and K* P_r = 1/√r with r = P_c/P_r, so

**(2.6)** Δθ_max = √r/(1 + √r) − 1/(1 + √r) = (√r − 1)/(√r + 1).

**(2.7) Optimum enthalpy.** From K* = (1/P°) exp[−ΔG*/(R T)]: ΔG* = R T ln(√(P_c P_r)/P°), so ΔH* = R T ln(√(P_c P_r)/P°) + T ΔS. For T = 298.15 K, P_c = 30 bar, P_r = 1.5 bar, ΔS = −8R: ΔH* = R T ln(6.708) − 8 R T = (4.716 − 19.83) kJ mol⁻¹ = −15.1 kJ mol⁻¹, which is the Bhatia–Myers result (benchmark B8).

**(2.8) General window.** For T_c ≠ T_r the optimum solves ∂Δθ/∂ε = 0 numerically; Δθ is unimodal in ε for all entropy models in the set (verified numerically on the grid used; a non-unimodal case would be flagged).

**(2.9) Heterogeneity theorem.** Let a landscape consist of a homogeneous population at ε₁ with site fraction (1 − φ) and a minority at ε₂ with fraction φ, for a fixed total site count. Then D ∝ (1 − φ)Δθ(ε₁) + φ Δθ(ε₂), and ∂D/∂φ ∝ Δθ(ε₂) − Δθ(ε₁). Hence adding stronger sites increases D if and only if Δθ(ε₂) > Δθ(ε₁). For ε₁ weaker than ε* (|ε₁| < |ε*|), because Δθ is unimodal there exists a unique threshold ε_c(ε₁; W) on the strong side of ε* with Δθ(ε_c) = Δθ(ε₁); strong sites help iff |ε₂| < |ε_c|. With mass accounting (Section 5), the condition becomes Δθ(ε₂) > Δθ(ε₁)(1 + δ_m) where δ_m is the fractional mass penalty per added site; the program computes ε_c with and without δ_m.

**(2.10) Feynman–Hibbs correction.** For an adsorbate of reduced mass μ in potential U(r), the effective classical potential to second order is U_FH(r) = U(r) + (ħ²/(24 μ k_B T)) ∇²U(r). For a pairwise Lennard-Jones H₂–C potential U = 4ε_LJ[(σ/r)¹² − (σ/r)⁶], ∇²U = 4ε_LJ[132 σ¹²/r¹⁴ − 30 σ⁶/r⁸]. Site energies are the minimum of the summed U_FH over carbon positions; the temperature dependence of ε is retained and evaluated at the window temperatures.

## 3. Atomic (dissociative) channel

**(3.1) Equilibrium.** For H₂(g) + 2 S ⇌ 2 H–S with per-H free energy change Δg = ε − T Δs (ε per H relative to ½ H₂):

θ²/(1 − θ)² = K P, K = (1/P°) exp[−2Δg/(R T)] ⇒ θ = √(K P)/(1 + √(K P)).

**(3.2)** For ε ≤ −0.3 eV per H and T ≤ 600 K, √(K P) ≫ 1 at all window pressures, so θ ≈ 1 − 1/√(K P); the thermodynamic swing is negligible and the channel is kinetically controlled.

**(3.3) Charging by atomic hydrogen.** With incident atomic flux Φ_H (atoms per site per second) and sticking coefficient s(T) = s₀ exp(−E_a†/(k_B T)), dθ/dt = s Φ_H (1 − θ) ⇒ θ(τ) = 1 − exp(−s Φ_H τ).

## 4. Release kinetics and admissibility

**(4.1) Polanyi–Wigner.** dθ/dt = −ν_n θⁿ exp[−E_d†/(k_B T)].

**(4.2) Isothermal solutions.** n = 1: θ(τ) = θ₀ e^{−kτ}, f = 1 − e^{−kτ}. n = 2: θ(τ) = θ₀/(1 + k θ₀ τ), f = k θ₀ τ/(1 + k θ₀ τ), with k = ν_n exp[−E_d†/(k_B T)].

**(4.3) Linear ramp.** T = T₀ + β t; the rate per unit temperature is r(T) = (ν_n/β) θⁿ exp[−E_d†/(k_B T)].

**(4.4) Redhead relations.** Peak temperature T_p satisfies, for n = 1: E_d†/(k_B T_p²) = (ν₁/β) exp[−E_d†/(k_B T_p)]; for n = 2: E_d†/(k_B T_p²) = (2 ν₂ θ_p/β) exp[−E_d†/(k_B T_p)] with θ_p the coverage at the peak. These are used as analytic checks of the ODE integrator.

**(4.5) Admissibility.** A_r = 1[f(T_r, τ_r) ≥ η_r], A_c = 1[θ(τ_c) ≥ η_c]. Smooth alternative: A = 1/(1 + exp[−κ ln(k τ/k₀τ₀)]) with κ = 4.

**(4.6) T₉₀.** The temperature at which f(T, τ_r) = 0.9 under isothermal release; for n = 1, T₉₀ = E_d†/(k_B ln(ν τ_r/ln 10)).

**(4.7) Multi-species basal plane.** Populations m (monomer), o (ortho dimer), p (para dimer), t (meta/transient): dm/dt = −k_mo m² ..., do/dt = ..., with transformation o → t → p and recombination from o directly and from p through t, each with its own barrier from the landscape; integrated with a stiff solver. The pathway set is fixed before benchmarking.

## 5. Mass accounting and the deliverable functional

**(5.1)** For class c with density x_c (defects per C) and sites s with multiplicity g_s and occupancy n_s, the hydrogen entities per carbon in state X are h_c(X) = x_c Σ_s g_s n_s θ_s(X).

**(5.2)** Sorbent mass per carbon in the charged state: m_tot = m_C + Σ_c x_c Δm_c + Σ_c h_c(charged) m_entity.

**(5.3)** Deliverable capacity (mass fraction): D = [Σ_c x_c Σ_s g_s n_s Δθ_s A_s m_entity] / m_tot, with m_entity = m_H₂ for the molecular channel and m_H for the atomic channel (two atomic sites release one H₂).

**(5.4) Bound.** For a single channel and window let d(ε) := N_max m_entity Δθ(ε; W) A(ε; W) / (m_C + N_max m_entity θ_c(ε; W) + Δm_min) be the deliverable mass fraction of a homogeneous landscape that places the channel's maximum site density N_max (sites per carbon) at ε, with Δm_min the smallest admissible lattice mass change (zero for basal sites), and let D_max(W) := sup_ε d(ε), attained at ε_D(W). Then D[ρ; W] ≤ D_max(W) for every landscape with Σ_c x_c Σ_s g_s n_s ≤ N_max. Proof: write D = Σ_s a_s / (m_C + Δm + Σ_s b_s) with a_s = w_s m_entity Δθ_s A_s and b_s = w_s m_entity θ_c,s for site weights w_s summing to at most N_max. For each site, d(ε_s) ≤ D_max gives N_max m_entity Δθ_s A_s ≤ D_max (m_C + Δm_min + N_max m_entity θ_c,s); multiplying by w_s/N_max and summing, Σ a_s ≤ D_max [(m_C + Δm_min) Σ w_s/N_max + Σ b_s] ≤ D_max (m_C + Δm + Σ b_s). For the atomic channel N_max = 1 (one H per carbon; graphane), giving the familiar 7.7 wt% ceiling times η. For the molecular channel N_max is set by the monolayer areal density on the maximum accessible area (both faces of graphene, 2630 m² g⁻¹) and is computed from the H₂–H₂ potential diameter. The earlier form of this bound, N_max m_entity η(W)/(m_C + N_max m_entity θ_c(ε*)) with η(W) := sup_ε Δθ A attained at the swing optimum ε*, evaluates d at ε* and therefore lies below sup_ε d by the small amount that the charged-hydrogen mass in the denominator shifts the optimum; a randomized test found a homogeneous landscape 0.002 eV weaker than ε* exceeding it by three parts in 10⁵ once the fugacity was evaluated with the reference equation of state, and the statement was corrected to the supremum of d. Because Δθ A ≤ η and θ_c ≥ 0, D_max ≤ N_max m_entity η/m_C still holds, and ε_D differs from ε* by about 10⁻³ eV.

**(5.5) Tightness.** Equality in (5.4) requires every site at ε_D and admissible (a homogeneous landscape at ε* comes within 3 × 10⁻⁴ of the bound for W1); the gap D_max − D[ρ] decomposes into an occupancy deficit (sites outside the window), an admissibility deficit (sites thermodynamically in the window but kinetically excluded), and a mass deficit (defect mass). The program reports the three deficits per class.

## 6. Electron knock-on displacement

**(6.1) McKinley–Feshbach cross section.** For incident electron kinetic energy E, β² = 1 − (m_e c²/(E + m_e c²))², maximum energy transfer T_max = 2E(E + 2 m_e c²)/(M c²) to a nucleus of mass M, charge Z, α = Z/137.04:

σ_d(E; E_d) = (π Z² e⁴ (1 − β²))/(m_e² c⁴ β⁴) · { (T_max/E_d − 1) − β² ln(T_max/E_d) + π α β [ 2(√(T_max/E_d) − 1) − ln(T_max/E_d) ] }, for T_max ≥ E_d; zero otherwise.

Numerically, π Z² e⁴/(m_e² c⁴) = π Z² r_e² with r_e = 2.8179 × 10⁻¹⁵ m, giving σ in barn when multiplied by 10²⁸.

**(6.2) Thermal smearing.** The target nucleus has an out-of-plane velocity v_z drawn from a Gaussian with variance ⟨v_z²⟩ set by the phonon spectrum (parameterized by an effective temperature or by the mean-square displacement). In the nucleus rest frame the electron energy and T_max are modified; the lab-frame energy transfer for a head-on collision becomes T′ = T_max + 2√(T_max E_v) cos χ + E_v, with E_v = ½ M v_z² and χ the angle between the recoil and the thermal velocity. The effective cross section is σ_eff(E) = ∫ σ_d(E; E_d − correction(v_z)) P(v_z) dv_z, implemented by adaptive quadrature. ⟨v_z²⟩ follows the Debye model, ⟨v_z²⟩ = (1/M) ∫₀^{ω_D} g(ω)(ħω/2) coth(ħω/2k_BT) dω with g(ω) = 3ω²/ω_D³, giving 3k_Bθ_D/(8M) at 0 K; θ_D = 1287 K (out-of-plane, ¹²C) and E_d = 22 eV are taken from the source without adjustment, so that benchmark B9 is a parameter-free validation of the smearing kinematics.

**(6.3) Angular threshold on nanotubes.** E_d(ϑ, d) = E_d,0 + ΔE(d) cos²ϑ-type dependence with ΔE(d) decreasing with diameter, following the published form; the parameters are configuration data with an interval.

**(6.4) Defect production rate.** dN_D/dt = σ_eff(E) · J/e · (N_C − N_D)/N_C per atom, with J the current density; fluence Φ = Jt/e.

## 7. Defect rate equations

**(7.1)** With x_V, x_I, x_DV, x_SW the per-carbon densities of monovacancies, interstitials (adatoms), divacancies, and Stone–Wales defects:

dx_V/dt = σ_V Φ̇ (1 − x_tot) − k_rec x_V x_I − 2 k_dv x_V²,
dx_I/dt = σ_V Φ̇ (1 − x_tot) − k_rec x_V x_I − k_loss x_I,
dx_DV/dt = k_dv x_V² + σ_DV Φ̇ x_V,
dx_SW/dt = σ_SW Φ̇ (1 − x_tot) − k_SW x_SW,

with k = ν exp(−E_m/(k_B T)) for each thermally activated process, E_m(V) ≈ 1.3 eV (migration) and other energies from the benchmark registry (B11), x_tot = x_V + 2x_DV + x_SW. Carbon conservation: carbons removed = x_V + 2 x_DV; conserved to solver tolerance (tested).

**(7.2) Amorphization bound.** With areal carbon density n_C = 38.2 nm⁻², the per-carbon defect density at mean inter-defect distance L_D is x = 1/(π L_D² n_C); the admissible set is x_tot ≤ x(L_D,min) with L_D,min ∈ [1, 2] nm.

## 8. Raman relations

**(8.1)** I_D/I_G(L_D) = C_A (r_A² − r_S²)/(r_A² − 2 r_S²) [ exp(−π r_S²/L_D²) − exp(−π (r_A² − r_S²)/L_D²) ] + C_S [ 1 − exp(−π r_S²/L_D²) ], with r_A, r_S, C_A(E_L), C_S taken from the cited sources and carried with uncertainty.

**(8.2)** I_D/I_D′ for a mixed population: first-order linear mixing by class fraction with class values from the cited source; limitation stated (ratios are not strictly additive).

## 9. Adsorption-energy-distribution inversion

**(9.1)** n_abs(T, P) = n_m Σ_j w_j θ(ε_j, T, P), Σ_j w_j = 1, on a grid ε_j.

**(9.2)** Excess–absolute: n_exc = n_abs − ρ_gas(T, P) V_pore.

**(9.3)** Likelihood: y_k ~ N(n_exc(T_k, P_k), s_k²), s_k² = (a_src n_exc)² + b_src² + (∂n/∂P · δP_dig)², with (a, b) per source type (tabular, digitized) estimated hierarchically and δP_dig the pressure digitization uncertainty.

**(9.4)** Prior: log w_j = u_j − logsumexp(u), u_{j+1} − 2u_j + u_{j−1} ~ N(0, λ²), λ ~ HalfNormal(λ₀); n_m ~ LogNormal; V_pore ~ LogNormal centered on reported pore volume (if any) else broad.

**(9.5)** Window occupancy: Ω(W) = Σ_j w_j 1[ε_j ∈ window(W)].

## 10. Hierarchical meta-analysis of published values

**(10.1)** y_i = μ_c(i) + β_func(i) + β_disp(i) + γ/N_cell,i + u_study(i) + e_i; μ_c ~ N(μ₀, τ_c²); β ~ N(0, σ_β²); u ~ N(0, τ_u²); e_i ~ N(0, s_i²); hyperpriors HalfNormal on τ, σ; identification by sum-to-zero constraints on β within each factor.

## 11. Multi-fidelity calibration

**(11.1)** E_T1 = a_c + b_c E_t + δ_c(x) + e, δ_c ~ GP(0, k(x, x′)), k Matérn-5/2 with ARD length scales; hyperparameters by marginal likelihood with weak priors; predictive distribution for a new configuration by the standard GP posterior; the Tier-0 offset o₀ ~ N(ô₀, s₀²) estimated from the available high-level anchors is added for the chemisorption and physisorption channels separately.

## 12. Fugacity

**(12.1)** Above 50 bar the Langmuir pressure is replaced by fugacity f = φ P with φ from the hydrogen equation of state (CoolProp, Leachman et al. formulation) or, as fallback, from the second virial coefficient: ln φ = B(T) P/(R T) with B(T) from the tabulated NIST correlation; agreement between the two to 2% at 100 bar and 298 K is tested.
