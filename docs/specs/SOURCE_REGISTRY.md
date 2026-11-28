# Public Source Registry

Every external input of the program is listed here. "Verification" records how the identifier and availability were confirmed (publisher landing page, arXiv record, database record, or a citing reference list in an open-access source). Values transcribed into `data/registry/published_values.yaml` carry a record ID referencing this registry.

## Data sources

### S1. NIST/ARPA-E Database of Novel and Emerging Adsorbent Materials (ISODB)

- Identity: NIST Data Resources for Adsorption Science and Technology, isotherm database.
- Persistent identifier: DOI 10.18434/T43882 (NIST data record FF429BC1787D8B3EE0431A570681E858236).
- Availability: public web application (https://adsorption.nist.gov/isodb) and official GitHub mirror https://github.com/NIST-ISODB/isodb-library (JSON, one file per isotherm, organized by DOI-derived directory).
- Acquisition: `git clone --depth 1` of the mirror, pinned by commit hash in `data/sources.lock`; parser `sitewindow.experimental.isodb`.
- Authentication: none.
- License: NIST open-data terms (https://www.nist.gov/open/license); derived tables redistributed with attribution and per-record DOIs; the raw clone is not redistributed.
- Quantities: single- and multi-component isotherms; fields include `adsorbent.name`, `adsorbates[].name`, `temperature` (K), `pressureUnits`, `adsorptionUnits`, `isotherm_data[]` (pressure, adsorption), `digitizer`, `tabular_data` (flag), `DOI`, `category`.
- Units: mixed (bar, kPa, MPa, atm; wt%, mmol g⁻¹, mol kg⁻¹, cm³(STP) g⁻¹, g L⁻¹, mg g⁻¹, …); normalized by `sitewindow.units`.
- Material characteristics: adsorbent names; BET area and pore volume only when present in the record; cross-reference to the NIST adsorbent registry by name.
- Suitability: at the pinned snapshot inspected during program design, 54 419 isotherm files; 5 461 single-component H₂ isotherms; 431 on carbon-like adsorbents by name rules (94 distinct names), including 152 at 77 K and 117 at 298 K and multi-temperature sets at 87, 100, 133, 177, 200, 293, 303, 318 K. Sufficient for corpus inference (C6).
- Limitations: adsorbent naming is free-text; excess/absolute status is not always explicit; digitized isotherms carry digitization error; surface areas are often absent.
- Overlap: none with other sources in the program.
- Verification: NIST program page and data.gov record (DOI, license); GitHub repository cloned and parsed during program design.

### S2. NIST Registry of Adsorbent Materials

- Identity: index of adsorbent materials with unique identifiers.
- Availability: https://adsorption.nist.gov (registry endpoint of the same portal). Used only for material-identity cross-reference; not essential.
- Verification: NIST program page describes the registry; programmatic endpoint usage is verified at execution and recorded in the source lock; absence of the endpoint degrades gracefully to name-based rules. At execution (`/isodb/api/material/<hashkey>.json`) a material record holds `hashkey`, `name`, `formula`, `synonyms` and `External_Resources` only: no surface area or pore volume, so the BET areas that benchmark B7 needs and the pore volumes of the excess-to-absolute conversion must be read in the source publications and are kept in `data/registry/material_metadata.yaml` with the place where each was read. The publishers' pages (Elsevier, ACS, Wiley, Springer, RSC) refuse automated clients, so the accessible text of most corpus publications is the abstract in the OpenAlex index; at this version the abstracts read give no sample-specific BET area, and the registry records instead which database records are model calculations rather than measurements (Fomkin et al. 2008; the model pore geometries of Thornton et al. 2012) and which generic "Carbon" records belong to papers with several carbons (Stadie et al. 2012; Kadono et al. 2003; Juan-Juan et al. 2010; Marco-Lozar et al. 2012).

### S3. arXiv open-access full texts (literature anchors and benchmark data)

- Availability: https://arxiv.org (public; PDF and abstracts).
- Acquisition: programmatic fetch of listed identifiers; text extraction; values recorded with page/table reference.
- License: arXiv non-exclusive license; facts and data values are used with citation.

### S4. Published-value registry (first-principles and experimental values used as anchors or benchmarks)

Stored in `data/registry/published_values.yaml`. Each record: `id`, `doi` or `arxiv`, `system`, `site`, `quantity`, `value`, `unit`, `uncertainty`, `method` (functional, dispersion, code, basis/cutoff, supercell, k-mesh, spin), `reference_state`, `extraction` (`tabular_text`, `digitized_figure`, `abstract_statement`), `source_id`.

First-principles chemisorption values enter only after the number was read in an open-access copy of the full text (arXiv version, a public institutional repository copy, or an abstract that states the number); the `where_read` field names the table, figure or paragraph. Values that a paper only cites from other work are not registered. Papers whose only open copies were unreachable (Hornekær et al. 2006; Šljivančanin et al. 2009, whose values the same authors reproduce in arXiv:1102.4984; Ivanovskaya et al. 2010; Ferro et al. 2002, 2003, 2008; Ma et al. 2009; Nikitin et al. 2008; Arellano et al. 2002; Lehtinen et al. 2004; Dobrota et al. 2016 on hydrogen beside epoxide and hydroxyl groups) are not represented, and no open-access first-principles value was found for hydrogen at the 555-777 divacancy, beside carbonyl or carboxyl groups, or on C60.

### S5. Electron knock-on cross sections for graphene (benchmark B9 data)

- Source: Meyer et al., Phys. Rev. Lett. 108, 196102 (2012); open-access preprint arXiv:1203.2372 (with supplementary information).
- Acquisition: measured cross sections versus beam energy extracted from the preprint (tabular where available; otherwise digitized with uncertainty).
- Verification: arXiv record and journal DOI 10.1103/PhysRevLett.108.196102 confirmed.

### S6. Raman defect-quantification parameters (benchmark B10 data)

- Sources: Cançado et al., Nano Lett. 11, 3190 (2011), DOI 10.1021/nl201432g, preprint arXiv:1105.0175 (equation parameters and corrected numerical factors in v2); Lucchese et al., Carbon 48, 1592 (2010), DOI 10.1016/j.carbon.2009.12.057; Eckmann et al., Nano Lett. 12, 3925 (2012), DOI 10.1021/nl300901a, preprint arXiv:1207.2058.
- Acquisition: parameters from the open-access preprints.
- Verification: arXiv records and DOIs confirmed.

### S7. DMC/DFT dataset for H on graphene

- Source: Dumi et al., J. Chem. Phys. 156, 144702 (2022), DOI 10.1063/5.0085982; preprint arXiv:2202.00194; associated data record DOI 10.18126/s1wc-tyac (Materials Data Facility).
- Acquisition: values from the open-access preprint (Table I); structures and inputs from the data record when needed.
- Verification: OSTI record, journal DOI, and data DOI confirmed.

### S8. CCSD(T)/CBS reference for H₂–coronene

- Source: Kocman, Jurečka, Dubecký, Otyepka, Cho, Kim, Phys. Chem. Chem. Phys. 17, 6423 (2015), DOI 10.1039/C4CP04354E; electronic supplementary information with geometries (https://www.rsc.org/suppdata/cp/c4/c4cp04354e/c4cp04354e1.pdf).
- Acquisition: reference binding energy from the abstract (−4.7 kJ mol⁻¹ per H₂ for coronene; −14.3 kJ mol⁻¹ for the B/Li-modified coronene) and geometries from the ESI.
- Verification: publisher landing page and ESI URL confirmed.

### S9. U.S. DOE onboard hydrogen storage technical targets (light-duty vehicles)

- Source: https://www.energy.gov/eere/fuelcells/doe-technical-targets-onboard-hydrogen-storage-light-duty-vehicles (table updated May 2017; targets include minimum delivery pressure 5 bar, maximum delivery pressure 12 bar, delivery temperature −40/85 °C, fill time 3–5 min, operating ambient −40/60 °C).
- Acquisition: transcribed with retrieval date.
- Verification: page content confirmed during program design.

## Literature references used in the manual (verified)

| Key | Reference | Identifier | Verification |
|---|---|---|---|
| BhatiaMyers2006 | S. K. Bhatia, A. L. Myers, "Optimum conditions for adsorptive storage," Langmuir 22, 1688–1700 (2006) | journal volume/pages confirmed in multiple citing sources and the authors' institutional copy | citing reference lists; author page |
| Dumi2022 | A. Dumi, S. Upadhyay, L. Bernasconi, H. Shin, A. Benali, K. D. Jordan, J. Chem. Phys. 156, 144702 (2022) | DOI 10.1063/5.0085982; arXiv:2202.00194 | OSTI; arXiv |
| Casolo2009 | S. Casolo, O. M. Løvvik, R. Martinazzo, G. F. Tantardini, "Understanding adsorption of hydrogen atoms on graphene," J. Chem. Phys. 130, 054704 (2009) | DOI 10.1063/1.3072333; arXiv:0808.1312 | publisher citation line in open copy; arXiv |
| ShaJackson2002 | X. Sha, B. Jackson, Surf. Sci. 496, 318 (2002) | — | citing reference lists (multiple) |
| Zecho2002 | T. Zecho, A. Güttler, X. Sha, B. Jackson, J. Küppers, J. Chem. Phys. 117, 8486 (2002) | — | citing reference lists (multiple) |
| Hornekaer2006 | L. Hornekær et al., Phys. Rev. Lett. 96, 156104 (2006) | DOI 10.1103/PhysRevLett.96.156104 | APS landing page |
| Hornekaer2006b | L. Hornekær et al., Phys. Rev. Lett. 97, 186102 (2006) | — | citing reference lists |
| Boukhvalov2008 | D. W. Boukhvalov, M. I. Katsnelson, "Chemical functionalization of graphene with defects," Nano Lett. 8, 4373 (2008) | arXiv:0807.3855 | arXiv; citing reference list |
| Nikitin2005 | A. Nikitin et al., "Hydrogenation of single-walled carbon nanotubes," Phys. Rev. Lett. 95, 225507 (2005) | arXiv:cond-mat/0510399 | arXiv; citing reference list |
| Nikitin2008 | A. Nikitin, X. Li, Z. Zhang, H. Ogasawara, H. Dai, A. Nilsson, "Hydrogen storage in carbon nanotubes through the formation of stable C–H bonds," Nano Lett. 8, 162–167 (2008) | DOI 10.1021/nl072325k | citation index record |
| Kocman2015 | M. Kocman et al., Phys. Chem. Chem. Phys. 17, 6423 (2015) | DOI 10.1039/C4CP04354E | publisher |
| Meyer2012 | J. C. Meyer et al., Phys. Rev. Lett. 108, 196102 (2012) | DOI 10.1103/PhysRevLett.108.196102; arXiv:1203.2372 | APS; arXiv |
| Krasheninnikov2010 | A. V. Krasheninnikov, K. Nordlund, J. Appl. Phys. 107, 071301 (2010) | DOI 10.1063/1.3318261 | journal page header in open copy |
| Jones2010 | J. D. Jones, K. K. Mahajan, W. H. Williams, P. A. Ecton, Y. Mo, J. M. Perez, "Formation of graphane and partially hydrogenated graphene by electron irradiation of adsorbates on graphene," Carbon 48, 2335–2340 (2010) | DOI 10.1016/j.carbon.2010.03.010 | citation index record |
| Lucchese2010 | M. M. Lucchese et al., Carbon 48, 1592–1597 (2010) | DOI 10.1016/j.carbon.2009.12.057 | citing reference list with DOI |
| Cancado2011 | L. G. Cançado et al., Nano Lett. 11, 3190–3196 (2011) | DOI 10.1021/nl201432g; arXiv:1105.0175 | arXiv; index record |
| Eckmann2012 | A. Eckmann et al., Nano Lett. 12, 3925–3930 (2012) | DOI 10.1021/nl300901a; arXiv:1207.2058 | index record; arXiv |
| Elias2009 | D. C. Elias et al., Science 323, 610 (2009) | — | citing reference lists (multiple) |
| Tozzini2013 | V. Tozzini, V. Pellegrini, "Prospects for hydrogen storage in graphene," Phys. Chem. Chem. Phys. 15, 80 (2013) | DOI 10.1039/C2CP42538F; arXiv:1207.5703 | publisher; arXiv |
| Tozzini2010 | V. Tozzini, V. Pellegrini, Phys. Rev. B 81, 113404 (2010) | — | reference list of Tozzini2013 |
| Chahine1994 | R. Chahine, T. K. Bose, Int. J. Hydrogen Energy 19, 161 (1994) | — | reference list of BhatiaMyers2006 |
| Dillon1997 | A. C. Dillon et al., Nature 386, 377 (1997) | — | citing reference lists |
| Knight2020 | E. W. Knight et al., "Determination of the enthalpy of adsorption of hydrogen in activated carbon at room temperature," Int. J. Hydrogen Energy (2020) | Elsevier PII S0360319920313495 | journal front matter in open copy |
| Langlois2025 | P. Langlois, C. Chilev, F. Lamari, "Reassessed ability of carbon-based physisorbing materials to keep pace with evolving practical targets for hydrogen storage," preprint (Dec 2025) with peer-reviewed version in C (MDPI) | DOI 10.20944/preprints202512.1196.v1; 10.3390/c12010009 | preprint server; index record |
| Sunnardianto2021 | G. K. Sunnardianto et al., Int. J. Hydrogen Energy 46, 5485–5494 (2021) | — | journal header in open copy |
| NIST-ISODB | NIST/ARPA-E Database of Novel and Emerging Adsorbent Materials | DOI 10.18434/T43882 | data.gov record |
| Lehtinen2004 | P. O. Lehtinen et al., Phys. Rev. Lett. 93, 187202 (2004) | — | reference list of Krasheninnikov & Banhart (Nat. Mater. 2007) |
| KrasheninnikovBanhart2007 | A. V. Krasheninnikov, F. Banhart, Nat. Mater. 6, 723 (2007) | DOI 10.1038/nmat1996 | doi.org landing |
| Gavardi2009 | E. Gavardi, H. M. Cuppen, L. Hornekær, "A kinetic Monte Carlo study of desorption of H₂ from graphite (0001)," Chem. Phys. Lett. (2009) | arXiv:0907.3361 | arXiv |
| BoukhvalovKatsnelsonLichtenstein2008 | D. W. Boukhvalov, M. I. Katsnelson, A. I. Lichtenstein, "Hydrogen on graphene: Electronic structure, total energy, structural distortions and magnetism from first-principles calculations," Phys. Rev. B 77, 035427 (2008) | DOI 10.1103/PhysRevB.77.035427; arXiv:0710.1971 | arXiv full text read |
| Koskinen2008 | P. Koskinen, S. Malola, H. Häkkinen, "Self-passivating edge reconstructions of graphene," Phys. Rev. Lett. 101, 115502 (2008) | DOI 10.1103/PhysRevLett.101.115502; arXiv:0802.2623 | arXiv full text read |
| Wassmann2008 | T. Wassmann, A. P. Seitsonen, A. M. Saitta, M. Lazzeri, F. Mauri, "Structure, stability, edge states, and aromaticity of graphene ribbons," Phys. Rev. Lett. 101, 096402 (2008) | DOI 10.1103/PhysRevLett.101.096402; arXiv:0808.2283 | arXiv full text read (edge energies per length; not a registry quantity) |
| Sljivancanin2011 | Ž. Šljivančanin, M. Andersen, L. Hornekær, B. Hammer, "Structure and stability of small H clusters on graphene," Phys. Rev. B 83, 205426 (2011) | arXiv:1102.4984 | arXiv full text read |
| Casartelli2014 | M. Casartelli, S. Casolo, G. F. Tantardini, R. Martinazzo, "Structure and stability of hydrogenated carbon atom vacancies in graphene," Carbon 77, 165 (2014) | arXiv:1312.7143 | arXiv full text read |
| Thomas2022 | J. A. Thomas, A. M. Asiri, N. D. Drummond, "Point defect formation energies in graphene from diffusion quantum Monte Carlo and density functional theory," Phys. Rev. B 105, 214106 (2022) | arXiv:2112.11539 | arXiv full text read |
| Duplock2004 | E. J. Duplock, M. Scheffler, P. J. D. Lindan, "Hallmark of perfect graphene," Phys. Rev. Lett. 92, 225502 (2004) | DOI 10.1103/PhysRevLett.92.225502 | public Max Planck Society repository copy read |
| Yildirim2001 | T. Yildirim, O. Gülseren, S. Ciraci, "Exohydrogenated single-wall carbon nanotubes," Phys. Rev. B 64, 075404 (2001) | DOI 10.1103/PhysRevB.64.075404; arXiv:cond-mat/0105244 | arXiv full text read |
| Gulseren2001 | O. Gülseren, T. Yildirim, S. Ciraci, "Tunable adsorption on carbon nanotubes," Phys. Rev. Lett. 87, 116802 (2001) | DOI 10.1103/PhysRevLett.87.116802; arXiv:cond-mat/0107485 | arXiv full text read |
| Gulseren2002 | O. Gülseren, T. Yildirim, S. Ciraci, "Effects of hydrogen adsorption on single-wall carbon nanotubes: Metallic hydrogen decoration," Phys. Rev. B 66, 121401 (2002) | DOI 10.1103/PhysRevB.66.121401; arXiv:cond-mat/0208510 | arXiv full text read |
| LiReichRobertson2005 | L. Li, S. Reich, J. Robertson, "Defect energies of graphite: Density-functional calculations," Phys. Rev. B 72, 184109 (2005) | DOI 10.1103/PhysRevB.72.184109 | abstract text from the OpenAlex metadata record; full text not open |

References without a DOI in this table are cited by journal, volume, and page as confirmed in open-access reference lists; the manuscript bibliography generator requires, for each entry, either a confirmed DOI or two independent confirmations of the journal citation, recorded in `data/registry/references.yaml`.
