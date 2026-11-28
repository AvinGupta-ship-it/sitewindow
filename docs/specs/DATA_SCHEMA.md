# Data Schema

All models are implemented as pydantic v2 classes in `sitewindow.schema` and exported as JSON Schema (`sitewindow schema export`). Canonical units are fixed per field; conversions happen only in `sitewindow.units`. Every record carries a `Provenance` block.

## Provenance

| Field | Type | Description |
|---|---|---|
| `kind` | enum | `experimental_observation`, `published_computation`, `computed_here`, `model_prediction`, `derived` |
| `source_id` | str | S# identifier or configuration hash |
| `identifier` | str | DOI, arXiv ID, ISODB filename, or cache key |
| `extraction` | enum or null | `tabular_text`, `digitized_figure`, `abstract_statement`, `database_record`, `computed` |
| `settings_hash` | str or null | SHA-256 of calculator/model settings |
| `package_version` | str | `sitewindow` version |
| `seed` | int or null | RNG seed |
| `created` | ISO 8601 | creation timestamp (provenance only; not used for organization) |

## Uncertainty

`UncertainScalar`: `value` (float), `sigma` (float ≥ 0), `distribution` (enum: `normal`, `lognormal`, `uniform`, `posterior_samples`), `samples_ref` (path to NetCDF variable when `posterior_samples`), `components` (optional dict: `numerical`, `digitization`, `method`, `measurement`).

## SiteClass

| Field | Type | Units | Description |
|---|---|---|---|
| `name` | str | — | canonical class name (`basal`, `monovacancy`, `divacancy_585`, `divacancy_555777`, `stone_wales`, `edge_zigzag`, `edge_armchair`, `cnt_armchair`, `cnt_zigzag`, `epoxide`, `hydroxyl`, `carbonyl_pair`, `carboxyl`, `mixed_vacancy_cnt`, `mixed_sw_cnt`, `pore_slit`, `cnt_interior`, `cnt_interstitial`, `cnt_groove`) |
| `channel` | enum | — | `molecular`, `atomic` |
| `recipe` | dict | — | deterministic construction parameters (supercell, removed atoms, rotations, groups, diameter, chirality, pore width) |
| `delta_mass_per_defect` | float | amu | lattice mass change per defect |
| `sites` | list[SiteSpec] | — | enumerated sites with multiplicity and max occupancy |
| `amorphization_x_max` | UncertainScalar | per C | admissible density bound |

`SiteSpec`: `site_id`, `carbon_indices` (list[int]), `sublattice` (enum A/B/null), `multiplicity` (int), `max_occupancy` (int), `descriptors` (dict: `pyramidalization_deg`, `coordination`, `n_neighbor_H`, `distance_to_defect_nm`, `curvature_inv_nm`).

## SiteConfiguration

`config_hash`, `site_class`, `site_id`, `n_H` (int), `structure` (extended XYZ string), `cell` (3×3), `pbc` (3 bool), `spin_multiplicity`.

## SiteEnergyRecord

| Field | Type | Units |
|---|---|---|
| `config_hash` | str | — |
| `tier` | enum | `T0`, `T1`, `T2`, `T3`, `potential` |
| `engine` | str | — |
| `quantity` | enum | `E_bind_H_half_H2`, `E_bind_H_atom`, `E_bind_H2`, `E_sequential`, `E_barrier_stick`, `E_barrier_recomb`, `E_formation_defect`, `E_migration` |
| `value` | UncertainScalar | eV |
| `reference_state` | enum | `half_H2`, `H_atom`, `H2_molecule`, `pristine_cell`, `adsorbed_state` (release barriers measured from the relaxed adsorbed initial state of the pathway) |
| `temperature_K` | float or null | K (for Feynman–Hibbs-corrected potentials) |
| `converged` | bool | — |
| `provenance` | Provenance | — |

## PublishedValue

`id` (PV-####), `doi`/`arxiv`, `system`, `site_class`, `site_id` (nullable), `quantity`, `value` (UncertainScalar), `unit`, `reference_state`, `method` (dict: `functional`, `dispersion`, `code`, `basis_or_cutoff`, `supercell`, `kmesh`, `spin`, `zpe`), `extraction`, `source_id`, `notes`, `compilation` (bool, optional: the entry summarizes a range of earlier studies and is not pooled), `exclude_from_meta_analysis` (str, optional: the reason a value is registered but left out of the pooled meta-analysis). Registry files may carry a `where_read` field naming the table, figure or paragraph of the open copy in which the value appears; `level` in `method` names a correlated wavefunction method (DMC, CCSD(T)) whose values form the Tier-0 anchors.

## Isotherm (from ISODB)

| Field | Type | Units |
|---|---|---|
| `isodb_filename` | str | — |
| `doi` | str | — |
| `material` | Material | — |
| `adsorbate` | str | `hydrogen` |
| `temperature` | float | K |
| `pressure` | ndarray | bar |
| `loading` | ndarray | mol kg⁻¹ |
| `loading_basis` | enum | `excess`, `absolute`, `unknown` |
| `original_units` | dict | as in source |
| `digitized` | bool | — |
| `tabular` | bool | — |
| `quality_flags` | list[str] | e.g. `nonmonotonic`, `duplicate_pressure`, `unit_ambiguous` |
| `provenance` | Provenance | — |

`Material`: `name`, `registry_id` (nullable), `material_class` (enum: `graphene`, `rGO`, `GO`, `SWCNT`, `MWCNT`, `activated_carbon`, `templated_carbon`, `carbon_other`, `non_carbon`), `bet_area_m2_g` (UncertainScalar, nullable), `pore_volume_cm3_g` (nullable), `oxygen_content_at` (nullable), `conflict_group` (nullable).

## DefectPopulation

`material`, `beam_energy_keV`, `fluence_e_nm2`, `temperature_K`, `densities` (dict class → UncertainScalar, per C), `x_total`, `admissible` (bool with reason), `raman` (dict: `ID_IG` by excitation, `ID_IDp`), `provenance`.

## OperatingWindow

`name`, `T_c_K`, `P_c_bar`, `T_r_K`, `P_r_bar`, `tau_c_s`, `tau_r_s`, `eta_c`, `eta_r`, `charging_mode` (enum: `H2_gas`, `atomic_H_source`), `fugacity` (bool).

## CapacityResult

`window`, `landscape_hash`, `population_hash` (nullable), `D_wt_percent` (UncertainScalar), `D_dischargeable_wt_percent` (nullable), `deficits` (dict: `occupancy`, `admissibility`, `mass`), `limiting_mechanism` (enum), `per_class` (dict), `provenance`.

## BenchmarkRecord

`benchmark_id`, `quantity`, `reference_value` (UncertainScalar), `predicted` (UncertainScalar), `metric`, `threshold`, `passed` (bool), `residual`, `notes`, `provenance`.

## Serialization

- Tables: Parquet (pyarrow) with embedded JSON schema version in file metadata (`sitewindow_schema_version`).
- Posteriors: NetCDF via ArviZ InferenceData; variable naming `<component>.<quantity>[class, site, draw]`.
- Configuration: YAML; validated to the pydantic models on load.
- Manifests: `MANIFEST.sha256` (file → hash) and `reproducibility_manifest.json` (package version, dependency lock hash, ISODB commit, seeds, settings hashes, hardware description, tier).

## Versioning and compatibility

`schema_version` follows semantic versioning; minor increments add optional fields only; readers accept any schema with the same major version and fill missing optional fields with nulls.
