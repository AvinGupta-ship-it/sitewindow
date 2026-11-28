"""Canonical data models (manual Section 6.4; ``docs/specs/DATA_SCHEMA.md``).

Every record carries a :class:`Provenance` block whose ``kind`` keeps experimental observations,
published computations, values computed by this program, model predictions and derived
quantities apart. Canonical units are fixed per field: energies in eV, temperatures in K,
pressures in bar (absolute), loadings in mol kg^-1, site densities per carbon atom, masses in
amu, lengths in nm unless the field name says otherwise.
"""

from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from typing import Any, TypeVar

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from sitewindow import __version__
from sitewindow.errors import SchemaError

SCHEMA_VERSION = "0.2.0"


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


# --------------------------------------------------------------------------- provenance


class ProvenanceKind(StrEnum):
    experimental_observation = "experimental_observation"
    published_computation = "published_computation"
    computed_here = "computed_here"
    model_prediction = "model_prediction"
    derived = "derived"


class Extraction(StrEnum):
    tabular_text = "tabular_text"
    digitized_figure = "digitized_figure"
    abstract_statement = "abstract_statement"
    database_record = "database_record"
    computed = "computed"


def _now() -> datetime:
    return datetime.now(UTC)


def _version() -> str:
    return __version__


class Provenance(_Model):
    """Origin of a record (DATA_SCHEMA, Provenance)."""

    kind: ProvenanceKind
    source_id: str
    identifier: str
    extraction: Extraction | None = None
    settings_hash: str | None = None
    package_version: str = Field(default_factory=_version)
    seed: int | None = None
    created: datetime = Field(default_factory=_now)


# --------------------------------------------------------------------------- uncertainty


class Distribution(StrEnum):
    normal = "normal"
    lognormal = "lognormal"
    uniform = "uniform"
    posterior_samples = "posterior_samples"


class UncertainScalar(_Model):
    """A scalar with uncertainty (DATA_SCHEMA, Uncertainty)."""

    value: float
    sigma: float = Field(default=0.0, ge=0.0)
    distribution: Distribution = Distribution.normal
    samples_ref: str | None = None
    components: dict[str, float] | None = None

    @model_validator(mode="after")
    def _samples_need_reference(self) -> UncertainScalar:
        if self.distribution is Distribution.posterior_samples and not self.samples_ref:
            raise ValueError("posterior_samples distributions require samples_ref")
        return self


# --------------------------------------------------------------------------- site classes


class Channel(StrEnum):
    molecular = "molecular"
    atomic = "atomic"


class Sublattice(StrEnum):
    A = "A"
    B = "B"


class SiteClassName(StrEnum):
    basal = "basal"
    monovacancy = "monovacancy"
    divacancy_585 = "divacancy_585"
    divacancy_555777 = "divacancy_555777"
    stone_wales = "stone_wales"
    edge_zigzag = "edge_zigzag"
    edge_armchair = "edge_armchair"
    cnt_armchair = "cnt_armchair"
    cnt_zigzag = "cnt_zigzag"
    epoxide = "epoxide"
    hydroxyl = "hydroxyl"
    carbonyl_pair = "carbonyl_pair"
    carboxyl = "carboxyl"
    mixed_vacancy_cnt = "mixed_vacancy_cnt"
    mixed_sw_cnt = "mixed_sw_cnt"
    pore_slit = "pore_slit"
    cnt_interior = "cnt_interior"
    cnt_interstitial = "cnt_interstitial"
    cnt_groove = "cnt_groove"


class SiteDescriptors(_Model):
    """Local descriptors of a site carbon (DATA_SCHEMA, SiteSpec.descriptors)."""

    pyramidalization_deg: float | None = None
    coordination: int | None = Field(default=None, ge=0)
    n_neighbor_H: int | None = Field(default=None, ge=0)
    distance_to_defect_nm: float | None = Field(default=None, ge=0.0)
    curvature_inv_nm: float | None = Field(default=None, ge=0.0)


class SiteSpec(_Model):
    site_id: str
    carbon_indices: list[int]
    sublattice: Sublattice | None = None
    multiplicity: int = Field(default=1, ge=1)
    max_occupancy: int = Field(default=1, ge=1)
    descriptors: SiteDescriptors = Field(default_factory=SiteDescriptors)


class SiteClass(_Model):
    name: SiteClassName
    channel: Channel
    recipe: dict[str, Any]
    delta_mass_per_defect_amu: float
    sites: list[SiteSpec]
    amorphization_x_max: UncertainScalar | None = None


class SiteConfiguration(_Model):
    config_hash: str
    site_class: SiteClassName
    site_id: str
    n_H: int = Field(ge=0)
    structure: str
    cell: list[list[float]]
    pbc: tuple[bool, bool, bool]
    spin_multiplicity: int | None = Field(default=None, ge=1)

    @model_validator(mode="after")
    def _cell_shape(self) -> SiteConfiguration:
        if len(self.cell) != 3 or any(len(row) != 3 for row in self.cell):
            raise ValueError("cell must be a 3x3 matrix")
        return self


# --------------------------------------------------------------------------- energies


class Tier(StrEnum):
    T0 = "T0"
    T1 = "T1"
    T2 = "T2"
    T3 = "T3"
    potential = "potential"


class Quantity(StrEnum):
    E_bind_H_half_H2 = "E_bind_H_half_H2"
    E_bind_H_atom = "E_bind_H_atom"
    E_bind_H2 = "E_bind_H2"
    E_sequential = "E_sequential"
    E_barrier_stick = "E_barrier_stick"
    E_barrier_recomb = "E_barrier_recomb"
    E_formation_defect = "E_formation_defect"
    E_migration = "E_migration"


class ReferenceState(StrEnum):
    half_H2 = "half_H2"
    H_atom = "H_atom"
    H2_molecule = "H2_molecule"
    pristine_cell = "pristine_cell"
    adsorbed_state = "adsorbed_state"


COMPUTED_TIERS = frozenset({Tier.T1, Tier.T2, Tier.T3, Tier.potential})


class SiteEnergyRecord(_Model):
    """One energy of one configuration at one fidelity tier (DATA_SCHEMA, SiteEnergyRecord).

    The optional fields beyond the specification (site class, site identifier, hydrogen count,
    descriptors, cell size) make the landscape table self-describing; the schema minor version
    permits optional additions only.
    """

    config_hash: str
    tier: Tier
    engine: str
    quantity: Quantity
    value: UncertainScalar
    reference_state: ReferenceState
    temperature_K: float | None = Field(default=None, gt=0.0)
    converged: bool
    provenance: Provenance
    site_class: SiteClassName | None = None
    site_id: str | None = None
    n_H: int | None = Field(default=None, ge=0)
    descriptors: SiteDescriptors | None = None
    n_atoms: int | None = Field(default=None, ge=1)
    supercell: str | None = None

    @model_validator(mode="after")
    def _computed_records_carry_settings(self) -> SiteEnergyRecord:
        if self.tier in COMPUTED_TIERS and not self.provenance.settings_hash:
            raise ValueError(f"tier {self.tier.value} records must carry provenance.settings_hash")
        if self.tier in COMPUTED_TIERS and self.provenance.kind is not ProvenanceKind.computed_here:
            raise ValueError("computed tiers require provenance.kind = computed_here")
        if abs(self.value.value) > 50.0:
            raise ValueError("energy outside the plausible range of +/-50 eV")
        return self


class PublishedValue(_Model):
    id: str = Field(pattern=r"^PV-\d{4}$")
    doi: str | None = None
    arxiv: str | None = None
    system: str
    site_class: SiteClassName | None = None
    site_id: str | None = None
    quantity: Quantity
    value: UncertainScalar
    unit: str
    reference_state: ReferenceState
    method: dict[str, Any] = Field(default_factory=dict)
    extraction: Extraction
    source_id: str
    notes: str | None = None
    compilation: bool = False
    """The entry summarizes a range of earlier studies rather than reporting one study's value."""
    exclude_from_meta_analysis: str | None = None
    """Reason the value is kept in the registry but left out of the pooled meta-analysis."""

    @model_validator(mode="after")
    def _identifier_present(self) -> PublishedValue:
        if not (self.doi or self.arxiv):
            raise ValueError("a published value needs a DOI or an arXiv identifier")
        return self


# --------------------------------------------------------------------------- corpus


class MaterialClass(StrEnum):
    graphene = "graphene"
    rGO = "rGO"
    GO = "GO"
    SWCNT = "SWCNT"
    MWCNT = "MWCNT"
    activated_carbon = "activated_carbon"
    templated_carbon = "templated_carbon"
    carbon_other = "carbon_other"
    non_carbon = "non_carbon"


class Material(_Model):
    name: str
    registry_id: str | None = None
    material_class: MaterialClass
    bet_area_m2_g: UncertainScalar | None = None
    pore_volume_cm3_g: float | None = Field(default=None, ge=0.0)
    oxygen_content_at: float | None = Field(default=None, ge=0.0, le=100.0)
    conflict_group: str | None = None


class LoadingBasis(StrEnum):
    excess = "excess"
    absolute = "absolute"
    unknown = "unknown"


class Isotherm(_Model):
    isodb_filename: str
    doi: str
    material: Material
    adsorbate: str = "hydrogen"
    temperature: float = Field(gt=0.0)
    pressure: list[float]
    loading: list[float]
    loading_basis: LoadingBasis = LoadingBasis.unknown
    original_units: dict[str, str] = Field(default_factory=dict)
    digitized: bool = False
    tabular: bool = False
    quality_flags: list[str] = Field(default_factory=list)
    provenance: Provenance

    @model_validator(mode="after")
    def _ranges(self) -> Isotherm:
        if len(self.pressure) != len(self.loading):
            raise ValueError("pressure and loading must have the same length")
        if len(self.pressure) < 3:
            raise ValueError("an isotherm needs at least three points")
        if any(p < 0.0 for p in self.pressure):
            raise ValueError("pressures must be non-negative")
        if any(n < 0.0 for n in self.loading) and "negative_loading" not in self.quality_flags:
            raise ValueError("negative loadings must be flagged with 'negative_loading'")
        return self


# --------------------------------------------------------------------------- populations


class DefectPopulation(_Model):
    material: str
    beam_energy_keV: float = Field(gt=0.0)
    fluence_e_nm2: float = Field(ge=0.0)
    temperature_K: float = Field(gt=0.0)
    densities: dict[str, UncertainScalar]
    x_total: UncertainScalar
    admissible: bool
    admissibility_reason: str | None = None
    raman: dict[str, Any] = Field(default_factory=dict)
    provenance: Provenance


# --------------------------------------------------------------------------- windows


class ChargingMode(StrEnum):
    H2_gas = "H2_gas"
    atomic_H_source = "atomic_H_source"


class OperatingWindow(_Model):
    """Operating window as configuration data (DATA_SCHEMA, OperatingWindow)."""

    name: str
    T_c_K: float = Field(gt=0.0)
    P_c_bar: float = Field(ge=0.0)
    T_r_K: float = Field(gt=0.0)
    P_r_bar: float = Field(ge=0.0)
    tau_c_s: float = Field(default=300.0, gt=0.0)
    tau_r_s: float = Field(default=1.0e3, gt=0.0)
    eta_c: float = Field(default=0.9, gt=0.0, lt=1.0)
    eta_r: float = Field(default=0.9, gt=0.0, lt=1.0)
    charging_mode: ChargingMode = ChargingMode.H2_gas
    fugacity: bool = True

    @model_validator(mode="after")
    def _gas_charging_needs_pressure(self) -> OperatingWindow:
        if self.charging_mode is ChargingMode.H2_gas and self.P_c_bar <= self.P_r_bar:
            raise ValueError("gas-phase charging requires P_c > P_r")
        return self


class LimitingMechanism(StrEnum):
    thermodynamic = "thermodynamic"
    kinetic_release = "kinetic_release"
    kinetic_charge = "kinetic_charge"
    mass_limited = "mass_limited"


class CapacityResult(_Model):
    window: str
    landscape_hash: str
    population_hash: str | None = None
    D_wt_percent: UncertainScalar
    D_dischargeable_wt_percent: UncertainScalar | None = None
    deficits: dict[str, float] = Field(default_factory=dict)
    limiting_mechanism: LimitingMechanism
    per_class: dict[str, Any] = Field(default_factory=dict)
    provenance: Provenance


class BenchmarkRecord(_Model):
    benchmark_id: str = Field(pattern=r"^B\d{1,2}$")
    quantity: str
    reference_value: UncertainScalar
    predicted: UncertainScalar
    metric: str
    threshold: float
    passed: bool
    residual: float
    notes: str | None = None
    provenance: Provenance


# --------------------------------------------------------------------------- helpers

M = TypeVar("M", bound=BaseModel)


def parse(model_type: type[M], data: Any) -> M:
    """Validate ``data`` against ``model_type``; raise :class:`SchemaError` on failure."""
    try:
        return model_type.model_validate(data)
    except ValidationError as e:
        raise SchemaError(str(e)) from e


def json_schema(model_type: type[BaseModel]) -> dict[str, Any]:
    """JSON Schema export of a model (``sitewindow schema export``)."""
    return model_type.model_json_schema()


__all__ = [
    "SCHEMA_VERSION",
    "ProvenanceKind",
    "Extraction",
    "Provenance",
    "Distribution",
    "UncertainScalar",
    "Channel",
    "Sublattice",
    "SiteClassName",
    "SiteDescriptors",
    "SiteSpec",
    "SiteClass",
    "SiteConfiguration",
    "Tier",
    "Quantity",
    "ReferenceState",
    "COMPUTED_TIERS",
    "SiteEnergyRecord",
    "PublishedValue",
    "MaterialClass",
    "Material",
    "LoadingBasis",
    "Isotherm",
    "DefectPopulation",
    "ChargingMode",
    "OperatingWindow",
    "LimitingMechanism",
    "CapacityResult",
    "BenchmarkRecord",
    "parse",
    "json_schema",
]
