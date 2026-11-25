"""Schema tests (manual Section 6.4; acceptance B6: records carry tier and settings hash)."""

import pytest

from sitewindow.errors import SchemaError
from sitewindow.schema import (
    Distribution,
    Isotherm,
    Material,
    OperatingWindow,
    Provenance,
    ProvenanceKind,
    Quantity,
    ReferenceState,
    SiteClassName,
    SiteEnergyRecord,
    Tier,
    UncertainScalar,
    frame_to_records,
    parse,
    records_to_frame,
)


def _prov(**kw):
    base = dict(
        kind=ProvenanceKind.computed_here, source_id="cfg", identifier="basal/graphene_4x4/monomer"
    )
    base.update(kw)
    return Provenance(**base)


def _record(**kw):
    base = dict(
        config_hash="abc",
        tier=Tier.T2,
        engine="tblite:GFN2-xTB",
        quantity=Quantity.E_bind_H_atom,
        value=UncertainScalar(value=-0.21, sigma=0.05),
        reference_state=ReferenceState.H_atom,
        converged=True,
        provenance=_prov(settings_hash="deadbeef"),
        site_class=SiteClassName.basal,
        site_id="top",
        n_H=1,
    )
    base.update(kw)
    return SiteEnergyRecord(**base)


def test_computed_records_require_settings_hash():
    with pytest.raises(ValueError):
        _record(provenance=_prov(settings_hash=None))


def test_computed_tiers_require_computed_here_kind():
    with pytest.raises(ValueError):
        _record(provenance=_prov(kind=ProvenanceKind.published_computation, settings_hash="x"))


def test_tier0_records_need_no_settings_hash():
    r = _record(
        tier=Tier.T0, engine="DMC", provenance=_prov(kind=ProvenanceKind.published_computation)
    )
    assert r.tier is Tier.T0


def test_energy_range_check():
    with pytest.raises(ValueError):
        _record(value=UncertainScalar(value=-120.0))


def test_posterior_samples_need_reference():
    with pytest.raises(ValueError):
        UncertainScalar(value=1.0, distribution=Distribution.posterior_samples)
    assert UncertainScalar(
        value=1.0, distribution=Distribution.posterior_samples, samples_ref="a.nc"
    )


def test_parse_raises_schema_error():
    with pytest.raises(SchemaError):
        parse(UncertainScalar, {"value": "not a number"})


def test_record_table_round_trip():
    recs = [
        _record(),
        _record(quantity=Quantity.E_sequential, reference_state=ReferenceState.half_H2, n_H=2),
    ]
    frame = records_to_frame(recs)
    back = frame_to_records(frame)
    assert [r.quantity for r in back] == [r.quantity for r in recs]
    assert back[0].provenance.settings_hash == "deadbeef"
    assert back[1].n_H == 2 and back[0].descriptors is not None


def test_isotherm_negative_loading_must_be_flagged():
    mat = Material(name="AX-21", material_class="activated_carbon")
    prov = Provenance(
        kind=ProvenanceKind.experimental_observation, source_id="S1", identifier="f.json"
    )
    with pytest.raises(ValueError):
        Isotherm(
            isodb_filename="f",
            doi="10.1/x",
            material=mat,
            temperature=77.0,
            pressure=[0, 1, 2],
            loading=[0.0, -0.1, 0.2],
            provenance=prov,
        )
    ok = Isotherm(
        isodb_filename="f",
        doi="10.1/x",
        material=mat,
        temperature=77.0,
        pressure=[0, 1, 2],
        loading=[0.0, -0.1, 0.2],
        quality_flags=["negative_loading"],
        provenance=prov,
    )
    assert ok.loading_basis.value == "unknown"


def test_operating_window_requires_pressure_swing_for_gas_charging():
    with pytest.raises(ValueError):
        OperatingWindow(name="bad", T_c_K=298.0, P_c_bar=5.0, T_r_K=298.0, P_r_bar=5.0)
    w = OperatingWindow(
        name="W4",
        T_c_K=300.0,
        P_c_bar=0.0,
        T_r_K=473.0,
        P_r_bar=5.0,
        charging_mode="atomic_H_source",
    )
    assert w.charging_mode.value == "atomic_H_source"
