"""Physisorption landscape (manual C2 physisorption set; experiment E5).

For every physisorption geometry of the catalogue (basal plane, vacancy-adjacent sheet, the
sites adjacent to the oxygen groups of the full tier, slit pores, nanotube interiors, bundle
grooves and interstitials) and every pair-potential parameter set, the strongest H2 site energy
is recorded classically and with the Feynman-Hibbs correction at the window temperatures.
Oxygen-adjacent sites use the oxygen pair parameters of ``physisorption.OXYGEN_TYPES`` and are
searched within 4 A of the group. Records carry tier ``potential`` and the parameter-set hash.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from sitewindow import __version__
from sitewindow.atomistic import physisorption as ph
from sitewindow.atomistic.descriptors import site_descriptors
from sitewindow.atomistic.parallel import is_master
from sitewindow.atomistic.structures import builders as b
from sitewindow.constants import EV_TO_KJ_PER_MOL
from sitewindow.io import sha256_file, write_table
from sitewindow.landscape.classes import catalogue
from sitewindow.schema import (
    Extraction,
    Provenance,
    ProvenanceKind,
    Quantity,
    ReferenceState,
    SiteClassName,
    SiteDescriptors,
    SiteEnergyRecord,
    Tier,
    UncertainScalar,
    records_to_frame,
)

TEMPERATURES_K: tuple[float | None, ...] = (None, 77.0, 233.15, 298.15)


@dataclass
class PhysisorptionSpec:
    tier: str = "verify"
    potentials: tuple[str, ...] = ("steele_buch", "steele_darkrim_levesque", "steele_buch_b3")
    temperatures_K: tuple[float | None, ...] = TEMPERATURES_K
    supercell: int = 4
    tube_length: int = 2
    recipes: dict[str, dict[str, Any]] = field(default_factory=dict)
    sigma_kJ_mol: float = 0.5


def _geometries(spec: PhysisorptionSpec) -> list[dict[str, Any]]:
    cat = catalogue()
    n = spec.supercell
    out: list[dict[str, Any]] = [
        {
            "site_class": SiteClassName.basal,
            "site_id": "physisorption_top",
            "name": f"graphene_{n}x{n}",
            "atoms": b.graphene_supercell(n),
            "kind": "sheet",
        },
    ]
    v, _ = b.monovacancy(n)
    out.append(
        {
            "site_class": SiteClassName.monovacancy,
            "site_id": "physisorption_adjacent",
            "name": f"monovacancy_{n}x{n}",
            "atoms": v,
            "kind": "sheet",
        }
    )
    if spec.tier == "full":
        for cls, build in (
            ("epoxide", b.epoxide),
            ("hydroxyl", b.hydroxyl),
            ("carbonyl_pair", b.carbonyl_pair),
            ("carboxyl", b.carboxyl),
        ):
            rec = cat[cls].recipe(spec.tier, spec.recipes.get(cls))
            m = int(rec["supercell"])
            a = build(m)
            out.append(
                {
                    "site_class": SiteClassName(cls),
                    "site_id": "physisorption_adjacent",
                    "name": f"{cls}_{m}x{m}",
                    "atoms": a,
                    "kind": "sheet_local",
                    "near": np.asarray(a.info["defect_center"], dtype=float),
                }
            )
    pore = cat["pore_slit"].recipe(spec.tier, spec.recipes.get("pore_slit"))
    for w in pore["widths_nm"]:
        out.append(
            {
                "site_class": SiteClassName.pore_slit,
                "site_id": f"slit_{w:.2f}nm",
                "name": f"slit_pore_{w:.2f}nm",
                "atoms": b.slit_pore(float(w), n=n),
                "kind": "pore",
                "width_nm": float(w),
            }
        )
    for cls, kind in (
        ("cnt_interior", "tube_interior"),
        ("cnt_interstitial", "interstitial"),
        ("cnt_groove", "groove"),
    ):
        if spec.tier not in cat[cls].recipes and spec.tier == "verify":
            continue
        rec = cat[cls].recipe(spec.tier, spec.recipes.get(cls))
        for d in rec["diameters_nm"]:
            nn, mm = b.chirality_for_diameter(float(d), "armchair")
            t = b.swcnt(nn, mm, length=spec.tube_length)
            atoms = t if kind == "tube_interior" else ph.bundle(t)
            out.append(
                {
                    "site_class": SiteClassName(cls),
                    "site_id": f"{kind}_{d:.2f}nm",
                    "name": f"swcnt_{nn}_{mm}_{kind}",
                    "atoms": atoms,
                    "kind": kind,
                    "diameter_nm": float(d),
                }
            )
    return out


def _site_energy(
    g: dict[str, Any], p: ph.PairPotential, T: float | None
) -> tuple[float, np.ndarray]:
    kind = g["kind"]
    if kind == "sheet":
        return ph.sheet_site_energy(g["atoms"], p, T)
    if kind == "sheet_local":
        return ph.sheet_site_energy(g["atoms"], p, T, near=g["near"])
    if kind == "pore":
        return ph.pore_site_energy(g["atoms"], p, T)
    if kind == "tube_interior":
        return ph.tube_site_energy(g["atoms"], p, True, T)
    return ph.bundle_site_energy(g["atoms"], p, kind, T)


def build_physisorption(spec: PhysisorptionSpec, out_dir: Path) -> dict[str, Any]:
    t0 = time.perf_counter()
    pots = ph.all_potentials()
    records: list[SiteEnergyRecord] = []
    rows: list[dict[str, Any]] = []
    for g in _geometries(spec):
        atoms = g["atoms"]
        carbons = b.carbon_indices(atoms)
        for pname in spec.potentials:
            p = pots[pname]
            for T in spec.temperatures_K:
                e_kj, q = _site_energy(g, p, T)
                nearest = (
                    int(
                        carbons[
                            np.argmin(
                                np.linalg.norm(ph.carbon_images(atoms)[: len(carbons)] - q, axis=1)
                            )
                        ]
                    )
                    if len(carbons)
                    else 0
                )
                desc = site_descriptors(atoms, nearest, atoms.info.get("defect_center"))
                if g["kind"] == "tube_interior":
                    desc["curvature_inv_nm"] = 1.0 / g["diameter_nm"]
                chash = f"{pname}:{g['name']}:{T}"
                records.append(
                    SiteEnergyRecord(
                        config_hash=chash,
                        tier=Tier.potential,
                        engine=f"LJ:{pname}",
                        quantity=Quantity.E_bind_H2,
                        value=UncertainScalar(
                            value=e_kj / EV_TO_KJ_PER_MOL,
                            sigma=spec.sigma_kJ_mol / EV_TO_KJ_PER_MOL,
                            components={"model_form": spec.sigma_kJ_mol / EV_TO_KJ_PER_MOL},
                        ),
                        reference_state=ReferenceState.H2_molecule,
                        temperature_K=T,
                        converged=True,
                        provenance=Provenance(
                            kind=ProvenanceKind.computed_here,
                            source_id=p.settings_hash,
                            identifier=chash,
                            extraction=Extraction.computed,
                            settings_hash=p.settings_hash,
                        ),
                        site_class=g["site_class"],
                        site_id=g["site_id"],
                        n_H=0,
                        descriptors=SiteDescriptors(**{k: v for k, v in desc.items()}),
                        n_atoms=len(atoms),
                        supercell=g["name"],
                    )
                )
                rows.append(
                    {
                        "geometry": g["name"],
                        "site_id": g["site_id"],
                        "potential": pname,
                        "T_K": T,
                        "E_kJ_mol": e_kj,
                        "position": q.round(3).tolist(),
                    }
                )
    frame = records_to_frame(records)
    out: dict[str, Any] = {"records": frame, "table": pd.DataFrame(rows)}
    if is_master():
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        p_rec = write_table(
            frame, out_dir / f"physisorption_{spec.tier}.parquet", {"tier": spec.tier}
        )
        manifest = {
            "tier": spec.tier,
            "package_version": __version__,
            "potentials": {
                k: {"epsilon_K": v.epsilon_K, "sigma_A": v.sigma_A, "source": v.source}
                for k, v in pots.items()
                if k in spec.potentials
            },
            "oxygen_types": {
                k: {"epsilon_K": v["epsilon_K"], "sigma_A": v["sigma_A"], "source": v["source"]}
                for k, v in ph.OXYGEN_TYPES.items()
            },
            "temperatures_K": [t for t in spec.temperatures_K],
            "n_records": int(len(frame)),
            "wall_s": time.perf_counter() - t0,
            "files": {p_rec.name: sha256_file(p_rec)},
        }
        (out_dir / f"physisorption_manifest_{spec.tier}.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        out["manifest"] = manifest
    return out


__all__ = ["PhysisorptionSpec", "build_physisorption", "TEMPERATURES_K"]
