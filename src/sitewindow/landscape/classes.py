"""Site-class catalogue: construction recipes, enumerated sites and hydrogen configurations.

Implements the minimum class set of manual C2 (table "Site classes enumerated") and the
``SiteClass`` extension point of Section 10.3: a class is a :class:`SiteClassDefinition` whose
``enumerate`` function returns the :class:`Configuration` objects (substrate, site, hydrogen
placements) to compute at a tier. Plugins register additional definitions through the entry
point group ``sitewindow.site_classes``.

Multiplicities count equivalent sites per defect; for the basal plane and curved sheets the
"defect" is one carbon atom, so the site density of the class is the number of carbons. Mass
changes are per defect relative to the pristine lattice.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from importlib.metadata import entry_points
from typing import Any

import numpy as np
from ase import Atoms

from sitewindow.atomistic.descriptors import site_descriptors
from sitewindow.atomistic.structures import builders as b
from sitewindow.constants import M_C_AMU, M_H_AMU, M_O_AMU
from sitewindow.schema import Channel, SiteClass, SiteClassName, SiteDescriptors, SiteSpec

VERIFY = "verify"
FULL = "full"
BOTH = (VERIFY, FULL)
FULL_ONLY = (FULL,)

Placement = tuple[int, str]


@dataclass(frozen=True)
class Configuration:
    """One hydrogenated configuration to compute: substrate plus hydrogen placements."""

    site_class: SiteClassName
    substrate_name: str
    substrate: Atoms
    site: SiteSpec
    name: str
    placements: tuple[Placement, ...]
    parent: str | None
    tiers: tuple[str, ...]
    pristine: Atoms | None = None
    pristine_name: str | None = None

    @property
    def n_H(self) -> int:
        return len(self.placements)

    @property
    def key(self) -> str:
        return f"{self.site_class.value}/{self.substrate_name}/{self.name}"


@dataclass(frozen=True)
class SiteClassDefinition:
    name: SiteClassName
    channel: Channel
    description: str
    delta_mass_per_defect_amu: float
    recipes: dict[str, dict[str, Any]]
    enumerate: Callable[[dict[str, Any], str], list[Configuration]]
    tiers: tuple[str, ...] = ("T1", "T2", "T3")

    def recipe(self, tier: str, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
        base = self.recipes.get(tier, self.recipes[FULL])
        return {**base, **(overrides or {})}

    def configurations(
        self, tier: str, overrides: dict[str, Any] | None = None
    ) -> list[Configuration]:
        return [c for c in self.enumerate(self.recipe(tier, overrides), tier) if tier in c.tiers]

    def site_class(self, tier: str, overrides: dict[str, Any] | None = None) -> SiteClass:
        """Schema record of the class: recipe, mass change and the unique sites enumerated."""
        sites: dict[str, SiteSpec] = {}
        for c in self.configurations(tier, overrides):
            sites.setdefault(c.site.site_id, c.site)
        return SiteClass(
            name=self.name,
            channel=self.channel,
            recipe=self.recipe(tier, overrides),
            delta_mass_per_defect_amu=self.delta_mass_per_defect_amu,
            sites=list(sites.values()),
        )


def _spec(
    site_id: str,
    indices: Sequence[int],
    atoms: Atoms,
    multiplicity: int,
    max_occupancy: int = 1,
    defect_center: np.ndarray | None = None,
) -> SiteSpec:
    d = site_descriptors(atoms, int(indices[0]), defect_center)
    return SiteSpec(
        site_id=site_id,
        carbon_indices=[int(k) for k in indices],
        multiplicity=multiplicity,
        max_occupancy=max_occupancy,
        descriptors=SiteDescriptors(**d),
    )


def _cfg(
    cls: SiteClassName,
    sub_name: str,
    sub: Atoms,
    site: SiteSpec,
    name: str,
    placements: Sequence[Placement],
    parent: str | None,
    tiers: tuple[str, ...],
    pristine: Atoms | None = None,
    pristine_name: str | None = None,
) -> Configuration:
    return Configuration(
        cls, sub_name, sub, site, name, tuple(placements), parent, tiers, pristine, pristine_name
    )


# --------------------------------------------------------------------------- flat sheets


def _basal(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    n = int(recipe["supercell"])
    g = b.graphene_supercell(n)
    c0 = b.central_atom_index(g)
    h = b.hexagon(g, c0)
    site = _spec("top", [c0], g, multiplicity=1)
    sub = f"graphene_{n}x{n}"
    cls = SiteClassName.basal
    up = "up"
    return [
        _cfg(cls, sub, g, site, "monomer", [(h[0], up)], None, BOTH),
        _cfg(cls, sub, g, site, "ortho", [(h[0], up), (h[1], up)], "monomer", BOTH),
        _cfg(cls, sub, g, site, "meta", [(h[0], up), (h[2], up)], "monomer", BOTH),
        _cfg(cls, sub, g, site, "para", [(h[0], up), (h[3], up)], "monomer", BOTH),
        _cfg(
            cls,
            sub,
            g,
            site,
            "ortho_alternating",
            [(h[0], up), (h[1], "down")],
            "monomer",
            FULL_ONLY,
        ),
        _cfg(
            cls,
            sub,
            g,
            site,
            "para_alternating",
            [(h[0], up), (h[3], "down")],
            "monomer",
            FULL_ONLY,
        ),
        _cfg(cls, sub, g, site, "trimer", [(h[0], up), (h[1], up), (h[3], up)], "ortho", FULL_ONLY),
        _cfg(
            cls,
            sub,
            g,
            site,
            "tetramer",
            [(h[0], up), (h[1], up), (h[3], up), (h[4], up)],
            "trimer",
            FULL_ONLY,
        ),
        _cfg(cls, sub, g, site, "graphane", b.graphane_placements(g), None, FULL_ONLY),
    ]


def _monovacancy(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    n = int(recipe["supercell"])
    v, _ = b.monovacancy(n)
    g = b.graphene_supercell(n)
    centre = np.asarray(v.info["defect_center"])
    nb = b.neighbor_map(v)
    dangling = sorted(k for k in b.carbon_indices(v) if len(nb[k]) == 2)
    adjacent = sorted({j for k in dangling for j in nb[k]} - set(dangling))
    s_d = _spec("dangling", dangling, v, 3, defect_center=centre)
    s_a = _spec("basal_adjacent", adjacent, v, 6, defect_center=centre)
    order: list[Placement] = [(k, "hole") for k in dangling] + [(k, "up") for k in adjacent]
    n_max = int(recipe.get("max_hydrogens", 6))
    out: list[Configuration] = []
    prev: str | None = None
    for i in range(1, n_max + 1):
        name = f"seq{i}"
        tiers = BOTH if i <= 3 else FULL_ONLY
        out.append(
            _cfg(
                SiteClassName.monovacancy,
                f"monovacancy_{n}x{n}",
                v,
                s_d if i <= 3 else s_a,
                name,
                order[:i],
                prev,
                tiers,
                g,
                f"graphene_{n}x{n}",
            )
        )
        prev = name
    return out


def _divacancy(kind: SiteClassName) -> Callable[[dict[str, Any], str], list[Configuration]]:
    def enumerate_(recipe: dict[str, Any], tier: str) -> list[Configuration]:
        n = int(recipe["supercell"])
        dv = b.divacancy_585(n) if kind is SiteClassName.divacancy_585 else b.divacancy_555777(n)
        g = b.graphene_supercell(n)
        centre = np.asarray(dv.info["defect_center"])
        recon = list(dv.info["reconstruction_carbons"])
        s_r = _spec("reconstruction", recon, dv, 4, defect_center=centre)
        sub = f"{kind.value}_{n}x{n}"
        out: list[Configuration] = []
        prev: str | None = None
        for i in range(1, 5):
            name = f"seq{i}"
            out.append(
                _cfg(
                    kind,
                    sub,
                    dv,
                    s_r,
                    name,
                    [(k, "up") for k in recon[:i]],
                    prev,
                    FULL_ONLY,
                    g,
                    f"graphene_{n}x{n}",
                )
            )
            prev = name
        if kind is SiteClassName.divacancy_555777:
            pair = list(dv.info["rotated_pair"])
            s_p = _spec("rotated_pair", pair, dv, 2, defect_center=centre)
            out.append(
                _cfg(
                    kind,
                    sub,
                    dv,
                    s_p,
                    "rotated_monomer",
                    [(pair[0], "up")],
                    None,
                    FULL_ONLY,
                    g,
                    f"graphene_{n}x{n}",
                )
            )
        return out

    return enumerate_


def _distance_shells(
    atoms: Atoms, centre: np.ndarray, n_shells: int, decimals: int = 1
) -> list[list[int]]:
    from ase.geometry import find_mic

    carbons = b.carbon_indices(atoms)
    _, dist = find_mic(atoms.positions[carbons] - centre, atoms.cell, atoms.pbc)
    keys = np.round(dist, decimals)
    shells: dict[float, list[int]] = {}
    for k, d in zip(carbons, keys, strict=True):
        shells.setdefault(float(d), []).append(int(k))
    return [sorted(shells[d]) for d in sorted(shells)[:n_shells]]


def _stone_wales(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    n = int(recipe["supercell"])
    sw = b.stone_wales(n)
    g = b.graphene_supercell(n)
    centre = np.asarray(sw.info["defect_center"])
    shells = _distance_shells(sw, centre, 4)
    names = ("rotated_bond", "first_neighbor", "second_neighbor_a", "second_neighbor_b")
    sub = f"stone_wales_{n}x{n}"
    out: list[Configuration] = []
    for shell, name in zip(shells, names, strict=True):
        site = _spec(name, shell, sw, len(shell), defect_center=centre)
        tiers = BOTH if name == "rotated_bond" else FULL_ONLY
        out.append(
            _cfg(
                SiteClassName.stone_wales,
                sub,
                sw,
                site,
                f"monomer_{name}",
                [(shell[0], "up")],
                None,
                tiers,
                g,
                f"graphene_{n}x{n}",
            )
        )
    pair = shells[0]
    s0 = _spec(names[0], pair, sw, len(pair), defect_center=centre)
    out.append(
        _cfg(
            SiteClassName.stone_wales,
            sub,
            sw,
            s0,
            "dimer_rotated_bond",
            [(pair[0], "up"), (pair[1], "up")],
            "monomer_rotated_bond",
            FULL_ONLY,
            g,
            f"graphene_{n}x{n}",
        )
    )
    return out


# --------------------------------------------------------------------------- edges


def _edge(kind: SiteClassName) -> Callable[[dict[str, Any], str], list[Configuration]]:
    edge = "zigzag" if kind is SiteClassName.edge_zigzag else "armchair"

    def enumerate_(recipe: dict[str, Any], tier: str) -> list[Configuration]:
        width, length = int(recipe["width"]), int(recipe["length"])
        r = b.nanoribbon(edge, width, length)
        ec = b.edge_carbons(r)
        sc = b.subedge_carbons(r)
        bare, ci = b.remove_edge_hydrogen(r, ec[0])
        s_edge = _spec(
            "edge_bare", [ci], bare, 1, defect_center=np.asarray(bare.info["defect_center"])
        )
        s_sub = _spec("subedge_basal", sc, r, 1)
        base = f"ribbon_{edge}_{width}x{length}"
        return [
            _cfg(
                kind, base + "_bare_edge", bare, s_edge, "edge_CH", [(ci, "edge")], None, FULL_ONLY
            ),
            _cfg(kind, base, r, s_sub, "subedge_monomer", [(sc[0], "up")], None, FULL_ONLY),
        ]

    return enumerate_


# --------------------------------------------------------------------------- nanotubes


def _cnt(kind: SiteClassName) -> Callable[[dict[str, Any], str], list[Configuration]]:
    chirality = "armchair" if kind is SiteClassName.cnt_armchair else "zigzag"

    def enumerate_(recipe: dict[str, Any], tier: str) -> list[Configuration]:
        length = int(recipe.get("length", 3))
        out: list[Configuration] = []
        for d in recipe["diameters_nm"]:
            n, m = b.chirality_for_diameter(float(d), chirality)
            t = b.swcnt(n, m, length=length)
            c0 = b.tube_middle_atom(t)
            h = b.hexagon(t, c0)
            sub = f"swcnt_{n}_{m}_x{length}"
            s_ext = _spec("exterior_top", [c0], t, 1)
            s_int = _spec("interior_top", [c0], t, 1)
            out += [
                _cfg(kind, sub, t, s_ext, "monomer_exterior", [(c0, "up")], None, BOTH),
                _cfg(
                    kind,
                    sub,
                    t,
                    s_ext,
                    "ortho_exterior",
                    [(c0, "up"), (h[1], "up")],
                    "monomer_exterior",
                    BOTH,
                ),
                _cfg(
                    kind,
                    sub,
                    t,
                    s_ext,
                    "para_exterior",
                    [(c0, "up"), (h[3], "up")],
                    "monomer_exterior",
                    BOTH,
                ),
                _cfg(kind, sub, t, s_int, "monomer_interior", [(c0, "down")], None, FULL_ONLY),
            ]
            z0 = t.positions[c0, 2]
            ring = sorted(int(k) for k in b.carbon_indices(t) if abs(t.positions[k, 2] - z0) < 0.3)
            if len(ring) >= 2:
                out.append(
                    _cfg(
                        kind,
                        sub,
                        t,
                        s_ext,
                        "stripe_exterior",
                        [(k, "up") for k in ring],
                        None,
                        FULL_ONLY,
                    )
                )
        return out

    return enumerate_


def _mixed_vacancy_cnt(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    out: list[Configuration] = []
    length = int(recipe.get("length", 3))
    for n, m in recipe["tubes"]:
        tv = b.swcnt_monovacancy(n, m, length)
        t = b.swcnt(n, m, length)
        nb = b.neighbor_map(tv)
        dangling = sorted(k for k in b.carbon_indices(tv) if len(nb[k]) == 2)
        centre = np.asarray(tv.info["defect_center"])
        site = _spec("dangling", dangling, tv, 3, defect_center=centre)
        sub = f"swcnt_{n}_{m}_x{length}_vacancy"
        out.append(
            _cfg(
                SiteClassName.mixed_vacancy_cnt,
                sub,
                tv,
                site,
                "seq1",
                [(dangling[0], "hole")],
                None,
                FULL_ONLY,
                t,
                f"swcnt_{n}_{m}_x{length}",
            )
        )
        out.append(
            _cfg(
                SiteClassName.mixed_vacancy_cnt,
                sub,
                tv,
                site,
                "seq2",
                [(dangling[0], "hole"), (dangling[1], "hole")],
                "seq1",
                FULL_ONLY,
                t,
                f"swcnt_{n}_{m}_x{length}",
            )
        )
    return out


def _mixed_sw_cnt(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    out: list[Configuration] = []
    length = int(recipe.get("length", 3))
    for n, m in recipe["tubes"]:
        ts = b.swcnt_stone_wales(n, m, length)
        t = b.swcnt(n, m, length)
        pair = list(ts.info["sw_pair"])
        centre = np.asarray(ts.info["defect_center"])
        site = _spec("rotated_bond", pair, ts, 2, defect_center=centre)
        sub = f"swcnt_{n}_{m}_x{length}_stone_wales"
        out.append(
            _cfg(
                SiteClassName.mixed_sw_cnt,
                sub,
                ts,
                site,
                "monomer",
                [(pair[0], "up")],
                None,
                FULL_ONLY,
                t,
                f"swcnt_{n}_{m}_x{length}",
            )
        )
        out.append(
            _cfg(
                SiteClassName.mixed_sw_cnt,
                sub,
                ts,
                site,
                "dimer",
                [(pair[0], "up"), (pair[1], "up")],
                "monomer",
                FULL_ONLY,
                t,
                f"swcnt_{n}_{m}_x{length}",
            )
        )
    return out


# --------------------------------------------------------------------------- oxygen groups


def _oxygen(kind: SiteClassName) -> Callable[[dict[str, Any], str], list[Configuration]]:
    def enumerate_(recipe: dict[str, Any], tier: str) -> list[Configuration]:
        n = int(recipe["supercell"])
        if kind is SiteClassName.epoxide:
            a = b.epoxide(n)
        elif kind is SiteClassName.hydroxyl:
            a = b.hydroxyl(n)
        elif kind is SiteClassName.carbonyl_pair:
            a = b.carbonyl_pair(n)
        else:
            a = b.carboxyl(n)
        centre = np.asarray(a.info["defect_center"])
        nb = b.neighbor_map(a)
        sub = f"{kind.value}_{n}x{n}"
        if kind is SiteClassName.carbonyl_pair:
            bare = list(a.info["bare_carbons"])
            site = _spec("hole_carbon", bare, a, 2, defect_center=centre)
            return [_cfg(kind, sub, a, site, "monomer_hole", [(bare[0], "hole")], None, FULL_ONLY)]
        func = set(a.info["functionalized_carbons"])
        adjacent = sorted({j for c in func for j in nb[c]} - func)
        site = _spec("adjacent", adjacent, a, len(adjacent), defect_center=centre)
        return [
            _cfg(kind, sub, a, site, "monomer_adjacent", [(adjacent[0], "up")], None, FULL_ONLY)
        ]

    return enumerate_


def _none(recipe: dict[str, Any], tier: str) -> list[Configuration]:
    return []


# --------------------------------------------------------------------------- catalogue


def _defs() -> list[SiteClassDefinition]:
    atomic = Channel.atomic
    molecular = Channel.molecular
    return [
        SiteClassDefinition(
            SiteClassName.basal,
            atomic,
            "Pristine graphene basal plane: monomer, ortho/meta/para dimers, clusters, graphane limit.",
            0.0,
            {VERIFY: {"supercell": 4}, FULL: {"supercell": 6}},
            _basal,
        ),
        SiteClassDefinition(
            SiteClassName.monovacancy,
            atomic,
            "Monovacancy: sequential hydrogenation of the dangling bonds, then adjacent basal carbons.",
            -M_C_AMU,
            {
                VERIFY: {"supercell": 4, "max_hydrogens": 3},
                FULL: {"supercell": 6, "max_hydrogens": 6},
            },
            _monovacancy,
        ),
        SiteClassDefinition(
            SiteClassName.divacancy_585,
            atomic,
            "5-8-5 divacancy: 1-4 H at the reconstruction carbons.",
            -2 * M_C_AMU,
            {FULL: {"supercell": 7}},
            _divacancy(SiteClassName.divacancy_585),
        ),
        SiteClassDefinition(
            SiteClassName.divacancy_555777,
            atomic,
            "555-777 divacancy: 1-4 H at the reconstruction carbons and the rotated pair.",
            -2 * M_C_AMU,
            {FULL: {"supercell": 7}},
            _divacancy(SiteClassName.divacancy_555777),
        ),
        SiteClassDefinition(
            SiteClassName.stone_wales,
            atomic,
            "Stone-Wales defect: H at the four inequivalent sites around the rotated bond.",
            0.0,
            {VERIFY: {"supercell": 4}, FULL: {"supercell": 6}},
            _stone_wales,
        ),
        SiteClassDefinition(
            SiteClassName.edge_zigzag,
            atomic,
            "Zigzag edge: edge C-H on a bare edge carbon and sub-edge basal site of the H-terminated ribbon.",
            M_H_AMU,
            {FULL: {"width": 6, "length": 4}},
            _edge(SiteClassName.edge_zigzag),
        ),
        SiteClassDefinition(
            SiteClassName.edge_armchair,
            atomic,
            "Armchair edge: edge C-H on a bare edge carbon and sub-edge basal site of the H-terminated ribbon.",
            M_H_AMU,
            {FULL: {"width": 10, "length": 3}},
            _edge(SiteClassName.edge_armchair),
        ),
        SiteClassDefinition(
            SiteClassName.cnt_armchair,
            atomic,
            "Armchair (n,n) nanotubes: exterior and interior monomers, dimers and a high-coverage stripe.",
            0.0,
            {
                VERIFY: {"diameters_nm": [0.68, 1.09, 1.36], "length": 3},
                FULL: {"diameters_nm": [0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5], "length": 3},
            },
            _cnt(SiteClassName.cnt_armchair),
        ),
        SiteClassDefinition(
            SiteClassName.cnt_zigzag,
            atomic,
            "Zigzag (n,0) nanotubes: exterior and interior monomers, dimers and a high-coverage stripe.",
            0.0,
            {FULL: {"diameters_nm": [0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5], "length": 2}},
            _cnt(SiteClassName.cnt_zigzag),
        ),
        SiteClassDefinition(
            SiteClassName.epoxide,
            atomic,
            "Epoxide on graphene: H on a carbon adjacent to the group.",
            M_O_AMU,
            {FULL: {"supercell": 6}},
            _oxygen(SiteClassName.epoxide),
        ),
        SiteClassDefinition(
            SiteClassName.hydroxyl,
            atomic,
            "Hydroxyl on graphene: H on a carbon adjacent to the group.",
            M_O_AMU + M_H_AMU,
            {FULL: {"supercell": 6}},
            _oxygen(SiteClassName.hydroxyl),
        ),
        SiteClassDefinition(
            SiteClassName.carbonyl_pair,
            atomic,
            "Carbonyl pair at a divacancy hole: H on a bare hole carbon.",
            2 * M_O_AMU - 2 * M_C_AMU,
            {FULL: {"supercell": 7}},
            _oxygen(SiteClassName.carbonyl_pair),
        ),
        SiteClassDefinition(
            SiteClassName.carboxyl,
            atomic,
            "Carboxyl on graphene: H on a carbon adjacent to the group.",
            M_C_AMU + 2 * M_O_AMU + M_H_AMU,
            {FULL: {"supercell": 6}},
            _oxygen(SiteClassName.carboxyl),
        ),
        SiteClassDefinition(
            SiteClassName.mixed_vacancy_cnt,
            atomic,
            "Monovacancy on (10,0) and (8,8) nanotubes: dangling-bond hydrogenation.",
            -M_C_AMU,
            {FULL: {"tubes": [[10, 0], [8, 8]], "length": 3}},
            _mixed_vacancy_cnt,
        ),
        SiteClassDefinition(
            SiteClassName.mixed_sw_cnt,
            atomic,
            "Stone-Wales defect on (10,0) and (8,8) nanotubes.",
            0.0,
            {FULL: {"tubes": [[10, 0], [8, 8]], "length": 3}},
            _mixed_sw_cnt,
        ),
        SiteClassDefinition(
            SiteClassName.pore_slit,
            molecular,
            "H2 physisorption in graphitic slit pores of width 0.6-1.2 nm (pairwise potential).",
            0.0,
            {
                VERIFY: {"widths_nm": [0.7, 1.0]},
                FULL: {"widths_nm": [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]},
            },
            _none,
            ("potential",),
        ),
        SiteClassDefinition(
            SiteClassName.cnt_interior,
            molecular,
            "H2 physisorption inside (n,n) nanotubes (pairwise potential).",
            0.0,
            {
                VERIFY: {"diameters_nm": [1.0, 1.4]},
                FULL: {"diameters_nm": [0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5]},
            },
            _none,
            ("potential",),
        ),
        SiteClassDefinition(
            SiteClassName.cnt_interstitial,
            molecular,
            "H2 physisorption in the interstitial channel of a nanotube bundle (pairwise potential).",
            0.0,
            {FULL: {"diameters_nm": [1.0, 1.4, 2.0]}},
            _none,
            ("potential",),
        ),
        SiteClassDefinition(
            SiteClassName.cnt_groove,
            molecular,
            "H2 physisorption in the groove between two tubes of a bundle (pairwise potential).",
            0.0,
            {FULL: {"diameters_nm": [1.0, 1.4, 2.0]}},
            _none,
            ("potential",),
        ),
    ]


_CATALOGUE: dict[str, SiteClassDefinition] | None = None


def catalogue(refresh: bool = False) -> dict[str, SiteClassDefinition]:
    """Built-in site classes plus plugins registered under ``sitewindow.site_classes``."""
    global _CATALOGUE
    if _CATALOGUE is None or refresh:
        out = {d.name.value: d for d in _defs()}
        for ep in entry_points(group="sitewindow.site_classes"):
            try:
                d = ep.load()
            except Exception:  # pragma: no cover - broken plugin
                continue
            if isinstance(d, SiteClassDefinition):
                out[d.name.value] = d
        _CATALOGUE = out
    return _CATALOGUE


def chemisorption_classes() -> list[str]:
    return [k for k, d in catalogue().items() if "T2" in d.tiers]


def physisorption_classes() -> list[str]:
    return [k for k, d in catalogue().items() if d.tiers == ("potential",)]


__all__ = [
    "VERIFY",
    "FULL",
    "Configuration",
    "SiteClassDefinition",
    "catalogue",
    "chemisorption_classes",
    "physisorption_classes",
]
