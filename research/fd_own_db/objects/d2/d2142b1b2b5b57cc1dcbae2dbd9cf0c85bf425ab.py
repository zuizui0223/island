"""Island x bird-species incidence from expert range maps, and nectarivore diet classes.

Implements the incidence rules frozen in docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md.
Geometries must already be in an equal-area projection; the runner script handles reprojection.
Only shapely, numpy, pandas and scipy are needed here.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse

PRESENCE_KEEP = frozenset({1, 2})  # extant, probably extant
ORIGIN_KEEP = frozenset({1, 2})  # native, reintroduced
SEASONAL_PRIMARY = frozenset({1, 2})  # resident, breeding
SEASONAL_WITH_NONBREEDING = frozenset({1, 2, 3, 4})  # sensitivity A

SPECIALIST_CUT = 50.0
OTHER, PARTIAL, SPECIALIST = 0, 1, 2
DIET_CLASS_NAMES = ("other", "partial", "specialist")


def normalise_binomial(name: str) -> str:
    """Lower-case 'genus species' with single spaces; drop authorities, subspecies and punctuation."""
    if not isinstance(name, str):
        raise TypeError(f"species name must be a string, got {type(name).__name__}")
    words = re.sub(r"[^a-z\- ]", " ", name.replace("_", " ").lower()).split()
    if len(words) < 2:
        raise ValueError(f"not a binomial: {name!r}")
    return f"{words[0]} {words[1]}"


def filter_range_records(
    records: pd.DataFrame, seasonal: frozenset[int] = SEASONAL_PRIMARY
) -> pd.DataFrame:
    """Keep range polygons by BirdLife PRESENCE / ORIGIN / SEASONAL codes."""
    missing = {"sci_name", "presence", "origin", "seasonal"} - set(records.columns)
    if missing:
        raise KeyError(f"range records missing columns: {sorted(missing)}")
    keep = (
        records["presence"].isin(PRESENCE_KEEP)
        & records["origin"].isin(ORIGIN_KEEP)
        & records["seasonal"].isin(seasonal)
    )
    return records.loc[keep].copy()


@dataclass(frozen=True)
class Incidence:
    """Island x species incidence. ``matrix[i, j]`` is True when species j occurs on island i."""

    island_ids: tuple
    species: tuple[str, ...]
    matrix: sparse.csr_matrix
    cover: sparse.csr_matrix  # covered fraction of island area, for every intersecting pair

    def at_threshold(self, min_cover: float) -> "Incidence":
        c = self.cover.tocoo()
        if min_cover <= 0:
            keep = c.data > 0
        else:
            keep = c.data >= min_cover - 1e-12
        m = sparse.csr_matrix(
            (np.ones(int(keep.sum()), dtype=bool), (c.row[keep], c.col[keep])), shape=c.shape
        )
        return Incidence(self.island_ids, self.species, m, self.cover)


def island_incidence(
    island_ids,
    island_geoms,
    range_species,
    range_geoms,
    min_cover: float = 0.5,
    island_buffer: float = 0.0,
) -> Incidence:
    """Fraction of each island covered by each species' retained range, and incidence at ``min_cover``.

    ``range_species`` and ``range_geoms`` are parallel; several polygons per species are unioned first.
    Areas are planar, so pass equal-area projected geometries (metres). ``island_buffer`` (same units)
    enlarges islands before intersection for the buffered sensitivity; the denominator stays the
    unbuffered island area.
    """
    from shapely import STRtree, area, buffer, intersection, union_all

    island_ids = tuple(island_ids)
    islands = np.asarray(island_geoms, dtype=object)
    if len(islands) != len(island_ids):
        raise ValueError("island_ids and island_geoms differ in length")
    if len(set(island_ids)) != len(island_ids):
        raise ValueError("duplicate island ids")
    island_area = area(islands)
    if np.any(~np.isfinite(island_area)) or np.any(island_area <= 0):
        raise ValueError("every island needs a positive finite area")
    probe = buffer(islands, island_buffer) if island_buffer > 0 else islands

    by_species: dict[str, list] = {}
    for s, g in zip(range_species, range_geoms, strict=True):
        by_species.setdefault(normalise_binomial(s), []).append(g)
    species = tuple(sorted(by_species))
    tree = STRtree(probe)
    rows, cols, vals = [], [], []
    for j, s in enumerate(species):
        rng = union_all(by_species[s])
        hits = tree.query(rng, predicate="intersects")
        if len(hits) == 0:
            continue
        frac = area(intersection(probe[hits], rng)) / island_area[hits]
        frac = np.minimum(frac, 1.0)
        ok = frac > 0
        rows.extend(hits[ok].tolist())
        cols.extend([j] * int(ok.sum()))
        vals.extend(frac[ok].tolist())
    cover = sparse.csr_matrix(
        (np.asarray(vals, float), (np.asarray(rows, int), np.asarray(cols, int))),
        shape=(len(island_ids), len(species)),
    )
    return Incidence(island_ids, species, sparse.csr_matrix(cover.shape, dtype=bool), cover).at_threshold(
        min_cover
    )


def island_incidence_streaming(island_ids, island_geoms, chunks, min_cover: float = 0.5, island_buffer: float = 0.0):
    """Same result as :func:`island_incidence`, reading range polygons chunk by chunk.

    ``chunks`` yields ``(species_names, geometries)``. Only the island-sized intersection pieces are
    kept per island x species pair, then unioned, so polygons of one species split across chunks are
    never double counted and the full range layer never sits in memory.
    """
    from shapely import STRtree, area, buffer, intersection, union_all

    island_ids = tuple(island_ids)
    islands = np.asarray(island_geoms, dtype=object)
    if len(set(island_ids)) != len(island_ids):
        raise ValueError("duplicate island ids")
    island_area = area(islands)
    if np.any(~np.isfinite(island_area)) or np.any(island_area <= 0):
        raise ValueError("every island needs a positive finite area")
    probe = buffer(islands, island_buffer) if island_buffer > 0 else islands
    tree = STRtree(probe)
    pieces: dict[tuple[int, str], list] = {}
    for names, geoms in chunks:
        geoms = np.asarray(geoms, dtype=object)
        for s, g in zip(names, geoms, strict=True):
            hits = tree.query(g, predicate="intersects")
            if len(hits) == 0:
                continue
            key = normalise_binomial(s)
            for h, piece in zip(hits, intersection(probe[hits], g)):
                if not piece.is_empty:
                    pieces.setdefault((int(h), key), []).append(piece)
    species = tuple(sorted({k for _, k in pieces}))
    col = {s: j for j, s in enumerate(species)}
    rows, cols, vals = [], [], []
    for (h, s), parts in pieces.items():
        frac = min(float(area(union_all(parts))) / island_area[h], 1.0)
        if frac > 0:
            rows.append(h)
            cols.append(col[s])
            vals.append(frac)
    cover = sparse.csr_matrix((vals, (rows, cols)), shape=(len(island_ids), len(species)))
    return Incidence(island_ids, species, sparse.csr_matrix(cover.shape, dtype=bool), cover).at_threshold(
        min_cover
    )


# Nectarivore clades for the clade-replication secondary analysis.
DREPANIDINI_GENERA = frozenset(
    {
        "akialoa", "akohekohe", "chlorodrepanis", "ciridops", "drepanis", "dysmorodrepanis", "hemignathus",
        "himatione", "loxops", "magumma", "melamprosops", "oreomystis", "palmeria", "paroreomyza",
        "pseudonestor", "viridonia", "vestiaria",
    }
)
CLADES = ("trochilidae", "meliphagidae", "nectariniidae", "drepanidini", "other")


def nectarivore_clade(family: str, genus: str) -> int:
    """Index into :data:`CLADES`."""
    f, g = family.strip().lower(), genus.strip().lower()
    if g in DREPANIDINI_GENERA:
        return CLADES.index("drepanidini")
    for c in ("trochilidae", "meliphagidae", "nectariniidae"):
        if f == c:
            return CLADES.index(c)
    return CLADES.index("other")


@dataclass(frozen=True)
class Crosswalk:
    mapping: dict[str, str]  # range-map name -> trait-table name
    unmatched: tuple[str, ...]
    via_synonym: int


def crosswalk_names(range_names, trait_names, synonyms: pd.DataFrame | None = None) -> Crosswalk:
    """Exact normalised match, then an explicit synonym table with columns ``from,to``. No fuzzy matching."""
    trait = {normalise_binomial(t) for t in trait_names}
    syn: dict[str, str] = {}
    if synonyms is not None and len(synonyms):
        if not {"from", "to"} <= set(synonyms.columns):
            raise KeyError("synonym table needs 'from' and 'to' columns")
        for a, b in zip(synonyms["from"], synonyms["to"]):
            a, b = normalise_binomial(a), normalise_binomial(b)
            if a in syn and syn[a] != b:
                raise ValueError(f"conflicting synonyms for {a}")
            syn[a] = b
    mapping, unmatched, n_syn = {}, [], 0
    for r in sorted({normalise_binomial(n) for n in range_names}):
        if r in trait:
            mapping[r] = r
        elif r in syn and syn[r] in trait:
            mapping[r] = syn[r]
            n_syn += 1
        else:
            unmatched.append(r)
    return Crosswalk(mapping, tuple(unmatched), n_syn)


def diet_class(diet_nect, specialist_cut: float = SPECIALIST_CUT) -> np.ndarray:
    """0 other (no nectar), 1 partial (0 < % < cut), 2 specialist (% >= cut)."""
    d = np.asarray(diet_nect, float)
    if np.any(~np.isfinite(d)) or np.any((d < 0) | (d > 100)):
        raise ValueError("Diet-Nect must be finite percentages in [0, 100]")
    return np.where(d >= specialist_cut, SPECIALIST, np.where(d > 0, PARTIAL, OTHER)).astype(int)


@dataclass(frozen=True)
class SpeciesTable:
    """Per-species labels aligned with ``Incidence.species`` after dropping unmatched names."""

    species: tuple[str, ...]
    diet: np.ndarray  # diet class codes
    genus: np.ndarray  # integer block ids
    family: np.ndarray  # integer block ids


def build_species_table(
    inc: Incidence,
    traits: pd.DataFrame,
    crosswalk: Crosswalk,
    max_unmatched_share: float = 0.05,
    specialist_cut: float = SPECIALIST_CUT,
    island_mask=None,
) -> tuple[Incidence, SpeciesTable, dict]:
    """Attach diet/genus/family to incidence species; fail closed on too many unmatched occurring species.

    ``traits`` needs columns ``scientific``, ``diet_nect``, ``family``. ``island_mask`` (bool per island)
    restricts the unmatched-share check to qualifying islands.
    """
    need = {"scientific", "diet_nect", "family"}
    if need - set(traits.columns):
        raise KeyError(f"trait table missing {sorted(need - set(traits.columns))}")
    t = traits.assign(key=traits["scientific"].map(normalise_binomial))
    if t["key"].duplicated().any():
        dup = t.loc[t["key"].duplicated(), "key"].iloc[0]
        raise ValueError(f"duplicate trait species after normalisation: {dup}")
    t = t.set_index("key")

    m = inc.matrix
    if island_mask is not None:
        m = m[np.asarray(island_mask, bool)]
    occurring = np.asarray(m.sum(axis=0)).ravel() > 0
    occ_names = [s for s, o in zip(inc.species, occurring) if o]
    unmatched_occ = [s for s in occ_names if s not in crosswalk.mapping]
    share = len(unmatched_occ) / max(len(occ_names), 1)
    report = {
        "species_in_ranges": len(inc.species),
        "species_occurring_on_qualifying_islands": len(occ_names),
        "unmatched_occurring": len(unmatched_occ),
        "unmatched_occurring_share": share,
        "matched_via_synonym": crosswalk.via_synonym,
    }
    if share > max_unmatched_share:
        raise ValueError(
            f"{share:.1%} of occurring species unmatched (limit {max_unmatched_share:.0%}); "
            "extend the synonym table rather than dropping them"
        )
    keep = [j for j, s in enumerate(inc.species) if s in crosswalk.mapping]
    names = tuple(inc.species[j] for j in keep)
    rows = t.loc[[crosswalk.mapping[s] for s in names]]
    genus_names = [s.split(" ")[0] for s in names]
    g_codes = {g: k for k, g in enumerate(sorted(set(genus_names)))}
    f_codes = {f: k for k, f in enumerate(sorted(set(rows["family"].astype(str))))}
    table = SpeciesTable(
        species=names,
        diet=diet_class(rows["diet_nect"].to_numpy(), specialist_cut),
        genus=np.array([g_codes[g] for g in genus_names], int),
        family=np.array([f_codes[f] for f in rows["family"].astype(str)], int),
    )
    sub = Incidence(inc.island_ids, names, inc.matrix[:, keep].tocsr(), inc.cover[:, keep].tocsr())
    report["species_kept"] = len(names)
    report["diet_class_counts"] = {
        n: int((table.diet == k).sum()) for k, n in enumerate(DIET_CLASS_NAMES)
    }
    return sub, table, report


def diet_counts(matrix: sparse.spmatrix, diet: np.ndarray) -> np.ndarray:
    """Island x 3 array of (other, partial, specialist) species counts."""
    onehot = sparse.csr_matrix(
        (np.ones(len(diet)), (np.arange(len(diet)), np.asarray(diet, int))), shape=(len(diet), 3)
    )
    return np.asarray((sparse.csr_matrix(matrix, dtype=float) @ onehot).todense())
