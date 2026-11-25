"""Source-disjoint, taxon-identity-matched island floral trait-pairing falsification.

One Direct species represents one fixed worldwide colour/form/symmetry tuple.
This is a *floristic repertoire* contrast, not an island-native evolutionary
adaptation or pollinator-guild shift.

Same-genus / same-candidate-family matching uses DIFFERENT accepted species
from island and candidate GIFT mainland source. Sampling is outcome blind.
The family candidate map MUST be independently strict CoL-XR + GBIF agreed.
Family matching is not equivalent to a dated phylogenetic tree.
"""
from __future__ import annotations
from collections import defaultdict
from hashlib import sha256
import numpy as np
import pandas as pd

from .regional_native_flora_fd import TRAITS, METRICS, exact_rarefied_flora_fd
from .regional_native_provenance_volume import cluster_equal_bootstrap

RANKS = ("genus", "family", "family_different_genus")
N_SPECIES = (5, 8)
REPEATS = 40
SUPPORT = (10, 5, 5)  # distinct islands, physical 10-degree blocks, GIFT sources


def stable_seed(*tokens: object) -> int:
    return int.from_bytes(sha256(":".join(map(str, tokens)).encode()).digest()[:8], "big")


def species_genus(name: str) -> str:
    words = name.strip().split()
    return words[0] if len(words) >= 2 and words[0] else ""


def _groups(names: set[str], *, rank: str,
            families: dict[str, str]) -> dict[str, tuple[str, ...]]:
    if rank not in RANKS:
        raise ValueError("Unfrozen taxon matching scope")
    groups = defaultdict(list)
    for name in sorted(names):
        key = species_genus(name) if rank == "genus" else families.get(name, "")
        if key:
            groups[key].append(name)
    return {k: tuple(v) for k, v in groups.items()}


def support_slots(local: set[str], donor: set[str], *,
                  rank: str, families: dict[str, str]) -> tuple[dict, dict, dict]:
    """Assess which taxon-ID matches are POSSIBLE, before trait scores."""
    if local & donor:
        raise ValueError("Source donor contains a focal island accepted plant name")
    a = _groups(local, rank=rank, families=families)
    b = _groups(donor, rank=rank, families=families)
    shared = {}
    for k in sorted(set(a) & set(b)):
        if rank == "family_different_genus":
            valid = {}
            for lhs in a[k]:
                possible = tuple(rhs for rhs in b[k]
                                 if species_genus(lhs) != species_genus(rhs))
                if possible:
                    valid[lhs] = possible
            if valid:
                shared[k] = valid
        else:
            shared[k] = True
    return a, b, shared


def _trait_score(names: list[str], reference: dict[str, tuple[str, str, str]],
                 n: int) -> np.ndarray:
    states = np.asarray([reference[x] for x in names], dtype=object)
    d = exact_rarefied_flora_fd(states, n=n)
    return np.asarray([d[k] for k in METRICS], dtype=float)


def sample_taxon_matched_pairing(
    local: set[str], donor: set[str], reference: dict[str, tuple[str,str,str]],
    *, rank: str, families: dict[str, str], n: int, repeats: int,
    seed: int,
) -> dict:
    """Means across repeated *source/target-matched*, one-species-per-group draws.

    Identity-free control draws the same N distinct genus/family IDs per SIDE,
    but lets the original IDs differ across sides. This is a separate target
    estimand on the EXACT same source-supported physical islands.
    """
    if n not in N_SPECIES or repeats < REPEATS:
        raise ValueError("Frozen n5/n8 or original within-island draw floor changed")
    if not (local | donor) <= set(reference):
        raise ValueError("Taxon input lacks original complete Direct triple states")
    a, b, common = support_slots(local, donor, rank=rank, families=families)
    if len(common) < n or len(a) < n or len(b) < n:
        raise ValueError("Insufficient source-matched botanical groups")
    rng = np.random.default_rng(seed)
    match, count_only = [], []
    keys = sorted(common)
    for _ in range(repeats):
        chosen = sorted(str(x) for x in rng.choice(keys, size=n, replace=False))
        left, right = [], []
        for k in chosen:
            if rank == "family_different_genus":
                lhs = str(rng.choice(sorted(common[k])))
                rhs = str(rng.choice(common[k][lhs]))
                if species_genus(lhs) == species_genus(rhs):
                    raise AssertionError("Different-genus family control failed")
            else:
                lhs = str(rng.choice(a[k]))
                rhs = str(rng.choice(b[k]))
            left.append(lhs); right.append(rhs)
        if len(set(left)) != n or len(set(right)) != n or set(left) & set(right):
            raise AssertionError("Matched plant comparison duplicated names")
        match.append(_trait_score(left, reference, n) -
                     _trait_score(right, reference, n))
        # Identity free: same count of distinct botanical rank labels per side.
        # Use all available local and donor genera/families; pre-select island
        # and donor groups without reading flower-phenotype outcomes.
        chosen_a = rng.choice(sorted(a), size=n, replace=False)
        chosen_b = rng.choice(sorted(b), size=n, replace=False)
        names_a = [str(rng.choice(a[str(x)])) for x in chosen_a]
        names_b = [str(rng.choice(b[str(x)])) for x in chosen_b]
        count_only.append(_trait_score(names_a, reference, n) -
                          _trait_score(names_b, reference, n))
    m = np.asarray(match).mean(axis=0)
    free = np.asarray(count_only).mean(axis=0)
    return {
        "matched": dict(zip(METRICS, map(float, m))),
        "same_taxon_richness_only": dict(zip(METRICS, map(float, free))),
        "n_identity_matchable_groups": len(common),
        "n_local_group_identities": len(a),
        "n_source_group_identities": len(b),
    }


def compare_native_source_taxon_pairing(
    direct: pd.DataFrame, classified: pd.DataFrame, geography: pd.DataFrame,
    near: dict[str,tuple[str,...]], candidate_flora: dict[str,set[str]],
    families: dict[str,str], *, repeats: int=REPEATS,
    seed: int=20261114, boots: int=1999
) -> dict:
    """Run protected original-data analysis, expose ONLY regional aggregates."""
    if len(direct) != 2234 or list(direct.columns) != list(TRAITS):
        raise ValueError("Original 2234 species-direct source morphology changed")
    if repeats < REPEATS or boots < 999:
        raise ValueError("Original matched taxon sensitivity resampling weakened")
    if not {"island_id","accepted_species","strict_source_native",
            "regional_native_compatible"}.issubset(classified):
        raise ValueError("Original WCVP flora provenance evidence incomplete")
    if classified.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island flora taxon duplicates")
    if not {"island_id","analysis_regime","spatial_block"}.issubset(geography):
        raise ValueError("Original geographic design missing")
    if geography.island_id.duplicated().any():
        raise ValueError("Original geographic island block duplicates")
    if len(families) < 1900 or any(x not in direct.index for x in families):
        raise ValueError("Strict candidate family source reference incomplete")
    reference = {
        str(name): tuple(map(str, values))
        for name, values in zip(direct.index, direct.to_numpy(dtype=str))
    }
    g = geography.set_index("island_id")
    ids = set(classified.island_id.unique())
    if len(ids) != 928 or set(near) != ids or not ids <= set(g.index):
        raise ValueError("Original WCVP/928-island source completeness gate changed")
    prepared = []
    outcomes = []
    for island, rows in classified.groupby("island_id", sort=True):
        island = str(island)
        all_names = set(rows.accepted_species)
        source_names = (set().union(*(candidate_flora.get(x, set()) for x in near[island]))
                        & set(reference)) - all_names
        if len(near[island])!=5 or len(set(near[island]))!=5:
            raise ValueError("Original source assignment not five independent candidate IDs")
        scope_names = {
            "all_recorded": all_names,
            "regional_WCVP_compatible": set(rows.loc[
                rows.regional_native_compatible, "accepted_species"]),
            "strict_source_native": set(rows.loc[
                rows.strict_source_native, "accepted_species"]),
        }
        for scope, names in scope_names.items():
            if len(names) < 30:
                continue
            for rank in RANKS:
                a, b, matched = support_slots(names, source_names, rank=rank,
                                              families=families)
                for n in N_SPECIES:
                    include = (len(matched) >= n and len(a) >= n and len(b) >= n)
                    private = dict(
                        island_id=island,
                        scope=scope,
                        regime=str(g.at[island,"analysis_regime"]),
                        spatial_block=str(g.at[island,"spatial_block"]),
                        source_signature=near[island][0],
                        matching_level=rank, n_each=int(n),
                        n_shared_ID_groups=len(matched),
                        source_comparable=bool(len(source_names)>=30),
                        matchable=bool(include),
                    )
                    prepared.append(private.copy())
                    if not include or len(source_names)<30:
                        continue
                    result = sample_taxon_matched_pairing(
                        names, source_names, reference, rank=rank,
                        families=families, n=n, repeats=repeats,
                        seed=stable_seed(seed,island,scope,rank,n))
                    for met in METRICS:
                        for k,prefix in (
                            ("matched", "matched_"),
                            ("same_taxon_richness_only","count_only_"),
                        ):
                            private[prefix+met] = result[k][met]
                        private["matched_minus_count_only_"+met] = (
                            private["matched_"+met] - private["count_only_"+met])
                    outcomes.append(private)
    support = pd.DataFrame(prepared)
    results = pd.DataFrame(outcomes)
    panels = []
    for (scope, regime, rank, n), sub in support.groupby(
            ["scope","regime","matching_level","n_each"],sort=True):
        usable = sub.loc[sub.matchable & sub.source_comparable]
        available = results.loc[
            results.scope.eq(scope) & results.regime.eq(regime) &
            results.matching_level.eq(rank) & results.n_each.eq(n)
        ]
        if len(usable) != len(available):
            raise AssertionError("Source-eligible matchable island score lost")
        n_blocks = usable.spatial_block.nunique()
        n_sources = usable.source_signature.nunique()
        enough = len(usable)>=SUPPORT[0] and n_blocks>=SUPPORT[1] and n_sources>=SUPPORT[2]
        rec = dict(scope=scope,regime=regime,
                   level=rank,n_matched_distinct_groups=n,
                   n_original_direct30_islands=len(sub),
                   n_taxon_matchable_source_comparable_islands=len(usable),
                   n_physical_10degree_blocks=int(n_blocks),
                   n_first_GIFT_source_signatures=int(n_sources),
                   n_matchable_taxon_groups_median=float(sub.n_shared_ID_groups.median()),
                   meets_pre_outcome_minimum_reporting_gate=bool(enough))
        if enough:
            for met in METRICS:
                for prefix in ("matched_","count_only_","matched_minus_count_only_"):
                    key=prefix+met
                    for weighting, cluster in (("block","spatial_block"),
                                               ("source","source_signature")):
                        rec[weighting+"_"+key] = cluster_equal_bootstrap(
                            available,key,cluster,draws=boots,
                            seed=stable_seed(seed,"boot",scope,regime,rank,n,key,weighting)
                        )
        panels.append(rec)
    if len(panels)<18:
        raise ValueError("Source WCVP 3 scope × native 4-region support ledger incomplete")
    return dict(
        status="ORIGINAL_ISLAND_NATIVE_TAXON_MATCHED_COLOUR_ARCHITECTURE_PAIRING",
        n_original_species_direct_phenotype_names=2234,
        n_original_island_floras=928,
        n_strict_source_native_eligible=87,
        n_WCVP_regionally_compatible_eligible=357,
        original_CoL_XR_exact_candidate_family_names=len(families),
        original_repeats_per_island=int(repeats),
        method=dict(
            n_matched_taxon_id_groups= list(N_SPECIES),
            matched_species_differ_between_island_and_native_GIFT_source=True,
            island_accepted_species_excluded_from_every_source=True,
            strict_original_species_direct_complete_states_only=True,
            all_original_island_species_used_for_source_exclusion=True,
            original_source_family_identifiers_are_unpromoted_CoLXR_candidates=True,
            WCVP_regionally_native_is_not_island_locally_native=True,
            opportunity_selected_before_flower_outcomes=True,
            equal_genus_or_family_richness_id_free_control=True,
            whole_spatial_blocks_and_GIFT_source_labels_are_not_evolutionary_replications=True,
            output_is_region_only_no_restricted_island_taxon_source_id=True,
            no_guild_decline_or_pollination_function_is_measured=True,
        ),
        regional_taxon_matched_panels=panels,
        raw_species_island_family_and_source_IDs_not_published=True,
    )
