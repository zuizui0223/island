"""Independent 10-trait direct-evidence opportunity for oceanic/continental endemics.

This is an ORIGINAL-source *support audit*, not a hypothesis test of adaptation.
Do not require colour+form+symmetry intersection before counting endemics.
Never infer endemicity from native origin; original endemic_status and
floristic_status must independently agree and be source-reviewed. WCVP
TDWG-L3 native-region corroboration alone is NOT an island-endemic label.

Potential GIFT nearest-mainland source congeners are NOT ancestors; absence
from a GIFT list is NOT failure to colonize nor evolutionary innovation.
Outputs contain only regional geological-type counts, no island/taxon IDs.
"""
from __future__ import annotations
from collections import defaultdict
import pandas as pd

TRAITS=(
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "self_incompatibility","autonomous_selfing_capacity",
    "mating_system","cleistogamy",
)
ORIGINS=("continental","oceanic")

def strict_original_direct_phenotypes(ledger:pd.DataFrame)->dict[str,set[str]]:
    needed={"accepted_species","trait_name","normalized_value","evidence_scope",
            "evidence_quality"}
    if not needed.issubset(ledger):
        raise ValueError("Original source Direct trait columns incomplete")
    df=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ].copy()
    if df.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Contradictory original Direct trait evidence")
    return {
        t:set(df.loc[df.trait_name.eq(t),"accepted_species"])
        for t in TRAITS
    }


def original_endemic_trait_opportunity(
    ledger:pd.DataFrame,raw_flora:pd.DataFrame,geology:pd.DataFrame,
    nearest:pd.DataFrame,source_flora:pd.DataFrame
)->dict:
    """Physical-island incidence, not duplicated 10-trait complete cases."""
    needed={"island_id","accepted_species","origin_status","endemic_status",
            "floristic_status"}
    if not needed.issubset(raw_flora):
        raise ValueError("Original independent endemic_status source missing")
    if raw_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island accepted plant repeats")
    if not raw_flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original plant native-status code unexpected")
    if not raw_flora.endemic_status.isin(("endemic","nonendemic","unresolved")).all():
        raise ValueError("Original endemic status code unexpected")
    if ((raw_flora.floristic_status.eq("endemic") &
         (~raw_flora.origin_status.eq("native") |
          ~raw_flora.endemic_status.eq("endemic"))).any() or
        (raw_flora.floristic_status.eq("native_nonendemic") &
         (~raw_flora.origin_status.eq("native") |
          ~raw_flora.endemic_status.eq("nonendemic"))).any()):
        raise ValueError("Native/endemic/floristic original source status conflict")
    if (not {"island_id","geological_origin","analysis_regime","spatial_block"}
        .issubset(geology) or geology.island_id.duplicated().any()):
        raise ValueError("Geologically qualified independent GSHHG island crosswalk missing")
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(nearest):
        raise ValueError("Original candidate nearest GIFT source assignments missing")
    if not {"entity_ID","work_species"}.issubset(source_flora):
        raise ValueError("Original GIFT source native checklist missing")
    sets=strict_original_direct_phenotypes(ledger)
    direct_global={t:len(names) for t,names in sets.items()}
    typed=geology.loc[geology.geological_origin.isin(ORIGINS),[
        "island_id","geological_origin","analysis_regime","spatial_block"
    ]]
    if typed.island_id.nunique()<4:
        raise ValueError("Independent geological type coverage not verified")
    names=raw_flora.loc[raw_flora.island_id.isin(typed.island_id)].merge(
        typed,on="island_id",how="inner",validate="many_to_one")
    if names.empty:
        raise ValueError("Geology has no original source flora")
    donor=nearest.loc[nearest.source_mode.eq("geo_k5")].copy()
    donor["source_rank"]=pd.to_numeric(donor.source_rank,errors="raise")
    donor=donor.loc[donor.source_rank.between(1,5)]
    if donor.duplicated(["island_id","source_rank"]).any():
        raise ValueError("GIFT nearest source rank repeated")
    five=donor.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    source_by=source_flora.drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()

    local_rows=[]
    for island,sub in names.groupby("island_id",sort=True):
        k=str(island)
        own=set(sub.accepted_species)
        candidates=five.get(k,())
        candidate_ok=(len(candidates)==5 and len(set(candidates))==5)
        nearby_all=(set().union(*(source_by.get(t,set()) for t in candidates))-own
                    if candidate_ok else set())
        for role,frame in (
            ("strict_endemic",sub.loc[sub.floristic_status.eq("endemic")]),
            ("strict_native_nonendemic",sub.loc[
                sub.floristic_status.eq("native_nonendemic")]),
        ):
            for trait in TRAITS:
                known=sorted(set(frame.accepted_species)&sets[trait])
                if candidate_ok:
                    same_genus_sources={w.split()[0] for w in
                        nearby_all&sets[trait] if len(w.split())>=2}
                    with_source=sum(1 for w in known
                        if len(w.split())>=2 and w.split()[0] in same_genus_sources)
                else:
                    with_source=None
                local_rows.append(dict(
                    geological_origin=sub.geological_origin.iloc[0],
                    regime=sub.analysis_regime.iloc[0],
                    spatial_block=sub.spatial_block.iloc[0],
                    role=role,trait=trait,
                    n_original_status_known_taxa=len(frame),
                    n_source_direct_phenotypes=len(known),
                    n_original_GIFT_geo_k5_source_assignments=5 if candidate_ok else 0,
                    n_with_direct_same_genus_GIFT_source=with_source,
                    n_candidates_with_frozen_donor=bool(candidate_ok),
                ))
    if not local_rows:
        raise ValueError("Original endemic/direct join empty")
    t=pd.DataFrame(local_rows)
    groups=[]
    for (origin,regime,trait,role),g in t.groupby(
        ["geological_origin","regime","trait","role"],sort=True):
        supported=g.loc[g.n_candidates_with_frozen_donor]
        n_endemic=int(g.n_original_status_known_taxa.sum())
        records=int(g.n_source_direct_phenotypes.sum())
        groups.append(dict(
            geological_origin=origin,regime=regime,trait=trait,
            original_status_tier=role,
            n_original_geologically_typed_islands=len(g),
            n_original_status_known_island_taxon_records=n_endemic,
            n_original_high_medium_Direct_island_taxon_trait_records=records,
            n_geologically_typed_islands_with_Direct_trait=int(
                g.n_source_direct_phenotypes.gt(0).sum()),
            n_islands_with_5_rank_GIFT_donors=int(len(supported)),
            n_Direct_taxon_rows_with_comparable_same_genus_source=(
                int(supported.n_with_direct_same_genus_GIFT_source.sum())
                if len(supported) else None),
            n_islands_with_Direct_source_congener=int(
                supported.n_with_direct_same_genus_GIFT_source.gt(0).sum()
                if len(supported) else 0),
            n_distinct_spatial_blocks_with_Direct=int(
                g.loc[g.n_source_direct_phenotypes.gt(0),"spatial_block"].nunique()),
            direct_status_coverage_fraction=(
                records/n_endemic if n_endemic>0 else None),
        ))
    return dict(
        status="INDEPENDENT_GEOLOGICAL_ENDEMIC_FLORAL_DIRECT_SINGLE_TRAIT_SUPPORT",
        n_original_geological_typed_islands=int(len(typed)),
        n_geologically_typed_islands_with_source_flora=int(names.island_id.nunique()),
        original_global_single_trait_Direct_species_counts=direct_global,
        trait_names=list(TRAITS),
        strict_original_status_counts_by_region_geology_trait=groups,
        measurement_only=True,
        no_single_taxon_names_island_source_IDs_or_study_coordinates_exported=True,
        no_actual_colonization_adaptive_selection_or_pollinator_loss_inferred=True,
    )
