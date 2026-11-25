"""Study-conditioned flower motif vs REAL European guild visitor/pollen records.

Primary botanical target: union of accepted plant species occurring in the
original 4,505 GIFT island floras. DoPI v2 is an EXTERNAL EUROPEAN visitor
literature source, NOT island-local visits or 4,505 independent guild surveys.

Each endpoint compares only plants POSITIVELY recorded with both guilds in
the SAME original study. It never interprets unreported guilds as absent.
Pollination grade 1 = confirmed; 2 = pollen carriage; 3 = visit only;
4 = nonpollinating interaction (excluded), unknown = not imputed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
import numpy as np

GUILDS=("Bombus","other_bee","butterfly_family",
        "Sphingidae_hawkmoth","Syrphidae_hoverfly")
PAIRS=(
    ("Sphingidae_hawkmoth","Bombus"),
    ("butterfly_family","Bombus"),
    ("Sphingidae_hawkmoth","butterfly_family"),
    ("Bombus","Syrphidae_hoverfly"),
    ("other_bee","Bombus"),
)
GRADES={
    "1":"confirmed",
    "2":"pollen_carriage",
    "3":"visit_only",
    "4":"nonpollinating_or_robbery",
}
SCOPES=("positive_visitor_or_pollen_123","pollen_evidence_12",
        "confirmed_pollination_1","direct_all_both_123")
PRIMARY=(
    ("Sphingidae_hawkmoth","Bombus","white_radial_tubular_deep"),
    ("butterfly_family","Bombus","yellow_orange_radial_deep"),
    ("Bombus","Syrphidae_hoverfly","blue_purple_zygomorphic_deep"),
    ("Bombus","Syrphidae_hoverfly","open_radial_shallow"),
)
MIN_STUDIES=10
MIN_DIRECT_PLANTS_EACH_GUILD=2


def _s(value)->str:
    return str(value or "").strip()


def botanical_genus(species:str)->str:
    bits=species.split()
    return bits[0] if len(bits)>1 else "_unknown_"+species


def _study(row:dict)->str:
    doi=_s(row.get("doi")).lower()
    # DOI may be an ordinary URL as per DoPI documentation, still a source key.
    if doi and doi not in {"na","n/a","unknown","none"}:
        return doi
    pub=_s(row.get("publication")).lower()
    if pub and pub not in {"na","n/a","unknown","none"}:
        return pub
    return ""


def screen_positive_dopi_rows(rows, *,columns:set[str],
                              island_plant_names:set[str],
                              motifs:dict, visitor_classifier)->dict:
    """Stream 676,900 source rows, retaining only positive study-plant-guild units.

    Returned 'indexed_units' contains private original accepted plant names
    and must NEVER be exported. Only summarize_study_guild_motifs output is
    eligible for an anonymous artifact.
    """
    need={"plant_species","pollinator_species","pollinator_family",
          "pollinator_order","doi","publication","pollination",
          "methodology","pollinator_survey","plant_survey"}
    if missing:=need-columns:
        raise ValueError("missing DoPI v2 original study/evidence taxonomy columns "+str(sorted(missing)))
    if not island_plant_names or not motifs:
        raise ValueError("original GIFT botanically accepted plant names and motifs missing")
    known_union=set().union(*(set(m["known"]) for m in motifs.values()))
    n=Counter()
    units=defaultdict(set)
    matching_all=set()
    grade_counts=Counter()
    source_studies=set()
    for row in rows:
        n["source_rows"]+=1
        sp=_s(row.get("plant_species"))
        if sp not in island_plant_names:
            continue
        n["rows_original_island_accepted_plant_match"]+=1
        matching_all.add(sp)
        if sp not in known_union:
            n["rows_original_island_plant_not_multiaxis_direct"]+=1
            continue
        n["rows_accepted_island_multiaxis_direct_plant_match"]+=1
        study=_study(row)
        if not study:
            n["rows_missing_study_identity"]+=1
            continue
        guild=visitor_classifier(row.get("pollinator_order"),
                                 row.get("pollinator_family"),
                                 row.get("pollinator_species"))
        if guild not in GUILDS:
            continue
        grade=_s(row.get("pollination"))
        grade_counts[grade]+=1
        if grade not in {"1","2","3"}:
            n["unsupported_evidence_or_nonpollinating_rows"]+=1
            continue
        source_studies.add(study)
        n["eligible_multiaxis_direct_positive_rows"]+=1
        units[(study,guild,SCOPES[0])].add(sp)
        if grade in {"1","2"}:
            units[(study,guild,SCOPES[1])].add(sp)
        if grade=="1":
            units[(study,guild,SCOPES[2])].add(sp)
        if (_s(row.get("methodology")).lower()=="do" and
            _s(row.get("pollinator_survey")).lower()=="all" and
            _s(row.get("plant_survey")).lower()=="all"):
            units[(study,guild,SCOPES[3])].add(sp)
    return dict(indexed_units=units,
                original_source_row_counts=dict(n),
                pollination_evidence_source_codes=dict(grade_counts),
                n_distinct_publication_keys=len(source_studies),
                n_original_island_accepted_plants_matching_DoPI=len(matching_all),
                n_source_study_guild_scope_cells=len(units),
                no_true_absence_inferred=True)


def _aggregate_studies(values: list[float], *,
                       seed: str,draws:int,
                       n_motif_hits_a:int,n_motif_hits_b:int,
                       n_motif_hits_either:int)->dict:
    v=np.asarray(values,float)
    if len(v)<MIN_STUDIES:
        return {"status":"NOT_IDENTIFIABLE_INDEPENDENT_STUDIES"}
    # All-zero combinations must NEVER masquerade as a precise zero
    # preference/syndrome result merely because known trait records exist.
    # A true exclusive guild motif may have zero among the competing guild.
    # Require observations of the motif in >=10 independent studies on
    # EITHER side; never require both guilds to contain that motif.
    if n_motif_hits_either<MIN_STUDIES:
        return {"status":"NOT_IDENTIFIABLE_MOTIF_POSITIVE_STUDY_SUPPORT"}
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    boot=v[rng.integers(len(v),size=(draws,len(v)))].mean(axis=1)
    return dict(status="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION",
                mean_guildA_minus_guildB=float(v.mean()),
                bootstrap95=[float(x) for x in np.quantile(boot,[.025,.975])])


def summarize_study_guild_motifs(index:dict, *,motifs:dict,draws:int=999)->dict:
    if draws<199:raise ValueError("whole-study bootstrap too short")
    units=index["indexed_units"]
    # Never export the study identifiers, plant species or exact interactions.
    studies=sorted({key[0] for key in units})
    panels=[]
    for scope in SCOPES:
        for ga,gb in PAIRS:
            paired_studies=[
                study for study in studies
                if (study,ga,scope) in units and (study,gb,scope) in units
            ]
            for motif_name,motif in motifs.items():
                known=set(motif["known"])
                positive=set(motif["joint"])
                values=[]
                shared_values=[]
                hits_a=hits_b=hits_either=0
                shared_hits_a=shared_hits_b=shared_hits_either=0
                unique_scored_a=set()
                unique_scored_b=set()
                n_common_species_studies=0
                for study in paired_studies:
                    a=units[(study,ga,scope)]&known
                    b=units[(study,gb,scope)]&known
                    if len(a)<MIN_DIRECT_PLANTS_EACH_GUILD or len(b)<MIN_DIRECT_PLANTS_EACH_GUILD:
                        continue
                    motif_a=bool(a&positive)
                    motif_b=bool(b&positive)
                    hits_a+=int(motif_a)
                    hits_b+=int(motif_b)
                    hits_either+=int(motif_a or motif_b)
                    values.append(len(a&positive)/len(a)-len(b&positive)/len(b))
                    unique_scored_a|=a
                    unique_scored_b|=b
                    if a&b:n_common_species_studies+=1
                    shared_genera={botanical_genus(x) for x in a}&{
                        botanical_genus(x) for x in b
                    }
                    if len(shared_genera)>=2:
                        aa={x for x in a if botanical_genus(x) in shared_genera}
                        bb={x for x in b if botanical_genus(x) in shared_genera}
                        if len(aa)>=2 and len(bb)>=2:
                            smotif_a=bool(aa&positive)
                            smotif_b=bool(bb&positive)
                            shared_hits_a+=int(smotif_a)
                            shared_hits_b+=int(smotif_b)
                            shared_hits_either+=int(smotif_a or smotif_b)
                            shared_values.append(
                                len(aa&positive)/len(aa)-len(bb&positive)/len(bb)
                            )
                cell=dict(guild_a=ga,guild_b=gb,motif=motif_name,
                          n_recorded_flower_axes=motif.get("axis_count",
                              motif.get("n_axes",len(motif.get("trait_axes",())))),
                          evidence_scope=scope,
                          predeclared_hypothesis=(ga,gb,motif_name) in PRIMARY,
                          n_publications_with_both_positive_guilds=len(paired_studies),
                          n_within_publication_trait_scored_contrasts=len(values),
                          n_publications_with_any_shared_accepted_plant=n_common_species_studies,
                          n_distinct_trait_scored_positive_plants_guild_a=len(unique_scored_a),
                          n_distinct_trait_scored_positive_plants_guild_b=len(unique_scored_b),
                          n_shared_genus_conditioned_publications=len(shared_values),
                          n_studies_with_positive_motif_guild_a=hits_a,
                          n_studies_with_positive_motif_guild_b=hits_b,
                          n_studies_with_positive_motif_either=hits_either,
                          n_shared_genus_studies_with_motif_guild_a=shared_hits_a,
                          n_shared_genus_studies_with_motif_guild_b=shared_hits_b,
                          n_shared_genus_studies_with_motif_either=shared_hits_either)
                cell["all_study"]= _aggregate_studies(values,
                    seed=f"{scope}|{ga}|{gb}|{motif_name}|all",draws=draws,
                    n_motif_hits_a=hits_a,n_motif_hits_b=hits_b,
                    n_motif_hits_either=hits_either)
                cell["shared_genus_study"]= _aggregate_studies(shared_values,
                    seed=f"{scope}|{ga}|{gb}|{motif_name}|genus",draws=draws,
                    n_motif_hits_a=shared_hits_a,n_motif_hits_b=shared_hits_b,
                    n_motif_hits_either=shared_hits_either)
                panels.append(cell)
    return dict(status="EXTERNAL_DOPI_V2_GIFT_REAL_MOTIF_WITHIN_PUBLICATION_GUILD_SCREEN",
                n_original_GIFT_botanical_islands=4505,
                n_GIFT_accepted_island_plant_species_matching_DoPI=index[
                    "n_original_island_accepted_plants_matching_DoPI"],
                n_original_DOPI_source_rows=index[
                    "original_source_row_counts"].get("source_rows",0),
                n_original_island_multiaxis_direct_positive_DOPI_rows=index[
                    "original_source_row_counts"].get("eligible_multiaxis_direct_positive_rows",0),
                n_unique_original_publication_keys=index["n_distinct_publication_keys"],
                n_indexed_study_guild_scope_cells=index["n_source_study_guild_scope_cells"],
                source_pollination_grade_counts=index[
                    "pollination_evidence_source_codes"],
                source_join_coverage=index["original_source_row_counts"],
                n_visitor_guild_pairs=len(PAIRS),
                n_motif_profiles=len(motifs),
                n_motifs_by_number_of_observed_axes=dict(sorted(
                    Counter(p.get("axis_count",p.get("n_axes",
                      len(p.get("trait_axes",())))) for p in motifs.values()).items()
                )),
                n_qualification_panels=len(panels),
                panels=panels,
                limits=dict(
                    original_all_4505_islands_are_botanical_population_not_guild_survey=True,
                    GIFT_morphology_species_wide_not_island_local=True,
                    no_true_guild_absence_inferred=True,
                    nonpollinating_grade4_excluded=True,
                    grade3_visitation_only_not_effective_pollination=True,
                    source_publications_not_proven_independent_of_GloBI=True,
                    same_study_positive_guilds_not_preference_or_fitness=True,
                    genus_conditioning_requires_two_shared_documented_genera=True,
                    minimum_positive_motif_source_studies_enforced=True,
                    all_zero_motif_per_guild_not_reported_as_no_association=True,
                    no_named_bird_pollination_calibration_from_European_DoPI=True,
                    no_syndrome_label_created_from_plant_morphology=True,
                    no_exact_accepted_plant_name_or_study_identity_exported=True))
