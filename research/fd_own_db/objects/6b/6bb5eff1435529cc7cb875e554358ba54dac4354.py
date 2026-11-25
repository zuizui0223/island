import pytest
from fd.cross_source_syndrome_concordance import (
    audit,exact_profile_label_benchmark,STRUCTURAL,COLOUR_FORM,
    DIRECT_SHA,FLORA_SHA,DOPI_SHA
)

A=[.109885844,-.013679663,.051724999,.044833249,.00990603,
   -.01397980,.03714125,.004126772]
B=[.02083574,-.03073092,.03312519,.00844642,.00140116,
   .01435465,-.04281590,-.00184126]


def fixtures():
    g=dict(n_full_original_botanical_islands=4505,
           original_Direct_species_trait_ledger_sha256=DIRECT_SHA,
           original_botanical_flora_sha256=FLORA_SHA,
           n_real_two_trait_pairs=11,panels=[])
    d=dict(n_original_GIFT_botanical_islands=4505,
           original_DoPI_v2_csv_sha256=DOPI_SHA,
           n_original_DOPI_source_rows=676900,n_motif_profiles=27,
           panels=[])
    for i,m in enumerate(STRUCTURAL+COLOUR_FORM):
        g["panels"].append(dict(region="tropical",plant_scope="all_recorded",
           documented_joint_sector=m,
           all__spatial_block=dict(status="DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY",
               n_islands=100,
               delta_geo5_pair_prevalence_block_equal=A[i],
               delta_geo5_pair_prevalence_bootstrap95=[A[i]-.001,A[i]+.001],
               delta_climate10_pair_prevalence_block_equal=A[i],
               delta_climate10_pair_prevalence_bootstrap95=[A[i]-.001,A[i]+.001])))
        d["panels"].append(dict(motif=m,guild_a="butterfly_family",
           guild_b="Bombus",evidence_scope="positive_visitor_or_pollen_123",
           n_within_publication_trait_scored_contrasts=60,
           n_studies_with_positive_motif_either=12,
           all_study=dict(status="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION",
               mean_guildA_minus_guildB=B[i],
               bootstrap95=[B[i]-.001,B[i]+.001]),
           shared_genus_study=dict(status="NOT_IDENTIFIABLE_MOTIF_POSITIVE_STUDY_SUPPORT")))
    return g,d


def test_source_audited_8_motif_results_not_picked_three_only():
    r=audit(*fixtures())
    q=next(x for x in r["summaries"] if x["region"]=="tropical"
           and x["scope"]=="all_recorded" and x["mode"]=="geo5")
    assert r["n_global_GIFT_botanical_island_floras"]==4505
    assert q["n_jointly_supported_motifs"]==8
    assert (q["structural_sign_matches"],q["structural_total"])==(5,5)
    assert (q["colour_form_sign_matches"],q["colour_form_total"])==(0,3)
    assert q["descriptive_profile_label_shuffle"]["sign_tail_fraction"]==pytest.approx(18/28)
    assert q["descriptive_profile_label_shuffle"]["abs_rank_tail_fraction"]==pytest.approx(.24305555555555555)
    assert r["limitations"]["overlapping_motifs_not_independent"]
    assert r["limitations"]["no_pollinator_shift_established"]


def test_wrong_flora_or_external_doi_source_fails_closed():
    g,d=fixtures()
    g["original_botanical_flora_sha256"]="wrong"
    with pytest.raises(ValueError,match="GIFT"):audit(g,d)
    g,d=fixtures()
    d["n_original_DOPI_source_rows"]=1
    with pytest.raises(ValueError,match="DoPI"):audit(g,d)


def test_exact_shuffled_labels_are_explicitly_not_ecological_significance():
    z=exact_profile_label_benchmark(A,B)
    assert z["matched_signs"]==5
    assert z["n_exact_label_permutations"]==40320
    assert z["overlapping_profiles_nonexchangeable_not_ecological_p"]
    with pytest.raises(ValueError):
        exact_profile_label_benchmark(A[:7],B[:7])


def test_unsupported_external_motif_is_not_silent_negative():
    g,d=fixtures()
    d["panels"][-1]["all_study"]={"status":"NOT_IDENTIFIABLE_MOTIF_POSITIVE_STUDY_SUPPORT"}
    z=audit(g,d)
    q=next(x for x in z["summaries"] if x["region"]=="tropical"
           and x["scope"]=="all_recorded" and x["mode"]=="geo5")
    assert q["n_jointly_supported_motifs"]==7
    assert q["colour_form_total"]==2
    assert "descriptive_profile_label_shuffle" not in q
