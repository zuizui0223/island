from collections import defaultdict
import pytest

from fd.dopi_gift_guild_motif_bridge import SCOPES
from fd.dopi_study_genus_equal import genus_equal_dopi_audit


def index_case(*,n_studies=12,variation=False,which_scope=SCOPES[0]):
    units=defaultdict(set)
    positive=set()
    for i in range(1,n_studies+1):
        study=f"paper{i}"
        g1=[f"Alpha p{j}" for j in range(9)]
        g2=[f"Beta n{j}" for j in range(9)]
        positive.update(g1)
        if not variation:
            # Raw 90% vs 10% positive because the *same two
            # genera* are weighted differently by the two guilds.
            A=set(g1)|{g2[0]}
            B={g1[0]}|set(g2)
        else:
            positive.update({f"Beta p{j}" for j in range(2)})
            A={g1[0],g1[1],g1[2],g2[0],g2[1],g2[2]}
            B={g1[0],g1[8],g2[0],g2[6],g2[7],g2[8]}
            positive.discard("Alpha p8")
            positive.discard("Beta n0")
            # Variation is observed because Alpha p8 is negative
            # while Alpha p0 etc are positive, and Beta also varies.
        units[(study,"Bombus",which_scope)].update(A)
        units[(study,"Syrphidae_hoverfly",which_scope)].update(B)
    known=set().union(*units.values())
    return dict(
        indexed_units=units,no_true_absence_inferred=True,
        n_distinct_publication_keys=n_studies,
        original_source_row_counts=dict(source_rows=n_studies*20),
        n_original_island_accepted_plants_matching_DoPI=len(known),
    ),dict(known=known,joint=positive&known,n_axes=1)


def one_panel(report):
    return next(q for q in report["panels"] if q["guild_a"]=="Bombus" and
                q["guild_b"]=="Syrphidae_hoverfly" and
                q["floral_motif"]=="deep_tube" and
                q["evidence_scope"]==SCOPES[0])


def test_equal_genus_weights_remove_simpson_style_genus_composition_confounding():
    idx,motif=index_case()
    out=genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=299)
    q=one_panel(out)
    assert out["n_original_GIFT_botanical_islands"]==4505
    assert q["n_eligible_publications_with_two_shared_genera"]==12
    assert q["genus_equal_all"]["status"]=="DESCRIPTIVE_POSITIVE_VISITOR_STUDY_EQUAL_GENUS_FIXED"
    assert q["genus_equal_all"]["raw_study_diff"]==pytest.approx(.8)
    assert q["genus_equal_all"]["shared_genus_pooled_diff"]==pytest.approx(.8)
    assert q["genus_equal_all"]["genus_equal_diff"]==pytest.approx(0)
    assert q["genus_equal_all"]["raw_minus_genus_equal"]==pytest.approx(.8)
    assert q["genus_equal_informative"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_STUDIES"
    assert out["interpretation"]["does_not_infer_undocumented_guild_absence"]


def test_informative_within_genus_split_is_retained_if_both_actual_states():
    idx,motif=index_case(variation=True)
    out=genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=299)
    q=one_panel(out)
    assert q["n_publications_with_variable_motif_shared_genus"]==12
    assert q["genus_equal_informative"]["status"]=="DESCRIPTIVE_POSITIVE_VISITOR_STUDY_EQUAL_GENUS_FIXED"
    assert q["genus_equal_informative"]["n_studies_with_variable_genus"]==12


def test_single_study_repeated_records_do_not_create_fake_study_replications():
    idx,motif=index_case(n_studies=1)
    out=genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=199)
    assert one_panel(out)["genus_equal_all"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_STUDIES"


def test_missing_competitor_study_guild_not_control_negative():
    idx,motif=index_case()
    idx["indexed_units"]={
        k:v for k,v in idx["indexed_units"].items() if k[1]=="Bombus"
    }
    z=genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=199)
    q=one_panel(z)
    assert q["n_eligible_publications_with_two_shared_genera"]==0
    assert "genus_equal_diff" not in q["genus_equal_all"]
    assert q["genus_equal_all"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_STUDIES"


def test_no_motif_positive_studies_is_missing_support_not_zero_effect():
    idx,motif=index_case()
    motif["joint"]=set()
    out=genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=299)
    q=one_panel(out)
    assert q["genus_equal_all"]["status"]=="NOT_IDENTIFIABLE_POSITIVE_MOTIF_STUDY_SUPPORT"
    assert "genus_equal_diff" not in q["genus_equal_all"]


def test_source_protocol_cannot_pretend_true_guild_zeros():
    idx,motif=index_case()
    idx["no_true_absence_inferred"]=False
    with pytest.raises(ValueError,match="positive source-only"):
        genus_equal_dopi_audit(idx,motifs={"deep_tube":motif},draws=199)
