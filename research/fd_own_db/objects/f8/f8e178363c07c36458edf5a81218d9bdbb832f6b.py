"""Synthetic-only congeneric source comparisons: no private original taxa."""
import numpy as np
import pandas as pd
import pytest

from fd.source_genus_same_island import (
    MORPHOLOGY, SCOPES, compare_island_source_genera,
    summarize_source_genus_frame,
)


def synthetic():
    full = {"Alpha island1", "Alpha island2", "Alpha island3",
            "Beta island1", "Beta island2", "Gamma island1"}
    g = {f"Alpha mainland{j}" for j in range(12)}
    g |= {f"Beta mainland{j}" for j in range(12)}
    c = {f"Alpha climate{j}" for j in range(14)}
    c |= {f"Beta climate{j}" for j in range(14)}
    states = {x: int(x.startswith("Alpha")) for x in full | g | c}
    # Retain known different accepted species; swap source genus traits to
    # test genuinely different reference vs identical island target names.
    for x in g:
        states[x] = 0
    for x in c:
        states[x] = int(x.startswith("Alpha"))
    return full, g, c, states


def test_shared_genus_source_donors_do_not_reuse_original_island_taxa():
    full, geo, clim, state = synthetic()
    result = compare_island_source_genera(
        full_island_names=full, scope_names=full, geo5_candidates=geo,
        climate10_candidates=clim, species_states=state
    )
    x = result["shared_genus_one_donor_each"]
    assert x["n_matched_island_direct"] == 5
    assert x["n_matched_genera"] == 2
    assert x["n_geo5_shared_genus_candidate_direct"] == 24
    assert x["n_climate10_shared_genus_candidate_direct"] == 28
    assert np.isclose(x["observed_island_matched_fraction"], .6)
    assert np.isclose(x["predicted_geo5_from_distinct_source_congeners"], 0)
    assert np.isclose(x["predicted_climate10_from_distinct_source_congeners"], .6)
    assert x["n_matched_island_total"] == 5
    assert result["n_scope_direct"] == 6  # Gamma lacks a mainland congener
    assert result["shared_genus_two_donors_each"]["n_matched_island_direct"] == 5


def test_island_accepted_names_never_used_as_own_sources():
    full, geo, clim, state = synthetic()
    with pytest.raises(ValueError, match="leakage"):
        compare_island_source_genera(
            full_island_names=full, scope_names=full,
            geo5_candidates=geo | {"Alpha island1"},
            climate10_candidates=clim, species_states=state
        )


def test_fixed_subset_must_be_actual_original_island_membership():
    full, geo, clim, state = synthetic()
    with pytest.raises(ValueError, match="scope"):
        compare_island_source_genera(
            full_island_names=full, scope_names=full | {"Alpha fake"},
            geo5_candidates=geo, climate10_candidates=clim,
            species_states=state
        )


def test_no_congener_donor_does_not_mean_no_phenotype():
    full, geo, clim, state = synthetic()
    x=compare_island_source_genera(
        full_island_names=full, scope_names=full,
        geo5_candidates=geo, climate10_candidates=set(),
        species_states=state
    )["shared_genus_one_donor_each"]
    assert x["n_matched_island_direct"] == 0
    assert np.isnan(x["observed_island_matched_fraction"])
    assert np.isnan(x["predicted_geo5_from_distinct_source_congeners"])


def example_frame(n=40, coverage_gap=False):
    full, geo, clim, state = synthetic()
    rows = []
    for i in range(n):
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                result=compare_island_source_genera(
                    full_island_names=full, scope_names=full,
                    geo5_candidates=geo, climate10_candidates=clim,
                    species_states=state
                )
                z=dict(island_id=f"island{i}", regime="tropical",
                       scope=scope, sector=sector,
                       spatial_block=f"b{i//2}",
                       source_signature=f"s{(i//2)%6}",
                       **{k:v for k,v in result.items()
                          if not k.startswith("shared_genus_")})
                for group in ("shared_genus_one_donor_each",
                              "shared_genus_two_donors_each"):
                    for key,value in result[group].items():
                        z[group+"__"+key]=value
                if coverage_gap:
                    z["shared_genus_one_donor_each__n_matched_island_total"]=50
                rows.append(z)
    return pd.DataFrame(rows)


def test_retains_full_4505_style_frame_and_estimates_independent_blocks():
    x=summarize_source_genus_frame(example_frame(), population=40, bootstrap=299)
    assert x["n_regime_scope_morphology_panels"]==100
    z=next(q for q in x["panels"] if q["region"]=="tropical" and
           q["plant_provenance_scope"]=="all_recorded" and
           q["floral_sector"]=="zygomorphic")
    assert z["n_original_islands_congeneric_both_sources"]==40
    assert z["n_original_islands_min2_donor_species_per_genus"]==40
    m=z["one_donor__full_matched_sample__spatial_block"]
    assert m["status"]=="ORIGINAL_CONGENERIC_SOURCE_DESCRIPTIVE_NOT_CAUSAL"
    assert np.isclose(m["observed_minus_geo5_congener_reference__mean"],.6)
    assert np.isclose(m["observed_minus_climate10_congener_reference__mean"],0)
    assert np.isclose(m["climate10_minus_geo5_congener_residual__mean"],-.6)


def test_original_documentation_imbalance_reduces_supported_panel():
    x=summarize_source_genus_frame(example_frame(40,coverage_gap=True),
                                   population=40,bootstrap=299)
    z=next(q for q in x["panels"] if q["region"]=="tropical" and
           q["plant_provenance_scope"]=="all_recorded" and
           q["floral_sector"]=="zygomorphic")
    assert z["n_original_islands_congeneric_both_sources"]==40
    assert z["n_original_islands_min2_donor_species_per_genus"]==40
    assert z["n_original_islands_native_and_source_status_not_inferred"]
    assert z["one_donor__n_coverage_gap_both_le10pp"]==0
    assert z["one_donor__matched_direct_documentation_gap_le10pp__spatial_block"][
        "status"]=="NOT_IDENTIFIABLE_CLUSTER_SUPPORT"


def test_missing_original_coverage_or_duplicate_roster_fails():
    x=example_frame()
    with pytest.raises(ValueError,match="duplicated"):
        summarize_source_genus_frame(pd.concat([x,x.iloc[:1]]),
                                     population=40,bootstrap=299)
    with pytest.raises(ValueError,match="must remain"):
        summarize_source_genus_frame(x,population=41,bootstrap=299)
