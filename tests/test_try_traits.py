import pandas as pd

from island_v2 import try_traits

BASE = {
    "DatasetID": "1",
    "Dataset": "Provider dataset",
    "SpeciesName": "Alpha beta",
    "AccSpeciesName": "Alpha beta",
    "ObservationID": "1",
    "ObsDataID": "a",
    "Reference": "Smith doi:10.1234/ABC",
}


def row(**kwargs):
    value = BASE.copy()
    value.update(kwargs)
    return value


def master():
    return pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "genus": "Alpha",
                "family": "Alphaceae",
            },
            {
                "accepted_species": "Gamma delta",
                "genus": "Gamma",
                "family": "Gammaceae",
            },
        ]
    )


def test_maps_three_strict_traits_and_corolla_sidecar():
    source = pd.DataFrame(
        [
            row(
                TraitID="207",
                TraitName="Flower color",
                DataID="583",
                DataName="Flower color",
                OrigValueStr="white",
            ),
            row(
                TraitID="2935",
                TraitName="Flower symmetry type",
                DataID="3793",
                DataName="Flower - Symmetry",
                OrigValueStr="radial",
                ObsDataID="b",
            ),
            row(
                TraitID="211",
                TraitName="Flower sexual self-incompatibility",
                DataID="589",
                DataName="Self-incompatibility",
                OrigValueStr="self-compatibel species",
                ObsDataID="c",
            ),
            row(
                TraitID="2936",
                TraitName="Flower corolla type",
                DataID="3794",
                DataName="Flower - Corolla",
                OrigValueStr="free",
                ObsDataID="d",
            ),
        ]
    )
    candidates, corolla, _, exclusions, _ = try_traits.build(source, master())
    assert set(candidates.trait_name) == {
        "flower_primary_color",
        "floral_symmetry",
        "self_incompatibility",
    }
    assert dict(zip(candidates.trait_name, candidates.standardized_value, strict=True)) == {
        "flower_primary_color": "white",
        "floral_symmetry": "actinomorphic",
        "self_incompatibility": "SC",
    }
    assert corolla[["trait_name", "standardized_value"]].iloc[0].tolist() == [
        "corolla_fusion",
        "free",
    ]
    assert exclusions.empty
    common = try_traits.common_evidence(candidates)
    assert set(common.source_group) == {"try"}
    assert dict(zip(common.trait_name, common.axis, strict=True)) == {
        "flower_primary_color": "flower_colour",
        "floral_symmetry": "floral_structural_complexity",
        "self_incompatibility": "reproductive_assurance",
    }


def test_multicolour_same_lineage_aggregates_not_duplicates():
    source = pd.DataFrame(
        [
            row(
                TraitID="207",
                TraitName="Flower color",
                DataID="583",
                DataName="Flower color",
                OrigValueStr="white",
            ),
            row(
                TraitID="207",
                TraitName="Flower color",
                DataID="583",
                DataName="Flower color",
                OrigValueStr="red",
                ObsDataID="b",
                ObservationID="2",
            ),
        ]
    )
    candidates, *_ = try_traits.build(source, master())
    assert len(candidates) == 1
    assert candidates.iloc[0].standardized_value == "multicolored_variable"
    assert candidates.iloc[0].source_lineage == "doi:10.1234/abc"


def test_symmetry_negative_boolean_does_not_infer_opposite():
    source = pd.DataFrame(
        [
            row(
                TraitID="2935",
                TraitName="Flower symmetry type",
                DataID="3987",
                DataName="Flower shape: zygomorphic",
                OrigValueStr="no",
            ),
            row(
                TraitID="2935",
                TraitName="Flower symmetry type",
                DataID="3988",
                DataName="Flower shape: actinomorphic",
                OrigValueStr="yes",
                ObsDataID="b",
            ),
        ]
    )
    candidates, _, _, exclusions, _ = try_traits.build(source, master())
    assert candidates.iloc[0].standardized_value == "actinomorphic"
    assert "negative_boolean_not_opposite_evidence" in set(exclusions.reason)


def test_si_rejects_non_species_and_inbreeding_but_accepts_mechanism():
    source = pd.DataFrame(
        [
            row(
                TraitID="211",
                TraitName="Flower sexual self-incompatibility",
                DataID="589",
                DataName="Self-incompatibility",
                OrigValueStr="self-compatibel genus",
            ),
            row(
                TraitID="211",
                TraitName="Flower sexual self-incompatibility",
                DataID="2045",
                DataName="Inbreeding",
                OrigValueStr="99",
                ObsDataID="b",
            ),
            row(
                TraitID="211",
                TraitName="Flower sexual self-incompatibility",
                DataID="592",
                DataName="Flower self-incompatibility mechanism",
                OrigValueStr="gametophytic self-incompatibility",
                ObsDataID="c",
                Reference="Other doi:10.9999/X",
            ),
        ]
    )
    candidates, _, _, exclusions, _ = try_traits.build(source, master())
    assert candidates[["trait_name", "standardized_value"]].iloc[0].tolist() == [
        "self_incompatibility",
        "SI",
    ]
    assert {
        "non_species_si_scope",
        "inbreeding_measure_not_self_incompatibility",
    } <= set(exclusions.reason)


def test_within_lineage_noncolour_conflict_is_excluded():
    source = pd.DataFrame(
        [
            row(
                TraitID="2935",
                TraitName="Flower symmetry type",
                DataID="3793",
                DataName="Flower - Symmetry",
                OrigValueStr="radial",
            ),
            row(
                TraitID="2935",
                TraitName="Flower symmetry type",
                DataID="3793",
                DataName="Flower - Symmetry",
                OrigValueStr="bilateral",
                ObsDataID="b",
            ),
        ]
    )
    candidates, _, _, exclusions, _ = try_traits.build(source, master())
    assert candidates.empty
    assert set(exclusions.reason) == {"within_original_source_lineage_conflict"}


def test_case_normalized_and_synonym_matches():
    source = pd.DataFrame(
        [
            row(
                TraitID="207",
                TraitName="Flower color",
                DataID="583",
                DataName="Flower color",
                OrigValueStr="yellow",
                AccSpeciesName="ALPHA BETA",
                SpeciesName="ALPHA BETA",
            ),
            row(
                TraitID="207",
                TraitName="Flower color",
                DataID="583",
                DataName="Flower color",
                OrigValueStr="blue",
                AccSpeciesName="Oldus nameus",
                SpeciesName="Oldus nameus",
                Reference="R2",
                ObservationID="2",
                ObsDataID="b",
            ),
        ]
    )
    candidates, *_ = try_traits.build(
        source,
        master(),
        {"oldus nameus": "Gamma delta"},
    )
    methods = dict(zip(candidates.accepted_species, candidates.name_match_method, strict=True))
    assert methods["Alpha beta"] == "accepted_name_exact"
    assert methods["Gamma delta"] == "exact_synonym_to_accepted"
