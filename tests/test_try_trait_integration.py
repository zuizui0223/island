import pandas as pd

from island_v2 import all_evidence_trait_audit as audit
from island_v2 import try_trait_integration as integration
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

ONTOLOGY = {
    "flower_primary_color": {
        "white",
        "red_pink",
        "yellow_orange",
        "blue_purple",
        "green_brown_inconspicuous",
        "multicolored_variable",
        "other_described",
    },
    "floral_symmetry": {"actinomorphic", "zygomorphic", "asymmetric"},
    "self_incompatibility": {"SI", "SC", "mixed_or_variable"},
}


def evidence(
    species: str,
    trait: str,
    value: str,
    lineage: str,
    *,
    source_group: str,
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "axis": audit.trait_axis(trait),
        "trait_name": trait,
        "normalized_value": value,
        "quality": "medium",
        "source_group": source_group,
        "source_provider": source_group,
        "source_url": "https://example.org/source",
        "source_record_id": f"{source_group}:{species}:{trait}",
        "source_citation": "Example original source",
        "source_excerpt": f"{trait}={value}",
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": lineage,
        "lineage_method": "test",
        "source_run_id": "test",
        "source_artifact": "test",
        "source_file": "test.csv",
        "acceptance_contract": "test",
    }


def empty_direct() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "accepted_species",
            "axis",
            "trait_name",
            "quality",
            "state_set",
            "normalized_value",
            "source_groups",
            "source_lineages",
            "classification",
            "genus",
        ]
    )


def empty_low() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "accepted_species",
            "genus",
            "axis",
            "trait_name",
            "state_set",
            "normalized_value",
            "source_lineage",
        ]
    )


def master() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"accepted_species": "Alpha beta", "genus": "Alpha"},
            {"accepted_species": "Beta gamma", "genus": "Beta"},
        ]
    )


def current_state(formal: pd.DataFrame):
    lineages, _ = audit.dedupe_direct_lineages(formal, ONTOLOGY)
    direct, _ = audit.resolve_direct_cells(lineages)
    if direct.empty:
        direct = empty_direct()
    else:
        direct["genus"] = direct["accepted_species"].str.split().str[0]
    low = empty_low()
    coverage = audit.species_axis_coverage(master(), direct, low)
    return direct, low, coverage


def test_try_new_direct_cell_updates_strict_axis_coverage():
    formal = pd.DataFrame(columns=EVIDENCE_COLUMNS)
    direct, low, coverage = current_state(formal)
    try_common = pd.DataFrame(
        [
            evidence(
                "Beta gamma",
                "floral_symmetry",
                "actinomorphic",
                "doi:10.1234/new",
                source_group="try",
            )
        ],
        columns=EVIDENCE_COLUMNS,
    )

    result = integration.integrate(
        try_common=try_common,
        formal_direct_evidence=formal,
        additional_common=[],
        current_direct=direct,
        current_low=low,
        current_coverage=coverage,
        master_genus_map=master(),
        ontology=ONTOLOGY,
    )

    updated = result["direct"].loc[
        result["direct"]["accepted_species"].eq("Beta gamma")
    ].iloc[0]
    assert updated["trait_name"] == "floral_symmetry"
    assert updated["normalized_value"] == "actinomorphic"
    cell = result["coverage"].loc[
        result["coverage"]["accepted_species"].eq("Beta gamma")
        & result["coverage"]["axis"].eq("floral_structural_complexity")
    ].iloc[0]
    assert cell["quality"] == "medium"
    assert result["summary"]["direct"]["added_species_trait"] == 1


def test_try_independent_si_conflict_removes_prior_resolved_cell():
    formal = pd.DataFrame(
        [
            evidence(
                "Alpha beta",
                "self_incompatibility",
                "SC",
                "doi:10.1111/prior",
                source_group="meyer_2026",
            )
        ],
        columns=EVIDENCE_COLUMNS,
    )
    direct, low, coverage = current_state(formal)
    assert len(direct) == 1
    try_common = pd.DataFrame(
        [
            evidence(
                "Alpha beta",
                "self_incompatibility",
                "SI",
                "doi:10.2222/try",
                source_group="try",
            )
        ],
        columns=EVIDENCE_COLUMNS,
    )

    result = integration.integrate(
        try_common=try_common,
        formal_direct_evidence=formal,
        additional_common=[],
        current_direct=direct,
        current_low=low,
        current_coverage=coverage,
        master_genus_map=master(),
        ontology=ONTOLOGY,
    )

    assert result["direct"].loc[
        result["direct"]["accepted_species"].eq("Alpha beta")
        & result["direct"]["trait_name"].eq("self_incompatibility")
    ].empty
    conflict = result["conflicts"].loc[
        result["conflicts"]["accepted_species"].eq("Alpha beta")
    ].iloc[0]
    assert conflict["classification"] == "unresolved_direct_conflict"
    assert conflict["resolution_status"] == "excluded"
    cell = result["coverage"].loc[
        result["coverage"]["accepted_species"].eq("Alpha beta")
        & result["coverage"]["axis"].eq("reproductive_assurance")
    ].iloc[0]
    assert cell["quality"] == ""
    assert result["summary"]["direct"]["removed_species_trait"] == 1
