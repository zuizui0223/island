from pathlib import Path

import pandas as pd

from analysis.prepare_near_rule_source_packages_20260809 import (
    _ncsu_value,
    _normalise_colour,
)
from island_v2.open_web_common import reviewed_source_package_evidence

PACKAGE = Path(
    "data/v2/staging/traits/direct_llm_pilot/20260809_near_rule_source_acquisition"
)


def test_colour_and_structured_ncsu_values_preserve_state_sets() -> None:
    assert _normalise_colour("White-pink") == "red_pink|white"
    assert _normalise_colour("green-golden yellow") == (
        "green_brown_inconspicuous|yellow_orange"
    )
    assert _ncsu_value("floral_form", "Flower Shape:: Radial; Tubular") == (
        "open_radial|tubular"
    )
    assert _ncsu_value(
        "inflorescence_display", "Flower Inflorescence:: Panicle; Umbel"
    ) == "raceme_spike_panicle|umbel_corymb"


def test_committed_near_rule_package_passes_common_gate() -> None:
    evidence = pd.read_csv(
        PACKAGE / "near_rule_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(PACKAGE / "near_rule_source_audit.csv", dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 417
    assert summary["accepted_correct"] == 417
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert len(evidence) == 1817
    # The exact autonomous-selfing record is retained, but its one-row trait
    # stratum remains below the production minimum of ten reviewed records.
    assert len(selected) == 1816
    assert evidence["source_record_id"].is_unique
    assert evidence["source_provider"].nunique() == 8
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
        "self_incompatibility",
        "tube_depth_class",
    }
    autonomous = scopes.loc[
        scopes["trait_name"].eq("autonomous_selfing_capacity")
    ].iloc[0]
    assert autonomous["reviewed"] == 1
    assert not autonomous["production_approved"]


def test_curated_primary_records_keep_reproductive_concepts_separate() -> None:
    evidence = pd.read_csv(
        PACKAGE / "near_rule_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    angraecum = evidence.loc[
        evidence["accepted_species"].eq("Angraecum cadetii")
        & evidence["source_provider"].eq("micheneau_etal_2010_primary_article")
    ]
    elaeocarpus = evidence.loc[
        evidence["accepted_species"].eq("Elaeocarpus angustifolius")
        & evidence["source_provider"].eq("flora_of_peninsular_india")
    ]

    assert set(zip(angraecum["trait_name"], angraecum["normalized_value"], strict=True)) == {
        ("self_incompatibility", "SC"),
        ("autonomous_selfing_capacity", "absent"),
    }
    assert angraecum["source_lineage"].eq("doi:10.1093/aob/mcp299").all()
    assert len(elaeocarpus) == 1
    assert elaeocarpus.iloc[0]["trait_name"] == "floral_form"
    assert elaeocarpus.iloc[0]["normalized_value"] == "bell_campanulate"


def test_flora_measurements_and_lineages_fail_closed() -> None:
    evidence = pd.read_csv(
        PACKAGE / "near_rule_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    plantnet = evidence.loc[evidence["source_provider"].eq("plantnet_nsw_flora")]
    size = plantnet.loc[plantnet["trait_name"].eq("flower_size_class")]
    ncsu = evidence.loc[evidence["source_provider"].eq("ncsu_extension_plant_toolbox")]
    bien = evidence.loc[evidence["source_provider"].eq("bien_trait_original_citation")]

    assert not size["source_excerpt"].str.contains(
        r"(?i)\b(?:pedicel|peduncle|cluster|cyme|inflorescence|spikelet|floret|flowered)\w*\b"
    ).any()
    assert not set(ncsu["trait_name"]).intersection({"flower_size_class"})
    assert ncsu["source_excerpt"].str.startswith("Flower ").all()
    assert bien["source_lineage"].str.startswith("origin:").all()
    assert evidence.loc[
        evidence["source_provider"].eq("sinnott_armstrong_2025_dryad"),
        "source_lineage",
    ].eq("doi:10.5061/dryad.r4xgxd2sc").all()
    assert evidence.loc[
        evidence["source_provider"].eq("goldberg_etal_2010_dryad"),
        "source_lineage",
    ].eq("doi:10.5061/dryad.1888").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type", "sex_system"}
    )
