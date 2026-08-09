from pathlib import Path

import pandas as pd

from analysis.prepare_near_rule_source_packages_20260809 import (
    _mating_system_from_outcrossing_rate,
    _ncsu_value,
    _normalise_colour,
    curated_primary_rows,
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


def test_outcrossing_rate_stays_a_mating_system_trait() -> None:
    assert _mating_system_from_outcrossing_rate(0.19) == "predominantly_selfing"
    assert _mating_system_from_outcrossing_rate(0.2) == "mixed_mating"
    assert _mating_system_from_outcrossing_rate(0.8) == "mixed_mating"
    assert _mating_system_from_outcrossing_rate(0.81) == "predominantly_outcrossing"
    assert _mating_system_from_outcrossing_rate(1.1) == ""


def test_committed_near_rule_package_passes_common_gate() -> None:
    evidence = pd.read_csv(
        PACKAGE / "near_rule_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(PACKAGE / "near_rule_source_audit.csv", dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 208
    assert summary["accepted_correct"] == 208
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert len(evidence) == 208
    # The exact autonomous-selfing record is retained, but its one-row trait
    # stratum remains below the production minimum of ten reviewed records.
    assert len(selected) == 206
    assert evidence["source_record_id"].is_unique
    assert evidence["source_provider"].nunique() == 3
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "mating_system",
        "self_incompatibility",
    }
    autonomous = scopes.loc[
        scopes["trait_name"].eq("autonomous_selfing_capacity")
    ].iloc[0]
    assert autonomous["reviewed"] == 1
    assert not autonomous["production_approved"]


def test_moeller_global_mating_rows_retain_rates_and_source_lineage() -> None:
    evidence = pd.read_csv(
        PACKAGE / "near_rule_incremental_evidence.csv.gz", dtype=str
    ).fillna("")
    moeller = evidence.loc[
        evidence["source_provider"].eq("moeller_etal_2017_dryad")
    ]

    assert len(moeller) == 206
    assert moeller["accepted_species"].nunique() == 161
    assert len(moeller.drop_duplicates(["accepted_species", "trait_name"])) == 173
    assert set(moeller["trait_name"]) == {"mating_system", "self_incompatibility"}
    assert moeller["source_url"].eq("https://doi.org/10.5061/dryad.577q1").all()
    assert moeller["source_excerpt"].str.contains(
        r"(?:mean\.tm:|si:)"
    ).all()
    assert moeller["source_lineage"].str.startswith(
        ("citation:", "doi:10.5061/dryad.577q1:row:")
    ).all()


def test_wave2_symmetry_records_are_exact_species_direct_evidence() -> None:
    master_names = {
        "Bulbophyllum singaporeanum",
        "Commelina erecta",
        "Drosera paradoxa",
        "Epipremnum aureum",
        "Pilea microphylla",
        "Scaevola taccada",
        "Tacca cristata",
        "Uvaria grandiflora",
        "Uvaria hirsuta",
        "Uvaria littoralis",
    }
    evidence = pd.DataFrame(curated_primary_rows(master_names, set()))
    symmetry = evidence.loc[
        evidence["source_run_id"].eq("manual-source-backed-web-20260809-wave2")
        & evidence["trait_name"].eq("floral_symmetry")
    ]

    assert len(symmetry) == 10
    assert symmetry["accepted_species"].nunique() == 10
    assert symmetry["source_lineage"].nunique() == 10
    assert set(symmetry["normalized_value"]) == {"actinomorphic", "zygomorphic"}
    assert symmetry["source_excerpt"].str.startswith("Flower Symmetry | ").sum() == 9
    assert symmetry.loc[
        symmetry["accepted_species"].eq("Tacca cristata"), "normalized_value"
    ].item() == "zygomorphic"


def test_curated_primary_records_keep_reproductive_concepts_separate() -> None:
    evidence = pd.DataFrame(
        curated_primary_rows(
            {"Angraecum cadetii", "Elaeocarpus angustifolius"}, set()
        )
    )
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
