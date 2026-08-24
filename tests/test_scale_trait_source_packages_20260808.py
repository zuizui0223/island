import inspect
from pathlib import Path

import pandas as pd

from analysis.acquire_ken_fern_reproductive_traits_20260808 import (
    parse_page as parse_ken_fern_page,
)
from analysis.acquire_pfaf_reproductive_traits_20260808 import (
    parse_page as parse_pfaf_page,
)
from analysis.prepare_scale_trait_source_packages_20260808 import (
    _baseflor_colour,
    _baseflor_inflorescence,
    _floraweb_selfing_trait,
    alien_plants_autogamy_rows,
    baseflor_rows,
)
from island_v2.open_web_common import reviewed_source_package_evidence

PACKAGE = Path("data/v2/staging/traits/direct_llm_pilot/20260808_scale_source_acquisition")


def test_baseflor_multistate_values_remain_trait_specific() -> None:
    assert _baseflor_colour("blanc, jaune, bleu") == ("blue_purple|white|yellow_orange")
    assert _baseflor_inflorescence("racème de capitules") == (
        "composite_display|raceme_spike_panicle"
    )
    assert _baseflor_inflorescence("fleur solitaire terminale") == "solitary"


def test_pollination_fields_do_not_substitute_for_autonomous_reproductive_output() -> None:
    assert _floraweb_selfing_trait("häufig Selbstbestäubung") == ("", "")
    assert _floraweb_selfing_trait("bei ausbleib. Fremdbest. Selbstbestäubung") == ("", "")
    assert _floraweb_selfing_trait("Kleistogamie") == ("cleistogamy", "facultative")
    assert 'trait="autonomous_selfing_capacity"' not in inspect.getsource(baseflor_rows)
    assert 'trait="autonomous_selfing_capacity"' not in inspect.getsource(
        alien_plants_autogamy_rows
    )


def test_committed_scale_package_passes_common_gate() -> None:
    evidence = pd.read_csv(PACKAGE / "scale_source_incremental_evidence.csv.gz", dtype=str).fillna(
        ""
    )
    audit = pd.read_csv(PACKAGE / "scale_source_manual_audit.csv", dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 1622
    assert summary["accepted_correct"] == 1622
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert len(evidence) == 13227
    assert len(selected) == 13221
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 12622
    assert evidence["source_record_id"].is_unique
    assert evidence["source_provider"].nunique() == 12
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "autonomous_selfing_capacity",
        "floral_form",
        "flower_primary_color",
        "inflorescence_display",
        "mating_system",
        "self_incompatibility",
        "tube_depth_class",
    }


def test_provider_redistribution_lineages_and_strict_axes_are_safe() -> None:
    evidence = pd.read_csv(PACKAGE / "scale_source_incremental_evidence.csv.gz", dtype=str).fillna(
        ""
    )
    nhm = evidence.loc[evidence["source_provider"].str.startswith("nhm_adept_flora_")]
    reused = nhm.loc[nhm["source_url"].str.contains("efloras.org", case=False)]

    assert not reused.empty
    assert reused["source_lineage"].str.startswith("url:http://efloras.org/").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type", "sex_system"}
    )
    assert not evidence["source_excerpt"].eq("").any()
    assert not evidence["source_url"].eq("").any()
    assert not evidence["source_lineage"].eq("").any()
    assert (
        not evidence.loc[evidence["source_provider"].eq("tree_of_sex_dryad"), "accepted_species"]
        .isin({"Thalictrum dioicum", "Thalictrum pubescens"})
        .any()
    )
    assert (
        not evidence.loc[evidence["source_provider"].eq("tree_of_sex_dryad"), "accepted_species"]
        .str.startswith("Phoenix ")
        .any()
    )
    package_species = set(evidence["accepted_species"])
    assert "Ephedra distachya" not in package_species
    assert "Ephedra major" not in package_species


def test_florabase_quotes_are_flower_clauses_not_search_snippets() -> None:
    evidence = pd.read_csv(PACKAGE / "scale_source_incremental_evidence.csv.gz", dtype=str).fillna(
        ""
    )
    florabase = evidence.loc[evidence["source_provider"].eq("florabase_western_australia")]

    assert len(florabase) == 694
    assert florabase["source_excerpt"].str.match(r"(?i)^Fl\.\s").all()
    assert florabase["source_url"].str.startswith("https://florabase.dbca.wa.gov.au/").all()


def test_floraweb_non_axis_traits_remain_supplemental() -> None:
    evidence = pd.read_csv(PACKAGE / "scale_source_incremental_evidence.csv.gz", dtype=str).fillna(
        ""
    )
    supplemental = pd.read_csv(
        PACKAGE / "scale_source_non_axis_candidates.csv.gz", dtype=str
    ).fillna("")
    floraweb = evidence.loc[evidence["source_provider"].eq("floraweb_biolflor_snapshot_2026_06_24")]

    assert len(supplemental) == 3655
    assert set(supplemental["trait_name"]) == {"pollen_vector_mode", "reward_type"}
    assert supplemental["axis"].eq("").all()
    assert not set(floraweb["trait_name"]).intersection({"pollen_vector_mode", "reward_type"})
    assert not floraweb["source_excerpt"].str.contains("Pseudokleistogamie").any()


def test_self_fertility_adapters_fail_closed_on_category_errors() -> None:
    def pfaf_html(statement: str) -> bytes:
        return (
            '<span id="ContentPlaceHolder1_lbldisplatinname">Testus example</span>'
            f'<span id="ContentPlaceHolder1_lblPhystatment">{statement}</span>'
        ).encode()

    assert (
        parse_pfaf_page(
            "Testus example", "https://example.test", pfaf_html("The plant is self-fertile.")
        )["status"]
        == "accepted_self_fertile_statement"
    )
    assert (
        parse_pfaf_page(
            "Testus example",
            "https://example.test",
            pfaf_html("The plant is self-fertile and reproduces apomictically."),
        )["status"]
        == "nonsexual_reproduction_category_error"
    )
    assert (
        parse_pfaf_page(
            "Testus example",
            "https://example.test",
            pfaf_html("The species is dioecious. The plant is self-fertile."),
        )["status"]
        == "sexual_system_internal_conflict"
    )

    def fern_html(cultivation: str) -> bytes:
        return (
            '<div class="latin_name"><h1>Testus example</h1></div>'
            '<div class="family"><h4>Testaceae</h4></div>'
            '<table class="tableProperties">'
            "<tr><td>Self-fertile</td><td>Yes</td></tr>"
            f"<tr><td>Cultivation Status</td><td>{cultivation}</td></tr>"
            "</table>"
        ).encode()

    assert (
        parse_ken_fern_page(
            "Testus example", "Testaceae", "https://example.test", fern_html("Wild")
        )["status"]
        == "accepted_wild_self_fertile"
    )
    assert (
        parse_ken_fern_page(
            "Testus example", "Testaceae", "https://example.test", fern_html("Cultivated")
        )["status"]
        == "cultivated_or_ornamental_only"
    )


def test_literature_synonym_is_not_relabelled_species_direct() -> None:
    evidence = pd.read_csv(PACKAGE / "scale_source_incremental_evidence.csv.gz", dtype=str).fillna(
        ""
    )
    fireweed = evidence.loc[evidence["accepted_species"].eq("Chamaenerion angustifolium")]

    assert len(fireweed) == 1
    assert fireweed.iloc[0]["name_match_method"] == "synonym_exact_two_backbone"
    assert fireweed.iloc[0]["evidence_scope"] == "synonym_direct"
