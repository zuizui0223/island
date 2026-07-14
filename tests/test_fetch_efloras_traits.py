from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace


SCRIPT = Path(__file__).parents[1] / "scripts" / "fetch_efloras_traits.py"
SPEC = importlib.util.spec_from_file_location("fetch_efloras_traits", SCRIPT)
assert SPEC and SPEC.loader
efloras = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = efloras
SPEC.loader.exec_module(efloras)


def listing_html(name: str, taxon_id: str, child: bool = True) -> str:
    child_link = (
        f'<a href="browse.aspx?flora_id=2&amp;start_taxon_id={taxon_id}">3</a>'
        if child
        else ""
    )
    return f"""
    <html><body><table><tr class="underline">
      <td>{taxon_id}</td>
      <td><a href="florataxon.aspx?flora_id=2&amp;taxon_id={taxon_id}">{name}</a></td>
      <td>{child_link}</td>
    </tr></table></body></html>
    """


def candidate_values(description: str) -> set[tuple[str, str]]:
    rows = efloras.trait_candidates(
        flora_id=2,
        flora_name="Flora of China",
        accepted_species="Acanthus ebracteatus",
        original_taxon_name="Acanthus ebracteatus",
        taxon_id="242300348",
        description=description,
        source_url="http://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=242300348",
        name_match_method="exact_accepted_name_match",
    )
    return {(row["trait_name"], row["candidate_value"]) for row in rows}


def test_binomial_key_strips_author_without_fuzzy_matching() -> None:
    assert efloras.binomial_key("Poa annua L.") == "poa annua"


def test_binomial_key_rejects_infraspecies_and_hybrids() -> None:
    assert efloras.binomial_key("Poa annua subsp. annua") is None
    assert efloras.binomial_key("Salix × sepulcralis") is None


def test_taxonomic_rank_keeps_species_infraspecies_and_hybrid_distinct() -> None:
    assert efloras.taxonomic_rank("Acanthus ebracteatus") == "species"
    assert efloras.taxonomic_rank("Acanthus ebracteatus var. alba") == "infraspecies"
    assert efloras.taxonomic_rank("Acanthus × example") == "hybrid"


def test_master_matcher_reports_all_conservative_match_classes() -> None:
    matcher = efloras.MasterMatcher(
        ["Acanthus ebracteatus", "Poa annua"],
        [("Dilivaria ebracteata", "Acanthus ebracteatus")],
    )
    assert matcher.match("Acanthus ebracteatus")["match_status"] == "exact_accepted_name_match"
    assert matcher.match("Poa annua L.")["match_status"] == "normalized_binomial_match"
    assert matcher.match("Dilivaria ebracteata")["match_status"] == "synonym_match"
    assert matcher.match("Unknown plant")["match_status"] == "unmatched"


def test_master_matcher_reports_ambiguous_normalized_binomial() -> None:
    matcher = efloras.MasterMatcher(["Poa annua", "Poa annua L."])
    assert matcher.match("Poa annua Smith")["match_status"] == "ambiguous"


def test_listing_parser_infers_family_from_flora_root() -> None:
    task = efloras.BrowseTask(url="http://example.test/browse", parent_rank="flora")
    row = efloras.parse_listing_page(listing_html("Acanthaceae", "10002"), task.url, 2, task)[0]
    assert row["listing_rank"] == "family"
    assert row["family"] == "Acanthaceae"
    assert row["browse_url"].endswith("start_taxon_id=10002")


def test_listing_parser_follows_family_to_genus() -> None:
    task = efloras.BrowseTask(
        url="http://example.test/browse",
        parent_rank="family",
        parent_taxon_id="10002",
        family="Acanthaceae",
    )
    row = efloras.parse_listing_page(listing_html("Acanthus", "100149"), task.url, 2, task)[0]
    assert row["listing_rank"] == "genus"
    assert row["genus"] == "Acanthus"
    assert row["parent_taxon_id"] == "10002"


def test_listing_parser_follows_genus_to_species() -> None:
    task = efloras.BrowseTask(
        url="http://example.test/browse",
        parent_rank="genus",
        parent_taxon_id="100149",
        family="Acanthaceae",
        genus="Acanthus",
    )
    row = efloras.parse_listing_page(
        listing_html("Acanthus ebracteatus", "242300348", child=False), task.url, 2, task
    )[0]
    assert row["listing_rank"] == "species"
    assert row["family"] == "Acanthaceae"
    assert row["genus"] == "Acanthus"


def test_treatment_parser_requires_bold_heading_and_description() -> None:
    parsed = efloras.parse_treatment(
        '<span id="lblTaxonDesc">1. <b>Acanthus ebracteatus</b> Vahl. '
        "Shrubs. Corolla white and tubular.</span>",
        "http://example.test/treatment",
    )
    assert parsed["taxon_heading"] == "Acanthus ebracteatus"
    assert parsed["taxonomic_rank"] == "species"
    assert parsed["treatment_confirmed"] is True


def test_treatment_parser_rejects_empty_description_and_non_species() -> None:
    empty = efloras.parse_treatment('<span id="lblTaxonDesc"><b>Acanthus ebracteatus</b></span>', "x")
    genus = efloras.parse_treatment(
        '<span id="lblTaxonDesc"><b>Acanthus</b> Genus description long enough to parse.</span>', "x"
    )
    assert empty["treatment_confirmed"] is False
    assert genus["taxonomic_rank"] == "genus"
    assert genus["treatment_confirmed"] is False


def test_treatment_parser_rejects_checklist_only_bibliography() -> None:
    checklist = efloras.parse_treatment(
        '<span id="lblTaxonDesc"><b>Tribulus terrestris</b> L., Sp. Pl. 387 (1753). '
        "Gokhur. Gardner 1653. 150 m; Almost Pantropical.</span>",
        "http://example.test/checklist",
    )
    assert checklist["taxonomic_rank"] == "species"
    assert checklist["narrative_chars"] > 20
    assert checklist["treatment_confirmed"] is False


def test_flower_colour_requires_a_floral_organ_in_same_fragment() -> None:
    values = candidate_values(
        "Leaves green; fruits red; seeds brown. Corolla white. Petals purple or violet."
    )
    assert ("flower_primary_color", "white") in values
    assert ("flower_primary_color", "purple|violet") in values
    assert all(value not in {"green", "red", "brown"} for trait, value in values if trait == "flower_primary_color")


def test_lip_and_labellum_are_valid_flower_colour_contexts() -> None:
    values = candidate_values("Lip yellow; labellum red and pinkish.")
    assert ("flower_primary_color", "yellow") in values
    assert ("flower_primary_color", "pink|red") in values


def test_colour_gate_rejects_indumentum_and_later_fruit_colour() -> None:
    values = candidate_values(
        "Sepals outside with dense white scales, inside sparsely white appressed puberulent; "
        "tepals white or pink, in fruit larger and green."
    )
    assert ("flower_primary_color", "white") not in values
    assert ("flower_primary_color", "green|pink|white") not in values
    assert ("flower_primary_color", "pink|white") in values


def test_form_symmetry_sex_and_self_incompatibility_are_unreviewed_candidates() -> None:
    description = (
        "Flowers bisexual and self-incompatible. Corolla campanulate and actinomorphic. "
        "Petals bilaterally symmetrical."
    )
    rows = efloras.trait_candidates(
        flora_id=2,
        flora_name="Flora of China",
        accepted_species="Example species",
        original_taxon_name="Example species",
        taxon_id="1",
        description=description,
        source_url="http://example.test/1",
        name_match_method="exact_accepted_name_match",
    )
    values = {(row["trait_name"], row["candidate_value"]) for row in rows}
    assert ("floral_form", "campanulate") in values
    assert ("floral_symmetry", "radial") in values
    assert ("floral_symmetry", "bilateral") in values
    assert ("sex_system", "bisexual") in values
    assert ("self_incompatibility", "self_incompatible") in values
    assert {row["evidence_status"] for row in rows} == {"source_backed_rule_extracted_unreviewed"}
    assert all(row["source_excerpt"] and row["source_url"] for row in rows)


def test_match_inventory_excludes_infraspecies_and_deduplicates_accepted_species() -> None:
    matcher = efloras.MasterMatcher(["Poa annua", "Poa annua var. annua"])
    taxa = [
        {
            "taxon_id": "1",
            "original_taxon_name": "Poa annua",
            "listing_rank": "species",
            "flora_id": "2",
            "family": "Poaceae",
            "genus": "Poa",
            "parent_taxon_id": "9",
            "treatment_url": "x",
            "browse_url": "",
            "discovery_url": "x",
        },
        {
            "taxon_id": "2",
            "original_taxon_name": "Poa annua var. annua",
            "listing_rank": "infraspecies",
            "flora_id": "2",
            "family": "Poaceae",
            "genus": "Poa",
            "parent_taxon_id": "1",
            "treatment_url": "y",
            "browse_url": "",
            "discovery_url": "y",
        },
    ]
    report, selected = efloras._match_inventory(taxa, 2, matcher)
    assert [row["taxon_id"] for row in selected] == ["1"]
    assert next(row for row in report if row["taxon_id"] == "2")["selection_status"] == "excluded_infraspecies"


def test_checkpoint_persistence_deduplicates_on_reload(tmp_path: Path) -> None:
    treatment = {field: "" for field in efloras.TREATMENT_FIELDS}
    treatment.update({"accepted_species": "A a", "flora_id": "2", "taxon_id": "1"})
    candidate = {field: "" for field in efloras.CANDIDATE_FIELDS}
    candidate.update(
        {
            "accepted_species": "A a",
            "flora_id": "2",
            "taxon_id": "1",
            "trait_name": "floral_form",
            "candidate_value": "tubular",
            "source_excerpt": "Corolla tubular.",
        }
    )
    efloras._persist_treatment_state(
        output_dir=tmp_path,
        treatments=[treatment],
        candidates=[candidate],
        errors=[],
        processed_taxon_ids={"1"},
        request_total=1,
    )
    checkpoint = json.loads((tmp_path / "treatment_checkpoint.json").read_text(encoding="utf-8"))
    assert checkpoint["processed_taxon_ids"] == ["1"]
    assert len(efloras.read_csv(tmp_path / "species_treatments.csv")) == 1


def test_acquisition_status_distinguishes_missing_species_depth_from_no_match() -> None:
    common = {
        "crawl_stop_reason": "queue_exhausted",
        "pending_selected": 0,
        "treatment_count": 0,
        "selected_count": 0,
    }
    assert efloras.acquisition_status(**common, species_count=0) == "no_species_taxa_discovered"
    assert efloras.acquisition_status(**common, species_count=12) == "no_master_matches"


def test_acquisition_status_does_not_hide_rejected_selected_records() -> None:
    common = {
        "crawl_stop_reason": "queue_exhausted",
        "pending_selected": 0,
        "selected_count": 20,
        "species_count": 50,
    }
    assert (
        efloras.acquisition_status(**common, treatment_count=19, rejected_count=1)
        == "partial_unconfirmed_treatments"
    )
    assert efloras.acquisition_status(**common, treatment_count=20, rejected_count=0) == "ok"


def test_treatment_circuit_breaker_defaults_to_three_exhausted_species_requests() -> None:
    args = efloras.build_parser().parse_args(["--flora-id", "11"])
    assert args.max_consecutive_treatment_failures == 3


def test_inventory_only_mode_is_explicit() -> None:
    assert efloras.build_parser().parse_args(["--flora-id", "2"]).inventory_only is False
    assert efloras.build_parser().parse_args(["--flora-id", "2", "--inventory-only"]).inventory_only is True


def test_inventory_coverage_reports_complete_reusable_inventory() -> None:
    taxa = [
        {"listing_rank": "genus", "original_taxon_name": "Poa", "match_status": "not_species"},
        {
            "listing_rank": "species",
            "original_taxon_name": "Poa annua",
            "match_status": "exact_accepted_name_match",
        },
        {"listing_rank": "species", "original_taxon_name": "Poa alpina", "match_status": "unmatched"},
    ]
    report = efloras.inventory_coverage(
        flora_id=2,
        flora_name="Flora of China",
        landing_url="http://example.test/flora/2",
        taxa=taxa,
        eligible_count=1,
        crawl_metrics={"inventory_complete": True, "inventory_stop_reason": "queue_exhausted"},
        client=SimpleNamespace(http_errors=0, errors=[]),
    )
    assert report["status"] == "inventory_ready"
    assert report["taxa_links_discovered"] == 3
    assert report["species_names_discovered"] == 2
    assert report["exact_master_matches"] == 1
    assert report["unmatched_names"] == 1


def test_combined_coverage_deduplicates_species_across_floras(tmp_path: Path) -> None:
    root = tmp_path / "input"
    for flora_id in (1, 2):
        target = root / f"flora-{flora_id}"
        treatment = {field: "" for field in efloras.TREATMENT_FIELDS}
        treatment.update({"accepted_species": "A a", "flora_id": str(flora_id), "taxon_id": str(flora_id)})
        candidate = {field: "" for field in efloras.CANDIDATE_FIELDS}
        candidate.update(
            {
                "accepted_species": "A a",
                "flora_id": str(flora_id),
                "taxon_id": str(flora_id),
                "trait_name": "flower_primary_color",
                "candidate_value": "white",
                "source_excerpt": "Flowers white.",
            }
        )
        efloras.write_csv(target / "species_treatments.csv", [treatment], efloras.TREATMENT_FIELDS)
        efloras.write_csv(target / "trait_candidates.csv", [candidate], efloras.CANDIDATE_FIELDS)
        efloras.write_json(target / "coverage_report.json", {"status": "ok", "flora_id": flora_id})
    report = efloras.combine_results(root, tmp_path / "combined")
    assert report["species_treatment_rows"] == 2
    assert report["unique_species_excluding_cross_flora_duplicates"] == 1
    assert report["species_repeated_across_floras"] == 1


def test_parser_exposes_resume_and_shard_inputs() -> None:
    parser = efloras.build_parser()
    args = parser.parse_args(["--flora-id", "2", "--shard-index", "1", "--shard-count", "3", "--no-resume"])
    assert isinstance(args, argparse.Namespace)
    assert args.shard_index == 1
    assert args.shard_count == 3
    assert args.resume is False
