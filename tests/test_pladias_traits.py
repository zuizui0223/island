from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from island_v2.integrated_trait_coverage import load_generic_bulk_candidates
from island_v2.pladias_traits import (
    CANDIDATE_COLUMNS,
    _normalise_cleistogamy,
    _normalise_colour,
    _normalise_floral_form,
    _normalise_inflorescence,
    _normalise_pollination,
    _normalise_reproduction,
    _parse_page,
    _stable_id,
    prepare_pladias,
)


def _html(species: str, family: str) -> str:
    return f"""
    <html><head><title>{species} - Pladias</title></head><body>
      <div class="container content"><h1>{species}</h1></div>
      <ol class="breadcrumb">
        <li><a href="/en/taxon/overview/{family}">{family}</a></li>
        <li><a href="/en/taxon/overview/{species.split()[0]}">{species.split()[0]}</a></li>
      </ol>
      <li class="featureDetail">
        <div><div><span class="font-italic">Shape of the sympetalous corolla or syntepalous perianth</span>: <b>
        tubular, campanulate</b></div><span data-target="#modal_feature_333">?</span></div>
        <div class="modal-body">Method. Data source and citation Corolla source, 2018.</div>
      </li>
      <li class="featureDetail">
        <div><div><span class="font-italic">Flower colour</span>: <b>
        yellow-white, violet</b></div><span data-target="#modal_feature_319">?</span></div>
        <div class="modal-body">Method. Data source and citation Colour source, 2019.</div>
      </li>
      <li class="featureDetail">
        <div><div><span class="font-italic">Flower symmetry</span>: <b>
        zygomorphic</b></div><span data-target="#modal_feature_9">?</span></div>
        <div class="modal-body">Method. Data source and citation Symmetry source, 2018.</div>
      </li>
      <li class="featureDetail">
        <div><div><span class="font-italic">Inflorescence type</span>: <b>
        flores solitarii, racemus</b></div><span data-target="#modal_feature_323">?</span></div>
        <div class="modal-body">Method. Data source and citation Display source, 2018.</div>
      </li>
      <li class="featureDetail">
        <div><div><span class="font-italic">Generative reproduction type</span>: <b>
        allogamy self-incompatibility, facultative allogamy</b></div>
        <span data-target="#modal_feature_145">?</span></div>
        <div class="modal-body">Method. Data source and citation Chrtek, 2018.</div>
      </li>
      <li class="featureDetail">
        <div><div><span class="font-italic">Pollination syndrome</span>: <b>
        insect-pollination</b></div><span data-target="#modal_feature_999">?</span></div>
      </li>
    </body></html>
    """


def _config(path: Path, taxon_list: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "source": {
                    "name": "pladias_czech_flora_2025",
                    "citation": "Pladias test citation",
                    "website": "https://pladias.cz",
                    "usage_terms_url": "https://pladias.cz/en/download/features",
                    "robots_url": "https://pladias.cz/robots.txt",
                    "taxon_list_url": "https://example.test/taxa.xls",
                    "expected_taxon_list_sha256": hashlib.sha256(
                        taxon_list.read_bytes()
                    ).hexdigest(),
                    "species_page_template": "https://pladias.cz/en/taxon/data/{species}",
                    "user_agent": "test-agent",
                    "accepted_feature_labels": [
                        "Flower colour",
                        "Flower symmetry",
                        "Inflorescence type",
                        "Shape of the sympetalous corolla or syntepalous perianth",
                        "Generative reproduction type",
                        "Pollination syndrome",
                    ],
                }
            }
        ),
        encoding="utf-8",
    )


def test_trait_specific_mapping_preserves_multistate_and_category_boundaries() -> None:
    assert json.loads(_normalise_colour("yellow-white, violet")) == [
        "blue_purple",
        "white",
        "yellow_orange",
    ]
    assert json.loads(_normalise_inflorescence("flores solitarii, racemus")) == [
        "raceme_spike_panicle",
        "solitary",
    ]
    assert json.loads(
        _normalise_inflorescence("corymbothyrsus ex anthodiis compositus")
    ) == ["composite_display", "umbel_corymb"]
    assert _normalise_inflorescence("flores solitarii, dichasium") == ""
    assert json.loads(_normalise_floral_form("tubular, campanulate")) == [
        "bell_campanulate",
        "tubular",
    ]
    assert _normalise_pollination("wind-pollination, insect-pollination, selfing") == (
        "mixed"
    )
    assert _normalise_pollination("selfing") == ""
    assert _normalise_pollination("wind-pollination, water-pollination") == ""
    assert _normalise_cleistogamy("cleistogamy") == "obligate"
    assert _normalise_cleistogamy("insect-pollination, selfing, cleistogamy") == (
        "facultative"
    )
    assert _normalise_cleistogamy("pseudocleistogamy") == ""
    assert _normalise_reproduction(
        "allogamy self-incompatibility, facultative allogamy"
    ) == [
        ("self_incompatibility", "SI"),
        ("mating_system", "predominantly_outcrossing"),
    ]
    assert _normalise_reproduction("facultative autogamy") == [
        ("mating_system", "predominantly_selfing")
    ]
    assert _normalise_reproduction("apomixis") == []
    assert _normalise_reproduction(
        "allogamy self-incompatibility, facultative apomixis"
    ) == [("self_incompatibility", "SI")]


def test_page_parser_retains_exact_quote_feature_citation_and_identity() -> None:
    parsed = _parse_page(_html("Alpha one", "FamA"))
    assert parsed["matched_page_name"] == "Alpha one"
    assert "FamA" in parsed["breadcrumbs"]
    colour = next(row for row in parsed["features"] if row["label"] == "Flower colour")
    assert colour["raw_value"] == "yellow-white, violet"
    assert colour["feature_id"] == "319"
    assert colour["feature_citation"] == "Colour source, 2019."
    assert colour["supporting_excerpt"] == "Flower colour: yellow-white, violet"


def test_prepare_uses_cache_exact_identity_and_emits_no_cross_trait_proxy(
    tmp_path: Path,
) -> None:
    taxon_list = tmp_path / "taxa.xlsx"
    pd.DataFrame(
        {"NameLatFinal": ["Alpha one", "Beta two", "Alpha one subsp. extra"]}
    ).to_excel(taxon_list, sheet_name="TaxonList", index=False)
    config = tmp_path / "config.yml"
    _config(config, taxon_list)
    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            ("Alpha one", "Alpha", "FamA"),
            ("Beta two", "Beta", "FamB"),
            ("Gamma three", "Gamma", "FamC"),
        ],
        columns=["accepted_species", "genus", "family"],
    ).to_csv(master, index=False)

    cache = tmp_path / "cache"
    cache.mkdir()
    for species, family in [("Alpha one", "FamA"), ("Beta two", "WrongFam")]:
        parsed = _parse_page(_html(species, family))
        record = {
            "accepted_species": species,
            "family": "FamA" if species == "Alpha one" else "FamB",
            "source_url": f"https://pladias.cz/en/taxon/data/{species.replace(' ', '%20')}",
            "http_status": 200,
            "matched_page_name": parsed["matched_page_name"],
            "page_family": family,
            "identity_status": (
                "accepted_name_and_family_exact"
                if species == "Alpha one"
                else "reject_family_conflict_or_missing"
            ),
            "page_title": parsed["page_title"],
            "breadcrumbs": parsed["breadcrumbs"],
            "features": parsed["features"],
            "retrieved_at_utc": "2026-08-08T00:00:00+00:00",
            "content_sha256": hashlib.sha256(_html(species, family).encode()).hexdigest(),
            "elapsed_seconds": 0.1,
            "cache_status": "fetched",
            "error": "",
        }
        (cache / f"{_stable_id(species)}.json").write_text(
            json.dumps(record), encoding="utf-8"
        )

    output = tmp_path / "output"
    manifest = prepare_pladias(
        taxon_list_xls=taxon_list,
        master_csv=master,
        output_dir=output,
        config_path=config,
        cache_dir=cache,
        workers=2,
    )
    candidates = pd.read_csv(output / "trait_candidates.csv.gz", dtype=str).fillna("")
    exclusions = pd.read_csv(output / "exclusions.csv.gz", dtype=str).fillna("")
    quality_audit = pd.read_csv(output / "quality_audit.csv.gz", dtype=str).fillna("")

    assert list(candidates.columns) == CANDIDATE_COLUMNS
    assert candidates["accepted_species"].eq("Alpha one").all()
    assert set(candidates["trait_name"]) == {
        "flower_primary_color",
        "floral_symmetry",
        "inflorescence_display",
        "floral_form",
        "self_incompatibility",
        "mating_system",
        "pollen_vector_mode",
    }
    assert "reward_type" not in set(candidates["trait_name"])
    assert candidates["confidence"].eq("high").all()
    assert candidates["source_lineage"].str.startswith("pladias:feature:").all()
    assert exclusions.loc[exclusions["accepted_species"].eq("Beta two"), "reason"].item() == (
        "reject_family_conflict_or_missing"
    )
    assert manifest["counts"]["pages"] == 2
    assert manifest["counts"]["candidate_species"] == 1
    assert manifest["quality_contract"]["pollen_vector_or_reward_projected_to_strict_axes"] is False
    assert manifest["quality_contract"]["genus_inference_emitted"] is False
    assert len(quality_audit) == 1
    assert quality_audit.loc[0, "accepted_correct"] == "True"
    assert manifest["quality_contract"]["audit"]["precision"] == 1.0
    assert manifest["quality_contract"]["audit"]["cultivar_contamination_rate"] == 0.0

    evidence, integration_audit = load_generic_bulk_candidates(
        output / "trait_candidates.csv.gz", {"sources": []}
    )
    assert set(integration_audit["reason"]) == {
        "accepted_by_strict_source_contract",
        "trait_outside_three_axes",
    }
    assert evidence["source_group"].eq("pladias").all()
    assert "pollen_vector_mode" not in set(evidence["trait_name"])
    assert evidence["source_lineage"].str.startswith("pladias:feature:").all()
