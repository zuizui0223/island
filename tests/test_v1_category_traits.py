from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.v1_category_traits import (
    CategoryTraitValidationError,
    OUTPUT_COLUMNS,
    derive_v1_categories,
    ingest_result_csv,
    prepare_prompt_batches,
    render_prompt,
    validate_result_table,
)


def _rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "species": "Campanula punctata",
                "flower_color": "white to pale pink",
                "flower_shape": "actinomorphic, open bell-shaped flowers",
                "pollination_guild": "bumblebees",
                "pollination_notes": "Visited mainly by long-tongued bumblebees.",
                "mating_system": "mixed_mating",
                "self_incompatibility": "likely_SC",
                "evidence_type": "field_study, review",
                "confidence": "high",
            },
            {
                "species": "Salvia splendens",
                "flower_color": "red",
                "flower_shape": "zygomorphic tubular flowers",
                "pollination_guild": "birds",
                "pollination_notes": "Bird-pollinated in much of its native range.",
                "mating_system": "obligate_outcrossing",
                "self_incompatibility": "likely_SI",
                "evidence_type": "flora",
                "confidence": "medium",
            },
            {
                "species": "Unknowna insularis",
                "flower_color": "unknown",
                "flower_shape": "unknown",
                "pollination_guild": "unknown",
                "pollination_notes": "unknown",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "inference",
                "confidence": "low",
            },
        ],
        columns=OUTPUT_COLUMNS,
    )


def test_render_prompt_preserves_exact_species_list() -> None:
    template = "Header\n{{SPECIES_LIST}}\nFooter\n"
    prompt = render_prompt(["Campanula punctata", "Salvia splendens"], template)
    assert prompt == (
        "Header\n- Campanula punctata\n- Salvia splendens\nFooter\n"
    )


def test_prepare_prompt_batches_writes_packets(tmp_path: Path) -> None:
    species_csv = tmp_path / "species.csv"
    pd.DataFrame(
        {
            "accepted_species": [
                "Campanula punctata",
                "Salvia splendens",
                "Campanula punctata",
                "Poa annua",
            ]
        }
    ).to_csv(species_csv, index=False)
    prompt_path = tmp_path / "prompt.txt"
    prompt_path.write_text(
        "Return CSV\n{{SPECIES_LIST}}\n",
        encoding="utf-8",
    )

    report = prepare_prompt_batches(
        species_csv=species_csv,
        output_dir=tmp_path / "packets",
        prompt_path=prompt_path,
        batch_size=2,
    )

    assert report["n_species_selected"] == 3
    assert report["n_batches"] == 2
    first_prompt = (tmp_path / "packets/batch_00001/prompt.txt").read_text(
        encoding="utf-8"
    )
    assert "- Campanula punctata" in first_prompt
    assert "- Salvia splendens" in first_prompt
    expected = pd.read_csv(
        tmp_path / "packets/batch_00002/expected_species.csv"
    )
    assert expected["species"].tolist() == ["Poa annua"]


def test_validate_and_derive_v1_categories() -> None:
    validated = validate_result_table(
        _rows(),
        expected_species=[
            "Campanula punctata",
            "Salvia splendens",
            "Unknowna insularis",
        ],
    )
    derived = derive_v1_categories(validated)

    campanula = derived.loc[
        derived["species"].eq("Campanula punctata")
    ].iloc[0]
    salvia = derived.loc[derived["species"].eq("Salvia splendens")].iloc[0]
    unknown = derived.loc[
        derived["species"].eq("Unknowna insularis")
    ].iloc[0]

    assert campanula["v1_color_binary"] == "conspicuous"
    assert campanula["v1_form_binary"] == "generalized"
    assert campanula["v1_self_compatibility_binary"] == "self_compatible"
    assert campanula["v1_animal_pollinated"] == "yes"
    assert campanula["v1_mating_binary"] == "mixed"
    assert campanula["evidence_type"] == "field_study|review"
    assert campanula["evidence_tier"] == "primary"
    assert bool(campanula["primary_color_eligible"])
    assert bool(campanula["primary_form_eligible"])
    assert bool(campanula["primary_self_eligible"])
    assert bool(campanula["primary_analysis_eligible"])

    assert salvia["v1_color_binary"] == "conspicuous"
    assert salvia["v1_form_binary"] == "specialized"
    assert salvia["v1_self_compatibility_binary"] == "self_incompatible"
    assert bool(salvia["primary_analysis_eligible"])

    assert unknown["v1_color_binary"] == "unknown"
    assert unknown["v1_form_binary"] == "unknown"
    assert unknown["evidence_tier"] == "inference_or_low"
    assert not bool(unknown["primary_color_eligible"])
    assert not bool(unknown["primary_form_eligible"])
    assert not bool(unknown["primary_self_eligible"])
    assert not bool(unknown["primary_analysis_eligible"])


def test_ingest_reports_nonunknown_coverage(tmp_path: Path) -> None:
    result_csv = tmp_path / "result.csv"
    expected_csv = tmp_path / "expected.csv"
    output_csv = tmp_path / "validated.csv"
    _rows().to_csv(result_csv, index=False)
    pd.DataFrame(
        {"species": _rows()["species"].tolist()}
    ).to_csv(expected_csv, index=False)

    report = ingest_result_csv(
        result_csv=result_csv,
        expected_species_csv=expected_csv,
        output_csv=output_csv,
    )

    assert output_csv.exists()
    assert report["n_species"] == 3
    assert report["coverage"] == {
        "flower_color_nonunknown": 2,
        "flower_shape_nonunknown": 2,
        "pollination_guild_nonunknown": 2,
        "mating_system_nonunknown": 2,
        "self_incompatibility_nonunknown": 2,
        "primary_color_eligible": 2,
        "primary_form_eligible": 2,
        "primary_self_eligible": 2,
        "primary_analysis_eligible": 2,
    }


def test_validation_rejects_invalid_enum_and_species_mismatch() -> None:
    invalid = _rows()
    invalid.loc[0, "pollination_guild"] = "beetles"
    with pytest.raises(CategoryTraitValidationError, match="invalid pollination_guild"):
        validate_result_table(invalid)

    with pytest.raises(CategoryTraitValidationError, match="species mismatch"):
        validate_result_table(
            _rows(),
            expected_species=["Campanula punctata", "Salvia splendens"],
        )


def test_unknown_information_requires_low_confidence() -> None:
    invalid = _rows().iloc[[2]].copy()
    invalid.loc[:, "confidence"] = "medium"
    with pytest.raises(CategoryTraitValidationError, match="must have low confidence"):
        validate_result_table(invalid)
