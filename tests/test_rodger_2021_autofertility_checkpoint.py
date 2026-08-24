from __future__ import annotations

import pandas as pd

from island_v2.rodger_2021_autofertility_checkpoint import (
    AUDIT_ROWS,
    AUTO_COLUMNS,
    build_checkpoint,
)


def epithet(index: int) -> str:
    return (
        "species"
        f"{chr(97 + index // (26 * 26))}"
        f"{chr(97 + (index // 26) % 26)}"
        f"{chr(97 + index % 26)}"
    )


def source_rows(n: int = 240) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    datasets = ("glopl", "konstanz", "stellenbosch")
    for index in range(n):
        species = f"Testus {epithet(index)}"
        row = {
            "source_row": str(index + 2),
            "pc.unique.number": str(index + 1),
            "source.dataset": datasets[index % len(datasets)],
            "study": f"Author {index % 40}, Journal, {2000 + index % 20}",
            "year.publication": str(2000 + index % 20),
            "year.study.estimate": str(1999 + index % 20),
            "taxon": species.replace(" ", "_"),
            "genus.species": species.replace(" ", "_"),
            "TPL.family": "Leguminosae",
            "alien": "native",
        }
        row.update({column: "" for column in AUTO_COLUMNS})
        row["auto.fs.x"] = "0" if index % 2 else "0.25"
        rows.append(row)
    return pd.DataFrame(rows)


def master(n: int = 240) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": [f"Testus {epithet(index)}" for index in range(n)],
            "family": ["Fabaceae"] * n,
        }
    )


def test_checkpoint_uses_only_pollinator_exclusion_and_skips_completed_pairs() -> None:
    baseline = pd.DataFrame(
        {
            "accepted_species": [f"Testus {epithet(5)}"],
            "trait_name": ["autonomous_selfing_capacity"],
        }
    )
    evidence, audit, selection = build_checkpoint(source_rows(), baseline, master())

    assert len(evidence) == 239
    assert len(audit) == AUDIT_ROWS
    assert set(evidence["trait_name"]) == {"autonomous_selfing_capacity"}
    assert set(evidence["axis"]) == {"reproductive_assurance"}
    assert set(evidence["quality"]) == {"high"}
    assert set(evidence["normalized_value"]) == {"autonomous", "absent"}
    assert f"Testus {epithet(5)}" not in set(evidence["accepted_species"])
    assert selection["selection_reason"].value_counts().to_dict() == {
        "selected": 239,
        "already_resolved_direct_pair": 1,
    }
    assert audit["accepted_correct"].eq("true").all()
    assert audit[["source_dataset", "normalized_value"]].drop_duplicates().shape[0] == 6


def test_checkpoint_rejects_family_conflicts_and_invalid_measurements() -> None:
    source = source_rows()
    source.loc[0, "TPL.family"] = "Wrongaceae"
    source.loc[1, "auto.fs.x"] = "not-a-number"
    evidence, _, selection = build_checkpoint(
        source,
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
        master(),
    )

    assert f"Testus {epithet(0)}" not in set(evidence["accepted_species"])
    assert f"Testus {epithet(1)}" not in set(evidence["accepted_species"])
    assert selection["selection_reason"].value_counts().to_dict() == {
        "selected": 238,
        "family_conflict": 1,
        "missing_or_invalid_pollinator_exclusion_measure": 1,
    }


def test_same_original_study_keeps_one_lineage_with_distinct_records() -> None:
    evidence, _, _ = build_checkpoint(
        source_rows(),
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
        master(),
    )
    shared = evidence.loc[evidence["source_lineage"].eq(evidence.iloc[0]["source_lineage"])]

    assert len(shared) > 1
    assert shared["source_record_id"].nunique() == len(shared)
    assert shared["lineage_method"].eq(
        "original_study_citation_not_database_redistributor"
    ).all()


def test_two_backbone_synonym_is_direct_but_one_backbone_row_is_rejected() -> None:
    source = source_rows()
    source.loc[0, "genus.species"] = "Oldus_speciesaaa"
    source.loc[0, "TPL.family"] = "Fabaceae"
    source.loc[1, "genus.species"] = "Oldus_speciesaab"
    source.loc[1, "TPL.family"] = "Fabaceae"
    synonyms = pd.DataFrame(
        {
            "submitted_name": [
                "Oldus speciesaaa",
                "Oldus speciesaab",
                "Unused oldname",
            ],
            "accepted_species": [
                f"Testus {epithet(0)}",
                f"Testus {epithet(1)}",
                f"Testus {epithet(2)}",
            ],
            "family": ["Fabaceae", "Fabaceae", "Fabaceae"],
            "gbif_species_key": ["1", "2", "3"],
            "gbif_confidence": [100, 90, 100],
            "gbif_match_type": ["EXACT", "EXACT", "EXACT"],
            "wfo_accepted_taxon_id": ["wfo-1", "wfo-2", "wfo-3"],
            "resolution_method": [
                "strict_exact_two_backbone_agreement",
                "strict_exact_two_backbone_agreement",
                "one_backbone_only",
            ],
        }
    )
    evidence, _, _ = build_checkpoint(
        source,
        pd.DataFrame(columns=["accepted_species", "trait_name"]),
        master(),
        synonyms,
    )
    row = evidence.loc[evidence["accepted_species"].eq(f"Testus {epithet(0)}")].iloc[0]

    assert row["evidence_scope"] == "synonym_direct"
    assert row["name_match_method"] == "synonym_exact_two_backbone"
    assert f"Testus {epithet(1)}" not in set(evidence["accepted_species"])
