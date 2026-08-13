from __future__ import annotations

import pandas as pd

from island_v2.cached_evidence_recovery_checkpoint import (
    _colour_is_supported,
    _whole_flower_size_is_supported,
    build_audit,
    build_individually_reviewed_batch,
    select_recoverable_evidence,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS


def _source_row(
    species: str,
    trait: str,
    value: str,
    record: str,
    *,
    provider: str = "official_flora",
    excerpt: str = "Flowers white.",
) -> dict[str, str]:
    row = {column: "" for column in EVIDENCE_COLUMNS}
    row.update(
        {
            "accepted_species": species,
            "axis": "flower_colour",
            "trait_name": trait,
            "normalized_value": value,
            "quality": "high",
            "source_group": "completed_acquisition",
            "source_provider": provider,
            "source_url": f"https://flora.example/{record}",
            "source_record_id": record,
            "source_citation": "Official flora treatment",
            "source_excerpt": excerpt,
            "evidence_scope": "species_direct",
            "name_match_method": "accepted_name_exact",
            "source_lineage": f"treatment:{record}",
            "lineage_method": "provider_treatment",
            "source_file": "pinned.csv.gz",
        }
    )
    return row


def _strict(species: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": species,
            "axis": ["flower_colour"] * len(species),
            "quality": [""] * len(species),
        }
    )


def test_selection_skips_completed_pairs_and_quarantines_conflicts() -> None:
    sources = pd.DataFrame(
        [
            _source_row("Alpha one", "flower_primary_color", "white", "a1"),
            _source_row("Alpha two", "flower_primary_color", "red", "a2"),
            _source_row("Alpha two", "flower_primary_color", "blue", "a3"),
            _source_row("Alpha three", "flower_primary_color", "white", "a4"),
        ]
    )
    prior = pd.DataFrame(
        {"accepted_species": ["Alpha three"], "trait_name": ["flower_primary_color"]}
    )

    selected = select_recoverable_evidence(
        _strict(["Alpha one", "Alpha two", "Alpha three"]), [sources], prior
    )

    assert list(selected["accepted_species"]) == ["Alpha one"]
    assert selected.iloc[0]["acceptance_contract"].endswith("_v1")


def test_audit_guarantees_ten_reviews_for_each_eligible_trait() -> None:
    rows = []
    for index in range(220):
        trait = "flower_primary_color" if index < 205 else "mating_system"
        excerpt = "Flowers white." if trait == "flower_primary_color" else "Mixed mating system."
        rows.append(_source_row(f"Species {index}", trait, "white", f"r{index}", excerpt=excerpt))
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)

    audit = build_audit(evidence)

    assert len(audit) == 200
    assert audit.groupby("trait_name").size().to_dict()["mating_system"] >= 10
    assert set(audit["candidate_id"]).issubset(set(evidence["source_record_id"]))


def test_individual_batch_is_exactly_one_hundred_and_strictly_direct() -> None:
    evidence = pd.DataFrame(
        [
            _source_row(
                f"Species {index}",
                "flower_primary_color",
                "white",
                f"r{index}",
            )
            for index in range(120)
        ],
        columns=EVIDENCE_COLUMNS,
    )

    curated = build_individually_reviewed_batch(evidence)

    assert len(curated) == 100
    assert set(curated["evidence_scope"]) == {"species_direct"}
    assert set(curated["name_match_method"]) == {"accepted_name_exact"}
    assert curated["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()


def test_colour_requires_each_state_on_a_flower_part() -> None:
    valid = pd.Series(
        {"normalized_value": "white|red_pink", "source_excerpt": "Petals white or pink."}
    )
    fruit_leakage = pd.Series(
        {
            "normalized_value": "green_brown_inconspicuous|red_pink",
            "source_excerpt": "Flowers green. Fruit bright red.",
        }
    )

    assert _colour_is_supported(valid)
    assert not _colour_is_supported(fruit_leakage)


def test_size_rejects_flower_part_and_hair_measurements() -> None:
    assert _whole_flower_size_is_supported("Flowers 8-10 mm in diameter.")
    assert _whole_flower_size_is_supported("Corolla 3 mm long; tube cylindrical.")
    assert not _whole_flower_size_is_supported(
        "Flowers rusty pubescent; hairs 0.1-0.2 mm long."
    )
    assert not _whole_flower_size_is_supported("Flowering stems 4-12 cm long.")
