from pathlib import Path

import pandas as pd

from island_v2.trait_pilot_eval import (
    build_stratified_pilot,
    evaluate_pilot,
    load_config,
)

CONFIG = load_config(Path("config/pilot_stratification.yml"))


def _candidates() -> pd.DataFrame:
    # One large stratum (many records) and several small ones. A record-count
    # selection would flood the pilot with the big stratum; stratified selection
    # must spread across strata instead.
    rows = []
    for i in range(20):
        rows.append(
            {
                "accepted_species": f"Bigfamilia sp{i:02d}",
                "native_status": "introduced",
                "growth_form": "herb",
                "family": "Asteraceae",
                "flower_conspicuousness": "conspicuous",
            }
        )
    extras = [
        ("Rareone woodya", "native", "woody", "Myrtaceae", "inconspicuous"),
        ("Rareone vineb", "native", "vine", "Convolvulaceae", "conspicuous"),
        ("Rareone epic", "native", "epiphyte", "Orchidaceae", "conspicuous"),
        ("Rareone windd", "native", "herb", "Poaceae", "wind_pollination_candidate"),
    ]
    for sp, ns, gf, fam, con in extras:
        rows.append(
            {
                "accepted_species": sp,
                "native_status": ns,
                "growth_form": gf,
                "family": fam,
                "flower_conspicuousness": con,
            }
        )
    return pd.DataFrame(rows)


def test_build_spreads_across_strata_not_record_count():
    pilot, summary = build_stratified_pilot(_candidates(), CONFIG)
    # All four rare strata must be represented, not drowned by the big one.
    assert summary["n_strata_occupied"] == 5
    species = set(pilot["accepted_species"])
    assert {"Rareone woodya", "Rareone vineb", "Rareone epic", "Rareone windd"} <= species
    # The big stratum should not exceed the whole pilot.
    n_big = sum(1 for s in species if s.startswith("Bigfamilia"))
    assert n_big < len(species)


def test_build_respects_target_max():
    # 500 species in one stratum, target_max caps the pilot.
    big = pd.DataFrame(
        {
            "accepted_species": [f"Sp{i:03d}" for i in range(500)],
            "native_status": "native",
            "growth_form": "herb",
            "family": "Fabaceae",
            "flower_conspicuousness": "conspicuous",
        }
    )
    pilot, summary = build_stratified_pilot(big, CONFIG)
    assert summary["n_selected"] == CONFIG["target_max"]
    assert len(pilot) == CONFIG["target_max"]


def test_build_is_deterministic():
    a, _ = build_stratified_pilot(_candidates(), CONFIG)
    b, _ = build_stratified_pilot(_candidates(), CONFIG)
    pd.testing.assert_frame_equal(a, b)


def _review() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": ["A", "A", "B", "B", "C", "C"],
            "literature_relevant": ["relevant", "irrelevant", "irrelevant", "relevant", "unknown", ""],
            "source_to_trait_relevant": ["yes", "no", "no", "yes", "", ""],
            "evidence_scope": ["species_direct", "", "genus_inference", "species_direct", "", ""],
            "trait_resolved": ["resolved", "", "unresolved", "resolved", "unresolved", ""],
            "hierarchical_transfer_correct": ["", "", "wrong", "", "correct", ""],
        }
    )


def test_evaluate_metrics_count_only_judgements():
    report = evaluate_pilot(_review(), CONFIG)
    m = report["metrics"]
    # 2 irrelevant of 4 relevant/irrelevant judgements
    assert m["irrelevant_literature_rate"]["rate"] == 0.5
    assert m["irrelevant_literature_rate"]["denominator"] == 4
    # 2 species_direct of 3 non-blank scopes
    assert m["species_direct_evidence_proportion"]["numerator"] == 2
    assert m["species_direct_evidence_proportion"]["denominator"] == 3
    # 2 unresolved of 4 resolved/unresolved judgements
    assert m["trait_specific_unresolved_rate"]["numerator"] == 2
    assert m["trait_specific_unresolved_rate"]["denominator"] == 4
    # 1 wrong of 2 transfer judgements
    assert m["wrong_hierarchical_transfer_rate"]["rate"] == 0.5


def test_evaluate_zero_denominator_is_null_not_fabricated():
    empty = pd.DataFrame({"accepted_species": ["A"], "literature_relevant": [""]})
    report = evaluate_pilot(empty, CONFIG)
    assert report["metrics"]["irrelevant_literature_rate"]["rate"] is None
    assert report["metrics"]["irrelevant_literature_rate"]["denominator"] == 0


def test_irrelevant_literature_rate_is_primary():
    assert "irrelevant_literature_rate" in CONFIG["primary_metrics"]
