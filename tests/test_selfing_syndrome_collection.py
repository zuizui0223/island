from pathlib import Path

import pandas as pd
import pytest

from island_v2.selfing_syndrome_collection import axis_status, load_config, read_fills


@pytest.fixture()
def config(tmp_path: Path):
    path = tmp_path / "config.yml"
    path.write_text(
        """
axis_definitions:
  colour_vividness:
    target_traits: [flower_primary_color]
  floral_complexity:
    target_traits: [floral_form, floral_symmetry]
  reproductive_assurance:
    target_traits: [self_incompatibility, autonomous_selfing_capacity, mating_system]
primary_completion:
  max_direct_rounds: 2
""",
        encoding="utf-8",
    )
    return load_config(path)


def rows(species, trait, tier, value):
    return {
        "accepted_species": species,
        "trait_name": trait,
        "fill_tier": tier,
        "filled_value": value,
    }


def test_three_direct_axes_complete(config):
    frame = pd.DataFrame(
        [
            rows("Alpha beta", "flower_primary_color", "species_direct", "red"),
            rows("Alpha beta", "floral_form", "species_direct", "tubular"),
            rows("Alpha beta", "self_incompatibility", "species_direct", "self_compatible"),
        ]
    )
    status = axis_status(frame, config, round_index=0).iloc[0]
    assert bool(status["all_three_axes_direct"])
    assert not bool(status["continue_direct_acquisition"])


def test_any_reproductive_trait_completes_reproductive_axis(config):
    frame = pd.DataFrame(
        [
            rows("Alpha beta", "flower_primary_color", "species_direct", "white"),
            rows("Alpha beta", "floral_symmetry", "synonym_direct", "radial"),
            rows("Alpha beta", "mating_system", "species_direct", "mixed"),
        ]
    )
    status = axis_status(frame, config, round_index=0).iloc[0]
    assert bool(status["reproductive_assurance_direct"])
    assert bool(status["all_three_axes_direct"])


def test_missing_reproductive_axis_is_highest_priority(config):
    frame = pd.DataFrame(
        [
            rows("Alpha beta", "flower_primary_color", "species_direct", "red"),
            rows("Alpha beta", "floral_form", "species_direct", "rotate"),
            rows("Alpha beta", "self_incompatibility", "unresolved_no_evidence", "unresolved"),
        ]
    )
    status = axis_status(frame, config, round_index=0).iloc[0]
    assert status["missing_direct_axes"] == "reproductive_assurance"
    assert status["collection_priority"] == 0
    assert bool(status["continue_direct_acquisition"])


def test_genus_then_family_resolve_but_do_not_count_as_direct(config):
    frame = pd.DataFrame(
        [
            rows("Alpha beta", "flower_primary_color", "species_direct", "red"),
            rows("Alpha beta", "floral_form", "genus_inference", "tubular"),
            rows("Alpha beta", "mating_system", "family_inference", "mixed"),
        ]
    )
    status = axis_status(frame, config, round_index=2).iloc[0]
    assert not bool(status["all_three_axes_direct"])
    assert bool(status["all_three_axes_resolved"])
    assert bool(status["ready_for_taxonomic_fallback"])
    assert not bool(status["continue_direct_acquisition"])


def test_missing_family_evidence_remains_unresolved(config):
    frame = pd.DataFrame(
        [
            rows("Alpha beta", "flower_primary_color", "species_direct", "red"),
            rows("Alpha beta", "floral_form", "genus_inference", "tubular"),
            rows("Alpha beta", "mating_system", "unresolved_no_evidence", "unresolved"),
        ]
    )
    status = axis_status(frame, config, round_index=2).iloc[0]
    assert status["unresolved_axes"] == "reproductive_assurance"
    assert not bool(status["all_three_axes_resolved"])


def test_global_fallback_is_rejected(tmp_path: Path):
    root = tmp_path / "fills" / "shards" / "shard_000"
    root.mkdir(parents=True)
    pd.DataFrame(
        [rows("Alpha beta", "flower_primary_color", "global_fallback", "red")]
    ).to_csv(root / "fills.csv.gz", index=False)
    with pytest.raises(Exception, match="global or unsupported fallback"):
        read_fills(tmp_path / "fills")
