"""The shared matching core, tested where the two lanes cannot see it.

Everything about *reading* a floral statement is exercised through the two lane
suites, which is where the sentences that motivated each rule live. What is
tested here is the machinery that makes sharing possible at all: the config
inheritance both lanes depend on, and the folding contract the offset map rests
on.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import typer

from island_v2 import floral_text_matching as matching

VOCABULARY = Path("config/floral_vocabulary.yml")


def _write(tmp_path: Path, name: str, body: str) -> Path:
    path = tmp_path / name
    path.write_text(body, encoding="utf-8")
    return path


def test_a_lane_config_inherits_the_shared_vocabulary(tmp_path: Path):
    _write(tmp_path, "base.yml", "colour_terms:\n  white: white\nwindow: 60\n")
    lane = _write(tmp_path, "lane.yml", "extends: base.yml\nlane_only: yes\n")

    config = matching.load_config(lane, {"colour_terms", "lane_only"})
    assert config["colour_terms"] == {"white": "white"}
    assert config["window"] == 60
    assert "extends" not in config


def test_a_lane_key_overrides_the_one_it_inherits(tmp_path: Path):
    # This is what lets the protologue lane keep a 60-character organ window
    # against the label lane's 40 without either owning a private vocabulary.
    _write(tmp_path, "base.yml", "window: 60\ncolour_terms:\n  white: white\n")
    lane = _write(tmp_path, "lane.yml", "extends: base.yml\nwindow: 40\n")

    assert matching.load_config(lane, {"window"})["window"] == 40


def test_the_base_is_resolved_next_to_the_config_that_names_it(tmp_path: Path):
    # Resolved relative to the lane config rather than the working directory,
    # so a job that runs from anywhere still finds it.
    nested = tmp_path / "config"
    nested.mkdir()
    _write(nested, "base.yml", "colour_terms:\n  white: white\n")
    lane = _write(nested, "lane.yml", "extends: base.yml\n")

    assert matching.load_config(lane, {"colour_terms"})["colour_terms"]


def test_a_config_missing_a_required_key_fails_closed(tmp_path: Path):
    lane = _write(tmp_path, "lane.yml", "colour_terms:\n  white: white\n")
    with pytest.raises(typer.BadParameter):
        matching.load_config(lane, {"colour_terms", "negation_markers"})


def test_a_config_that_is_not_a_mapping_fails_closed(tmp_path: Path):
    lane = _write(tmp_path, "lane.yml", "- one\n- two\n")
    with pytest.raises(typer.BadParameter):
        matching.load_config(lane, set())


def test_the_shared_vocabulary_stands_alone():
    # It is data, not a lane: nothing in it may depend on being merged into one.
    config = matching.load_config(VOCABULARY, {"colour_terms", "floral_organ_terms"})
    assert set(matching.expand_colour_terms(config).values()) <= {
        "white",
        "green_brown_inconspicuous",
        "yellow_orange",
        "red_pink",
        "blue_purple",
    }
    assert matching.plain_colour_values(config) == {
        "white",
        "green_brown_inconspicuous",
    }


def test_every_vocabulary_term_survives_folding_intact():
    """A term that folds to something else can never match itself.

    The dakuten bug was exactly this: `ピンク色` was stored in the vocabulary and
    folded to `ヒンク色`, so the ledger recorded a term that appears nowhere on
    the page. Terms are written folded by convention, and this is what enforces
    it.
    """
    config = matching.load_config(VOCABULARY, {"colour_terms"})
    mangled = [
        term
        for term in config["colour_terms"]
        if matching.fold(str(term))[0] != str(term)
    ]
    assert not mangled, mangled


def test_folding_is_injective_across_the_vocabulary():
    # Two terms folding together would make one of them unreachable and the
    # other ambiguous. Checked over the expanded vocabulary, Latin included.
    config = matching.load_config(VOCABULARY, {"colour_terms"})
    collisions = {}
    for term, value in config["colour_terms"].items():
        folded = matching.fold(str(term))[0]
        collisions.setdefault(folded, set()).add(str(value))
    disagreeing = {k: v for k, v in collisions.items() if len(v) > 1}
    assert not disagreeing, disagreeing


def test_the_offset_map_is_one_entry_per_folded_character():
    for text in ("Blüten weiß", "花はピンク色", "ﾋﾟﾝｸ", "цветки белые", "corolla alba"):
        folded, origin = matching.fold(text)
        assert len(folded) == len(origin), text
        assert all(0 <= index < len(text) for index in origin), text
        assert origin == sorted(origin), text
