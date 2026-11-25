"""End-to-end, hand-calculated checks for the floral-profile analysis."""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest


def test_cli_counts_species_profiles_and_removes_shared_species(tmp_path):
    """Catch counting records as species or removing shared profiles instead of species."""
    profiles = pd.DataFrame(
        {
            "accepted_species": ["Alpha a", "Alpha b", "Alpha c", "Beta d", "Alpha", "Alpha sp."],
            "flower_primary_color": ["white", "white", "red_pink", "yellow_orange", "blue_purple", "white"],
            "floral_form": ["open_radial", "open_radial", "tubular", "bell_campanulate", "spurred", "spurred"],
        }
    )
    profiles.to_csv(tmp_path / "profiles.csv", index=False)
    pd.DataFrame(
        {
            "island_id": ["one", "one", "one", "two", "two", "excluded"],
            "accepted_species": ["Alpha a", "Alpha b", "Alpha c", "Alpha a", "Beta d", "Alpha a"],
        }
    ).to_csv(tmp_path / "flora.csv", index=False)
    pd.DataFrame(
        {
            "island_id": ["one", "two"],
            "analysis_regime": ["test_region", "test_region"],
            "island_latitude": [10.0, 20.0],
            "island_longitude": [120.0, 140.0],
        }
    ).to_csv(tmp_path / "geography.csv", index=False)
    output = tmp_path / "output"
    result = subprocess.run(
        [
            sys.executable, "-m", "island_v2.floral_space",
            "--profiles", str(tmp_path / "profiles.csv"),
            "--flora", str(tmp_path / "flora.csv"),
            "--geography", str(tmp_path / "geography.csv"),
            "--output", str(output),
            "--rarefaction-size", "2", "--min-species", "1", "--pair-count", "1",
        ],
        capture_output=True, text=True,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    assert result.returncode == 0, result.stderr
    metrics = pd.read_csv(output / "island_metrics.csv").set_index("island_id")
    assert set(metrics.index) == {"one", "two"}
    assert metrics.loc["one", "n_species"] == 3
    assert metrics.loc["one", "n_profiles"] == 2
    assert metrics.loc["one", "singleton_fraction"] == pytest.approx(0.5)
    assert metrics.loc["one", "rarefied_profiles"] == pytest.approx(5 / 3)
    assert metrics.loc["one", "expected_profiles_genus_pool"] == pytest.approx(2)
    assert metrics.loc["two", "expected_profiles_genus_pool"] == pytest.approx(2)
    pairs = pd.read_csv(output / "pairwise_comparisons.csv")
    assert len(pairs) == 1
    assert pairs.loc[0, "species_jaccard"] == pytest.approx(1 / 4)
    assert pairs.loc[0, "profile_jaccard"] == pytest.approx(1 / 3)
    assert pairs.loc[0, "profile_jaccard_without_shared_species"] == 0
    summary = json.loads((output / "summary.json").read_text())
    assert summary["n_profiles"] == 3
    assert summary["n_complete_species"] == 4
    assert summary["geography_excluded_flora_rows"] == 1
    assert summary["pcoa"]["correction"] == "sqrt"
    assert summary["n_profile_population_species"] == 4
    assert summary["n_excluded_uninomial_names"] == 1
    assert summary["n_excluded_unresolved_names"] == 1


def test_cli_regional_genus_expectation_is_not_global_pool(tmp_path):
    """Catch replacing genus-conditioned occupancy with unconditioned richness."""
    pd.DataFrame(
        {
            "accepted_species": ["Alpha a", "Alpha b", "Beta c", "Beta d"],
            "flower_primary_color": ["white", "red_pink", "white", "yellow_orange"],
            "floral_form": ["open_radial"] * 4,
        }
    ).to_csv(tmp_path / "profiles.csv", index=False)
    pd.DataFrame(
        {
            "island_id": ["target", "target", "pool", "pool", "pool", "pool"],
            "accepted_species": ["Alpha a", "Beta c", "Alpha a", "Alpha b", "Beta c", "Beta d"],
        }
    ).to_csv(tmp_path / "flora.csv", index=False)
    pd.DataFrame(
        {"island_id": ["target", "pool"], "analysis_regime": ["r", "r"]}
    ).to_csv(tmp_path / "geography.csv", index=False)
    result = subprocess.run(
        [
            sys.executable, "-m", "island_v2.floral_space",
            "--profiles", str(tmp_path / "profiles.csv"),
            "--flora", str(tmp_path / "flora.csv"),
            "--geography", str(tmp_path / "geography.csv"),
            "--output", str(tmp_path / "out"),
            "--min-species", "1", "--rarefaction-size", "2", "--pair-count", "0",
        ], capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"},
    )
    assert result.returncode == 0, result.stderr
    metrics = pd.read_csv(tmp_path / "out/island_metrics.csv").set_index("island_id")
    assert metrics.loc["target", "n_profiles"] == 1
    assert metrics.loc["target", "expected_profiles_regional_pool"] == pytest.approx(11 / 6)
    assert metrics.loc["target", "expected_profiles_genus_pool"] == pytest.approx(7 / 4)
    assert metrics.loc["target", "expected_singleton_profiles_genus_pool"] == pytest.approx(3 / 2)


def test_common_reference_space_does_not_inflate_scope_profile_count(tmp_path):
    """A common All atlas must not label unobserved types as Direct observations."""
    from island_v2.floral_space import analyze

    reference = pd.DataFrame({
        "accepted_species": ["Alpha a", "Alpha b", "Beta c"],
        "flower_primary_color": ["white", "red_pink", "yellow_orange"],
        "floral_form": ["open_radial", "tubular", "bell_campanulate"],
    })
    flora = pd.DataFrame({"island_id": ["one", "one"], "accepted_species": ["Alpha a", "Alpha b"]})
    geography = pd.DataFrame({"island_id": ["one"], "analysis_regime": ["r"]})
    result = analyze(reference.iloc[:2], flora, geography, tmp_path,
                     space_profiles=reference, min_species=1, rarefaction_size=2, pair_count=0)
    assert result["n_profiles"] == 2
    assert result["n_reference_profiles"] == 3
