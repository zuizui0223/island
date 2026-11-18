"""Matched membership selection and a real, reproducible pair-analysis CLI."""
from __future__ import annotations

from collections import Counter
import hashlib
import importlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

from island_v2.floral_space_null import REGIMES, SCOPES


def module():
    return importlib.import_module("island_v2.floral_space_pair_analysis")


def membership_design(members, *, regime=REGIMES[0], reverse=False):
    names = sorted(set().union(*map(set, members.values())))
    islands = sorted(members, reverse=reverse)
    rows = [(i, names.index(taxon)) for i, island in enumerate(islands)
            for taxon in sorted(members[island], reverse=reverse)]
    if reverse:
        rows.reverse()
    return SimpleNamespace(analysis_regime=regime, island_ids=tuple(islands),
                           member_island_indices=np.array([r[0] for r in rows]),
                           member_taxon_indices=np.array([r[1] for r in rows]))


def nonnested_designs(reverse=False):
    return (membership_design({"A": {"s", "x", "y"}, "B": {"s", "x", "z"},
                               "C": {"q", "r", "t"}, "D": {"v"}}, reverse=reverse),
            membership_design({"A": {"a", "b", "c"}, "B": {"a", "d", "e"},
                               "C": {"a", "d", "f"}, "E": {"v"}}, reverse=reverse))


def test_support_uses_both_nonnested_memberships_and_matched_eligibility():
    pairs, audit = module().select_matched_pairs(*nonnested_designs(), min_unique=2)
    assert pairs[["island_a", "island_b"]].values.tolist() == [["A", "C"]]
    assert audit["n_eligible_islands_all"] == audit["n_eligible_islands_direct"] == 4
    assert audit["n_matched_eligible_islands"] == 3
    assert audit["n_candidate_pairs"] == 3
    assert audit["n_supported_pairs_all"] == audit["n_supported_pairs_direct"] == 2
    assert audit["n_joint_supported_pairs"] == 1
    assert audit["n_excluded_pairs"] == 2
    assert audit["n_islands_with_supported_pairs"] == audit["n_sampled_islands"] == 2
    assert audit["n_sampled_pairs"] == 1


def test_zero_supported_pairs_have_valid_empty_schema_and_audit():
    all_design, direct = nonnested_designs()
    pairs, audit = module().select_matched_pairs(all_design, direct, min_unique=4)
    assert pairs.columns.tolist() == ["island_a", "island_b"]
    assert pairs.empty and audit["n_candidate_pairs"] == 3
    assert audit["n_joint_supported_pairs"] == audit["n_sampled_pairs"] == 0
    assert audit["n_islands_with_supported_pairs"] == audit["n_sampled_islands"] == 0


def test_no_matched_islands_keeps_zero_candidate_counts():
    pairs, audit = module().select_matched_pairs(
        membership_design({"A": {"a"}}), membership_design({"B": {"b"}}), min_unique=1)
    assert pairs.empty
    assert audit["n_matched_eligible_islands"] == audit["n_candidate_pairs"] == 0


def test_selection_is_blind_to_profiles_and_input_row_order():
    first, audit1 = module().select_matched_pairs(*nonnested_designs(), min_unique=2, seed=17)
    second, audit2 = module().select_matched_pairs(*nonnested_designs(reverse=True),
                                                 min_unique=2, seed=17)
    # These designs deliberately expose no profile values: selection can only
    # use incidence and counts. Renumbered member rows cannot change the cohort.
    pd.testing.assert_frame_equal(first, second)
    assert audit1 == audit2


def test_seeded_reservoir_caps_pairs_without_duplicates_or_order_bias():
    members = {f"I{i:02}": {f"taxon{i}"} for i in range(12)}
    design = membership_design(members)
    first, audit = module().select_matched_pairs(design, design, min_unique=1, max_pairs=7, seed=123)
    second, _ = module().select_matched_pairs(membership_design(members, reverse=True), design,
                                            min_unique=1, max_pairs=7, seed=123)
    different, _ = module().select_matched_pairs(design, design, min_unique=1, max_pairs=7, seed=124)
    pd.testing.assert_frame_equal(first, second)
    assert not first.equals(different)
    assert len(first) == 7 and not first.duplicated().any()
    assert list(first.itertuples(index=False, name=None)) == sorted(first.itertuples(index=False, name=None))
    assert (first.island_a < first.island_b).all()
    assert audit["n_joint_supported_pairs"] == audit["n_candidate_pairs"] == 66
    assert audit["n_islands_with_supported_pairs"] == 12


def test_reservoir_gives_every_qualifying_pair_comparable_inclusion_frequency():
    design = membership_design({str(i): {str(i)} for i in range(4)})
    frequencies = Counter()
    for seed in range(600):
        pairs, _ = module().select_matched_pairs(design, design, min_unique=1, max_pairs=2, seed=seed)
        frequencies.update(pairs.itertuples(index=False, name=None))
    assert set(frequencies) == set(itertools.combinations(map(str, range(4)), 2))
    assert all(150 < count < 250 for count in frequencies.values())


@pytest.mark.parametrize("keyword,value", [("min_unique", 0), ("min_unique", True),
                                           ("max_pairs", 0), ("max_pairs", 1.5)])
def test_invalid_support_or_cap_fails_explicitly(keyword, value):
    with pytest.raises(ValueError, match="positive integer"):
        module().select_matched_pairs(*nonnested_designs(), **{keyword: value})


def test_different_regimes_fail_before_pair_selection():
    first, second = nonnested_designs()
    second.analysis_regime = REGIMES[1]
    with pytest.raises(ValueError, match="same regime"):
        module().select_matched_pairs(first, second)


def fixture_inputs(tmp_path):
    results = tmp_path / "source"
    results.mkdir()
    names = [f"Alpha {x}" for x in "abcdef"]
    memberships = [
        {"A": "abc", "B": "ade", "C": "bef", "pool": "f"},
        {"A": "abc", "B": "abc"},
        {"A": "ab"},
        {"A": "abc", "B": "ade", "pool": "f"},
    ]
    flora, geography = [], []
    for region, members in zip(REGIMES, memberships):
        for island, letters in members.items():
            key = f"{region}:{island}"
            geography.append((key, region))
            flora.extend((key, f"Alpha {letter}") for letter in letters)
    flora_path, geography_path = tmp_path / "flora.csv.gz", tmp_path / "geography.csv"
    pd.DataFrame(flora, columns=["island_id", "accepted_species"]).to_csv(flora_path, index=False)
    pd.DataFrame(geography, columns=["island_id", "analysis_regime"]).to_csv(geography_path, index=False)
    for scope in SCOPES:
        directory = results / scope
        directory.mkdir()
        profile_ids = [0, 1, 1, 2, 2, 3] if scope.startswith("direct") else [0, 0, 1, 1, 2, 2]
        profile_lookup = dict(zip("abcdef", profile_ids))
        pd.DataFrame({"accepted_species": names, "genus": ["Alpha"] * 6,
                      "profile_id": profile_ids}).to_csv(directory / "species_profiles.csv.gz", index=False)
        metrics = []
        for region, members in zip(REGIMES, memberships):
            pool_letters = set("".join(members.values()))
            pool_counts = Counter(profile_lookup[x] for x in pool_letters)
            for island, letters in members.items():
                k, n = len(letters), len(pool_letters)
                expected = sum(1 - math.comb(n - count, k) / math.comb(n, k)
                               if n - count >= k else 1 for count in pool_counts.values())
                metrics.append({"island_id": f"{region}:{island}", "analysis_regime": region,
                                "n_species": k, "n_profiles": len({profile_lookup[x] for x in letters}),
                                "expected_profiles_genus_pool": expected})
        pd.DataFrame(metrics).to_csv(directory / "island_metrics.csv", index=False)
        (directory / "input_manifest.json").write_text(json.dumps({"fixture": scope}))
    (results / "provenance_manifest.json").write_text(json.dumps({"status": "fixture", "run_id": 42}))
    return results, flora_path, geography_path


def run_cli(results, flora, geography, output, **overrides):
    options = dict(permutations=11, min_species=3, min_unique=2, max_pairs=2, seed=123)
    options.update(overrides)
    return subprocess.run([
        sys.executable, "-m", "island_v2.floral_space_pair_analysis",
        "--results", str(results), "--flora", str(flora), "--geography", str(geography),
        "--output", str(output),
        *[value for key, number in options.items() for value in ("--" + key.replace("_", "-"), str(number))],
    ], text=True, capture_output=True, env={**os.environ, "PYTHONPATH": "src"})


def test_cli_outputs_matched_pairs_fixed_family_and_reconstructible_draws(tmp_path):
    results, flora, geography = fixture_inputs(tmp_path)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output)
    assert process.returncode == 0, process.stderr
    assert sorted(p.name for p in output.iterdir()) == sorted([
        "manifest.json", "pair_cohort_audit.csv", "pairwise_observed_and_null.csv.gz",
        "permutation_summary.csv", "permutation_null_draws.csv.gz"])
    summary = pd.read_csv(output / "permutation_summary.csv")
    draws = pd.read_csv(output / "permutation_null_draws.csv.gz")
    pairs = pd.read_csv(output / "pairwise_observed_and_null.csv.gz")
    audit = pd.read_csv(output / "pair_cohort_audit.csv")
    assert list(summary[["evidence_scope", "analysis_regime"]].itertuples(index=False, name=None)) == list(itertools.product(SCOPES, REGIMES))
    assert len(audit) == 8
    assert len(pairs) == 4 * 3 and len(draws) == 8 * 11
    assert set(audit.encoding) == {"full_state_sets", "single_state"}
    for _, group in pairs.groupby(["encoding", "analysis_regime"]):
        scope_pairs = [list(g[["island_a", "island_b"]].itertuples(index=False, name=None))
                       for _, g in group.groupby("evidence_scope")]
        assert len(scope_pairs) == 2 and scope_pairs[0] == scope_pairs[1]
    from island_v2.floral_space_null import finite_two_sided_p
    for row in summary.itertuples():
        part = draws.loc[draws.evidence_scope.eq(row.evidence_scope) & draws.analysis_regime.eq(row.analysis_regime)]
        if row.status == "no_eligible_pairs":
            assert row.n_pairs == row.n_islands == row.n_permutations == 0
            assert pd.isna(row.p_two_sided) and part.empty
        else:
            assert len(part) == row.n_permutations == 11
            inference = finite_two_sided_p(row.observed_median_jaccard, part.median_profile_jaccard)
            for column, value in inference.items():
                assert getattr(row, column) == pytest.approx(value)
            assert row.n_taxa_pool == 6  # Includes the below-threshold pool island.
    safe_p = summary.p_two_sided.fillna(1)
    valid = summary.p_two_sided.notna()
    np.testing.assert_allclose(summary.loc[valid, "q_bh_16"], false_discovery_control(safe_p, method="bh")[valid])
    np.testing.assert_allclose(summary.loc[valid, "q_by_16"], false_discovery_control(safe_p, method="by")[valid])
    assert summary.loc[~valid, ["q_bh_16", "q_by_16"]].isna().all().all()
    assert "accepted_species" not in pairs.columns
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["fixed_comparison_family_size"] == 16
    assert manifest["parameters"] == {"permutations": 11, "min_species": 3, "min_unique": 2, "max_pairs": 2, "seed": 123}
    assert len(manifest["seed_streams"]["pair_selection"]) == 8
    assert len(manifest["seed_streams"]["permutation"]) == 16
    assert manifest["seed_streams"]["pair_selection"][0]["entropy"] == [123, 17, 0, 0]
    assert manifest["seed_streams"]["permutation"][-1]["entropy"] == [123, 29, 3, 3]
    assert len(manifest["parent_commit"]) == 40
    assert manifest["software"]["python"]
    assert {"numpy", "pandas", "scipy"} <= set(manifest["software"]["libraries"])
    assert set(manifest["code"]) == {"floral_space_pairs", "floral_space_pair_analysis"}
    assert manifest["source_result_provenance"]
    for record in [*manifest["code"].values(), *manifest["outputs"].values()]:
        assert hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest() == record["sha256"]
    assert manifest["status"] == "complete"


def test_shuffled_sources_produce_identical_tables_and_gzip_headers(tmp_path):
    results, flora, geography = fixture_inputs(tmp_path)
    first, second = tmp_path / "first", tmp_path / "second"
    assert run_cli(results, flora, geography, first).returncode == 0
    for path in [flora, geography, *results.glob("*/species_profiles.csv.gz"),
                 *results.glob("*/island_metrics.csv")]:
        frame = pd.read_csv(path)
        frame.sample(frac=1, random_state=42).to_csv(path, index=False)
    process = run_cli(results, flora, geography, second)
    assert process.returncode == 0, process.stderr
    for name in ["permutation_summary.csv", "pair_cohort_audit.csv",
                 "pairwise_observed_and_null.csv.gz", "permutation_null_draws.csv.gz"]:
        assert (first / name).read_bytes() == (second / name).read_bytes()
        if name.endswith(".gz"):
            assert (first / name).read_bytes()[4:8] == b"\0\0\0\0"


@pytest.mark.parametrize("option,value", [("permutations", 0), ("min_species", 0),
                                          ("min_unique", 0), ("max_pairs", 0), ("seed", -1)])
def test_cli_rejects_invalid_options_before_writing_output(tmp_path, option, value):
    results, flora, geography = fixture_inputs(tmp_path)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output, **{option: value})
    assert process.returncode != 0
    assert "positive integer" in process.stderr or "nonnegative integer" in process.stderr
    assert not output.exists()


def test_cli_retains_saved_metric_consistency_gate(tmp_path):
    results, flora, geography = fixture_inputs(tmp_path)
    path = results / "direct" / "island_metrics.csv"
    metrics = pd.read_csv(path)
    metrics.loc[0, "expected_profiles_genus_pool"] = 100
    metrics.to_csv(path, index=False)
    process = run_cli(results, flora, geography, tmp_path / "output")
    assert process.returncode != 0
    assert "expectation" in process.stderr


def test_suite_matches_scope_eligibility_within_each_encoding(tmp_path):
    results, flora, geography = fixture_inputs(tmp_path)
    path = results / "direct" / "island_metrics.csv"
    metrics = pd.read_csv(path)
    metrics = metrics.loc[~metrics.island_id.eq(f"{REGIMES[0]}:C")]
    metrics.to_csv(path, index=False)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output)
    assert process.returncode == 0, process.stderr
    summary = pd.read_csv(output / "permutation_summary.csv")
    first = summary.loc[summary.analysis_regime.eq(REGIMES[0])].set_index("evidence_scope")
    assert first.n_pairs.to_dict() == {"all": 1, "direct": 1, "all_single_state": 2,
                                     "direct_single_state": 2}
    audit = pd.read_csv(output / "pair_cohort_audit.csv")
    cohort = audit.loc[audit.analysis_regime.eq(REGIMES[0])].set_index("encoding")
    assert cohort.loc["full_state_sets", "n_matched_eligible_islands"] == 2
    assert cohort.loc["single_state", "n_matched_eligible_islands"] == 3


def test_all_empty_suite_keeps_complete_schema_and_undefined_inference(tmp_path):
    results, flora, geography = fixture_inputs(tmp_path)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output, min_unique=4)
    assert process.returncode == 0, process.stderr
    summary = pd.read_csv(output / "permutation_summary.csv")
    assert len(summary) == 16 and set(summary.status) == {"no_eligible_pairs"}
    assert (summary[["n_pairs", "n_islands", "n_permutations"]] == 0).all().all()
    assert summary[["p_lower", "p_upper", "p_two_sided", "q_bh_16", "q_by_16",
                    "null_median", "observed_median_jaccard", "observed_to_null_median"]].isna().all().all()
    pairs = pd.read_csv(output / "pairwise_observed_and_null.csv.gz")
    draws = pd.read_csv(output / "permutation_null_draws.csv.gz")
    assert pairs.empty and {"island_a", "island_b", "profile_jaccard", "null_sd_jaccard"} <= set(pairs)
    assert draws.empty and {"permutation", "median_profile_jaccard"} <= set(draws)
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["actual_counts"]["nonempty_comparisons"] == 0
    assert manifest["actual_counts"]["summary_slots"] == 16


def test_suite_uses_declared_permutation_stream_and_hashes_each_input(tmp_path):
    from island_v2.floral_space_null import build_membership, prepare_region
    from island_v2.floral_space_pairs import prepare_pair_design, run_pair_comparison
    results, flora, geography = fixture_inputs(tmp_path)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output)
    assert process.returncode == 0, process.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    records = [manifest["sources"]["flora"], manifest["sources"]["geography"],
               *manifest["source_result_provenance"].values()]
    records.extend(record for scope in SCOPES for record in manifest["sources"][scope].values())
    assert len(manifest["source_result_provenance"]) == 5
    for record in records:
        assert hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest() == record["sha256"]
    profiles = pd.read_csv(results / "all" / "species_profiles.csv.gz")
    metrics = pd.read_csv(results / "all" / "island_metrics.csv")
    design = prepare_region(build_membership(profiles, pd.read_csv(flora), pd.read_csv(geography)),
                            metrics, REGIMES[0], min_species=3)
    pairs = pd.read_csv(output / "pairwise_observed_and_null.csv.gz")
    chosen = pairs.loc[pairs.evidence_scope.eq("all") & pairs.analysis_regime.eq(REGIMES[0]),
                       ["island_a", "island_b"]]
    _, expected, _ = run_pair_comparison(prepare_pair_design(design, chosen, min_unique=2),
                                       permutations=11, seed=np.random.SeedSequence([123, 29, 0, 0]))
    draws = pd.read_csv(output / "permutation_null_draws.csv.gz")
    actual = draws.loc[draws.evidence_scope.eq("all") & draws.analysis_regime.eq(REGIMES[0]),
                       expected.columns].reset_index(drop=True)
    pd.testing.assert_frame_equal(actual, expected)


def test_cli_bh_and_by_match_scipy_with_nontrivial_adjustments(tmp_path):
    """Matching rare concentrated sets produce small p, exposing BH/BY errors."""
    results = tmp_path / "source"
    results.mkdir()
    names = [f"Alpha taxon{i:02}" for i in range(24)]
    profiles = pd.DataFrame({"accepted_species": names, "genus": ["Alpha"] * 24,
                             "profile_id": [0] * 6 + list(range(1, 19))})
    flora_rows = [("A", name) for name in names[:3]] + [("B", name) for name in names[3:6]]
    flora_rows += [(f"pool{i:02}", name) for i, name in enumerate(names[6:])]
    flora, geography = tmp_path / "flora.csv", tmp_path / "geography.csv"
    pd.DataFrame(flora_rows, columns=["island_id", "accepted_species"]).to_csv(flora, index=False)
    pd.DataFrame({"island_id": sorted({island for island, _ in flora_rows}),
                  "analysis_regime": REGIMES[0]}).to_csv(geography, index=False)
    expectation = 1 - math.comb(18, 3) / math.comb(24, 3) + 18 * 3 / 24
    metrics = pd.DataFrame({"island_id": ["A", "B"], "analysis_regime": REGIMES[0],
                            "n_species": 3, "n_profiles": 1,
                            "expected_profiles_genus_pool": expectation})
    for scope in SCOPES:
        directory = results / scope
        directory.mkdir()
        profiles.to_csv(directory / "species_profiles.csv.gz", index=False)
        metrics.to_csv(directory / "island_metrics.csv", index=False)
    output = tmp_path / "output"
    process = run_cli(results, flora, geography, output, permutations=199)
    assert process.returncode == 0, process.stderr
    summary = pd.read_csv(output / "permutation_summary.csv")
    valid = summary.p_two_sided.notna()
    safe = summary.p_two_sided.fillna(1)
    assert valid.sum() == 4
    assert (summary.loc[valid, "q_bh_16"] < summary.loc[valid, "q_by_16"]).all()
    assert (summary.loc[valid, "q_by_16"] < 1).all()
    np.testing.assert_allclose(summary.loc[valid, "q_bh_16"], false_discovery_control(safe, method="bh")[valid])
    np.testing.assert_allclose(summary.loc[valid, "q_by_16"], false_discovery_control(safe, method="by")[valid])
