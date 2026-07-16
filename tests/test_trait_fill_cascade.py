import hashlib
import json

import pandas as pd
import pytest

from island_v2 import trait_fill_cascade as cascade


def _config(master_path) -> dict:
    return {
        "master_taxa_csv": str(master_path),
        "species_column": "accepted_species",
        "genus_column": "genus",
        "family_column": "family",
        "target_traits": ["flower_primary_color"],
        "min_genus_support": 1,
        "min_family_support": 3,
        "evidence_sources": {},
        "tier_labels": {
            "species_direct": {"evidence_scope": "species_direct", "confidence": "direct_reported"},
            "synonym_direct": {"evidence_scope": "synonym_direct", "confidence": "synonym_reported"},
            "genus_inference": {"evidence_scope": "genus_inference", "confidence": "genus_inference_low"},
            "family_inference": {"evidence_scope": "family_inference", "confidence": "family_inference_very_low"},
            "global_fallback": {"evidence_scope": "global_fallback", "confidence": "global_fallback_negligible"},
        },
        "benchmark_sample_size": 10,
    }


def _master() -> pd.DataFrame:
    # Genus Aaa: one direct species + one gap (genus inference).
    # Genus Bbb (family Fff): no direct, but family Fff has >=3 direct -> family inference.
    # Genus Zzz (family Ggg): nothing in family -> global fallback.
    rows = [
        ("Aaa one", "Aaa", "Fam1"),
        ("Aaa two", "Aaa", "Fam1"),   # gap -> genus inference from Aaa one
        ("Bbb gap", "Bbb", "Fff"),    # gap, genus Bbb has no direct -> family Fff inference
        ("Fff a", "Ccc", "Fff"),
        ("Fff b", "Ddd", "Fff"),
        ("Fff c", "Eee", "Fff"),
        ("Zzz gap", "Zzz", "Ggg"),    # nothing in genus/family -> global fallback
    ]
    return pd.DataFrame(rows, columns=["accepted_species", "genus", "family"])


def _evidence() -> pd.DataFrame:
    rows = [
        ("Aaa one", "flower_primary_color", "red", 5.0),
        ("Fff a", "flower_primary_color", "white", 1.0),
        ("Fff b", "flower_primary_color", "white", 1.0),
        ("Fff c", "flower_primary_color", "blue", 1.0),  # family modal = white (2 vs 1)
    ]
    return pd.DataFrame(rows, columns=["accepted_species", "trait_name", "value", "weight"])


def test_cascade_tiers_and_modal(tmp_path):
    master = _master()
    config = _config(tmp_path / "m.csv")
    fills = cascade.build_fills(master, _evidence(), config).set_index("accepted_species")

    assert fills.loc["Aaa one", "fill_tier"] == "species_direct"
    assert fills.loc["Aaa one", "filled_value"] == "red"
    assert fills.loc["Aaa one", "reported_value"] == "red"
    assert fills.loc["Aaa one", "inferred_value"] == ""
    # genus Aaa has one direct (red) -> genus inference for the gap
    assert fills.loc["Aaa two", "fill_tier"] == "genus_inference"
    assert fills.loc["Aaa two", "filled_value"] == "red"
    assert fills.loc["Aaa two", "reported_value"] == ""
    assert fills.loc["Aaa two", "inferred_value"] == "red"
    # genus Bbb empty, family Fff has 3 direct -> family inference, modal white
    assert fills.loc["Bbb gap", "fill_tier"] == "family_inference"
    assert fills.loc["Bbb gap", "filled_value"] == "white"
    # nothing in genus/family -> global fallback; global prior is one species,
    # one vote, so two white species beat one red and one blue species.
    assert fills.loc["Zzz gap", "fill_tier"] == "global_fallback"
    assert fills.loc["Zzz gap", "filled_value"] == "white"

    # The selected value is never duplicated across reported and inferred channels.
    assert (fills["reported_value"].ne("") ^ fills["inferred_value"].ne("")).all()

    # every master species got a fill: unknown driven to zero
    assert len(fills) == len(master)

    # inference tiers retain a distribution; direct does not
    assert fills.loc["Aaa two", "value_distribution"] != ""
    assert fills.loc["Aaa one", "value_distribution"] == ""
    dist = json.loads(fills.loc["Bbb gap", "value_distribution"])
    assert dist["white"] == 2.0 and dist["blue"] == 1.0


def test_family_threshold_blocks_thin_support(tmp_path):
    master = _master()
    config = _config(tmp_path / "m.csv")
    config["min_family_support"] = 5  # Fff has only 3 direct -> below threshold
    fills = cascade.build_fills(master, _evidence(), config).set_index("accepted_species")
    # Bbb gap can no longer use family Fff (3 < 5) -> falls through to global fallback
    assert fills.loc["Bbb gap", "fill_tier"] == "global_fallback"


def test_equal_weight_species_conflict_is_not_promoted_as_reported(tmp_path):
    master = pd.DataFrame(
        [("Conflict one", "Conflict", "Fam"), ("Anchor one", "Anchor", "Other")],
        columns=["accepted_species", "genus", "family"],
    )
    evidence = pd.DataFrame(
        [
            ("Conflict one", "self_incompatibility", "SI", 1.0),
            ("Conflict one", "self_incompatibility", "SC", 1.0),
            ("Anchor one", "self_incompatibility", "SI", 1.0),
        ],
        columns=["accepted_species", "trait_name", "value", "weight"],
    )
    config = _config(tmp_path / "m.csv")
    config["target_traits"] = ["self_incompatibility"]
    fills = cascade.build_fills(master, evidence, config).set_index("accepted_species")

    conflict = fills.loc["Conflict one"]
    assert conflict["reported_value"] == ""
    assert conflict["inferred_value"] == "SI"
    assert conflict["fill_tier"] == "global_fallback"
    assert conflict["direct_conflict"] == "true"
    assert json.loads(conflict["direct_conflict_distribution"]) == {"SC": 1.0, "SI": 1.0}


def test_tied_global_prior_uses_reproducible_distribution_draw(tmp_path):
    master = pd.DataFrame(
        [
            ("Open one", "Open", "Openaceae"),
            ("Deep one", "Deep", "Deepaceae"),
            ("Gap one", "Gap", "Gapaceae"),
        ],
        columns=["accepted_species", "genus", "family"],
    )
    evidence = pd.DataFrame(
        [
            ("Open one", "tube_depth_class", "absent_or_open", 1.0),
            ("Deep one", "tube_depth_class", "deep", 1.0),
        ],
        columns=["accepted_species", "trait_name", "value", "weight"],
    )
    config = _config(tmp_path / "m.csv")
    config["target_traits"] = ["tube_depth_class"]

    first = cascade.build_fills(master, evidence, config).set_index("accepted_species")
    second = cascade.build_fills(master, evidence, config).set_index("accepted_species")
    gap = first.loc["Gap one"]

    assert gap["fill_tier"] == "global_fallback"
    assert gap["reported_value"] == ""
    assert gap["inferred_value"] in {"absent_or_open", "deep"}
    assert gap["filled_value"] == second.loc["Gap one", "filled_value"]
    assert json.loads(gap["value_distribution"]) == {
        "absent_or_open": 1.0,
        "deep": 1.0,
    }


def test_tied_genus_prior_stays_at_genus_tier(tmp_path):
    master = pd.DataFrame(
        [
            ("Aaa one", "Aaa", "Fff"),
            ("Aaa two", "Aaa", "Fff"),
            ("Aaa gap", "Aaa", "Fff"),
            ("Other one", "Other", "Otheraceae"),
        ],
        columns=["accepted_species", "genus", "family"],
    )
    evidence = pd.DataFrame(
        [
            ("Aaa one", "floral_symmetry", "actinomorphic", 1.0),
            ("Aaa two", "floral_symmetry", "zygomorphic", 1.0),
            ("Other one", "floral_symmetry", "actinomorphic", 1.0),
        ],
        columns=["accepted_species", "trait_name", "value", "weight"],
    )
    config = _config(tmp_path / "m.csv")
    config["target_traits"] = ["floral_symmetry"]
    fills = cascade.build_fills(master, evidence, config).set_index("accepted_species")

    gap = fills.loc["Aaa gap"]
    assert gap["fill_tier"] == "genus_inference"
    assert gap["filled_value"] in {"actinomorphic", "zygomorphic"}
    assert json.loads(gap["value_distribution"]) == {
        "actinomorphic": 1.0,
        "zygomorphic": 1.0,
    }


def test_candidate_long_reads_standardized_values_and_excludes_inference(tmp_path):
    candidates = tmp_path / "trait_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "self_incompatibility",
                "standardized_value": "SC",
                "candidate_kind": "source_backed",
                "evidence_scope": "species_direct",
                "source_url": "https://example.org/alpha-one",
                "evidence_excerpt": "Alpha one is self-compatible.",
            },
            {
                "accepted_species": "Alpha two",
                "trait_name": "self_incompatibility",
                "standardized_value": "SC",
                "candidate_kind": "hierarchical_inference",
                "evidence_scope": "genus_inference",
                "source_url": "https://example.org/alpha-two",
                "evidence_excerpt": "Genus-level inference.",
            },
            {
                "accepted_species": "Alpha three",
                "trait_name": "self_incompatibility",
                "standardized_value": "SI",
                "candidate_kind": "source_backed",
                "evidence_scope": "species_indirect",
                "source_url": "https://example.org/alpha-three",
                "evidence_excerpt": "Indirect statement.",
            },
        ]
    ).to_csv(candidates, index=False)
    loaded = cascade._evidence_candidate_long({"glob": str(candidates)})
    assert loaded.to_dict("records") == [
        {
            "accepted_species": "Alpha one",
            "trait_name": "self_incompatibility",
            "value": "SC",
            "weight": 1.0,
        }
    ]


def test_candidate_long_reads_provisional_values_only_after_bundle_validation(tmp_path):
    candidates = tmp_path / "v2_llm_reported_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "sex_system",
                "provisional_candidate_value": "dioecious",
                "evidence_scope": "species_direct",
                "evidence_quote": "Plants dioecious.",
                "extraction_method": "provider_neutral_llm_from_frozen_source_text",
                "source_url": "https://example.org/alpha-one",
            }
        ]
    ).to_csv(candidates, index=False)

    loaded = cascade._evidence_candidate_long(
        {"glob": str(candidates), "require_candidate_kind": False}
    )

    assert loaded.to_dict("records") == [
        {
            "accepted_species": "Alpha one",
            "trait_name": "sex_system",
            "value": "dioecious",
            "weight": 1.0,
        }
    ]


def test_validated_llm_adapter_rejects_self_consistent_but_source_free_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    candidate = bundle / "v2_llm_reported_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "sex_system",
                "provisional_candidate_value": "dioecious",
                "evidence_scope": "species_direct",
            }
        ]
    ).to_csv(candidate, index=False)
    files = {
        "packet_input_sha256": bundle / "packet_input.json",
        "packet_manifest_sha256": bundle / "packet_manifest.json",
        "prompt_sha256": bundle / "prompt.txt",
        "result_jsonl_sha256": bundle / "result.jsonl",
    }
    for path in files.values():
        path.write_text("locked\n", encoding="utf-8")
    pd.DataFrame(columns=["reason"]).to_csv(bundle / "llm_rejected_claims.csv", index=False)

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest = {
        "validation_status": "success",
        "accepted_csv_filename": candidate.name,
        "accepted_csv_sha256": sha(candidate),
        "accepted_csv_row_count": 1,
        "n_rejected_claims": 0,
        **{key: sha(path) for key, path in files.items()},
    }
    (bundle / "llm_ingest_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    # Matching file hashes alone are insufficient: the canonical validator must
    # also be able to bind every row and exact quote to a frozen source chunk.
    with pytest.raises(Exception):
        cascade._evidence_validated_llm_bundle({"glob": str(candidate)})


def test_allmaster_adapter_keeps_only_source_backed_and_maps_legacy_fields(tmp_path):
    evidence = tmp_path / "all_species_trait_evidence.csv.gz"
    pd.DataFrame(
        [
            {
                "species": "Alpha one",
                "field": "flower_color",
                "value": "red/pink",
                "source_backed": "True",
                "evidence_type": "flora",
            },
            {
                "species": "Alpha two",
                "field": "self_incompatibility",
                "value": "likely_SC",
                "source_backed": "False",
                "evidence_type": "inference",
            },
            {
                "species": "Alpha three",
                "field": "flower_shape",
                "value": "actinomorphic / radially symmetric",
                "source_backed": "yes",
                "evidence_type": "flora",
            },
            {
                "species": "Alpha four",
                "field": "pollination_guild",
                "value": "bees",
                "source_backed": "True",
                "evidence_type": "inference",
            },
        ]
    ).to_csv(evidence, index=False, compression="gzip")

    loaded = cascade._evidence_allmaster_long({"glob": str(evidence)})

    assert loaded.to_dict("records") == [
        {
            "accepted_species": "Alpha one",
            "trait_name": "flower_primary_color",
            "value": "red_pink",
            "weight": 1.0,
        },
        {
            "accepted_species": "Alpha three",
            "trait_name": "floral_symmetry",
            "value": "actinomorphic",
            "weight": 1.0,
        },
    ]


def test_direct_value_normalization_is_ontology_safe_and_conservative():
    assert cascade._normalize_direct_value("flower_primary_color", "pink|white") == (
        "multicolored_variable"
    )
    assert cascade._normalize_direct_value("floral_form", "campanulate") == (
        "bell_campanulate"
    )
    assert cascade._normalize_direct_value("floral_form", "campanulate|tubular") == (
        "other_described"
    )
    assert cascade._normalize_direct_value("floral_symmetry", "radial") == "actinomorphic"
    assert cascade._normalize_direct_value("floral_symmetry", "bilateral|radial") == ""
    assert cascade._normalize_direct_value("self_incompatibility", "self_compatible") == "SC"
    assert cascade._normalize_direct_value("sex_system", "dioecious|unisexual") == "dioecious"
    assert cascade._normalize_direct_value("sex_system", "unisexual") == ""


def test_coverage_summary_zero_unknown(tmp_path):
    master = _master()
    config = _config(tmp_path / "m.csv")
    fills = cascade.build_fills(master, _evidence(), config)
    summary = cascade.build_coverage_summary(fills, config, len(master))
    trait = summary["by_trait"]["flower_primary_color"]
    assert trait["n_filled"] == len(master)
    assert trait["n_unknown_remaining"] == 0
    assert trait["n_species_direct"] == 4  # Aaa one + Fff a/b/c have direct evidence


def test_stable_shard_is_deterministic_and_bounded():
    for name in ("Aaa one", "Zzz gap", "Rosa canina"):
        assert cascade.stable_shard(name, 16) == cascade.stable_shard(name, 16)
        assert 0 <= cascade.stable_shard(name, 16) < 16


def test_model_fingerprint_changes_with_evidence(tmp_path):
    config = _config(tmp_path / "m.csv")
    base = cascade.model_fingerprint(_evidence(), config)
    assert base == cascade.model_fingerprint(_evidence(), config)  # stable
    more = pd.concat(
        [_evidence(), pd.DataFrame([("New sp", "flower_primary_color", "green", 1.0)],
                                   columns=["accepted_species", "trait_name", "value", "weight"])],
        ignore_index=True,
    )
    assert cascade.model_fingerprint(more, config) != base


def test_model_and_shard_fingerprints_change_with_taxonomy(tmp_path):
    config = _config(tmp_path / "m.csv")
    master = _master()
    changed = master.copy()
    changed.loc[changed["accepted_species"].eq("Aaa two"), "family"] = "ChangedFam"

    assert cascade.model_fingerprint(_evidence(), config, master) != (
        cascade.model_fingerprint(_evidence(), config, changed)
    )
    assert cascade.shard_species_fingerprint(master) != cascade.shard_species_fingerprint(
        changed
    )


def test_sharded_fills_partition_matches_whole_universe(tmp_path):
    master = _master()
    config = _config(tmp_path / "m.csv")
    evidence = _evidence()

    whole = cascade.build_fills(master, evidence, config)
    whole_summary = cascade.build_coverage_summary(whole, config, len(master))

    # Apply the same model per shard and aggregate the compact summaries.
    model = cascade.build_model(master, evidence, config)
    shard_count = 8
    assignments = master["accepted_species"].map(lambda s: cascade.stable_shard(s, shard_count))
    shard_summaries = []
    all_species = set()
    for shard_index in range(shard_count):
        shard_master = master.loc[assignments.eq(shard_index)].reset_index(drop=True)
        shard_fills = cascade.fill_species_frame(shard_master, model, config)
        all_species.update(shard_fills["accepted_species"].tolist())
        shard_summaries.append(cascade.shard_summary_from_fills(shard_fills, config))

    agg = cascade.aggregate_shard_summaries(shard_summaries, config, len(master))
    # Every eligible species landed in exactly one shard, and the aggregate equals
    # the single-pass coverage.
    assert all_species == set(master["accepted_species"])
    assert agg["fills_by_tier"] == whole_summary["fills_by_tier"]
    assert agg["by_trait"]["flower_primary_color"] == whole_summary["by_trait"]["flower_primary_color"]
