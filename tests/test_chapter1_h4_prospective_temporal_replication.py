from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h4_prospective_temporal_replication import (
    PLACEHOLDER_COMMIT,
    _verify_unblinding_lock,
    build_support_lock,
    load_config,
    match_metadata_to_traits,
    prepare_frozen_trait_states,
    run_replication,
    validate_metadata_outcome_blind,
)


CONFIG_PATH = Path("config/chapter1_h4_prospective_temporal_replication_v1.yml")


def _config() -> dict:
    with CONFIG_PATH.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _metadata(n: int = 40) -> pd.DataFrame:
    rows = []
    for i in range(n):
        rows.append(
            {
                "experiment_key": f"e{i:03d}",
                "publication_id": f"p{i // 2:03d}",
                "doi": f"10.0000/{i // 2:03d}",
                "publication_date": f"{2016 + (i % 10)}-06-01",
                "study_key": f"p{i // 2:03d}",
                "site_key": f"site{i:03d}",
                "accepted_species": f"Species {i:03d}",
                "analysis_regime": "tropical" if (i // 2) % 2 else "northern_midlatitude",
                "z_distance": (i - n / 2) / 10.0,
            }
        )
    return pd.DataFrame(rows)


def _traits(n: int = 40) -> pd.DataFrame:
    rows = []
    for i in range(n):
        state = i % 2
        rows.append(
            {
                "accepted_species": f"Species {i:03d}",
                "autonomous_selfing": state,
                "generalized_form": state,
                "actinomorphic_symmetry": state,
                "shallow_open_tube": state,
            }
        )
    return pd.DataFrame(rows)


def test_preflight_rejects_outcome_columns() -> None:
    metadata = _metadata()
    metadata["PL_Effect_Size"] = 0.2
    with pytest.raises(
        Exception,
        match="outcome columns are forbidden before support lock",
    ):
        validate_metadata_outcome_blind(metadata, _config())


def test_temporal_holdout_rejects_legacy_glopl_year() -> None:
    metadata = _metadata()
    metadata.loc[0, "publication_date"] = "2015-12-31"
    with pytest.raises(Exception, match="outside frozen temporal holdout"):
        validate_metadata_outcome_blind(metadata, _config())


def test_accessibility_score_keeps_all_three_frozen_components() -> None:
    traits = _traits(3)
    traits.loc[0, "generalized_form"] = 1
    traits.loc[0, "actinomorphic_symmetry"] = 0
    traits.loc[0, "shallow_open_tube"] = 1
    prepared = prepare_frozen_trait_states(traits, _config())
    row = prepared.loc[prepared["accepted_species"].eq("Species 000")].iloc[0]
    assert row["accessibility_component_count"] == 3
    assert row["accessibility_generalization_score"] == pytest.approx(2 / 3)


def test_support_lock_is_outcome_blind_and_requires_commit(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.csv"
    traits_path = tmp_path / "traits.csv"
    config_path = tmp_path / "config.yml"
    metadata = validate_metadata_outcome_blind(_metadata(), _config())
    traits = prepare_frozen_trait_states(_traits(), _config())
    matched = match_metadata_to_traits(metadata, traits)

    metadata.to_csv(metadata_path, index=False)
    _traits().to_csv(traits_path, index=False)
    config_path.write_text(CONFIG_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    lock = build_support_lock(
        metadata_path,
        traits_path,
        config_path,
        metadata,
        traits,
        matched,
        _config(),
    )
    assert lock["outcome_status"] == "unopened"
    assert lock["frozen_commit"] == PLACEHOLDER_COMMIT
    assert set(lock["evaluable_primary_hypotheses"]) == {
        "H4a_reproductive_assurance",
        "H4b_accessibility_generalization",
    }


def test_analysis_refuses_uncommitted_support_lock(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.csv"
    traits_path = tmp_path / "traits.csv"
    config_path = tmp_path / "config.yml"
    metadata = validate_metadata_outcome_blind(_metadata(), _config())
    traits = prepare_frozen_trait_states(_traits(), _config())
    matched = match_metadata_to_traits(metadata, traits)

    metadata.to_csv(metadata_path, index=False)
    _traits().to_csv(traits_path, index=False)
    config_path.write_text(CONFIG_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    lock = build_support_lock(
        metadata_path,
        traits_path,
        config_path,
        metadata,
        traits,
        matched,
        _config(),
    )
    with pytest.raises(Exception, match="must be committed"):
        _verify_unblinding_lock(
            lock,
            metadata_path=metadata_path,
            trait_states_path=traits_path,
            config_path=config_path,
        )


def test_primary_replication_recovers_negative_effects() -> None:
    metadata = validate_metadata_outcome_blind(_metadata(80), _config())
    traits = prepare_frozen_trait_states(_traits(80), _config())
    matched = match_metadata_to_traits(metadata, traits)

    joined = matched.copy()
    joined["PL_Effect_Size_Type1"] = "FS"
    joined["PL_Effect_Size_Type2"] = "Sup"
    joined["Constant_added"] = "False"
    joined["Level_of_Supplementation"] = "whole_plant"
    joined["PL_Effect_Size"] = (
        0.25 * joined["z_distance"]
        - 0.55 * joined["autonomous_selfing"]
        - 0.30 * joined["accessibility_generalization_score"]
        + joined.index.to_series().mod(5) * 0.002
    )

    lock = {
        "support_counts": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": True},
        }
    }
    result = run_replication(joined, lock, _config())
    assert result["primary_results"]["H4a_reproductive_assurance"]["supported"]
    assert result["primary_results"]["H4b_accessibility_generalization"]["supported"]
    assert result["n_supported_primary_hypotheses"] == 2


def test_contract_declares_two_coprimary_alpha_0025() -> None:
    config = load_config(CONFIG_PATH)
    assert (
        config["primary_hypotheses"]["H4a_reproductive_assurance"]["one_sided_alpha"]
        == 0.025
    )
    assert (
        config["primary_hypotheses"]["H4b_accessibility_generalization"][
            "one_sided_alpha"
        ]
        == 0.025
    )
    assert (
        config["primary_hypotheses"]["H4b_accessibility_generalization"]["components"][
            "shallow_open_tube"
        ]
        == 1.0
    )
