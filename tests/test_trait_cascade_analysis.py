from pathlib import Path

import pandas as pd
import pytest

from island_v2 import trait_cascade_analysis as module


def test_load_shards_requires_complete_128(tmp_path: Path) -> None:
    root = tmp_path / "fill_cascade"
    shard = root / "shards" / "shard_00000"
    shard.mkdir(parents=True)
    pd.DataFrame(
        {
            "accepted_species": ["Species alpha"],
            "trait_name": ["flower_primary_color"],
            "filled_value": ["white"],
            "fill_tier": ["species_direct"],
            "evidence_scope": ["species_direct"],
            "confidence": ["direct_reported"],
        }
    ).to_csv(shard / "fills.csv.gz", index=False, compression="gzip")

    with pytest.raises(Exception, match="expected 128 fill shards"):
        module._load_shards(root)


def test_analysis_tier_contract_is_ordered() -> None:
    config = module._load_tiers(Path("config/trait_cascade_analysis_tiers.yml"))
    tiers = config["analysis_tiers"]
    assert tiers["primary"]["allowed_fill_tiers"] == ["species_direct"]
    assert "global_fallback" not in tiers["broad"]["allowed_fill_tiers"]
    assert "global_fallback" not in tiers["sensitivity_all"]["allowed_fill_tiers"]
    assert "unresolved_no_evidence" not in tiers["sensitivity_all"]["allowed_fill_tiers"]


def test_analysis_tier_contract_rejects_unresolved_or_global(tmp_path: Path) -> None:
    config = tmp_path / "tiers.yml"
    config.write_text(
        "analysis_tiers:\n  bad:\n    allowed_fill_tiers: [species_direct, global_fallback]\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="non-analysable fills"):
        module._load_tiers(config)


def test_historical_artifact_without_analysis_flag_is_fail_closed() -> None:
    frame = pd.DataFrame(
        {
            "fill_tier": [
                "species_direct",
                "genus_inference",
                "family_inference",
                "global_fallback",
                "unresolved_no_evidence",
                "llm_inference",
                "future_unknown_tier",
            ]
        }
    )

    assert module._analysis_eligibility_mask(frame).tolist() == [
        True,
        True,
        True,
        False,
        False,
        False,
        False,
    ]


def test_analysis_tier_contract_rejects_unknown_or_llm_tiers(tmp_path: Path) -> None:
    config = tmp_path / "tiers.yml"
    config.write_text(
        "analysis_tiers:\n  bad:\n    allowed_fill_tiers: [species_direct, llm_inference]\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="non-analysable fills"):
        module._load_tiers(config)
