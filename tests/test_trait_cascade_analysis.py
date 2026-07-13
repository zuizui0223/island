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
    assert "global_fallback" in tiers["sensitivity_all"]["allowed_fill_tiers"]
