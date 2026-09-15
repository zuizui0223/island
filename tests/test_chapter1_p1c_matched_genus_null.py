import json
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_p1c_matched_genus_null import aggregate_shards


def test_aggregate_uses_frozen_schedule(tmp_path: Path) -> None:
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    config = {
        "contract": "chapter1_p1_assembly_depth_defense_v1",
        "p1c_matched_complexity_genus_null": {
            "permutations": 4,
            "execution_shards": 2,
            "permutations_per_shard": 2,
            "minimum_valid_permutations": 3,
            "permutation_seed_rule": "master_seed_plus_zero_based_permutation_id",
            "random_seed": 10,
            "primary_pass_rule": {"max_randomization_p_value": 0.5},
        },
    }
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    for shard, ids in ((0, [0, 1]), (1, [2, 3])):
        frame = pd.DataFrame(
            {
                "permutation_id": ids,
                "seed": [10 + x for x in ids],
                "valid": [True, True],
                "primary_median_conditional_attenuation": [0.1, 0.2],
            }
        )
        frame.to_csv(shard_dir / f"p1c_null_shard_{shard:02d}.csv.gz", index=False)
        (shard_dir / f"p1c_null_shard_{shard:02d}_manifest.json").write_text(
            json.dumps({"observed_primary_statistic": 0.3}), encoding="utf-8"
        )
    result = aggregate_shards(
        shard_dir=shard_dir,
        p1_config_path=config_path,
        output_dir=tmp_path / "out",
    )
    assert result["valid_permutations"] == 4
    assert result["null_exceedances_greater_equal_observed"] == 0
    assert result["one_sided_randomization_p_value"] == 0.2
    assert result["strong_assembly_depth_defense_pass"] is True
