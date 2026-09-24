import os
from pathlib import Path

src = Path(os.environ["ISLAND_CORRECTION_SOURCES"])
out = Path(os.environ["ISLAND_CORRECTION_OUTPUT"])
data = Path(os.environ["ISLAND_CORRECTION_GEOMETRY"])
out.mkdir(parents=True, exist_ok=True)
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2 import chapter1_v14_h2_decomposition as h2

cfg = yaml.safe_load((src / "island--chapter1_v14_h2_decomposition.yml").read_text(encoding="utf8"))
cov = pd.read_csv(
    src / "progressive-input/fixed/isolation/results/purpose_shortest_island_data.csv"
)
checks = []
for scope, evidence in [("all", "all_analysis_eligible"), ("direct", "direct_only")]:
    counts = pd.read_csv(src / f"v14-artifact/{scope}/h1/all_data_probability_counts.csv.gz")
    scores = pd.read_csv(src / f"all-data-primary/syndrome/{scope}/island_syndrome_scores.csv.gz")
    result, summary = h2.run_decomposition(counts, scores, cov, cfg, evidence_scope=evidence)
    frozen = pd.read_csv(src / f"v14-artifact/{scope}/h2/h2_decomposition_models.csv")
    merged = result.merge(
        frozen, on=["context", "response"], suffixes=("_replay", "_frozen"), validate="one_to_one"
    )
    assert len(merged) == 12
    assert (merged.n_islands_replay == merged.n_islands_frozen).all()
    for col in ["distance_estimate", "distance_se", "distance_p", "primary_H2b_q"]:
        x = merged[col + "_replay"]
        y = merged[col + "_frozen"]
        assert np.allclose(x, y, atol=1e-4, rtol=1e-4, equal_nan=True), (scope, col)
    result.to_csv(out / f"h2_original_reproduced_{scope}.csv", index=False)
    checks.append(
        {
            "scope": scope,
            "rows": len(result),
            "max_coefficient_difference": float(
                abs(merged.distance_estimate_replay - merged.distance_estimate_frozen).max()
            ),
            "n_matched": True,
        }
    )
    print(checks[-1], flush=True)
(out / "h2_original_reproduction.json").write_text(json.dumps(checks, indent=2), encoding="utf8")
