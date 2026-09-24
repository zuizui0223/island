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

from island_v2 import chapter1_v13_colour_architecture_audit as arch
from island_v2 import chapter1_v13_raw_colour_audit as color
from island_v2 import chapter1_v13_raw_colour_coupling_audit as coupling

cfg = yaml.safe_load(
    (src / "island--chapter1_all_data_biogeographic_pattern.yml").read_text(encoding="utf8")
)
cfg["strata"] = ["all_observed"]
cfg["colours"] = list(color.FOCAL_COLOURS)
oldcov = pd.read_csv(
    src / "progressive-input/fixed/isolation/results/purpose_shortest_island_data.csv"
)
newcov = pd.read_csv(out / "corrected_geography_covariates.csv")
checks = []
for scope in ["all", "direct"]:
    scores = pd.read_csv(src / f"all-data-primary/syndrome/{scope}/island_syndrome_scores.csv.gz")
    base = src / f"raw-coupling/results/{scope}"
    tasks = [
        (
            "colour",
            "raw_colour_counts.csv.gz",
            "raw_colour_model_results.csv",
            color.fit_raw_colour_models,
            ["colour"],
        ),
        (
            "colour",
            "raw_colour_counts.csv.gz",
            "raw_colour_joint_omnibus.csv",
            color.fit_within_context_joint_omnibus,
            [],
        ),
        (
            "architecture",
            "raw_colour_architecture_counts.csv.gz",
            "raw_colour_architecture_model_results.csv",
            arch.fit_colour_architecture_models,
            ["combination"],
        ),
        (
            "coupling",
            "raw_colour_conditioned_architecture_counts.csv.gz",
            "raw_colour_conditioned_architecture_model_results.csv",
            coupling.fit_colour_conditioned_architecture_models,
            ["combination"],
        ),
    ]
    for folder, countfile, resultfile, fn, extra in tasks:
        print("Reproduce then correct", scope, resultfile, flush=True)
        counts = pd.read_csv(base / folder / countfile)
        frozen = pd.read_csv(base / folder / resultfile)
        original = fn(counts, scores, oldcov, cfg)
        keys = ["stratum", "context", "support_tier", "model", *extra]
        merged = original.merge(
            frozen, on=keys, suffixes=("_replay", "_frozen"), validate="one_to_one"
        )
        assert len(merged) == len(frozen)
        cols = (
            ["distance_estimate", "distance_se", "distance_p", "distance_q"]
            if "distance_estimate" in frozen
            else ["wald_chisq", "p_value"]
        )
        for col in cols:
            assert np.allclose(
                merged[col + "_replay"],
                merged[col + "_frozen"],
                atol=1e-6,
                rtol=1e-5,
                equal_nan=True,
            ), (scope, resultfile, col)
        new = fn(counts, scores, newcov, cfg)
        dest = out / scope / "raw_patterns"
        dest.mkdir(parents=True, exist_ok=True)
        new.to_csv(dest / resultfile, index=False)
        new.merge(
            frozen,
            on=keys,
            suffixes=("_corrected", "_original"),
            how="outer",
            validate="one_to_one",
        ).to_csv(dest / ("comparison_" + resultfile), index=False)
        checks.append(
            {"scope": scope, "result": resultfile, "original_reproduced": True, "rows": len(new)}
        )
(out / "raw_pattern_refit_completion.json").write_text(
    json.dumps(checks, indent=2), encoding="utf8"
)
print("All raw pattern refits complete.", flush=True)
