import os
from pathlib import Path

src = Path(os.environ["ISLAND_CORRECTION_SOURCES"])
out = Path(os.environ["ISLAND_CORRECTION_OUTPUT"])
data = Path(os.environ["ISLAND_CORRECTION_GEOMETRY"])
out.mkdir(parents=True, exist_ok=True)
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2 import chapter1_all_data_probability as h1
from island_v2 import chapter1_v14_h2_decomposition as h2

check = json.loads((data / "spherical_geometry_validation.json").read_text(encoding="utf8"))
assert check["status"] == "pass"
dist = pd.read_csv(data / "gshhg_spherical_distances_all.csv")
assert len(dist) == 8264 and dist.island_id.nunique() == 8264
assert (
    dist.spherical_coast_distance_km.gt(0).all()
    and np.isfinite(dist.spherical_coast_distance_km).all()
)
cov = pd.read_csv(
    src / "progressive-input/fixed/isolation/results/purpose_shortest_island_data.csv"
)
bad = "gshhg_2.3.7_h_8b13189234b949ee1ff6"
assert set(cov.island_id) - set(dist.island_id) == {bad}
cov = cov.merge(
    dist[["island_id", "spherical_coast_distance_km"]], on="island_id", validate="one_to_one"
)
cov["original_projected_distance_km"] = cov.distance_to_continent_km
cov["distance_to_continent_km"] = cov.spherical_coast_distance_km
cov["log_distance_to_continent_km"] = np.log1p(cov.spherical_coast_distance_km)
cov.to_csv(out / "corrected_geography_covariates.csv", index=False)
receipt = {
    "status": "geometry_gate_passed_before_trait_fit",
    "metric": "Minimum great-circle arc-to-arc distance on a mean-radius sphere (6371.0088 km), GSHHG 2.3.7 high resolution; not WGS84 ellipsoid",
    "retained_universe": 8264,
    "excluded_mainland_component": bad,
    "distance_zero_count": 0,
    "validation": "analytic boundary minimization, candidate pruning compared with exhaustive evaluation; independently optimized 24 cases; no unresolved boundary intersections in retained cohort",
    "distance_csv_sha256": hashlib.sha256(
        (data / "gshhg_spherical_distances_all.csv").read_bytes()
    ).hexdigest(),
    "limitations": "Spherical Earth approximation and underlying shoreline measurement accuracy remain. Post-hoc correction; no claim of a new confirmatory analysis.",
}
(out / "corrected_geometry_gate.json").write_text(json.dumps(receipt, indent=2), encoding="utf8")
c1 = yaml.safe_load(
    (src / "island--chapter1_v14_all_data_probability.yml").read_text(encoding="utf8")
)
c2 = yaml.safe_load((src / "island--chapter1_v14_h2_decomposition.yml").read_text(encoding="utf8"))
for scope, evidence in [("all", "all_analysis_eligible"), ("direct", "direct_only")]:
    dest = out / scope
    dest.mkdir(exist_ok=True)
    counts = pd.read_csv(src / f"v14-artifact/{scope}/h1/all_data_probability_counts.csv.gz")
    scores = pd.read_csv(src / f"all-data-primary/syndrome/{scope}/island_syndrome_scores.csv.gz")
    print("Starting corrected H2", scope, flush=True)
    result, summary = h2.run_decomposition(counts, scores, cov, c2, evidence_scope=evidence)
    result.to_csv(dest / "h2_decomposition_models.csv", index=False)
    (dest / "h2_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf8")
    print("Starting corrected H1", scope, flush=True)
    tables = h1.run_probability_analysis(counts, cov, c1)
    for table, name in zip(
        tables,
        [
            "beta_binomial_within_slopes",
            "beta_binomial_between_slopes",
            "beta_binomial_within_omnibus",
            "beta_binomial_between_omnibus",
        ],
    ):
        table.to_csv(dest / f"{name}.csv", index=False)
    print("Completed corrected H1/H2", scope, flush=True)
(out / "refit_completion.json").write_text(
    json.dumps(
        {
            "status": "completed",
            "scopes": ["all", "direct"],
            "H1_strata": c1["strata"],
            "model_specs_unchanged": True,
            "no_results_selected_out": True,
        },
        indent=2,
    ),
    encoding="utf8",
)
