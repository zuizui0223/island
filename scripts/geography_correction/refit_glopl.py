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

from island_v2 import (
    chapter1_h5_glopl_global_distance as g,
)
from island_v2 import (
    chapter1_v14_h4_exact_h2_score_bridge as h,
)

cfg = yaml.safe_load(
    (src / "island--chapter1_h5_glopl_global_distance_v1.yml").read_text(encoding="utf8")
)
hcfg = yaml.safe_load(
    (src / "island--chapter1_v14_h4_exact_h2_score_bridge.yml").read_text(encoding="utf8")
)
distance = (
    pd.read_csv(out / "glopl_corrected_site_distances.csv")
    .set_index("site_key")
    .spherical_distance_km
)
rows = pd.read_csv(src / "glopl-global/out/ELIGIBLE_EFFECT_ROWS.csv.gz")
lock = json.loads((src / "glopl-global/out/RESULT.json").read_text(encoding="utf8"))
scores = pd.read_csv(src / "all-data-primary/syndrome/direct/species_syndrome_concordance.csv.gz")
ra = pd.read_csv(src / "h4-ra/out/MATCHED_EFFECT_ROWS.csv.gz")
arch = pd.read_csv(src / "h4-arch/out/MATCHED_EFFECT_ROWS.csv.gz")
original_h4, _ = h.run_exact_h2_bridge(scores, ra, arch, hcfg)
frozen = pd.read_csv(src / "v14-artifact/h4-exact/h4_exact_h2_score_bridge_results.csv")
match = original_h4.merge(
    frozen, on=["family", "analysis"], suffixes=("_replay", "_frozen"), validate="one_to_one"
)
for col in ["estimate", "se", "two_sided_p"]:
    assert np.allclose(match[col + "_replay"], match[col + "_frozen"], atol=1e-10, rtol=1e-8)
assert len(match) == 6
original_h4.to_csv(out / "h4_exact_original_reproduced.csv", index=False)
results = {}
for mode in ["original", "corrected"]:
    r = rows.copy()
    if mode == "corrected":
        r["log1p_distance_to_major_continent_km"] = np.log1p(r.site_key.map(distance))
        assert r.log1p_distance_to_major_continent_km.notna().all()
    cells = g.aggregate_measurement_cells(r, cfg)
    cells, mean, sd = g._standardize_distance(cells)
    result = g.fit_global_gradient(cells, cfg)
    if mode == "original":
        assert abs(result["distance_slope"] - lock["global_gradient"]["distance_slope"]) < 1e-12
        assert (
            abs(result["distance_slope_se"] - lock["global_gradient"]["distance_slope_se"]) < 1e-12
        )
    sensitivities = {}
    for name, mask in {
        "supplemental_only": r.PL_Effect_Size_Type2.astype(str).eq("Sup"),
        "no_zero_constant": ~g._truthy(r.Constant_added),
    }.items():
        sc = g.aggregate_measurement_cells(r[mask], cfg)
        sc = g._reuse_distance_scaling(sc, mean, sd)
        sensitivities[name] = g.fit_global_gradient(sc, cfg)
        if mode == "original":
            assert (
                abs(
                    sensitivities[name]["distance_slope"]
                    - lock["sensitivities"][name]["global"]["distance_slope"]
                )
                < 1e-12
            )
    results[mode] = {
        "global_gradient": result,
        "sensitivities": sensitivities,
        "distance_standardization": {"mean": mean, "sd": sd},
    }
    if mode == "corrected":
        cells.to_csv(out / "h3_corrected_measurement_cells.csv.gz", index=False)
        r.to_csv(out / "h3_corrected_effect_rows.csv.gz", index=False)
        for frame in [ra, arch]:
            frame["distance_to_major_continent_km"] = frame.site_key.map(distance)
            assert frame.distance_to_major_continent_km.notna().all()
            frame["log1p_distance_to_major_continent_km"] = np.log1p(
                frame.distance_to_major_continent_km
            )
            frame["z_distance"] = (frame.log1p_distance_to_major_continent_km - mean) / sd
        corrected_h4, summary = h.run_exact_h2_bridge(scores, ra, arch, hcfg)
        corrected_h4.to_csv(out / "h4_exact_corrected.csv", index=False)
        corrected_h4.merge(
            frozen,
            on=["family", "analysis"],
            suffixes=("_corrected", "_original"),
            validate="one_to_one",
        ).to_csv(out / "h4_exact_comparison.csv", index=False)
        ra.to_csv(out / "h4_ra_corrected_effect_rows.csv.gz", index=False)
        arch.to_csv(out / "h4_arch_corrected_effect_rows.csv.gz", index=False)
        results["H4_corrected_manifest"] = summary
(out / "h3_original_corrected_comparison.json").write_text(
    json.dumps(results, indent=2), encoding="utf8"
)
print(json.dumps(results, indent=2))
