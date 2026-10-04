"""Restore original poster coefficient inference; add WCVP under identical rules."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_all_data_probability import _fit_within, _prepare, build_broad_counts
from island_v2.chapter1_wcvp_native_compatibility import build_regional_native_compatible_flora


def original_poster_columns(frame):
    out = frame.copy()
    out["ci_low"] = out.geography_slope_log_odds - 1.96 * out.cluster_robust_se
    out["ci_high"] = out.geography_slope_log_odds + 1.96 * out.cluster_robust_se
    out["nominal_supported"] = out.p_value.lt(0.05)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, required=True)
    a = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    out = root / "results/h1_poster_original_wcvp_20261004"
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((root / "config/chapter1_h1_traitwise_20261004.yml").read_text())
    cfg.update(strata=["all_native"], minimum_outcomes_per_vector=2, max_iter=1500)
    flora_path = (
        a.inputs / "sources/progressive-input/fixed/canonical/input/chapter1_status_flora.csv.gz"
    )
    mapping_path = a.inputs / "traitwise-inputs/status/resolved/island_tdwg_l3_mapping.csv"
    ranges_path = (
        a.inputs / "traitwise-inputs/wcvp/ch1-wcvp/ranges/wcvp_native_range_summary.csv.gz"
    )
    flora = pd.read_csv(flora_path, dtype={"island_id": str})
    native, coverage = build_regional_native_compatible_flora(
        flora, pd.read_csv(ranges_path), pd.read_csv(mapping_path, dtype={"island_id": str})
    )
    covpath = root / "results/geography_20260924/corrected_geography_covariates.csv"
    cov = pd.read_csv(covpath, dtype={"island_id": str})
    frames = []
    support = []
    inputs = [flora_path, mapping_path, ranges_path, covpath]
    for scope in ["all", "direct"]:
        primary = root / f"results/geography_20260924/{scope}/beta_binomial_within_slopes.csv"
        inputs.append(primary)
        original = pd.read_csv(primary).query("stratum=='all_observed'").copy()
        original = original_poster_columns(original)
        original.insert(0, "flora_scope", "broad")
        original.insert(0, "evidence_scope", scope)
        frames.append(original)
        auditpath = (
            a.inputs
            / f"sources/progressive-input/snapshot/chapter1_trait_state_audit_{scope}.csv.gz"
        )
        inputs.append(auditpath)
        counts = build_broad_counts(native, pd.read_csv(auditpath), cfg)
        prepared = _prepare(counts, cov, cfg)
        for context in cfg["contexts"]:
            slopes, _ = _fit_within(
                prepared, stratum="all_native", context_value=context, threshold=50, config=cfg
            )
            slopes = original_poster_columns(slopes)
            slopes.insert(0, "flora_scope", "wcvp")
            slopes.insert(0, "evidence_scope", scope)
            frames.append(slopes)
            for trait, part in prepared.loc[prepared[cfg["context_column"]].eq(context)].groupby(
                "outcome"
            ):
                support.append(
                    {
                        "evidence_scope": scope,
                        "context": context,
                        "outcome": trait,
                        "n_islands": part.island_id.nunique(),
                        "n_spatial_blocks": part[cfg["cluster_column"]].nunique(),
                        "n_island_species_records": int(part.trials.sum()),
                    }
                )
            print(scope, context, "complete", flush=True)
    result = pd.concat(frames, ignore_index=True)
    assert len(result) == 112
    result.to_csv(out / "traitwise_results.csv", index=False)
    pd.DataFrame(support).to_csv(out / "wcvp_support.csv", index=False)
    receipt = {
        "rule": "Original poster: cluster sandwich covariance from original seven-trait stack; two-sided normal p; estimate +/- 1.96 SE; unadjusted individual p<.05. No Holm or t substitution.",
        "primary": "Original corrected tables copied without refitting or changing their p/SE/coefficients.",
        "wcvp_coverage": coverage,
        "unique_native_species_before_trait_filter": int(
            native.loc[native.origin_status.eq("native"), "accepted_species"].nunique()
        ),
        "all_optimizers_converged": bool(result.optimizer_success.all()),
        "inputs": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    }
    (out / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
