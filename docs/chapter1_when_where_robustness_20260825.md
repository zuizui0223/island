# Chapter 1 when / where robustness — 2026-08-25

## Purpose

This checkpoint tests whether the canonical Chapter 1 when/where result depends on the original latitude grouping, log-distance functional form, zero-distance islands, or inclusion of weakly oceanic islands.

The canonical headline is:

> isolation-associated floral/reproductive filtering is detectable in northern mid-latitude and tropical island floras, persists in native non-endemic assemblages, and the northern-versus-tropical multivariate isolation-response vectors differ.

These robustness scenarios were declared before inspecting their outputs. No pollinator variable enters any scenario.

## Run provenance

- workflow run: `32838910957`
- artifact: `chapter1-when-where-main-32838910957`
- artifact id: `9559815910`
- artifact digest: `sha256:4a0f9a6f1b8e9b753582ce65b1e899bdd21525464aecbaa1fe87daa0f540a8e9`
- direct-trait input run: `32702160934`
- floristic-status input run: `32559322028`
- isolation/context input run: `29228212586`
- genus-fixed null draws: `1000`

## Scenarios

Five scenarios were evaluated for both `all_native` and `native_nonendemic` strata:

1. `canonical_log_all` — original context grouping and log-distance;
2. `core_latitudes_log_all` — buffered core latitude groups (`|lat|<20` tropical, 30–60° northern mid-latitude, >=70° northern high-latitude, 30–60° southern extratropical; transition bands unresolved);
3. `canonical_sqrt_all` — square-root rather than log-distance;
4. `canonical_log_positive_distance` — remove zero-distance islands;
5. `canonical_log_oceanic_50km` — retain only islands at least 50 km from a major continent.

For each scenario the headline is replicated only when all three confirmatory tests are supported:

- northern-mid-latitude WHERE vector != 0;
- tropical WHERE vector != 0;
- northern-mid-latitude vector != tropical vector.

A scenario that cannot meet the declared >=50-island response-support threshold is classified as not testable, not as biological contradiction.

## Result

All **10 of 10** scenario × status-stratum opportunities were confirmatorily testable in this run, and all **10 of 10** replicated the complete headline result.

| scenario | stratum | north WHERE q | tropical WHERE q | north vs tropical q | common responses | headline |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| canonical log, all islands | all native | 1.76e-39 | 2.52e-27 | 2.35e-08 | 17 | replicated |
| canonical log, all islands | native non-endemic | 7.43e-39 | 7.43e-39 | 7.13e-07 | 17 | replicated |
| buffered core latitudes | all native | 4.21e-22 | 5.17e-27 | 2.57e-09 | 16 | replicated |
| buffered core latitudes | native non-endemic | 2.81e-28 | 2.81e-34 | 1.42e-07 | 16 | replicated |
| sqrt distance | all native | 3.82e-61 | 2.99e-43 | 8.40e-16 | 17 | replicated |
| sqrt distance | native non-endemic | 9.09e-60 | 2.14e-30 | 2.24e-11 | 17 | replicated |
| positive-distance only | all native | 4.10e-74 | 2.95e-33 | 2.31e-05 | 16 | replicated |
| positive-distance only | native non-endemic | 1.62e-54 | 1.02e-31 | 1.68e-03 | 16 | replicated |
| oceanic >=50 km | all native | 4.38e-125 | 4.38e-25 | 1.20e-19 | 16 | replicated |
| oceanic >=50 km | native non-endemic | 5.05e-147 | 1.04e-34 | 2.62e-15 | 16 | replicated |

The weakest replicated between-region result is the positive-distance-only `native_nonendemic` comparison (`q=0.00168`), still well below the declared 0.05 threshold.

## Interpretation

The current north/tropical when/where conclusion is not an artifact of:

- the exact canonical latitude cut points;
- a log-distance rather than square-root isolation transform;
- zero-distance islands; or
- inclusion of islands closer than 50 km to a continent.

This materially strengthens the inference that isolation-associated filtering is detectable in both northern mid-latitude and tropical island floras and that the two response vectors differ.

## Remaining statistical robustness issue

The primary WHERE tests use cluster-robust joint Wald tests with:

- northern mid-latitude: 36 spatial clusters for a 21-dimensional vector;
- tropical: 46 spatial clusters for a 17-dimensional vector;
- northern-vs-tropical comparison: 82 spatial clusters for a 17-dimensional contrast vector.

The sensitivity results above address geographic definitions and isolation specification, but they do not by themselves eliminate finite-cluster/high-dimensional Wald concerns. A block-level resampling or cluster-jackknife sensitivity is therefore the next statistical robustness check before treating the p-values as fully submission-ready.

## Claim boundary

This robustness analysis changes confidence in **when/where**, not mechanism. It provides no evidence that Bombus, birds, Lepidoptera, or another pollinator guild caused either regional response vector.
