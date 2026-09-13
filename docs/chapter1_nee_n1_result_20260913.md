# Chapter 1 NEE N1 canonical result — 2026-09-13

## Status

The frozen prospective N1 gate failed. Per the pre-outcome contract, H5a stops before N2. No reduced-channel rescue, channel substitution, Search rerun, or post-hoc threshold change is permitted.

## Upstream qualification

Qualification run `34750153702` admitted four confirmatory channels: Bombus, non-Bombus bees, flower-visiting birds, and Diptera. Lepidoptera was not qualified because it had zero disrupted islands under the frozen exhaustive non-detection rule.

## Canonical N1 run

- workflow run: `34750351009`
- artifact: `10315671249`
- artifact SHA256: `01d938a1f42b7f64634232a2fcb76c68a60183824684c5f135f567107a4c5d85`
- primary rows: 6,256
- islands represented: 3,150
- spatial blocks: 225
- raw source regions: 217; model source regions: 145

The frozen joint cluster-robust Wald test of isolation × channel terms was not significant: `W = 1.6187`, `df = 3`, `p = 0.65516`. The heterogeneous-slope model improved log likelihood over the common-slope model by `2.1161`, but this does not satisfy the predeclared N1 gate without significant global heterogeneity.

Descriptive standardized isolation slopes were all negative: Bombus `-1.0823`, Diptera `-0.7915`, flower-visiting birds `-0.4926`, and non-Bombus bees `-0.1597` log-odds per SD. These are descriptive only because the global heterogeneity gate failed.

## Deletion robustness

All 225 leave-one-spatial-block-out fits and all 217 leave-one-source-region-out fits completed. Every deletion fit retained a global heterogeneity `p >= 0.05`; the minimum p-values were `0.1271` and `0.1255`, respectively. The stored `reversal_fraction = 1.0` therefore means 100% of deletion fits were non-significant, not that slope signs reversed.

## Required action

`N1_pass = false`; `failure_action = stop_before_N2_and_keep_frozen_Chapter1`; `claim_ceiling = N1_not_promoted`. N2 must remain closed for this H5a prospective chain.

The exact machine-readable evidence is frozen in `config/chapter1_nee_n1_result_lock.json`. The earlier exact-island Search transport limitation also remains in force: inference is limited to transport-evaluable exact island geometries, with large/complex-island Search errors retained as unresolved rather than rescued.
