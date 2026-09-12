# Chapter 1 NEE N1 execution wrapper

## Role

This wrapper executes the already frozen N1 model. It does not alter the N1 estimand,
channel gate, model formula, global heterogeneity test, robustness criteria, or pass rule.

The execution order is:

```text
validated geo_k5 source + island observation
        -> retained / disrupted projection
        -> confirmatory channel qualification
        -> canonical Chapter 1 covariates
        -> frozen N1 common and heterogeneous retention models
        -> joint channel x isolation Wald test
        -> leave-one-spatial-block-out
        -> leave-one-raw-source-region-out
        -> N1 gate receipt
```

If fewer than three channels reach the frozen confirmatory input tier, the wrapper writes
an `N1_not_evaluable_for_NEE_gate` receipt and does not fit a rescue model.

## Canonical covariate freeze

N1 reuses the same geography/covariate artifact consumed by the canonical Chapter 1
progressive analysis:

- workflow run: `29228212586`
- artifact: `purpose-shortest-distance-regime-29228212586`
- artifact ID: `8270544465`
- artifact ZIP SHA-256: `695f35b97bae07e81b05deab537dd73fa687b2d99e9efdb1cb3babd2fa12dfb6`
- file: `results/purpose_shortest_island_data.csv`

The table has 8,265 unique islands and 291 frozen 10-degree spatial blocks. Climate
PC1-PC4 each have 46 missing island rows; those are handled by the already frozen N1
complete-case rule.

The historical artifact uses `log_distance_to_continent_km` and
`log_island_area_km2`. The runner verifies row-by-row that these are exactly
`log1p(distance_to_continent_km)` and `log1p(area_km2)` before renaming them to the N1
contract names `log1p_distance_to_continent_km` and `log_area`. A mismatch is a hard
stop. No new geography, distance metric, area transform, climate PCA, or spatial block
is estimated by this wrapper.

## Required robustness

A qualified N1 fit always runs both frozen deletion analyses:

1. leave one `spatial_block` out;
2. leave one raw `source_region_id` out.

The raw source-region deletion is used rather than the internally pooled model fixed-
effect label. The frozen gate requires reversal fractions no greater than 20% for both.

## Outputs

- `N1_model_support.csv`
- `N1_coefficients.csv`
- `N1_channel_slopes.csv`
- `N1_global_heterogeneity_test.json`
- `N1_model_comparison.json`
- `N1_leave_one_spatial_block_out.csv`
- `N1_leave_one_source_region_out.csv`
- `N1_gate_receipt.json`
- `N1_run_metadata.json`

If N1 fails, the receipt's required action is
`stop_before_N2_and_keep_frozen_Chapter1`. Only an N1 pass may return `N2_may_open`.

## Leakage boundary

The runner accepts projected pollinator-channel states, the outcome-blind qualification
receipt, and frozen geography/covariates only. It does not read floral traits, plant
assemblage outcomes, genus-entry outcomes, or N2 evidence.
