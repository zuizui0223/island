# Chapter 1 all-data probability result — 2026-09-15

Status: **exploratory expansion completed; do not replace the frozen v10 submission result yet**.

## Execution provenance

- branch: `ch1-all-data-primary`
- workflow: `Run Chapter 1 all-data probability analysis`
- successful run: `34961775336`
- artifact: `chapter1-all-data-probability-34961775336`
- artifact ID: `10394245237`
- artifact digest: `sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b`
- frozen upstream Chapter 1 progressive artifact: run `34232450884`, artifact `chapter1-progressive-analysis-34232450884`

Both validation and analysis jobs completed successfully. All beta-binomial optimizers used in the reported fitted cells converged.

## Model hierarchy

Primary model:

- beta-binomial logit regression on island-level trait counts;
- six predeclared atomic outcomes;
- island-specific denominator = number of species with a resolved value for that outcome;
- no arbitrary minimum species-per-island cutoff;
- distance, island area and climate PC1--PC4 in the mean model;
- spatial-block cluster-robust sandwich covariance;
- direct North--Tropical difference tested by a joint Wald test of distance x context coefficients.

Evidence hierarchy:

1. `all_analysis_eligible` (High + Medium + trait-specific Validated Low) = primary;
2. `direct_only` (High + Medium) = evidence-quality sensitivity.

Flora hierarchy:

1. `all_observed` = broad observed-flora composition;
2. `all_native` = status-resolved sensitivity;
3. `native_nonendemic` = stricter status-resolved sensitivity.

## Expanded all-observed result

### All-analysis evidence

Within-context six-outcome vector support:

- northern mid-latitude: **2,164 islands**, 86 spatial blocks, joint p = **1.06e-6**;
- northern high latitude: **409 islands**, 63 blocks, p = **1.10e-6**;
- tropical: **1,462 islands**, 110 blocks, p = **1.56e-5**;
- southern extratropical: **299 islands**, 41 blocks, p = **2.09e-13**.

Formal North--Tropical interaction-vector test:

- **3,626 unique islands**;
- 187 spatial blocks;
- chi-square = **23.1011**, df = 6;
- **p = 0.0007633**.

Thus the broad observed-flora North--Tropical response difference is supported under the beta-binomial model.

### High/Medium direct-evidence sensitivity

Within-context vector support remains in all four broad contexts. The formal North--Tropical test uses **3,569 unique islands**, 187 spatial blocks and remains supported:

- chi-square = **18.8072**, df = 6;
- **p = 0.004502**.

Therefore the broad observed-flora result is not created solely by Validated-Low trait evidence.

## Floristic-status result changes the interpretation

The broad result does **not** transfer cleanly to status-resolved native flora.

All-analysis evidence, North--Tropical test:

- `all_observed`: p = **0.000763**;
- `all_native`: p = **0.2491**;
- `native_nonendemic`: p = **0.05097**.

Direct-only sensitivity:

- `all_observed`: p = **0.004502**;
- `all_native`: p = **0.3018**;
- `native_nonendemic`: p = **0.2079**.

The difference is not explained simply by the smaller number of status-resolved islands. Holding the island set fixed at the native-supported frame gives:

### All-analysis evidence, same island frame

- all observed species on those islands: **375 islands**, 82 blocks, **p = 0.0005607**;
- native species only on those same islands: **375 islands**, 82 blocks, **p = 0.2491**.

### Direct-only evidence, same island frame

- all observed species: **p = 0.0003959**;
- native species only: **p = 0.3018**.

Therefore the status dependence is not merely a sample-size or island-selection effect. Which species within the same islands are admitted by floristic status materially changes the inferred North--Tropical response difference.

## Unresolved-status diagnostic

On the same status-supported island frame, unresolved-status species alone remain informative:

- all-analysis evidence: unresolved-only North--Tropical joint p = **0.01914**;
- direct-only evidence: unresolved-only p = **0.03044**.

Introduced-only data do not meet the frozen 50-island-per-outcome support rule and are not testable as a six-outcome vector. The combined non-native-status diagnostic is supported for all-analysis evidence (p = **0.01939**) but is just above 0.05 with direct evidence (p = **0.05711**).

These diagnostics do **not** establish that unresolved records are introduced, nor that introduction causes the broad branch. They show that the expanded all-observed result is partly carried by records whose native status is unresolved and therefore cannot currently be promoted as a native island-assembly result.

## Decision

The beta-binomial expansion is useful and should be retained because it establishes three things:

1. the available global data support a broad contemporary observed-flora response pattern across thousands of islands;
2. that pattern survives restriction from all-analysis evidence to High/Medium direct trait evidence;
3. the native-flora inference is the present bottleneck, not raw trait coverage or sample size alone.

However, the broad result must currently be described as **observed island-flora composition**, not as defended native colonisation filtering or in-situ floral evolution.

The frozen v10 Chapter 1 result therefore remains canonical until the floristic-status problem is addressed. The next gate is to improve or explicitly model status uncertainty, then rerun the beta-binomial hierarchy before deciding whether the expanded route can replace the frozen native-status analysis.
