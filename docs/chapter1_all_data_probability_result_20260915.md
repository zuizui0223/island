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

## Model-form sensitivity: the probability result is not a beta-binomial artefact

A matched sensitivity was run with the repository's existing grouped-binomial logit model, keeping the **same six atomic outcomes, same island sets, same covariates and same spatial-block robust inference** as the beta-binomial analysis.

Provenance:

- workflow: `Run Chapter 1 all-data model-form sensitivity`;
- run: **34963020944**;
- artifact: `chapter1-all-data-model-form-34963020944`;
- artifact ID: **10393863668**;
- digest: `sha256:dbd5b20b09659404470c5e706449fe90bb870878aed3fc3716037b7fedcd6507`.

Formal all-observed North--Tropical six-outcome test:

- all-analysis evidence: **3,626 islands**, 187 blocks, chi-square = **21.6768**, df = 6, **p = 0.001385**;
- direct-only evidence: **3,569 islands**, 187 blocks, chi-square = **18.8079**, df = 6, **p = 0.004501**.

These closely reproduce the beta-binomial decision. The broad all-observed difference is therefore **not created by the beta-binomial likelihood**.

## Response-definition audit: six atomic outcomes and the two composite axes are not interchangeable

The all-observed branch was also rerun through the existing two-axis composite route (`accessibility_generalization` + `reproductive_assurance`) using the same final trait snapshot.

Provenance:

- workflow: `Run Chapter 1 all-data primary candidate`;
- run: **34962433803**;
- artifact: `chapter1-all-data-primary-34962433803`;
- artifact ID: **10393997800**;
- digest: `sha256:b6295be7d4bce20f7ca3dddeb4a58106a647635023c6e7ab8fda025108470ff1`.

For the **two composite axes**, the all-observed North--Tropical vector difference is not supported:

- all-analysis evidence: **3,253 islands**, 184 blocks, p = **0.73998**, BH q = **0.88797**;
- direct-only evidence: **3,173 islands**, 183 blocks, p = **0.06142**, BH q = **0.13094**.

By contrast, the matched six-atomic probability vector remains supported under both beta-binomial and grouped-binomial models.

This isolates the source of the discrepancy: **response compression, not probability-model choice**. Collapsing multiple floral and reproductive traits into two composite scores discards a detectable part of the North--Tropical multivariate difference.

The all-observed atomic difference is concentrated in a subset of components rather than a clean reversal of one two-dimensional syndrome. In the beta-binomial all-analysis comparison, the larger Tropical-minus-Northern distance responses are strongest for actinomorphic symmetry, selfing mating system and autonomous selfing, with generalized form near the single-axis boundary; shallow/open tube contributes little. Direct-only evidence gives the clearest differences for generalized form and autonomous selfing.

Therefore the expanded route should not be described as evidence that the old two-axis syndrome simply becomes stronger with more islands. It supports a more specific statement:

> **Across thousands of observed island floras, isolation-associated floral and reproductive composition differs among biogeographic contexts at the multivariate atomic-trait level, but that heterogeneity is partly lost when compressed into the current two composite axes.**

## Decision

The all-data probability expansion is useful and should be retained because it establishes four things:

1. the available global data support a broad contemporary observed-flora response pattern across thousands of islands;
2. that pattern survives restriction from all-analysis evidence to High/Medium direct trait evidence;
3. the six-outcome result survives replacement of beta-binomial by the matched grouped-binomial model, so model family is not driving the conclusion;
4. the native-flora inference and the choice of response representation are the present bottlenecks, not raw trait coverage or sample size alone.

However, the broad result must currently be described as **observed island-flora composition**, not as defended native colonisation filtering or in-situ floral evolution. It also should not be collapsed back into the old two-axis syndrome without acknowledging that the compression removes supported multivariate structure.

The frozen v10 Chapter 1 result therefore remains canonical until the floristic-status problem and response-representation decision are addressed. The next gate is to improve or explicitly model status uncertainty and decide whether the paper should promote the six-atomic response geometry or retain the two-axis compression as a secondary summary.
