# Chapter 1 all-data probability result — 2026-09-15

Status: **exploratory expansion completed; integrated decision is now in `docs/chapter1_all_data_route_decision_20260915.md`. Frozen v11 remains the canonical native-assembly submission surface.**

## Primary analysis

The expanded route uses a beta-binomial logit model on six atomic island-level trait counts with no arbitrary minimum species-per-island cutoff. The mean model contains distance, island area and climate PC1--PC4, with spatial-block cluster-robust covariance. `all_analysis_eligible` (High + Medium + trait-specific Validated Low) is primary and `direct_only` (High + Medium) is the evidence sensitivity.

Broad observed-flora North--Tropical result:

- all-analysis: **3,626 islands / 187 blocks**, p = **0.0007633**;
- direct-only: **3,569 islands / 187 blocks**, p = **0.004502**.

Matched grouped-binomial sensitivity gives p = **0.001385** and **0.004501**, respectively, so the result is not created by the beta-binomial likelihood.

## Response representation

The historical species-level two-axis concordance route is weak on the expanded all-observed flora:

- all-analysis: p = **0.73998**;
- direct-only: p = **0.06142**.

However, projecting the six fitted interaction coefficients onto two predeclared family-level contrasts retains support:

- all-analysis: **p = 0.0016467**;
- direct-only: **p = 0.0120921**.

The two-family subspace retains **78.22%** and **64.82%** of Euclidean signal, respectively; discarded within-family contrasts are not supported. Relaxing the historical per-species `minimum_informative_traits` from 2 to 1 does not rescue the old composite score. A trait-first island-level family refit is also weak. Therefore the primary expanded estimand should remain the **six-atomic multivariate slope vector**, with family summaries derived as coefficient-level contrasts rather than as new species-level concordance scores.

## Floristic-status boundary

Source-backed native restriction does not retain the global six-atomic contrast:

- all-native: p = **0.2491** all-analysis / **0.3018** direct;
- native-nonendemic: p = **0.05097** all-analysis / **0.2079** direct.

This is not simply loss of islands: on the same 375 status-supported islands, all observed species retain p = **0.0005607** / **0.0003959**, whereas native species only give p = **0.2491** / **0.3018**.

An outcome-free status-resolution audit shows that origin-status resolution changes with distance within contexts, but the North--Tropical difference in that status-resolution slope is absent (**p = 0.849**). Thus status uncertainty is not reducible to a simple differential missingness gradient between the two focal contexts.

## WCVP regional-native compatibility expansion

A hash-pinned 2026 WCVP native-range sensitivity was added. Exact accepted WCVP matches exist for **101,101** of the fixed island-flora species, and **100,887** have an extant non-doubtful native TDWG-L3 range.

Unresolved source-status records were upgraded only when the focal island mapped unambiguously to a TDWG-L3 region included in that species' WCVP native range. This is regional-native compatibility, not exact island-native proof.

Coverage gain:

- **362,789** unresolved rows upgraded;
- **2,281** islands receive at least one upgrade;
- resulting regional-native-compatible flora: **513,320 rows / 2,372 islands**.

North--Tropical six-atomic test after this expansion:

- all-analysis: **1,848 analysis islands / 178 blocks**, p = **0.07954**;
- direct: **1,827 islands / 178 blocks**, p = **0.17472**.

So greatly expanding native-compatible coverage does not restore the broad observed-flora contrast.

## WCVP partition diagnostic

The broad result is concentrated outside the regional-native-compatible partition.

All-analysis:

- source-backed native: p = **0.24915**;
- WCVP-compatible unresolved: p = **0.06562**;
- WCVP-incompatible unresolved: p = **0.002010**;
- WCVP-unclassifiable unresolved: p = **0.01368**;
- source-introduced + WCVP-incompatible unresolved: p = **0.001102**.

Direct-only:

- source-backed native: p = **0.30184**;
- WCVP-compatible unresolved: p = **0.35251**;
- WCVP-incompatible unresolved: p = **0.003065**;
- WCVP-unclassifiable unresolved: p = **0.002687**;
- source-introduced + WCVP-incompatible unresolved: p = **0.009484**.

These results do **not** establish that every WCVP-incompatible or unclassifiable record is introduced. They establish that the thousands-island observational contrast cannot currently be promoted as a native island-assembly result.

## Decision

The all-data route is retained as the broad descriptive layer:

> **Across thousands of contemporary observed island floras, isolation-associated floral and reproductive composition differs among biogeographic contexts at the six-atomic-trait level.**

It survives direct-only evidence and model-form sensitivity.

For biological inference, however, the current global native-status gate is not passed. Frozen v11 therefore remains the defended native-assembly core, with the Palearctic genus-structured result retained as the strongest lineage-assembly evidence.

See `docs/chapter1_all_data_route_decision_20260915.md` for the full claim hierarchy and poster/paper consequences.
