# Chapter 1 matched hierarchical-depth result — 2026-09-13

## Decision

The matched re-expression of the frozen final PR142 dataset supports a **bounded hierarchical-depth signal**, not a universal taxonomic-depth rule.

The same four-component floral-architecture response family was compared before and after source-matched genus adjustment in northern-midlatitude and tropical contexts. The formal statistic was the multivariate context difference in the distance slope of `beyond_genus - pre_genus`, with spatial-block cluster-robust covariance.

The result is strongest in **native non-endemic assemblages**: all four frozen source definitions pass BH-FDR in both all-analysis and direct-only evidence. The result is weaker when island endemics are pooled into `all_native`: all-analysis supports 0/4 source modes after FDR and direct-only supports 3/4.

Therefore the supported generalisation is bounded to widespread native assemblages. It does not establish a universal hierarchy across all island floras.

## Frozen input and execution

- final PR142 workflow run: `34232450884`
- final PR142 artifact: `10058653212`
- final PR142 digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- matched-depth audit run: `34763034251`
- audit artifact: `10319522647`
- audit artifact digest: `sha256:2c5fb8ac5dcce158ede148eb580eedf5494e1c181ab07eb38fda4a62270e017d`
- synthetic tests: 3 passed
- Ruff: passed
- workflow: green

The first workflow attempt failed only because a Typer single-command CLI was invoked with an extra `run` token. No biological output was opened by that failed invocation. The workflow call was corrected without changing inputs, model, thresholds or scientific rules; the next run completed successfully.

## Matched multivariate test

| evidence scope | floristic stratum | source mode | chi-square | df | p | q | supported |
|---|---|---|---:|---:|---:|---:|---|
| all-analysis | all-native | geo_k5 | 5.602 | 3 | 0.1327 | 0.1529 | no |
| all-analysis | all-native | geo_k10 | 6.279 | 3 | 0.0988 | 0.1529 | no |
| all-analysis | all-native | geo_k20 | 5.272 | 3 | 0.1529 | 0.1529 | no |
| all-analysis | all-native | geo50_climate10 | 8.472 | 3 | 0.0372 | 0.1488 | no |
| all-analysis | native-nonendemic | geo_k5 | 13.623 | 3 | 0.00347 | 0.00462 | yes |
| all-analysis | native-nonendemic | geo_k10 | 14.424 | 3 | 0.00238 | 0.00462 | yes |
| all-analysis | native-nonendemic | geo_k20 | 9.662 | 3 | 0.0217 | 0.0217 | yes |
| all-analysis | native-nonendemic | geo50_climate10 | 13.681 | 3 | 0.00337 | 0.00462 | yes |
| direct-only | all-native | geo_k5 | 11.909 | 3 | 0.00770 | 0.0103 | yes |
| direct-only | all-native | geo_k10 | 12.116 | 3 | 0.00700 | 0.0103 | yes |
| direct-only | all-native | geo_k20 | 7.416 | 3 | 0.0598 | 0.0598 | no |
| direct-only | all-native | geo50_climate10 | 14.576 | 3 | 0.00222 | 0.00887 | yes |
| direct-only | native-nonendemic | geo_k5 | 11.274 | 3 | 0.0103 | 0.0138 | yes |
| direct-only | native-nonendemic | geo_k10 | 12.900 | 3 | 0.00486 | 0.00972 | yes |
| direct-only | native-nonendemic | geo_k20 | 7.850 | 3 | 0.0492 | 0.0492 | yes |
| direct-only | native-nonendemic | geo50_climate10 | 13.142 | 3 | 0.00434 | 0.00972 | yes |

The robust covariance of the four architecture components has rank three, so the joint Wald statistic has 3 df. This is expected from the strong shared architecture structure already identified by V4; it should not be rewritten as four independent pollination-syndrome tests.

## Endemicity does not supply a confirmatory time axis

The latest status-stratified support audit shows that tropical endemic assemblages remain below the frozen confirmatory support threshold for the focal broad responses:

- generalized form: 30 islands;
- plain colour: 46 islands;
- self compatibility: 29 islands.

The tropical endemic WHEN/WHERE omnibus is unsupported (`p = 0.204481`). Endemicity is therefore retained only as a pilot distribution-history / evolutionary-opportunity stratum. It is not used as a temporal proxy or as a causal separation of assembly from in-situ evolution.

## Interpretation

The correct new statement is:

> **Within widespread native non-endemic island assemblages, the taxonomic attenuation of isolation-associated floral architecture differs between northern-midlatitude and tropical contexts.**

More concretely, source-matched genus adjustment changes the isolation-associated architecture vector differently in the two regions, and that context-by-taxonomic-stage difference is reproduced across all four source definitions and both evidence scopes in native non-endemics.

The stronger universal statement is not supported:

> taxonomic response depth differs among all island floras.

Pooling endemics into `all_native` weakens the result, especially in all-analysis evidence. That attenuation is part of the result and must remain visible in the manuscript.

## Publication consequence

This result creates a defensible route to an Ecology Letters first submission **only with a bounded general claim**. The paper should not be sold as a new universal island rule. Its broader contribution is that ecological filtering need not be expressed at the same taxonomic depth across biogeographic contexts, demonstrated here in widespread native island assemblages under matched source controls.

Global Ecology and Biogeography remains the strong fallback because the island-specific biogeographic branching and source/genus assembly result is already sufficient for that scope.

## Claim boundary

This audit does not identify:

- pollinator identity;
- historical pollinator loss;
- functional replacement;
- effective pollination service;
- in-situ evolution;
- a temporal sequence;
- a causal mechanism for the taxonomic-depth contrast.

It is a matched secondary re-expression of a previously frozen dataset. It should be reported as such rather than described as a newly preregistered outcome-blind primary test.
