# Chapter 1 figure and table plan

Every display must be generated from the artifact identified by
`config/chapter1_submission_freeze_lock.json`. Do not manually copy values from development
checkpoints when a locked CSV exists.

## Main figures

### Figure 1. Hypothesis and observable boundary

Content:

```text
distance / source accessibility x area / capacity
|-- source-pool lineage entry and loading
`-- latent pollination-channel filter
     |-- reproductive assurance
     `-- attraction / access
```

Visually separate observed plant/geographic variables from unmeasured pollinator processes.
Bombus, bird, and Lepidoptera should appear only as non-exhaustive Discussion examples.

Source: conceptual, constrained by `docs/manuscript_submission_contract.md`.

### Figure 2. Global WHEN / WHERE response vectors

Panels:

- within-context atomic vector support by operational context and floristic stratum;
- direct northern-midlatitude versus tropical vector difference;
- formal Palearctic, Nearctic, and Neotropical branch summary;
- sample/support annotations showing unresolved regions.

Locked sources:

- `frozen/frozen/input/canonical/when_where/when_where_within_context_omnibus.csv`;
- `frozen/frozen/input/canonical/when_where/when_where_between_context_omnibus.csv`;
- `frozen/frozen/input/input/source/current/primary/observed_within_outcome_slopes.csv`;
- embedded realm sensitivity outputs.

### Figure 3. Source-matched lineage entry versus loading

Panels:

- broad Palearctic entry slopes across four source modes and two strata;
- broad loading increments on the same scale;
- direct north-tropical and Palearctic-Neotropical entry contrasts;
- exact-SI entry/loading comparison as a bounded sensitivity.

Locked sources:

- `frozen/result/lineage_representation_context_slopes.csv`;
- `frozen/result/lineage_representation_between_context.csv`;
- `frozen/result/lineage_representation_island_scores.csv.gz`.

### Figure 4. Continuous area moderation

Panels:

- Palearctic broad genus-entry distance slopes at area z = -1, 0, +1;
- source-mode x stratum interaction coefficients;
- Neotropical reproductive-assurance bounded result;
- northern-midlatitude versus tropical accessibility interaction difference.

Locked sources:

- `result/area_moderation_coefficients.csv`;
- `result/area_moderation_within_context.csv`;
- `result/area_moderation_between_coefficients.csv`;
- `result/area_moderation_between_context.csv`.

Avoid labelling z = -1 and z = +1 as discrete small/large island categories. They are
illustrative conditional slopes from a continuous model.

## Main tables

### Table 1. Data support and estimands

Report full-universe attrition, response-specific support rule, floristic strata, geographic
contexts, and the conditional observed-flora estimand.

### Table 2. Frozen primary tests

Report within-context vector tests, direct between-context contrasts, retained response
counts, islands, spatial blocks, q-values, and persistence classification. Do not replace
the direct contrast with separate significance labels.

### Table 3. Source-pool and area decomposition

Report Palearctic entry/loading counts, exact-SI boundary, area interaction states,
Neotropical sample size, and unsupported direct mechanisms.

## Supplementary figures and tables

- Figure S1: observation and trait-resolution coverage;
- Figure S2: information-weight sensitivity;
- Figure S3: distance functional forms;
- Figure S4: leave-one-block-out influence;
- Figure S5: source-mode x IPW matrix;
- Figure S6: exact-SI genus-fixed null;
- Figure S7: genus influence deletion;
- Figure S8: evidence/template sensitivities;
- Table S1: trait ontology and evidence tiers;
- Table S2: response-specific island and cluster support;
- Table S3: all source-matching scenarios;
- Table S4: all area-moderation coefficients and hierarchical states;
- Table S5: preserved negative, pilot, unresolved, and not-evaluable outcomes.

## Visual QA gate

Before submission:

1. regenerate every display from the locked artifact in a clean environment;
2. confirm all panel counts and headline values against the freeze verification report;
3. inspect labels at final journal dimensions;
4. use colour-blind-safe palettes and redundant shape/line encodings;
5. include island and spatial-block support in captions;
6. verify that no caption crosses the causal claim ceiling.
