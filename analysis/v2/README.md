# Chapter 1 analysis execution surface

## Current canonical model contract

The current Chapter 1 hypothesis contrast is **coherent island syndrome vs component-specific floral reorganization**.

Current frozen evidence favors **component-specific reorganization** and does not support a confirmatory claim of biogeographic slope heterogeneity.

Canonical scientific definitions live in:

- `config/analysis_models.yml`
- `docs/chapter1_frozen_result_20260825.md`
- `docs/v2_pollination_regime_framework.md`
- `docs/manuscript_submission_contract.md`
- `THESIS_CHAPTER_POSITIONING.md`

The canonical model ladder is:

```text
M0 context-aware baseline
M1 common isolation component effects          [primary]
M2 isolation × context joint heterogeneity     [secondary boundary test]
M3 genus-fixed status / lineage guardrail      [required]
M4 atomic category-preserving decomposition    [primary interpretation]
```

Pollination-syndrome comparisons are performed only after M0–M4 results and formal heterogeneity classifications are frozen.

## Current canonical implementation

The active implementation is now on `main`:

- `src/island_v2/chapter1_context_input.py`
  - builds Bombus-free status-aware trait composition;
  - expands colour/form multistates to atomic category presences;
  - retains `SC|SI` once as `mixed_or_variable`;
  - fails closed on conflicting direct trait records.
- `src/island_v2/chapter1_context_analysis.py`
  - fits nested M0/M1/M2 grouped-binomial models;
  - uses cluster-robust covariance;
  - computes context-specific simple slopes;
  - computes the formal joint Wald test for all `isolation × context` interactions.
- `src/island_v2/genus_fixed_trait_null.py`
  - runs the vectorized genus-composition-preserving M3 null.
- `src/island_v2/status_stratified_lineage_analysis.py`
  - tests broad residual outcomes by status/context.
- `src/island_v2/chapter1_trait_vector_freeze.py`
  - applies declared FDR rules;
  - distinguishes within-context slope signals from supported regional heterogeneity;
  - freezes the component vector before syndrome interpretation.

The canonical workflow is `.github/workflows/run-chapter1-context-main.yml`.

## Frozen run

Canonical joint-Wald run: `32833362756`.

Key frozen counts:

- fitted atomic categories: 43
- context-specific slopes: 95
- FDR-supported within-context slopes: 17
- confirmatory-count-supported within-context slopes: 16, all northern-midlatitude
- FDR-supported joint contingency categories: 1
- confirmatory joint contingency categories: 0

See `docs/chapter1_frozen_result_20260825.md` for the scientific interpretation.

## Interpretation rule

Do not equate:

```text
significant in northern_midlatitude
+
nonsignificant in tropical
```

with a regional difference.

Regional heterogeneity requires the formal joint interaction test to survive the declared multiplicity correction and adequate context support.

The present manuscript-level result is therefore:

```text
isolation-associated component-specific floral reorganization
+
strong northern-midlatitude within-context signal concentration
+
no confirmatory regional heterogeneity
+
no confirmatory broad generalized/plain/SC syndrome after M3
```

## Legacy script in this directory

`run_bayesian_m0_m4_main.R` is a **historical Bombus-primary analysis**. It remains for reproducibility and sensitivity work but is not the current Chapter 1 manuscript workflow.

The associated GitHub Actions workflow remains manual-only legacy status.

Do not cite its old M0–M4 labels as the current model ladder.

## Protected boundaries

The canonical workflow must not:

- require `bombus_deficit` to enter the primary analysis universe;
- infer pollinator guild from flower colour/form;
- treat opportunistic absence as historical pollinator loss;
- treat one-region significance as evidence of interaction;
- collapse all categories to one syndrome index before decomposition;
- omit status/lineage safeguards from manuscript-level inference;
- call the current northern trait vector a demonstrated Bombus-loss syndrome;
- force bird/Lepidoptera counter-syndromes where formal regional contrasts are unsupported.
