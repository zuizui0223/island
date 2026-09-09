# Chapter 1 paper pipeline

This is the shortest reproducibility map from the database to the paper.

## Pipeline at a glance

```text
[1] island + flora + source/status infrastructure
                    +
[2] provenance-preserving 3-axis trait database
                    |
                    v
[3] frozen trait snapshot
                    |
                    v
[4] H1-H5 scientific contract
                    |
                    v
[5] progressive Chapter 1 analysis
                    |
                    v
[6] V1-V5 validation / explanation-gap layers
                    |
                    v
[7] frozen interpretation + v7 manuscript
```

## 1. Database inputs

See `docs/DATABASE_BUILD.md`.

Paper-level fixed quantities:

- 8,265-island universe;
- 106,295 analysis-applicable species;
- 318,885 species-axis denominator;
- final snapshot Run `34191508045` / `source-scale-batch-integration-34191508045`;
- 222,688 resolved cells.

## 2. Freeze the trait snapshot

The progressive workflow validates the incoming species-axis ledger before fitting the paper analysis. Snapshot rules include fixed denominator, required columns, evidence quality, provenance retention, and explicit handling of audited evidence revision.

Implementation:

- `src/island_v2/chapter1_trait_snapshot.py`
- `src/island_v2/chapter1_trait_resolution_mnar.py`

## 3. Apply the frozen H1-H5 contract

Canonical config:

- `config/chapter1_progressive_analysis.yml`

Hypotheses are evaluated in order:

### H1 — Universal syndrome rival

Test whether the same coherent floral/reproductive response vector appears across contexts. Direct between-context heterogeneity falsifies the universal form.

### H2 — Biogeographic branching

Test within-context multivariate responses and direct between-context vector differences, including native-nonendemic persistence and observation/source sensitivities.

### H3 — Source / lineage assembly

Decompose the observed response through taxonomic depth and source-matched lineage representation. The final Palearctic result is retained through family adjustment but disappears after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`).

Implementation includes:

- `src/island_v2/chapter1_taxonomic_depth_decomposition.py`
- `src/island_v2/chapter1_pr138_lineage_representation_bridge.py`

### H4 — Area / capacity moderation

Test continuous distance x continuous area as a modifier and explicitly audit whether heteroskedastic measurement alone can generate apparent moderation.

Implementation:

- `src/island_v2/chapter1_area_capacity_moderation.py`
- `src/island_v2/chapter1_area_support_artifact.py`

### H5 — Channel-gated residual mechanism

Mechanistic promotion requires independent source-channel exposure, retained/disrupted/structurally-absent contrast, effort-aware or functional channel measurement, and added explanatory value beyond source/lineage/observation controls.

The final Chapter 1 plant database does **not** satisfy that evidence chain, so H5 is not identified.

## 4. Secondary floral-architecture concordance

`large_bee_like`, `butterfly_like`, and `bird_like` are post-primary plant architecture scores. They are decomposed to quantify their shared architecture component and are never treated as realised pollinator identities.

Implementation:

- `src/island_v2/chapter1_pollination_architecture_factor.py`

## 5. Mandatory robustness / validation layers

The canonical contract requires, among other checks:

- direct-only versus all analysis-eligible evidence;
- equal/capped island information;
- trait-resolution selection adjustment / MNAR stress;
- source-pool sensitivity;
- alternative distance transforms;
- leave-one-spatial-block-out checks;
- genus-composition guardrail;
- continuous area moderation;
- climate-overlap diagnostic.

Current validation modules include:

- `src/island_v2/chapter1_climate_overlap_validation.py`
- `src/island_v2/chapter1_trait_resolution_mnar.py`

The frozen explanation-gap sequence is summarized by the final checkpoints:

- `docs/chapter1_v1_climate_overlap_checkpoint.md`
- `docs/chapter1_v2_taxonomic_depth_checkpoint.md`
- `docs/chapter1_v3_area_support_checkpoint.md`
- `docs/chapter1_v4_architecture_decomposition_checkpoint.md`
- `docs/chapter1_v5_mnar_tipping_point_checkpoint.md`

These V1-V5 labels are validation stages of the **current Chapter 1 pipeline** and are unrelated to the historical repository `legacy/v1/`.

## 6. Canonical execution

Single paper-level workflow:

- `.github/workflows/run-chapter1-progressive-trait-analysis.yml`

Final execution:

- Run `34232450884`;
- artifact `chapter1-progressive-analysis-34232450884`;
- artifact ID `10058653212`;
- SHA-256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

## 7. Frozen paper surface

Read in this order:

1. `docs/chapter1_submission_freeze_20260909.md`
2. `docs/chapter1_submission_hypothesis_framework_20260909.md`
3. `docs/chapter1_figure1_hypothesis_tree_spec_20260909.md`
4. `docs/chapter1_h3_h5_causal_hierarchy_20260909.md`
5. `docs/chapter1_literature_positioning_20260909.md`
6. `docs/chapter1_manuscript_full_v7_submission_order_20260909.md`

## 8. Legacy v1 boundary

`legacy/v1/` is a frozen historical analysis, not a validation stage of the current pipeline. It remains intact and independently inspectable.

Tree recorded before this cleanup:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

No current Chapter 1 claim should silently mix v1 outputs into the v2 estimand. If a v1 result is discussed historically, it must be labelled as v1 provenance.
