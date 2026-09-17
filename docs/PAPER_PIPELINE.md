# Chapter 1 v13 paper pipeline

This document is the shortest reproducibility map for the **current v13 global-only H1–H4 submission**. Pre-v13 hypothesis ladders and defended submission surfaces are historical provenance and are archived from the active paper interface under `legacy/chapter1-pre-v13/`.

## Pipeline at a glance

```text
island geography + flora/status infrastructure
                    +
provenance-preserving plant trait evidence
                    |
                    v
frozen 8,265-island / 106,295-species database
                    |
                    +------------------------------+
                    |                              |
                    v                              v
six-atomic plant isolation responses      full-global GloPL distance analysis
                    |                              |
                    +---------------+--------------+
                                    |
                                    v
                         v13 H1–H4 synthesis
                                    |
                    +---------------+---------------+
                    |                               |
                    v                               v
      unified v13 result lock           post-hoc functional bridge lock
                    |                               |
                    +---------------+---------------+
                                    |
                                    v
                     manuscript + 3 main figures
                                    |
                                    v
                         v13 submission audit
```

## 1. Frozen database

Paper-level fixed quantities:

- **8,265** islands;
- **106,295** analysis-applicable plant species;
- **318,885** possible species-by-raw-axis cells;
- source run **34191508045**;
- resolved cells **222,688 / 318,885 = 69.83%**.

The raw database axes are flower colour, floral structural complexity and reproductive assurance. They are evidence-coverage axes; the v13 biological responses are derived frozen floral/reproductive states.

Database build details: `docs/DATABASE_BUILD.md`.

## 2. H1 — recurrent global island syndrome

The plant-side parent evidence uses six identically oriented responses:

1. generalized floral form;
2. actinomorphic symmetry;
3. shallow/open tube;
4. self-compatibility;
5. selfing mating system;
6. autonomous selfing.

All four predeclared geographic strata retain supported multivariate isolation responses in both evidence scopes. The arithmetic mean of the six classic-oriented slopes is positive in every stratum. The strata are used as replication of recurrence, not as a current between-stratum comparison.

Canonical synthesis: `config/chapter1_v13_unified_island_syndrome_result_lock.json`.

## 3. H2 — global pollination constraint

The independent GloPL layer contains **2,969 experiments, 1,248 sites and 919 publications**. The standardized global distance coefficient is `+0.07937` (one-sided positive `p=0.01772`).

This is direct evidence that pollen limitation increases with geographic separation in the sampled global experiment database. It does not identify which component of pollination service generated the limitation.

## 4. H3 — two plant response pathways

The six plant responses are summarized into two biological families:

- **reproductive assurance:** self-compatibility, selfing mating system, autonomous selfing;
- **floral accessibility/generalization:** generalized form, actinomorphy, shallow/open tube.

Both family means are positive across all four geographic replication strata in both evidence scopes. Frozen conditional analyses support treating them as partially separable components rather than one obligatory causal sequence.

## 5. H4 — post-hoc functional triangulation

Exact-species GloPL matches ask whether frozen trait states are associated with current pollen limitation after adjustment for distance, geography and measurement structure.

Canonical run:

- run **35141624253**;
- artifact **10465048981**;
- digest `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`.

Autonomous selfing is the strongest bridge (`beta=-0.44672`, two-sided `p=2.60e-8`) and remains negative in within-study checks. Actinomorphy and generalized floral form are directionally concordant but weaker. This layer is explicitly post-hoc and cannot establish historical selection or mediation.

Canonical lock: `config/chapter1_v13_functional_bridge_result_lock.json`.

## 6. Canonical paper outputs

Read in this order:

1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
3. `docs/chapter1_v13_submission_figure_sync_20260917.md`
4. `config/chapter1_v13_functional_bridge_result_lock.json`
5. `docs/chapter1_unified_hypothesis_20260917.md`
6. `docs/chapter1_submission_freeze_v13_20260917.md`

The figure contract contains **three main figures**. Evidence roles are kept explicit: frozen primary/confirmatory parent results, frozen sensitivity/robustness, descriptive synthesis and post-hoc functional triangulation are not interchangeable.

## 7. Active audit contract

Current validation:

- `tests/test_chapter1_v13_submission_lock.py`
- `tests/test_chapter1_v13_submission_surface.py`
- `tests/test_chapter1_v13_functional_bridge.py`
- `.github/workflows/audit-chapter1-v13-submission.yml`

The workflow runs on changes targeting `main`, executes the repository suite with only the two documented main-baseline CLI failures temporarily deselected, lints v13 code/tests and guards historical frozen provenance from accidental modification.

## 8. Claim ceiling

Allowed conclusion: the three evidence layers converge on a recurrent global floral/reproductive island syndrome that is aligned with increasing experimental pollen limitation and expressed through partially separable reproductive-assurance and floral-accessibility strategies.

Not identified by v13:

- historical causal mediation from pollen limitation to trait evolution;
- a global decline in pollinator abundance with isolation;
- loss of a named pollinator guild;
- equal response vectors across geographic strata;
- within-lineage evolutionary change as opposed to assemblage composition.

## 9. Chapter 1 / Chapter 2 handoff

Chapter 1 establishes recurrence, the pollen-limitation gradient, and functional compatibility. Chapter 2 (`izu-core`) can test the mechanistic sequence within a resolved biological system:

`interaction state -> effective service -> reproductive outcome -> phenotype`

This is where causal ordering among pollination service, pollen limitation, reproductive assurance and floral accessibility can be tested prospectively.

## 10. Legacy boundary

The prior publication-facing README, paper-pipeline narrative and retired pre-v13 submission-surface audit are preserved under `legacy/chapter1-pre-v13/`.

Historical analysis implementations and frozen parent result files remain in their original code/config locations when v13 reproducibility depends on them. Their presence is provenance, not an active competing paper claim.

`legacy/v1/` remains a separate frozen historical analysis.
