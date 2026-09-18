# Chapter 1 v14 paper pipeline

The canonical v14 surface changes the **question order and H1/H2 analysis surface**
while leaving frozen v13 GloPL and functional-bridge estimates intact.

    H1  island syndrome
        reproductive assurance + plain colour + accessible/generalized structure
          |
          v
    H2  floral decomposition
        selfing route vs selfing-adjusted pollinator-facing route
        + raw colour x raw form/tube concordance
          |
          v
    H3  independent pollen-limitation gradient (GloPL)
          |
          v
    H4  tiered functional bridge
        individual-trait discovery
        + literal H2 species-score discovery
        + atomic reconstruction sensitivity
        + prospective validation layers

H1 is now a seven-response beta-binomial multivariate test. H2 uses selfing_core as
the reproductive mechanism covariate and tests plain_colour and
generalized_accessible conditionally. Pollination-syndrome concordance is carried by
the frozen raw five-colour and raw colour-conditioned form/tube analyses; weighted
bee/butterfly/bird scores are historical secondary summaries only. H4 first bridges
the literal Direct-only H2 species scores to GloPL; atomic family reconstruction is a
separate post-hoc sensitivity.

See:
- docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md
- config/chapter1_v14_canonical_result_lock.json
- docs/chapter1_submission_freeze_v14_20260918.md
- config/chapter1_v14_hypothesis_architecture.yml
- config/chapter1_v14_all_data_probability.yml
- config/chapter1_v14_h2_decomposition.yml
- config/chapter1_v14_h4_evidence_hierarchy.yml
- config/chapter1_v14_h4_exact_h2_score_bridge.yml
- config/chapter1_v14_h4_exact_h2_score_bridge_preflight_result_lock.json
- config/chapter1_v14_h4_family_bridge.yml
- config/chapter1_v14_h4_family_bridge_preflight_result_lock.json
- docs/chapter1_v14_hypothesis_reorder_20260918.md
- docs/chapter1_v14_preflight_results_20260918.md

The v13 pipeline below remains immutable parent provenance.

---

## Frozen v13 parent paper pipeline

This parent section is the shortest reproducibility map for the frozen **v13 global-only H1–H4 submission**. Pre-v13 hypothesis ladders, mechanism branches and defended submission surfaces are historical provenance and are archived/indexed from the active paper interface under `legacy/chapter1-pre-v13/`.

## Pipeline at a glance

```text
island geography + flora/status infrastructure
                    +
provenance-preserving plant trait evidence
                    |
                    v
frozen 8,265-island / 106,295-species database
                    |
          +---------+----------------------+----------------------+
          |                                |                      |
          v                                v                      v
six-atomic plant isolation        raw colour + coupling     full-global GloPL
responses (frozen H1 core)        extension (H3 display)    distance analysis
          |                                |                      |
          +----------------+---------------+----------------------+
                           |
                           v
                upgraded v13 H1–H4 synthesis
                           |
          +----------------+----------------+
          |                                 |
          v                                 v
unified v13 result lock          post-hoc functional bridge lock
          |                                 |
          +----------------+----------------+
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

The raw database axes are flower colour, floral structural complexity and reproductive assurance. They are evidence-coverage axes; the v13 biological responses are derived frozen floral/reproductive states plus the reproduced raw-colour extension.

Database build details: `docs/DATABASE_BUILD.md`.

## 2. H1 — recurrent global functional island syndrome

The plant-side parent evidence uses six identically oriented responses:

1. generalized floral form;
2. actinomorphic symmetry;
3. shallow/open tube;
4. self-compatibility;
5. selfing mating system;
6. autonomous selfing.

All four predeclared geographic strata retain supported multivariate isolation responses in both evidence scopes. The arithmetic mean of the six classic-oriented slopes is positive in every stratum. The strata are used as replication of recurrence, not as a current between-stratum comparison.

These six responses remain the frozen confirmatory H1 core. The later raw-colour analyses do not become a seventh H1 atomic trait.

Canonical synthesis: `config/chapter1_v13_unified_island_syndrome_result_lock.json`.

## 3. H2 — global pollination constraint

The independent GloPL layer contains **2,969 experiments, 1,248 sites and 919 publications**. The standardized global distance coefficient is `+0.07937` (one-sided positive `p=0.01772`).

This is direct evidence that pollen limitation increases with geographic separation in the sampled global experiment database. It does not identify which component of pollination service generated the limitation.

## 4. H3 — three partially separable plant-response modules

### H3a Reproductive assurance

The frozen recurrent core includes self-compatibility, selfing mating system and autonomous selfing. The family mean is positive across all four geographic replication strata in both evidence scopes.

### H3b Floral accessibility/generalization

The frozen recurrent core also includes generalized form, actinomorphy and shallow/open tube. The family mean is positive across all four strata in both evidence scopes, although individual atomic traits do not move uniformly.

### H3c Pollinator-facing display reorganization

PR #233 added a reproduced raw-trait extension using five reported flower-colour states plus raw floral form and tube-depth architecture. The extension asks two stronger questions than a single syndrome score:

1. does raw colour composition change with isolation after `selfing_core` adjustment?;
2. among species retaining a focal colour, does the conditional probability of a prespecified architecture change with isolation?

Main reproduced outcomes:

- northern mid-latitude: robust `red_pink` decline in both evidence scopes;
- northern high latitude: replicated decline of `blue_purple` coupling to butterfly-associated form and intermediate/deep large-bee-associated tube architecture;
- tropical Direct evidence: `yellow_orange × deep butterfly-associated tube` increases with isolation;
- southern extratropical: mixed yellow/orange restructuring rather than a coherent named syndrome.

The biological synthesis is therefore **recurrent but non-uniform response**: reproductive assurance and accessibility/generalization form a recurrent global functional core, while pollinator-facing display composition and colour–architecture coupling can be reorganized in context-dependent ways.

Reproduced audits:

- `docs/chapter1_v13_raw_colour_audit_20260917.md`
- `docs/chapter1_v13_raw_colour_coupling_audit_20260917.md`

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
6. `docs/chapter1_v13_raw_colour_audit_20260917.md`
7. `docs/chapter1_v13_raw_colour_coupling_audit_20260917.md`
8. `docs/chapter1_submission_freeze_v13_20260917.md`

The figure contract still contains **three main figures**. Evidence roles remain explicit: frozen primary/confirmatory parent results, frozen sensitivity/robustness, reproduced extended H3 evidence, descriptive synthesis and post-hoc functional triangulation are not interchangeable.

## 7. Active audit contract

Current validation:

- `tests/test_chapter1_v13_submission_lock.py`
- `tests/test_chapter1_v13_submission_surface.py`
- `tests/test_chapter1_v13_functional_bridge.py`
- `tests/test_chapter1_v13_raw_colour_audit.py`
- `tests/test_chapter1_v13_colour_architecture_audit.py`
- `tests/test_chapter1_v13_raw_colour_coupling_audit.py`
- `tests/test_chapter1_v13_legacy_boundary.py`
- `tests/test_v13_active_surface_legacy_boundary.py`
- `.github/workflows/audit-chapter1-v13-submission.yml`
- `.github/workflows/run-chapter1-v13-raw-colour-audit.yml`

The workflows execute the repository suite with the documented exclusions where applicable, lint v13 code/tests, guard the active/legacy boundary and preserve the frozen parent result-lock provenance required by v13.

## 8. Claim ceiling

Allowed conclusion: the evidence layers converge on a recurrent global **functional** floral/reproductive island syndrome aligned with increasing experimental pollen limitation. The recurrent core is expressed through partially separable reproductive-assurance and floral-accessibility strategies, while raw pollinator-facing display contains additional context-dependent isolation responses not reducible to measured reproductive assurance.

Not identified by v13:

- historical causal mediation from pollen limitation to trait evolution;
- a global decline in pollinator abundance with isolation;
- loss of a named pollinator guild;
- flower colour as a unique identifier of realized pollinator identity;
- reduced attraction investment, pigment concentration, animal-visual contrast or UV signalling as directly measured quantities;
- equal detailed phenotype vectors across geographic strata;
- within-lineage evolutionary change as opposed to assemblage composition.

## 9. Chapter 1 / Chapter 2 handoff

Chapter 1 establishes two levels simultaneously:

1. **what recurs:** a functional island-syndrome core plus a global pollen-limitation gradient;
2. **what does not have to recur identically:** the detailed pollinator-facing display response.

Chapter 2 (`izu-core`) takes the next question: why can a broad ecological constraint recur while response trajectories and determinant rankings differ among systems? Its transportability framework treats receiving functional state, interaction regime and nonlinear response geometry as reasons that one ranking or response pathway need not transport unchanged across ecological contexts.

Programme-level bridge:

`Ch1: recurrent functional core + context-dependent display`

`-> Ch2: state/context/response geometry -> non-uniform response branches and determinant rankings`

This is a conceptual handoff, not a claim that Chapter 2 theory causally validates the Chapter 1 geographic associations.

## 10. Legacy boundary

The active package no longer exposes retired pre-v13 publication or mechanism implementations. Representative historical files are retained under `legacy/chapter1-pre-v13/`, while the complete retired implementation remains exactly recoverable from the pre-cleanup Git commit recorded in `legacy/chapter1-pre-v13/ARCHIVE_MANIFEST.md`.

The active `config/` directory deliberately retains the frozen parent result locks required by `config/chapter1_v13_unified_island_syndrome_result_lock.json`. They are result provenance, not competing active hypotheses.

`legacy/v1/` remains a separate frozen historical analysis.