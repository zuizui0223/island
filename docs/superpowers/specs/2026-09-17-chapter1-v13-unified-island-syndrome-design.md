# Chapter 1 v13 unified island-syndrome paper — design

Date: 2026-09-17
Status: approved design; implementation not yet started
Base: `main` at the post-PR #229 Chapter 1 integration state
Target branch: `ch1-v13-unified-island-syndrome`

## 1. Goal

Replace the current submission-facing biological story with a new versioned v13 surface while preserving v11/v12 as immutable provenance.

The v13 paper asks whether geographic isolation is associated with a recurrent global floral/reproductive island syndrome and whether independent experimental pollen-limitation evidence provides a functional bridge to two partially independent response pathways:

1. reproductive assurance;
2. generalized / pollinator-accessible floral architecture.

Taxonomic assembly is retained as the level at which the common syndrome can be realized, not as the definition of whether the global syndrome exists.

## 2. Core synthesis to freeze

The paper-level hypothesis is:

> Island isolation is associated globally with stronger pollen limitation and with a recurrent floral/reproductive island syndrome. The syndrome contains two partially independent response pathways: greater reproductive assurance and greater floral accessibility/generalization. Experimental pollen-limitation data provide a functional bridge, strongest for autonomous selfing, while taxonomic assembly determines how the common syndrome is realized across flora layers and regions.

This is a triangulation claim, not a causal mediation claim.

## 3. Hypothesis architecture

### H1 — Global island syndrome

Question: Is there a positive global common component in the isolation response toward the predeclared classic island-syndrome direction?

Primary evidence:
- all-observed contemporary flora;
- all-analysis-eligible trait scope primary;
- direct-only sensitivity;
- four analysis regimes retained simultaneously.

The target direction combines increased reproductive assurance and increased generalized/accessibility architecture. Regional coefficients may differ in magnitude; exact equality is not required.

Regional heterogeneity is secondary modulation, not the headline hypothesis.

### H2 — Global pollination-service constraint

Question: Does experimental pollen limitation increase with geographic separation from major continental landmasses?

Frozen parent evidence:
- GloPL full-global distance model;
- 2,969 experiments;
- 1,248 sites;
- 919 publications;
- global distance estimate approximately +0.07937;
- one-sided positive p approximately 0.01772;
- post-hoc offshore-only shape evidence remains descriptive rather than confirmatory.

Interpretation ceiling: this is an isolation-associated reproductive-service constraint, not proof of declining pollinator abundance or visitation.

### H3 — Dual response pathways

Question: Does the global plant-side syndrome contain two partially independent positive isolation responses?

H3a reproductive-assurance branch:
- self compatibility;
- selfing mating system;
- autonomous selfing;
- strict `selfing_core` where available.

H3b generalized-accessibility branch:
- generalized floral form;
- actinomorphic symmetry;
- shallow/open tube where supported;
- predeclared generalized-accessible composite where available.

The global common effect and regional deviations must be reported separately.

### H4 — Functional bridge

Question: Are the predeclared response traits associated with lower experimentally measured pollen limitation under the frozen exact-species GloPL overlap?

This is a post-hoc functional triangulation layer introduced after the failed frozen Route A/B interaction tests. It must never be relabelled confirmatory.

Primary post-hoc functional contrasts to reproduce and freeze without recoding:
- reproductive assurance: self compatibility, autonomous selfing; selfing mating system remains non-evaluable if support remains inadequate;
- floral architecture: generalized form, actinomorphic symmetry; shallow/open tube remains non-evaluable if support remains inadequate.

Primary model form:

`PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects`

with:
- exact species matching only;
- the same frozen trait recodes as existing Route A/B;
- publication total weight = 1;
- publication-cluster robust covariance;
- the same GloPL measurement fixed effects;
- no support-threshold relaxation;
- no new synonym or genus fallback.

Required sensitivities:
1. supplemental-only;
2. no-zero-constant;
3. within-publication contrast where both states occur;
4. within-publication × site contrast where both states occur and estimability is adequate.

Observed preliminary targets that must be independently reproduced before locking include:
- autonomous selfing: strong negative PL association, approximately beta -0.44 globally;
- autonomous selfing within publication: approximately beta -0.32;
- autonomous selfing within publication × site: approximately beta -0.26;
- generalized form: negative PL association, approximately beta -0.19;
- actinomorphic symmetry: negative PL association, approximately beta -0.38;
- self compatibility: negative but currently imprecise.

Any discrepancy must be reported and the reproduced values, not these preliminary values, become canonical.

Interpretation:
- autonomous selfing can be treated as direct functional evidence that reproductive assurance reduces current pollen limitation if the frozen reanalysis and within-study/site sensitivities reproduce;
- generalized/actinomorphic architecture can be treated as functional compatibility evidence only if global and sensitivity directions reproduce;
- this layer does not establish that pollen limitation historically selected the trait.

### H5 — Taxonomic realization

Question: At what taxonomic depth is the common syndrome represented?

Retain the two existing layers:

1. Broad all-observed H3A:
   - source-free taxonomic residualization;
   - the broad response remains after genus residualization;
   - conclusion: broad contemporary context dependence is not erased by raw genus composition.

2. Defended native Palearctic H3B:
   - source-matched family/genus decomposition;
   - 4/4 -> 4/4 -> 0/4 support trajectory;
   - matched pseudo-genus null p approximately 0.02899;
   - conclusion: defended native Palearctic response is strongly genus-structured.

These are different flora layers and must not be forced into one causal mediation chain.

## 4. Evidence hierarchy

### Primary paper spine

1. Global plant syndrome.
2. Global experimental pollen-limitation gradient.
3. Dual plant pathways.
4. Post-hoc functional bridge.
5. Taxonomic realization and native-status boundary.

### Secondary modifiers

- biogeographic differences in the magnitude/composition of the common syndrome;
- area moderation;
- observation-bias robustness;
- response geometry.

### Supplementary / non-promoted evidence

- GloBI source-genus channel breadth;
- exact-island individual pollinator channels;
- identity-aware two-route tests;
- named large-bee/butterfly/bird template identity claims.

GloBI is retained only as sampling-effort-sensitive interaction-structure context and must not carry the mechanism argument.

## 5. Versioning and provenance

Do not edit or delete v11/v12 result locks or manuscripts.

Create a new v13 submission surface consisting of:

1. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
2. `config/chapter1_v13_functional_bridge_v1.yml`
3. `config/chapter1_v13_functional_bridge_result_lock.json`
4. `src/island_v2/chapter1_v13_functional_bridge.py`
5. tests for model construction, weighting, exact matching, support gates, and within-study/site contrasts
6. a GitHub Actions workflow that reproduces the v13 functional bridge from pinned sources
7. `docs/chapter1_unified_hypothesis_20260917.md`
8. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
9. `docs/chapter1_v13_submission_figure_sync_20260917.md`
10. submission-surface updates in README and `docs/PAPER_PIPELINE.md`

The v13 lock must point back to all reused canonical v11/v12 locks and their run/artifact/digest provenance rather than copying undocumented numbers.

## 6. Manuscript structure

### Working title

`A global floral island syndrome emerges under increasing pollination limitation through reproductive assurance and floral generalization`

Title is working copy and may be shortened after the complete draft is assembled.

### Abstract logic

1. Isolation can constrain pollination service, but global floral island syndromes are usually inferred without independent service evidence.
2. Across thousands of island floras, a common classic-syndrome direction recurs across all four broad regions, with regional modulation in magnitude and composition.
3. Independent GloPL experiments show increasing pollen limitation with geographic separation.
4. Plant response separates into reproductive-assurance and generalized-accessibility branches rather than one obligatory selfing-to-simplification chain.
5. Exact-species GloPL triangulation tests whether those traits are functionally associated with current PL; autonomous selfing is expected to be the strongest bridge if reproduction succeeds.
6. Taxonomic decomposition shows that the common broad response and the defended native Palearctic response are realized differently across flora layers.
7. Conclusion: a recurrent global island syndrome is compatible with a shared pollination-service constraint expressed through dual pathways and taxonomic assembly, without claiming historical causal mediation.

### Results order

1. Global common syndrome component.
2. Four-region modulation around the common component.
3. Global pollen-limitation gradient.
4. Reproductive-assurance branch.
5. Generalized-accessibility branch and independence from selfing core.
6. Post-hoc functional bridge to current PL.
7. Taxonomic realization: broad H3A versus native Palearctic H3B.
8. Robustness boundaries: native-status transfer, observation bias, area, geometry, and non-promoted GloBI/channel evidence.

## 7. Figure architecture

### Figure 1 — Unified hypothesis and evidence ladder

Isolation -> pollination-service constraint -> two partially independent plant pathways -> assemblage syndrome, with taxonomic realization shown as an orthogonal layer rather than a downstream causal mediator.

### Figure 2 — Global island syndrome

Four-region slope vectors plus a pooled/common component for reproductive assurance and generalized accessibility. Regional deviations are shown but not framed as the main discovery.

### Figure 3 — Experimental pollen limitation and functional bridge

A. global GloPL distance slope;
B. offshore-only descriptive shape;
C. trait-state PL differences for reproductive assurance;
D. trait-state PL differences for floral architecture, with within-study/site sensitivity markers.

### Figure 4 — Taxonomic realization and claim boundary

A. broad all-observed response after source-free family/genus residualization;
B. defended native Palearctic family-to-genus attenuation;
C. native-status transfer failure / WCVP compatibility boundary;
D. causal-claim ceiling.

## 8. Claim ceiling

Allowed if reproduced:
- geographic isolation is associated with a recurrent global floral/reproductive island-syndrome direction;
- experimental pollen limitation increases with geographic isolation globally;
- reproductive assurance and generalized accessibility form partially independent plant-side response pathways;
- autonomous selfing is associated with lower current experimental pollen limitation under exact-species matched GloPL data, including within-study/site sensitivity if reproduced;
- generalized / actinomorphic architecture is functionally compatible with lower current PL if the frozen post-hoc reanalysis reproduces;
- broad contemporary and defended native responses differ in taxonomic realization.

Prohibited:
- pollen limitation causally evolved the observed traits;
- pollinator abundance or visitation globally declines with isolation;
- GloBI proves pollinator-community restructuring;
- autonomous selfing mediates the geographic isolation effect;
- generalized architecture necessarily evolved because specialist pollinators disappeared;
- Palearctic genus structure explains the global all-observed response;
- all observed flora equals native assembly;
- post-hoc functional triangulation is confirmatory;
- lack of Route A/B interaction support means pollination is unimportant.

## 9. Verification requirements before v13 promotion

1. Reproduce every reused canonical value from its existing result lock or pinned artifact.
2. Add tests before implementing new v13 functional-bridge model code.
3. Run the new workflow from pinned GloPL and trait inputs.
4. Verify artifact digest and write a result lock only from a successful run.
5. Re-run the exact v13 paper-level lock validator after manuscript generation.
6. Ensure v11/v12 files are unchanged byte-for-byte.
7. Ensure README and `PAPER_PIPELINE.md` point to v13 only after all v13 validation passes.
8. Open a PR against `main`; do not merge automatically without a separate integration decision.

## 10. Acceptance criterion

v13 is ready for review when a fresh clone can recover, from committed locks and pinned artifacts, the complete evidence chain:

`isolation -> global syndrome`

`isolation -> experimental pollen limitation`

`global syndrome -> reproductive assurance + generalized accessibility`

`functional traits <-> current experimental pollen limitation`

`syndrome -> flora-layer-specific taxonomic realization`

with every arrow labelled as confirmatory, post-hoc triangulation, descriptive, or unresolved and with no causal claim exceeding the evidence class.
