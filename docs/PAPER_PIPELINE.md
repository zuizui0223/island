# Chapter 1 paper pipeline

This is the shortest reproducibility map from the database to the current v8 paper.

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
[6] V1-V6 validation / observation-bias layers
                    |
                    v
[7] calibrated response-geometry audit
                    |
                    v
[8] prospectively gated H5 mechanism tests
                    |
                    v
[9] court-style evidence ledger
                    |
                    v
[10] canonical v8 manuscript
```

## 1. Database inputs

See `docs/DATABASE_BUILD.md`.

Paper-level fixed quantities:

- 8,265-island universe;
- 106,295 analysis-applicable species;
- 318,885 species-axis denominator;
- final snapshot Run `34191508045` / `source-scale-batch-integration-34191508045`;
- 222,688 resolved cells (69.83%).

The raw database axes are flower colour, floral structural complexity and reproductive assurance. They are coverage axes, not three scalar H2 responses.

## 2. Freeze the trait snapshot

The progressive workflow validates the incoming species-axis ledger before fitting the paper analysis. Snapshot rules include fixed denominator, required columns, evidence quality, provenance retention, and explicit handling of audited evidence revision.

Implementation:

- `src/island_v2/chapter1_trait_snapshot.py`
- `src/island_v2/chapter1_trait_resolution_mnar.py`

## 3. Apply the frozen H1-H5 contract

Canonical config:

- `config/chapter1_progressive_analysis.yml`

### H1 — universal-syndrome rival

Test whether the same coherent floral/reproductive response vector appears across contexts. Direct between-context heterogeneity is required to reject the universal form.

### H2 — biogeographic branching

Test within-context multivariate responses and direct between-context vector differences. The primary response is pollinator-name-free accessibility/generalization + reproductive assurance. Named pollinator-like floral templates are secondary architecture summaries only.

### H3 — source / lineage assembly

Decompose the supported response through taxonomic depth and source-matched lineage representation. The final Palearctic result is retained at the observed stage and after family adjustment but fails the predeclared vector gate after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`). Genus adjustment attenuates roughly 78.8–85.9% of the observed vector, with most attenuation concentrated at family -> genus.

Implementation includes:

- `src/island_v2/chapter1_taxonomic_depth_decomposition.py`
- `src/island_v2/chapter1_pr138_lineage_representation_bridge.py`

### H4 — area / capacity moderation

Test continuous distance x continuous area as a modifier and audit whether heteroskedastic measurement can generate apparent moderation. No small/large island cutoff is introduced after outcome inspection.

Implementation:

- `src/island_v2/chapter1_area_capacity_moderation.py`
- `src/island_v2/chapter1_area_support_artifact.py`

### H5 — independent mechanism gate

A pollinator mechanism may be promoted only with independent pollinator-side information. Floral architecture cannot itself identify historical channel loss or replacement.

The primary Chapter 1 plant pattern can therefore be strong while H5 remains unidentified.

## 4. Primary execution

Single paper-level workflow:

- `.github/workflows/run-chapter1-progressive-trait-analysis.yml`

Final primary execution:

- Run `34232450884`;
- artifact `chapter1-progressive-analysis-34232450884`;
- artifact ID `10058653212`;
- SHA-256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

## 5. Validation and observation-bias layers

The explanation-gap sequence is no longer V1-V5 only. It now includes V6 as a separate species-list problem.

### V1 — climate/common-support validation

Tests whether regional contrasts can be transported over common measured-climate support. Failure to establish climate-independent categorical causation narrows interpretation without erasing H2.

### V2 — taxonomic-depth decomposition

Localizes the strongest Palearctic response to the family-to-genus assembly transition.

### V3 — area-support falsification

Prevents suggestive distance x area patterns from being promoted as founder/capacity mechanisms without passing heteroskedastic and support safeguards.

### V4 — architecture decomposition

Shows that named large-bee-like, butterfly-like and bird-like templates are dominated by a shared plant-architecture factor (~86.9% / 86.4%), so they are not independent visitor identities.

### V5 — trait-resolution MNAR tipping

Tests state-dependent trait missingness. The broad Palearctic primary response survives the finite predeclared MNAR grid, while some reproductive details and context contrasts remain more bounded.

### V6 — species-detection / list-completeness tipping

Tests a distinct problem: species present on islands but missing from compiled lists.

Canonical V6:

- Run `34800498716`;
- artifact `10331282464`;
- baseline OR_D=1 reproduces the frozen H2 outputs;
- Palearctic accessibility survives 99/100 baseline-supported surfaces and all 80/80 remote-under-survey scenarios in the biologically concerning direction;
- tropical accessibility is more sensitive (40/75 baseline-supported surfaces tip);
- the North–Tropical vector contrast is more robust than the tropical single-axis result.

V6 is a tipping analysis, not an occupancy estimate and not a claim that true GBIF completeness is known.

## 6. Calibrated response geometry

Geometry candidates were fixed as flat, cline, step, hinge and reversal. A naive AICc selector failed qualification because clustered true clines were often falsely labelled nonlinear. The calibrated V2 geometry layer therefore estimated a cell-specific critical nonlinear evidence threshold under monotonic truths and validated it on an independent seed.

Qualification showed that true steps of the target magnitude are detectable in almost all design cells, while hinges are not reliably identifiable. The observed opening then produced 12/12 `monotonic_or_unresolved` broad response cells: zero promoted step, hinge or reversal shapes.

Interpretation: no common global assemblage breakpoint is supported under the frozen contract. This does not rule out local-system or lineage-specific thresholds.

## 7. Prospectively gated H5 extensions

These analyses were added after the primary claim architecture had been frozen. They are allowed to strengthen or constrain H5, not rescue a failed mechanism by moving thresholds.

### N1 — independent channel heterogeneity

- joint isolation x channel Wald `W=1.6187`, df=3, p=0.65516;
- not promoted;
- simulation indicates limited power for modest true heterogeneity, so the result is non-identification rather than evidence of channel equality.

### Source-side GloBI breadth

Effort-matched sampled source interaction breadth promoted 0/4 context x stratum cells. Sampled partner breadth is therefore not elevated as the missing H3 mechanism.

### H5c — independent biotic vs wind specificity

External GIFT `pollen_vector_mode` provides a negative control independent of the response trait architecture.

Qualification:

- Run `34803493307`;
- 5,771 unambiguous mode species;
- support-qualified 4/8 cells;
- fully power-qualified 1/8 cell;
- only direct-only x native-nonendemic x Palearctic was allowed to open.

Observed:

- Run `34803837463`;
- `distance x biotic` = +0.06495;
- 95% CI [-0.09030, 0.22020];
- p=0.41221;
- classification: `no_pollination_mode_specificity_support`.

This does not prove pollinators do not matter. It prevents pollination mode from being claimed as an independently supported generator of the Palearctic accessibility gradient.

### H5d — distributed lineage-threshold identifiability

Outcome-closed simulations compare a distributed-threshold generator with heterogeneous smooth-cline generators on the realized assemblage design.

- Run `34803545574`;
- qualified 0/8 cells;
- classification accuracy about 0.733–0.778;
- false distributed-threshold selection under smooth clines 0.19–0.255;
- observed genus-level threshold distributions remain closed.

This is the formal reason the Chapter 1 paper should not infer lineage thresholds from smooth global gradients.

## 8. Evidence ledger and claim hierarchy

The manuscript-facing evidence table is:

- `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`

It deliberately places favourable and adverse evidence together. The strongest claim surviving cross-examination is:

> geographic isolation is associated with reproducible but biogeographically contingent floral/reproductive assemblage change, and the strongest Palearctic syndrome is principally expressed at the family-to-genus assembly transition.

The global upstream pollinator mechanism remains unidentified.

## 9. Canonical paper surface

Read in this order:

1. `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`
2. `docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md`
3. `docs/chapter1_literature_positioning_20260909.md`
4. `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`

Previous v7 and 2026-09-09 framing documents remain historical development surfaces and git provenance; they no longer define the canonical narrative.

## 10. Chapter 1 / Chapter 2 handoff

Chapter 1 identifies response direction, component decoupling, assembly depth and failure of common global mechanisms. It cannot identify the local causal chain.

Chapter 2 / `izu-core` should resolve:

`interaction state -> effective service -> reproductive outcome -> phenotype`

within one biological system and test cline, step, shared-breakpoint and channel-specific response geometries prospectively.

The two chapters therefore answer different levels of the same question:

- **Chapter 1:** where, which components, and at what assembly depth?
- **Chapter 2:** how and why functionally within a resolved system?

## 11. Legacy v1 boundary

`legacy/v1/` is a frozen historical analysis, not a validation stage of the current pipeline. Preserved tree:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

No current Chapter 1 claim should silently mix historical v1 outputs into the v2 estimand.