# Chapter 1 paper pipeline

This is the shortest reproducibility map from the database to the current **v13 unified global island-syndrome paper**, with v11/v12 retained as historical defended provenance.

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
[9] P0 immutable claim reconciliation
                    |
                    v
[10] P1 assembly-inference defense
     same support -> matched genus null -> paired spatial uncertainty
                    |
                    v
[11] P2 component-support defense
     same islands -> paired blocks -> same species denominator
                    |
                    v
[12] P3 joint observation-bias defense
     V5 reproduction -> V6 reproduction -> joint surface -> partial identification
                    |
                    v
[13] historical island-first v11 defended manuscript
                    |
                    v
[14] v13 global syndrome + GloPL + functional bridge synthesis
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

Historical H3 decomposes the supported Palearctic response through observed, family-adjusted and source-matched genus-adjusted stages. The predeclared support ladder is `4/4 -> 4/4 -> 0/4`; descriptive total genus attenuation is roughly 78.8–85.9% of the observed vector.

Historical H3 alone is **not** the final v9 assembly claim. P1 later tests whether this genus sensitivity could be caused by sample loss or arbitrary fine grouping and quantifies spatial uncertainty in the family→genus increment.

Implementation includes:

- `src/island_v2/chapter1_taxonomic_depth_decomposition.py`
- `src/island_v2/chapter1_pr138_lineage_representation_bridge.py`

### H4 — area / capacity moderation

Test continuous distance x continuous area as a modifier and audit whether heteroskedastic measurement can generate apparent moderation. No small/large island cutoff is introduced after outcome inspection.

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

### V1 — climate/common-support validation

Tests whether regional contrasts can be transported over common measured-climate support. Failure to establish climate-independent categorical causation narrows interpretation without erasing H2.

### V2 — taxonomic-depth decomposition

Historical localization layer only. Its descriptive family→genus contrast is now interpreted through the P1 safeguards below rather than used by itself as a precise taxonomic breakpoint.

### V3 — area-support falsification

Prevents suggestive distance x area patterns from being promoted as founder/capacity mechanisms without passing heteroskedastic and support safeguards.

### V4 — architecture decomposition

Shows that named large-bee-like, butterfly-like and bird-like templates are dominated by a shared plant-architecture factor (~86.9% / 86.4%), so they are not independent visitor identities.

### V5 — trait-resolution MNAR tipping

Tests state-dependent trait missingness. The broad Palearctic primary response survives the finite predeclared MNAR grid, while some reproductive details and context contrasts remain more bounded.

### V6 — species-detection / list-completeness tipping

Canonical V6:

- Run `34800498716`;
- artifact `10331282464`;
- baseline OR_D=1 reproduces frozen H2;
- Palearctic accessibility survives 99/100 baseline-supported surfaces and all 80/80 remote-under-survey scenarios in the biologically concerning direction;
- tropical accessibility is more sensitive (40/75 baseline-supported surfaces tip);
- North–Tropical vector contrast is more robust than the tropical single axis.

V6 is a tipping analysis, not an occupancy estimate and not a claim that true GBIF completeness is known.

## 6. Calibrated response geometry

Geometry candidates were flat, cline, step, hinge and reversal. A naive AICc selector failed qualification because clustered true clines were often falsely labelled nonlinear. The calibrated geometry layer estimated cell-specific critical nonlinear evidence thresholds and validated them on independent seeds.

Observed opening: **12/12 `monotonic_or_unresolved`; zero promoted nonlinear shapes**.

Interpretation: no common global assemblage breakpoint is supported. Local-system or lineage-specific thresholds remain possible.

## 7. Prospectively gated H5 extensions

### N1 — independent channel heterogeneity

- joint isolation x channel Wald `W=1.6187`, df=3, p=0.65516;
- not promoted;
- limited power for modest true heterogeneity, so non-identification rather than channel equality.

### Source-side GloBI breadth

Effort-matched sampled source interaction breadth promoted 0/4 context x stratum cells.

### H5c — independent biotic vs wind specificity

External GIFT `pollen_vector_mode` provides a negative control independent of response-trait architecture.

- qualification Run `34803493307`;
- fully qualified 1/8 cell;
- observed Run `34803837463`;
- `distance x biotic` = +0.06495;
- 95% CI [-0.09030, 0.22020];
- p=0.41221;
- `no_pollination_mode_specificity_support`.

### H5d — distributed lineage-threshold identifiability

- Run `34803545574`;
- qualified 0/8 cells;
- classification accuracy about 0.733–0.778;
- false distributed-threshold selection under smooth clines 0.19–0.255;
- observed genus-level threshold distributions remain closed.

## 8. P0 — immutable claim reconciliation

Design and ledger:

- `docs/chapter1_nee_existing_data_design_island_first_20260915.md`
- `docs/chapter1_p0_claim_ledger_20260915.md`

P0 verifies exact run/artifact/table/estimand provenance before any P1 extension. Missing evidence is treated as unavailable, not negative.

## 9. P1 — defend the assembly inference

Canonical integrated decision:

- `config/chapter1_p1_final_decision_result_lock.json`
- `docs/chapter1_p1_final_decision_20260915.md`

### P1a — exact paired support

- Run `34935183075`;
- artifact `10382749051`;
- same focal island–species observations and same `n_species` weights across observed/family/genus stages;
- stage-specific sample loss is not a viable attenuation explanation.

### P1c — matched-complexity pseudo-genus null

Species are randomized **within family** while preserving the exact real genus-count and genus group-size multiset.

- source permutation Run `34936193944`;
- 40 shards x 50 = **2,000/2,000 valid permutations**;
- aggregate-only Run `34941827774` reuses those exact permutations;
- artifact `10385820775`, digest `861d18fe...`;
- true-genus median conditional attenuation `0.7206615`;
- null median `0.2243375`;
- 57/2,000 null permutations ≥ observed;
- one-sided randomization **p=0.0289855**;
- verdict: `true_genus_exceeds_matched_complexity_null`.

This rejects the simple “any equally fine grouping absorbs the signal” explanation.

### P1d — paired spatial-block uncertainty

- Run `34939113182`;
- artifact `10383754545`;
- 2,000 paired block bootstrap draws;
- total genus attenuation remains large;
- **0/8** direct-only primary profiles have the 95% interval for the additional family→genus attenuation wholly above zero.

Integrated inference:

> **The Palearctic floral-island response is strongly genus-structured beyond matched grouping complexity, while the exact incremental family-to-genus attenuation remains spatially imprecise.**

P1 supports taxonomic localization, not a precisely estimated family→genus breakpoint and not a causal assembly mechanism.


## 10. P2 — defend the component/context contrast

Canonical result:

- `config/chapter1_p2_component_nonconcordance_result_lock.json`
- `docs/chapter1_p2_component_nonconcordance_result_20260915.md`

Formal direct H2 contrast:

- `northern_midlatitude` versus `tropical` within `analysis_regime`;
- common-island direct-only NNE vector difference `p=0.000947`;
- 348 common islands, 78 spatial blocks;
- paired-block determinant interval includes zero, so strong non-collinearity is not established;
- common-species sensitivity uses 853 direct-only co-observed species and retains the NNE joint vector difference (`p=0.00517`);
- frozen Palearctic–Neotropical direct tests remain unsupported (`p=0.078–0.396`).

P2 therefore strengthens the same-layer joint branching claim while prohibiting the cross-layer shorthand “Palearctic versus tropical” as a formal direct H2 test.


## 11. P3 — jointly bound observation bias

Canonical result:

- `config/chapter1_p3_joint_observation_bias_result_lock.json`
- `docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md`
- Run `34949880409` / artifact `10389197309`.

Execution order was fail-closed: V5 reproduced first, V6 reproduced second, both were reconciled to the frozen PR142 baseline, and only then was the joint surface opened. The finite domain contains 1,575 V5×V6 assumption surfaces per evidence scope and never increases regression precision for hypothetical species.

Primary direct-only native-nonendemic results:

- North–Tropical vector difference: `1541/1575` robust cells;
- Palearctic accessibility: `1575/1575` robust cells;
- tropical accessibility: `1269/1575` robust cells.

The deterministic 48-corner partial-identification envelope preserves the positive Palearctic accessibility sign, but support for the formal North–Tropical vector is not identified across every corner and tropical accessibility crosses zero. The resulting labels are **observation-robust core**, **finite-domain robust / partially identified**, and **observation-fragile**, respectively.

Grid fractions are assumption-domain coverage, not probabilities. P3 does not estimate true species-list completeness or arbitrary-MNAR latent truth.

## 12. Locked figures

- Figure 1: `config/chapter1_v8_figure1_result_lock.json`
- **Figure 2: `config/chapter1_v10_figure2_p2_result_lock.json`**
- **Figure 3: `config/chapter1_v9_figure3_p1_defense_result_lock.json`**
- **Figure 4: `config/chapter1_v11_figure4_result_lock.json`**

New Figure 3:

- Run `34942779753`;
- artifact `10386555245`;
- digest `sha256:93957ed4e12193fd2630aa3c274a99f004af4c8e40a6f63985e46a77c4cdc239`;
- Panel A: frozen attenuation trajectories;
- Panel B: true genus versus 2,000 matched pseudo-genus permutations;
- Panel C: paired spatial-block CIs showing the imprecise increment;
- Panel D: final claim boundary.

This figure deliberately displays favourable and adverse evidence together.

## 13. Canonical paper surface

The current publication-facing surface is **v13 global-only**, with an **H1–H4** scientific spine. Historical v11/v12 surfaces remain immutable provenance.

Read in this order:

1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
3. `docs/chapter1_v13_submission_figure_sync_20260917.md`
4. `config/chapter1_v13_functional_bridge_result_lock.json`
5. `docs/chapter1_unified_hypothesis_20260917.md`
6. `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md` — historical defended manuscript

The current publication-facing hierarchy is:

- H1: recurrent global classic-island direction across four geographic replication strata;
- H2: full-global experimental pollen-limitation gradient;
- H3: partially separable reproductive-assurance and floral-accessibility pathways;
- H4: explicitly post-hoc exact-species functional triangulation.

The four geographic strata are used to demonstrate recurrence of the global direction. Between-stratum differences are not part of the current submission spine.

Canonical v13 functional-bridge run:

- run `35141624253`;
- artifact `10465048981`;
- digest `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`.

The frozen Route A/B distance-by-trait moderation failures remain negative results. v13 does not relabel them as support. GloBI remains supplementary sampling-sensitive evidence.

Previous v8/v9/v10/v11/v12 surfaces remain historical provenance and continue to document the chronology of claim defense.

## 14. Chapter 1 / Chapter 2 handoff

Chapter 1 identifies response direction, component decoupling, genus structuring and failure of common global mechanisms. It cannot identify the local causal chain.

Chapter 2 / `izu-core` should resolve:

`interaction state -> effective service -> reproductive outcome -> phenotype`

within one biological system and test cline, step, shared-breakpoint and channel-specific response geometries prospectively.

- **Chapter 1:** where, which components, and at what lineage-assembly level is the syndrome represented?
- **Chapter 2:** how and why functionally within a resolved system?

## 15. Legacy v1 boundary

`legacy/v1/` is a frozen historical analysis, not a validation stage of the current pipeline. Preserved tree:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

No current Chapter 1 claim should silently mix historical v1 outputs into the v2 estimand.
