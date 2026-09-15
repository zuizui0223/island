# Island — Chapter 1 paper repository

This repository contains the **current Chapter 1 / global island-flora analysis**. The paper-facing tree is organized around one traceable line:

```text
source data
  -> reproducible Chapter 1 database
  -> frozen H1-H5 hypothesis contract
  -> canonical progressive analysis
  -> prospective robustness / falsification extensions
  -> P0 provenance + P1 assembly defense + P2 component-support defense
  -> P3 joint observation-bias / partial-identification defense
  -> court-style evidence ledger
  -> locked Figures 1-4
  -> canonical island-first v11 manuscript
```

Development experiments and superseded analyses remain recoverable in git history. **`legacy/v1/` is a separately frozen historical analysis and is preserved intact.**

## 1. Canonical analysis database

The paper uses a fixed universe of **8,265 islands** and a fixed trait denominator of **106,295 analysis-applicable plant species x 3 raw measurement axes = 318,885 species-axis cells**.

The three strict raw axes are:

- flower colour;
- floral structural complexity;
- reproductive assurance.

These database axes are not three arbitrary fitted scalar responses. The primary Chapter 1 response is the pollinator-name-free pair `accessibility_generalization` + `reproductive_assurance`; broader atomic colour/structure/self-compatibility contrasts and named floral templates are secondary diagnostics.

Core data layers are:

- GSHHG-derived island geography / area / distance infrastructure;
- GBIF records as an opportunistic sample of realised island floras, never a complete census or absence layer;
- frozen floristic source/status infrastructure for native/non-endemic/endemic and source-pool comparisons;
- provenance-preserving species-level trait evidence.

Trait evidence precedence is **species-direct High/Medium first, then trait-specific Validated Low only where direct evidence is absent**. Family inference and global fallback are prohibited in the canonical Chapter 1 contract.

Database construction is documented in [`docs/DATABASE_BUILD.md`](docs/DATABASE_BUILD.md). The current reproducible trait-build line is the recovered checkpoint plus source-scale integration merged through PR #151. Source manifests are retained under `data/v2/staging/traits/source_batches/`, with promotion/reporting rules in [`docs/chapter1_trait_coverage_contract.md`](docs/chapter1_trait_coverage_contract.md).

Final trait snapshot used by the paper:

- source run: **34191508045**;
- artifact: `source-scale-batch-integration-34191508045`;
- resolved cells: **222,688 / 318,885 = 69.83%**;
- reproductive-assurance cells: **48,497 / 106,295 = 45.63%**.

The exact paper database is identified as **Database 1.0** through a versioned, SHA-locked manifest. See [`DATABASE_RELEASE.md`](DATABASE_RELEASE.md) for the Zenodo release route and the database-version entry point. Future Database 2.0+ snapshots can replace the database manifest while retaining the same scientific contract.

## 2. Frozen hypothesis ladder

Canonical primary config:

- [`config/chapter1_progressive_analysis.yml`](config/chapter1_progressive_analysis.yml)
- contract: `chapter1_progressive_analysis_v1`

Hypothesis order:

1. **H1 — Universal-syndrome rival:** does one floral/reproductive island syndrome hold across contexts?
2. **H2 — Biogeographic branching:** do response vectors differ among contexts?
3. **H3 — Source/lineage assembly:** how strongly is a supported response represented by taxonomic composition and does a robust beyond-genus response remain?
4. **H4 — Area/capacity moderation:** does continuous island area modify isolation-associated filtering strongly enough to promote a mechanism?
5. **H5 — Independent mechanism gate:** a named pollinator mechanism requires independent evidence and cannot be inferred from floral phenotype.

The claim ceiling is part of the contract. Floral architecture cannot be used to infer historical pollinator loss, replacement, exact source ancestry, or in-situ evolution.

## 3. Canonical primary analysis

Single paper-level workflow:

- [`.github/workflows/run-chapter1-progressive-trait-analysis.yml`](.github/workflows/run-chapter1-progressive-trait-analysis.yml)

Primary final execution:

- run: **34232450884**;
- artifact: `chapter1-progressive-analysis-34232450884`;
- artifact ID: **10058653212**;
- digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

The primary workflow fixes H1-H4, taxonomic depth, climate/source safeguards, trait-MNAR sensitivity, floral-architecture decomposition and the original H5 claim ceiling.

## 4. Prospective robustness and falsification extensions

The following analyses were added only under outcome-blind or prospectively frozen gates. They do **not** reopen failed mechanisms or relax prior thresholds.

### V6 — species-list / detection tipping analysis

- canonical run: **34800498716**;
- artifact: **10331282464**;
- baseline reproduces frozen H2 to numerical precision;
- Palearctic accessibility survives **99/100** baseline-supported bias surfaces and **80/80** scenarios with distance-dependent under-survey in the biologically concerning direction;
- tropical accessibility is more sensitive (**40/75** baseline-supported surfaces reach a tipping event);
- this is a robustness/tipping analysis, not an estimate of true GBIF completeness.

### Calibrated response geometry

Naive breakpoint selection was rejected after spatial-block simulations exposed high false nonlinear selection. Calibrated gates then showed that a step of the target magnitude would be identifiable in nearly all design cells. When observed outcomes were opened, **0/12** broad atomic response cells promoted a nonlinear shape. The result rules out a common global assemblage breakpoint under the tested geometry contract, not local or lineage-specific thresholds.

### H5c — independent biotic vs wind specificity

External GIFT `pollen_vector_mode` was used as a negative-control mechanism test and was never inferred from floral phenotype.

- qualification run: **34803493307**;
- independent unambiguous mode species: **5,771**;
- fully qualified cells: **1/8**;
- that single cell was opened once in run **34803837463**;
- direct-only native-nonendemic Palearctic `distance x biotic` interaction = **+0.06495**, 95% CI **[-0.09030, 0.22020]**, p=**0.41221**;
- verdict: **no independent pollination-mode specificity support**.

### H5d — distributed-threshold identifiability

Outcome-closed simulation asked whether lineage-specific thresholds can be distinguished from heterogeneous smooth clines in the realized assemblage design.

- run: **34803545574**;
- qualified cells: **0/8**;
- classification accuracy roughly **0.733–0.778**;
- false distributed-threshold selection under smooth clines: **0.19–0.255**;
- observed genus threshold distributions remain closed.

These results are not mechanism “failures” to be rescued. They define the present global claim boundary.

## 5. P0/P1/P2 defense of the island-syndrome inference

The NEE-oriented existing-data development plan starts from the island question rather than from a generic macroecological theory:

- [`docs/chapter1_nee_existing_data_design_island_first_20260915.md`](docs/chapter1_nee_existing_data_design_island_first_20260915.md)

### P0 — immutable claim reconciliation

- [`docs/chapter1_p0_claim_ledger_20260915.md`](docs/chapter1_p0_claim_ledger_20260915.md)
- headline result provenance was reconciled before new P1 modelling.

### P1a — exact paired support

- run **34935183075**;
- artifact **10382749051**;
- observed/family/genus stages use the same focal island–species observations and information weights;
- stage-specific sample loss is not a viable explanation for the frozen attenuation.

### P1c — true genus versus matched-complexity pseudo-genus

Within every family, species were randomly repartitioned while preserving the exact real genus count and genus group-size multiset.

- frozen permutation run: **34936193944**;
- **2,000/2,000** valid permutations;
- final aggregate run: **34941827774**;
- artifact: **10385820775**;
- true-genus conditional attenuation statistic: **0.72066**;
- matched pseudo-genus null median: **0.22434**;
- null ≥ observed: **57/2,000**;
- one-sided randomization **p=0.02899**;
- verdict: **true genus boundaries contain attenuation-relevant structure beyond arbitrary fine grouping of identical complexity**.

### P1d — paired spatial-block uncertainty

- run **34939113182**;
- artifact **10383754545**;
- 2,000 paired spatial-block bootstrap draws;
- total genus attenuation remains large;
- **0/8** direct-only primary profiles have a 95% interval for the *additional family→genus attenuation* entirely above zero.

Integrated decision:

> **The Palearctic floral-island response is strongly genus-structured, but the exact incremental family-to-genus attenuation is spatially imprecise.**

Canonical lock:

- [`config/chapter1_p1_final_decision_result_lock.json`](config/chapter1_p1_final_decision_result_lock.json)

### P2 — common-support defense of the H2 component contrast

The formal direct H2 comparison is `northern_midlatitude` versus `tropical` within the `analysis_regime` layer. `Palearctic` remains a separate within-context realm result and the focal H3/P1 genus-structure system.

- run **34945548775**;
- artifact **10387261614**;
- direct-only native-nonendemic common-island North–Tropical vector difference: **p=0.000947**;
- 2,000 paired spatial-block draws: determinant CI includes zero, so strong vector non-collinearity is **not established**;
- common-species sensitivity: **853** direct-only co-observed species, native-nonendemic joint vector difference **p=0.00517**;
- frozen Palearctic–Neotropical direct tests remain unsupported (`p=0.078–0.396`).

Canonical lock:

- [`config/chapter1_p2_component_nonconcordance_result_lock.json`](config/chapter1_p2_component_nonconcordance_result_lock.json)


### P3 — joint observation-bias and partial-identification defense

V5 trait-resolution MNAR and V6 species-list detection were first reproduced separately on the pinned PR142 input. Only then were the two processes crossed under a prospectively frozen joint contract.

- canonical P3 run: **34949880409**;
- artifact: **10389197309**;
- finite joint surfaces: **1,575 per evidence scope**;
- primary direct-only native-nonendemic North–Tropical vector: **1,541/1,575** robust cells;
- Palearctic accessibility, native non-endemics: **1,575/1,575** robust in both evidence scopes;
- tropical accessibility, native non-endemics: **1,161/1,575** all-analysis and **1,269/1,575** direct-only;
- deterministic partial-identification envelope preserves the positive Palearctic accessibility sign but does not identify formal North–Tropical support across every corner; tropical accessibility crosses zero.

Grid-cell fractions describe the declared sensitivity domain and are **not probabilities**. P3 estimates neither true flora completeness nor an arbitrary-MNAR latent truth.

Canonical lock:

- [`config/chapter1_p3_joint_observation_bias_result_lock.json`](config/chapter1_p3_joint_observation_bias_result_lock.json)

## 6. Current submission surface

Read these first, in this order:

1. [`docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md`](docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md) — current v11 submission state and claim ceiling;
2. [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md) — canonical island-first v11 manuscript;
3. [`docs/chapter1_v11_submission_figure_sync_20260915.md`](docs/chapter1_v11_submission_figure_sync_20260915.md) — final v11 panel mapping and figure-reference contract;
4. [`docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md`](docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md) — P3 finite-grid and partial-identification result;
5. [`docs/chapter1_p2_component_nonconcordance_result_20260915.md`](docs/chapter1_p2_component_nonconcordance_result_20260915.md) — P2 common-support and claim-boundary result;
6. [`docs/chapter1_p1_final_decision_20260915.md`](docs/chapter1_p1_final_decision_20260915.md) — integrated P1a/P1c/P1d decision;
7. [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`](docs/chapter1_court_evidence_and_theory_synthesis_20260914.md) — historical court-style evidence ledger;
8. [`docs/chapter1_literature_positioning_20260909.md`](docs/chapter1_literature_positioning_20260909.md) — frozen literature-positioning note.

The **single canonical working manuscript** is now:

- [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md)

v10 remains available as the previous P1/P2-defended surface. v11 additionally localizes joint observation robustness through P3 and must be used for current quoting.

### Locked main figures

- **Figure 1:** `config/chapter1_v8_figure1_result_lock.json` — hierarchy-of-assembly inference map;
- **Figure 2:** `config/chapter1_v10_figure2_p2_result_lock.json` — formal same-layer North–Tropical contrast plus P2 common-support defense;
- **Figure 3:** `config/chapter1_v9_figure3_p1_defense_result_lock.json` — genus-specific matched-null defense plus paired uncertainty;
- **Figure 4:** `config/chapter1_v11_figure4_result_lock.json` — joint observation-bias robustness, partial identification, and retained mechanism boundaries.

### Current paper-level verdict

- **H1:** one universal floral/reproductive island syndrome is not recovered;
- **H2/P2:** the formal same-layer North–Tropical joint response difference survives common-island and common-species restrictions; strong geometric non-collinearity is not precisely established;
- **H3/P1:** the Palearctic response is strongly structured by true genus composition; true genera outperform arbitrary within-family partitions of matched complexity (`p=0.02899`), but the exact additional family→genus attenuation is not precisely estimated across spatial blocks;
- **H4:** area remains a measurement-sensitive modifier rather than a promoted mechanism;
- **geometry:** no common nonlinear assemblage threshold is promoted;
- **P3:** the Palearctic accessibility branch is the observation-robust core; the formal North–Tropical vector is highly finite-domain robust but only partially identified, while tropical accessibility is observation-fragile;
- **H5:** a global pollinator-specific mechanism is not identified by N1, sampled source breadth, independent biotic-vs-wind specificity, or distributed-threshold tests.

Publication-facing concept:

> **A floral island syndrome can be a genus-structured assemblage syndrome: a visible community-level trait pattern can be carried by non-random lineage composition rather than one repeated organismal response.**

`assembly depth` remains useful as a localization concept, but v11 neither presents the family→genus increment as a perfectly sharp taxonomic breakpoint nor treats Palearctic and tropical as labels of one formal direct contrast.

`large_bee_like`, `butterfly_like`, and `bird_like` remain secondary **floral-architecture concordance scores**, not pollinator classifiers.

## 7. Chapter 1 / Chapter 2 division of labour

Chapter 1 answers:

- where isolation-associated responses occur;
- which floral/reproductive components move;
- whether directions differ among contexts;
- how strongly the strongest syndrome is represented in non-random genus composition;
- which simple global explanations survive explicit falsification.

Chapter 2 (`izu-core`) is the mechanistic-resolution layer. It can measure the within-system chain:

`interaction state -> effective service -> reproductive outcome -> phenotype`

and test cline versus threshold-like response geometry directly. The H5d non-identifiability result is the reason this local resolution is necessary; Chapter 1 should not retrofit a threshold mechanism from global assemblage averages.

## 8. Repository map

```text
README.md
  <- start here
docs/DATABASE_BUILD.md
  <- how the analysis database was built
DATABASE_RELEASE.md
  <- versioned DB / Zenodo entry point
docs/PAPER_PIPELINE.md
  <- database -> H1-H5 -> robustness/falsification -> P0/P1/P2/P3 -> v11
config/chapter1_progressive_analysis.yml
config/chapter1_database_versions/
config/chapter1_p1_final_decision_result_lock.json
config/chapter1_p2_component_nonconcordance_result_lock.json
config/chapter1_p3_joint_observation_bias_result_lock.json
config/chapter1_v8_figure1_result_lock.json
config/chapter1_v10_figure2_p2_result_lock.json
config/chapter1_v9_figure3_p1_defense_result_lock.json
config/chapter1_v11_figure4_result_lock.json
.github/workflows/run-chapter1-progressive-trait-analysis.yml
src/island_v2/
data/v2/
docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md
docs/chapter1_v11_submission_figure_sync_20260915.md
docs/chapter1_p2_component_nonconcordance_result_20260915.md
docs/chapter1_p1_final_decision_20260915.md
docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md
legacy/v1/
```

## 9. Legacy v1 preservation

`legacy/v1/` is retained as the complete historical v1 analysis. The preserved tree is:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

No current Chapter 1 claim should silently mix historical v1 outputs into the present v2 estimand.
