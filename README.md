# Island — Chapter 1 paper repository

This repository contains the **current Chapter 1 / global island-flora analysis**. The paper-facing tree is organized around one traceable line:

```text
source data
  -> reproducible Chapter 1 database
  -> frozen H1-H5 hypothesis contract
  -> canonical progressive analysis
  -> prospective robustness / falsification extensions
  -> court-style evidence ledger
  -> canonical v8 manuscript
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
3. **H3 — Source/lineage assembly:** at what taxonomic depth is a supported response represented?
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

## 5. Current evidence synthesis and canonical manuscript

Read these first:

- [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`](docs/chapter1_court_evidence_and_theory_synthesis_20260914.md) — supported claims, adverse evidence, defenses and claim ceilings;
- [`docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md`](docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md) — new Figure 1 inference map;
- [`docs/chapter1_literature_positioning_20260909.md`](docs/chapter1_literature_positioning_20260909.md) — frozen literature-positioning note.

The **single canonical working manuscript** is now:

- [`docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`](docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md)

v7 remains available in git history and as the previous submission-order draft; it is no longer the canonical narrative surface.

### Current paper-level verdict

- **H1:** one universal floral/reproductive island syndrome is not recovered;
- **H2:** source-separation responses branch among biogeographic contexts and trait components can decouple;
- **H3:** the strongest Palearctic response is concentrated at the family-to-genus assembly transition; genus adjustment attenuates roughly **78.8–85.9%** of the observed vector;
- **H4:** area remains a measurement-sensitive modifier rather than a promoted mechanism;
- **geometry:** no common nonlinear assemblage threshold is promoted;
- **V5/V6:** the Palearctic core survives strong trait-missingness and specified species-detection challenges, while tropical accessibility is less robust;
- **H5:** a global pollinator-specific mechanism is not identified by N1, sampled source breadth, independent biotic-vs-wind specificity, or distributed-threshold tests.

Publication-facing concept:

> **A floral island syndrome can be a hierarchically assembled phenotypic syndrome: a visible community-level trait pattern whose direction and assembly depth depend on biogeographic context.**

`large_bee_like`, `butterfly_like`, and `bird_like` remain secondary **floral-architecture concordance scores**, not pollinator classifiers.

## 6. Chapter 1 / Chapter 2 division of labour

Chapter 1 answers:

- where isolation-associated responses occur;
- which floral/reproductive components move;
- whether directions differ among contexts;
- at what taxonomic depth the strongest syndrome is expressed;
- which simple global explanations survive explicit falsification.

Chapter 2 (`izu-core`) is the mechanistic-resolution layer. It can measure the within-system chain:

`interaction state -> effective service -> reproductive outcome -> phenotype`

and test cline versus threshold-like response geometry directly. The H5d non-identifiability result is the reason this local resolution is necessary; Chapter 1 should not retrofit a threshold mechanism from global assemblage averages.

## 7. Repository map

```text
README.md
  <- start here
docs/DATABASE_BUILD.md
  <- how the analysis database was built
DATABASE_RELEASE.md
  <- versioned DB / Zenodo entry point
docs/PAPER_PIPELINE.md
  <- database -> H1-H5 -> robustness/falsification -> v8
config/chapter1_progressive_analysis.yml
config/chapter1_database_versions/
.github/workflows/run-chapter1-progressive-trait-analysis.yml
src/island_v2/
data/v2/
docs/chapter1_court_evidence_and_theory_synthesis_20260914.md
docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md
docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md
legacy/v1/
```

## 8. Legacy v1 preservation

`legacy/v1/` is retained as the complete historical v1 analysis. The preserved tree is:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

No current Chapter 1 claim should silently mix historical v1 outputs into the present v2 estimand.