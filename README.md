# Island — Chapter 1 paper repository

This repository contains the **current Chapter 1 / global island-flora analysis**. The paper-facing working tree is intentionally organized around one traceable line:

```text
source data
  -> reproducible Chapter 1 database
  -> frozen H1-H5 hypothesis contract
  -> canonical progressive analysis
  -> frozen result + canonical manuscript
```

Development experiments and superseded analyses remain recoverable in git history. **`legacy/v1/` is a separately frozen historical analysis and is not modified by the v2 paper cleanup.**

## 1. Canonical database

The final Chapter 1 analysis uses a fixed universe of **8,265 islands** and a fixed trait denominator of **106,295 analysis-applicable plant species x 3 axes = 318,885 species-axis cells**.

The three analysis axes are:

- flower colour;
- floral structural complexity;
- reproductive assurance.

Core data layers are:

- island geography / area / distance infrastructure derived from the frozen GSHHG-based island universe;
- GBIF records used as an opportunistic sample of realised island floras, never as a complete census;
- frozen floristic source/status infrastructure used for native/non-endemic/endemic and source-pool comparisons;
- species-level trait evidence with provenance and quality labels.

Trait evidence precedence is **species-direct High/Medium first, then trait-specific Validated Low only where direct evidence is absent**. Family inference and global fallback are prohibited in the canonical Chapter 1 contract.

The current database-construction route is documented in [`docs/DATABASE_BUILD.md`](docs/DATABASE_BUILD.md). The authoritative trait-build line is the recovered checkpoint and source-scale integration merged through PR #151, with source manifests retained under `data/v2/staging/traits/source_batches/` and the evidence contract in `docs/chapter1_trait_coverage_contract.md`.

Final trait snapshot used by the paper:

- source run: **34191508045**;
- artifact: `source-scale-batch-integration-34191508045`;
- resolved cells: **222,688 / 318,885 = 69.83%**;
- reproductive-assurance cells: **48,497 / 106,295 = 45.63%**.

Coverage is not itself an analysis gate. Missingness and trait-resolution processes are audited explicitly.

## 2. Frozen hypotheses

There is one canonical scientific contract:

- [`config/chapter1_progressive_analysis.yml`](config/chapter1_progressive_analysis.yml)
- contract name: `chapter1_progressive_analysis_v1`

The hypothesis order is fixed:

1. **H1 — Universal-syndrome rival:** does one floral/reproductive island syndrome hold across contexts?
2. **H2 — Biogeographic branching:** do response vectors differ among contexts?
3. **H3 — Source/lineage assembly:** how much of the response is represented by source-available lineage entry/loading rather than a uniform within-lineage shift?
4. **H4 — Area/capacity moderation:** does continuous island area modify the isolation-associated filtering?
5. **H5 — Channel-gated residual mechanism:** only independent pollinator-channel evidence can promote a mechanistic channel-loss claim.

The claim ceiling is part of the same config. Floral architecture cannot be used to infer historical pollinator loss, replacement, source ancestry, or in-situ evolution.

## 3. Canonical analysis pipeline

The single paper-level workflow is:

- [`.github/workflows/run-chapter1-progressive-trait-analysis.yml`](.github/workflows/run-chapter1-progressive-trait-analysis.yml)

It snapshots the current trait ledger and reruns the same predeclared contract. The current implementation uses the Chapter 1 modules for:

- trait snapshot validation;
- trait-resolution / MNAR sensitivity;
- climate-overlap validation;
- pollination-associated architecture factorization;
- family/genus taxonomic-depth decomposition;
- area/capacity moderation;
- source-matched lineage representation.

Detailed routing is in [`docs/PAPER_PIPELINE.md`](docs/PAPER_PIPELINE.md).

Final paper analysis:

- workflow run: **34232450884**;
- artifact: `chapter1-progressive-analysis-34232450884`;
- artifact ID: **10058653212**;
- digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

## 4. Frozen result and manuscript

The final interpretation and provenance are fixed in:

- [`docs/chapter1_submission_freeze_20260909.md`](docs/chapter1_submission_freeze_20260909.md)
- [`docs/chapter1_submission_hypothesis_framework_20260909.md`](docs/chapter1_submission_hypothesis_framework_20260909.md)
- [`docs/chapter1_figure1_hypothesis_tree_spec_20260909.md`](docs/chapter1_figure1_hypothesis_tree_spec_20260909.md)
- [`docs/chapter1_literature_positioning_20260909.md`](docs/chapter1_literature_positioning_20260909.md)
- [`docs/chapter1_h3_h5_causal_hierarchy_20260909.md`](docs/chapter1_h3_h5_causal_hierarchy_20260909.md)

The **single canonical manuscript** is:

- [`docs/chapter1_manuscript_full_v7_submission_order_20260909.md`](docs/chapter1_manuscript_full_v7_submission_order_20260909.md)

Headline result:

- H1: one universal floral/reproductive island syndrome is not recovered;
- H2: source-separation responses branch among biogeographic contexts;
- H3: the broad Palearctic response is represented at genus-level lineage assembly beyond family composition;
- H4: area remains a measurement-sensitive modifier rather than a promoted causal result;
- H5: pollination-channel causation is not identified from the Chapter 1 plant database.

`large_bee_like`, `butterfly_like`, and `bird_like` are secondary floral-architecture concordance scores, **not pollinator classifiers**.

## 5. Repository map

```text
README.md                              <- start here
docs/DATABASE_BUILD.md                 <- how the analysis database was built
docs/PAPER_PIPELINE.md                 <- database -> H1-H5 -> analysis -> manuscript
config/chapter1_progressive_analysis.yml
.github/workflows/run-chapter1-progressive-trait-analysis.yml
src/island_v2/                         <- analysis/data-construction implementation
data/v2/                               <- current v2 data layers and provenance
docs/chapter1_*20260909.md             <- final paper framing and freeze
legacy/v1/                             <- frozen historical v1, preserved intact
```

Older PR-specific workflows, superseded August result notes, and baseline manuscript copies have been removed from the current working tree to prevent them from appearing canonical. They remain available in git history.

## 6. Legacy v1 preservation

`legacy/v1/` is intentionally retained as the complete historical v1 analysis, including its R code, data objects, metadata, trait files, artifacts, and migration documentation.

At the start of the paper-tree cleanup, the `legacy/v1` git tree was:

`8febaeb4e77f1c595f34dd672c95e5926fa58b0a`

The v2 cleanup must not alter this tree. v1 is historical provenance; it is not silently mixed into the Chapter 1 v2 estimand.
