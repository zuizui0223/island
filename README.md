# Island floral syndrome — v2

This repository supports **Chapter 1 / macroecological pattern-and-boundary analysis** of island floral and reproductive composition.

The current Chapter 1 result no longer supports a headline claim of a universal or confirmatorily biogeographically contingent floral island syndrome. The main contrast is now:

> **Does island isolation generate one coherent floral/reproductive syndrome, or does it reorganize particular floral components independently?**

The current frozen evidence favors **component-specific floral reorganization**:

> **Isolation is associated with changes in particular floral categories, while broad `generalized + plain + SC` syndrome summaries weaken or disappear after category preservation and lineage control. Strong within-context signals are concentrated in northern mid-latitude islands, but formal joint interaction tests do not yet establish confirmatory regional differences in the isolation effect.**

See [`docs/chapter1_frozen_result_20260825.md`](docs/chapter1_frozen_result_20260825.md) for the locked result, run IDs, FDR rules and claim boundary, [`docs/chapter1_pollination_syndrome_concordance_20260825.md`](docs/chapter1_pollination_syndrome_concordance_20260825.md) for the post-freeze Discussion-only literature audit, and [`THESIS_CHAPTER_POSITIONING.md`](THESIS_CHAPTER_POSITIONING.md) for the dissertation architecture.

## Manuscript status

The repository contains substantial development history. For submission, a result is canonical only when it follows the rules in:

- [`docs/manuscript_submission_contract.md`](docs/manuscript_submission_contract.md) — submission surface, evidence-tier rules, model safeguards, and archival requirements;
- [`docs/chapter1_frozen_result_20260825.md`](docs/chapter1_frozen_result_20260825.md) — current scientific result checkpoint;
- [`docs/chapter1_pollination_syndrome_concordance_20260825.md`](docs/chapter1_pollination_syndrome_concordance_20260825.md) — Discussion-only Bombus / bird / Lepidoptera concordance and mismatch audit;
- [`docs/v2_pollination_regime_framework.md`](docs/v2_pollination_regime_framework.md) — current scientific framework;
- [`docs/v2_channel_architecture.md`](docs/v2_channel_architecture.md) — measurement and interpretation guardrails;
- [`config/bombus_applicability.yml`](config/bombus_applicability.yml) — retained Bombus diagnostic provenance rules.

Historical pilots, scouting workflows, alternative model variants, and acquisition experiments are development records, not manuscript methods unless explicitly promoted into the submission contract.

## Primary reproducible data path

```text
frozen exact island universe
→ exact point-in-polygon flora assignment
→ locked direct trait evidence
→ floristic-status / endemicity support
→ lineage-preserving safeguard
→ locked geographic/environmental covariates
→ context-aware baseline
→ common isolation effect
→ formal isolation × context joint test
→ category-preserving trait decomposition
→ genus-fixed broad-outcome guardrail
→ result freeze
→ pollination-syndrome concordance or mismatch in Discussion
```

The primary analysis does **not** require classifying every island as Bombus retained / lost / absent.

## Current result in one paragraph

The canonical joint-Wald run fitted 43 atomic trait categories. Seventeen within-context isolation slopes survived FDR; 16 met the count component of confirmatory support and all 16 occurred in the northern-midlatitude context. However, only one atomic category showed FDR-supported joint `isolation × context` heterogeneity, and that comparison used only 32 versus 30 islands, so it remains pilot-level. Meanwhile the genus-composition-preserving M3 analysis did not retain confirmatory northern or tropical isolation slopes for the broad `generalized_form`, `plain_colour`, or `self_compatibility` outcomes. The manuscript therefore treats the main result as **trait-component-specific reorganization rather than a coherent classic syndrome or a confirmed region-specific syndrome**.

## Evidence tiers

Complete fill is not complete evidence. Trait resolution remains visible in every analysis.

- **Confirmatory:** direct source-backed species evidence.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback.

The `sensitivity_all` layer may be useful for stress-testing conclusions, but it is not itself confirmatory evidence.

## Bombus and other pollinator-syndrome interpretation boundary

Existing Bombus applicability, environmental-niche and occurrence-diagnostic products remain in the repository for provenance, exploratory diagnostics and sensitivity analyses. They are **not primary predictors for Chapter 1**.

The retained environmental-niche estimator measures climatic-environmental compatibility. It is not realized occurrence probability, abundance, visitation rate, pollination service, or evidence of historical loss. Opportunistic non-detection is not treated as proof of absence.

The frozen northern-midlatitude trait vector may be compared with Bombus / large- or long-tongued-bee literature only for **partial concordance or mismatch**. Because the broad `generalized + plain + SC` signal is not recovered after the M3 lineage guardrail, the result must not be described as a demonstrated Bombus-loss syndrome.

The post-freeze concordance audit finds only broad conceptual agreement with the idea that island pollinator functional change can reorganize floral architecture. It does **not** recover the Izu-style SC shift or a simple loss-of-specialization axis. Likewise, the present data do not support a bird- or Lepidoptera-mediated tropical/southern counter-syndrome.

## Thesis handoff

```text
Chapter 1 — island
isolation
→ component-specific floral reorganization
→ test whether effects are shared or genuinely context dependent
→ candidate ecological concordance / mismatch

Chapter 2 — izu-core
which interaction states make particular components branch,
propagate, or buffer?

Chapter 3 — shimahotarubukuro
focal lineage
→ realized multidimensional floral divergence
```

Chapter 1 identifies **what parts of the floral/reproductive phenotype reorganize under isolation and where signals are detectable**. It deliberately leaves the question “which interaction mechanism makes those components respond?” for the stricter mechanistic framework of Chapter 2.

## Repository layout

- `src/island_v2/` — reusable v2 data and analysis utilities
- `analysis/v2/` — statistical analysis scripts and execution notes
- `config/` — frozen contracts, ontology, and artifact locks
- `data/v2/` — external/staging/curated/template data layers
- `docs/` — scientific design, data policy, methods, and reproducibility notes
- `.github/workflows/` — active validation, materialization, and analysis workflows
- `legacy/v1/` — frozen v1 provenance only

## Reproducibility rule before submission

GitHub Actions artifacts are temporary and are not a permanent supplement. The manuscript release must archive all critical inputs and outputs durably, record checksums, identify one canonical workflow per main analysis, and report attrition from the frozen 8,265-island universe to every fitted model.

GBIF request catchments are retrieval devices only. Final occurrence assignment is always against the original exact island polygons.
