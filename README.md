# Island floral syndrome — v2

This repository supports **Chapter 1 / macroecological pattern-and-boundary analysis** of island floral and reproductive composition.

The working hypothesis is **biogeographically contingent, not universal**:

> **Island isolation does not generate one universal floral and reproductive syndrome. The magnitude and direction of isolation-associated trait filtering can differ among biogeographic contexts after floristic status, lineage/source-pool structure, climate, island geography and observation support are represented.**

The primary Chapter 1 question is therefore:

> **When and where do different floral and reproductive island syndromes emerge?**

Pollinator identity is not used as the primary causal answer to that question. Bombus-, bird-, Lepidoptera- and other pollination-syndrome expectations are used primarily as **post-result ecological concordance / discussion frameworks** unless independent mechanistic evidence exists.

See [`THESIS_CHAPTER_POSITIONING.md`](THESIS_CHAPTER_POSITIONING.md) for the dissertation architecture.

## Manuscript status

The repository contains substantial development history. For submission, a result is canonical only when it follows the rules in:

- [`docs/manuscript_submission_contract.md`](docs/manuscript_submission_contract.md) — submission surface, evidence-tier rules, model safeguards, and archival requirements;
- [`docs/v2_pollination_regime_framework.md`](docs/v2_pollination_regime_framework.md) — current context-dependent scientific framework;
- [`docs/v2_channel_architecture.md`](docs/v2_channel_architecture.md) — measurement and interpretation guardrails;
- [`config/bombus_applicability.yml`](config/bombus_applicability.yml) — retained Bombus diagnostic provenance rules.

Historical pilots, scouting workflows, alternative model variants, and acquisition experiments are development records, not manuscript methods unless explicitly promoted into the submission contract.

## Primary reproducible data path

```text
frozen exact island universe
→ exact point-in-polygon flora assignment
→ locked trait evidence with explicit evidence tier
→ floristic-status / endemicity support
→ lineage / source-pool-aware expectation
→ locked geographic/environmental covariates
→ attrition and model-support audit
→ global isolation baseline
→ isolation × biogeographic-context analysis
→ category-preserving trait decomposition
→ pollination-syndrome concordance in Discussion
→ declared sensitivity analyses
```

The primary analysis does **not** require classifying every island as Bombus retained / lost / absent.

## Evidence tiers

Complete fill is not complete evidence. Trait resolution must remain visible in every analysis.

- **Confirmatory:** direct source-backed species evidence.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback.

The `sensitivity_all` layer may be useful for stress-testing conclusions, but it is not itself confirmatory evidence.

## Bombus interpretation boundary

Existing Bombus applicability, environmental-niche and occurrence-diagnostic products remain in the repository for provenance, exploratory diagnostics and sensitivity analyses.

They are **not required to define the primary Chapter 1 model**.

The retained environmental-niche estimator measures climatic-environmental compatibility. It is not realized occurrence probability, abundance, visitation rate, pollination service, or evidence of historical loss. Opportunistic non-detection is not treated as proof of absence.

If a northern-temperate trait vector resembles patterns reported in systems with reduced Bombus / large- or long-tongued-bee function, that resemblance may be discussed as **concordance with a candidate pollination mechanism**, not as causal identification.

The same rule applies symmetrically to tropical or southern patterns that resemble bird-, Lepidoptera-, hawkmoth- or other pollination syndromes.

## Thesis handoff

```text
Chapter 1 — island
isolation × biogeographic context
→ regional floral / reproductive trait syndromes
→ candidate ecological interpretations

Chapter 2 — izu-core
candidate interaction mechanisms
→ branching / propagation / buffering

Chapter 3 — shimahotarubukuro
focal lineage
→ realized multidimensional floral divergence
```

Chapter 1 identifies the **pattern and boundary conditions**. It deliberately leaves the question “which pollination mechanism actually generated the pattern?” for the stricter mechanistic framework of Chapter 2.

## Repository layout

- `src/island_v2/` — reusable v2 data and analysis utilities
- `analysis/v2/` — statistical analysis scripts
- `config/` — frozen contracts, ontology, and artifact locks
- `data/v2/` — external/staging/curated/template data layers
- `docs/` — scientific design, data policy, methods, and reproducibility notes
- `.github/workflows/` — active validation, materialization, and analysis workflows
- `legacy/v1/` — frozen v1 provenance only

## Reproducibility rule before submission

GitHub Actions artifacts are temporary and are not a permanent supplement. The manuscript release must archive all critical inputs and outputs durably, record checksums, identify one canonical workflow per main analysis, and report attrition from the frozen 8,265-island universe to every fitted model.

GBIF request catchments are retrieval devices only. Final occurrence assignment is always against the original exact island polygons.
