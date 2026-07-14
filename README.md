# Island floral syndrome — v2

This repository supports a global comparative analysis of island floral and reproductive composition. The working hypothesis is **conditional**, not universal: island isolation may be associated with different floral outcomes in different pollination regimes, and a Bombus-channel hypothesis is evaluated only where that channel is biologically interpretable.

The analysis does **not** treat absence of a Bombus record as pollinator loss and does not claim that cross-sectional island associations prove causal floral evolution.

## Manuscript status

The repository contains substantial development history. For submission, a result is canonical only when it follows the rules in:

- [`docs/manuscript_submission_contract.md`](docs/manuscript_submission_contract.md) — submission surface, evidence-tier rules, model safeguards, and archival requirements;
- [`docs/v2_pollination_regime_framework.md`](docs/v2_pollination_regime_framework.md) — scientific scope and conditional pollination-regime framework;
- [`docs/v2_channel_architecture.md`](docs/v2_channel_architecture.md) — channel measurement architecture;
- [`config/bombus_applicability.yml`](config/bombus_applicability.yml) — outcome-blind applicability rules.

Historical pilots, scouting workflows, alternative model variants, and acquisition experiments are development records, not manuscript methods unless explicitly promoted into the submission contract.

## Reproducible data path

```text
frozen exact island universe
→ exact point-in-polygon flora assignment
→ locked trait evidence with explicit evidence tier
→ locked geographic/environmental covariates
→ Bombus environmental/occurrence diagnostics with provenance
→ attrition and model-support audit
→ global composition analysis
→ northern-midlatitude conditional Bombus analysis
→ declared sensitivity analyses
```

## Evidence tiers

Complete fill is not complete evidence. Trait resolution must remain visible in every analysis.

- **Confirmatory:** direct source-backed species evidence.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback.

The `sensitivity_all` layer may be useful for stress-testing conclusions, but it is not itself confirmatory evidence.

## Bombus interpretation boundary

The retained environmental-niche estimator measures climatic-environmental compatibility. It is not realized occurrence probability, source-pool membership, abundance, visitation rate, pollination service, or evidence of historical loss.

A global all-species environmental maximum is diagnostic only. Any analysis-ready Bombus predictor must carry explicit provenance and satisfy the downstream semantic guardrails.

## Repository layout

- `src/island_v2/` — reusable v2 data and analysis utilities
- `analysis/v2/` — statistical analysis scripts
- `config/` — frozen contracts, ontology, and artifact locks
- `data/v2/` — external/staging/curated/template data layers
- `docs/` — scientific design, data policy, methods, and reproducibility notes
- `.github/workflows/` — active validation, materialization, and analysis workflows
- `legacy/v1/` — frozen v1 provenance only

## Reproducibility rule before submission

GitHub Actions artifacts are temporary and are not a permanent supplement. The manuscript release must archive all critical inputs and outputs durably, record checksums, identify one canonical workflow per main analysis, and report the attrition from the frozen 8,265-island universe to every fitted model.

GBIF request catchments are retrieval devices only. Final occurrence assignment is always against the original exact island polygons.
