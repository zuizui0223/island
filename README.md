# Island floral syndrome  -  v2

This repository is organized around a reproducible v2 pipeline for a
**conditional** island floral-syndrome analysis. It does not assume that islands
share one universal Bombus-driven white-flower or inconspicuous-flower syndrome.

The framework first asks how isolation, source pools, establishment, observation
processes, and pollination-function regimes shape whole-flora pollen-vector and
reproductive-assurance composition. It then tests a Bombus-channel hypothesis
only where Bombus is biologically applicable to the island's predeclared
source-region context.

```text
frozen exact island universe
-> GBIF download blocks and campaign ledger
-> archive collection and exact point-in-polygon assignment
-> flora and Bombus observation-process diagnostics
-> bulk taxonomic / trait evidence tables
-> coverage and attrition audit
-> conditional island-level inference
```

The governing design documents are:

- [`docs/v2_pollination_regime_framework.md`](docs/v2_pollination_regime_framework.md)  -  scientific scope, island strata, preregistration, models, and continuation rules;
- [`docs/v2_channel_architecture.md`](docs/v2_channel_architecture.md)  -  compact channel-audit architecture;
- [`docs/v2_trait_acquisition_quantitative_plan.md`](docs/v2_trait_acquisition_quantitative_plan.md)  -  category-first trait acquisition with optional quantitative enrichment;
- [`docs/public_web_trait_shard_campaign.md`](docs/public_web_trait_shard_campaign.md)  -  128-shard public-web extraction, evidence packets, checkpoints, and resume;
- [`config/bombus_applicability.yml`](config/bombus_applicability.yml)  -  outcome-blind rule for the biological applicability of the Bombus hypothesis;
- [`config/bombus_observation_diagnostics.yml`](config/bombus_observation_diagnostics.yml)  -  effort-aware Bombus detection and non-detection policy.

## Active code

- `src/island_v2/`  -  v2 Python package
- `config/`  -  current data-acquisition and inference configuration
- `data/v2/`  -  external, staging, curated, and template data layers
- `docs/`  -  v2 architecture, data policy, and operational notes
- `.github/workflows/`  -  validation, GBIF submission/polling/collection, and manual trait acquisition

## Global trait acquisition

The expected species master is too large for one-by-one web or LLM retrieval to
be the global coverage path. The primary route is a downloaded, public bulk trait
source joined to `island_taxa.csv` with a source-specific codebook profile.
The 128-shard public-web campaign is a resumable gap-recovery layer after those
bulk and machine candidates are applied; it is not a replacement for bulk joins.

`ingest_bulk_trait_source` is the manual GitHub Actions workflow for this route.
It needs a direct source-download URL and **no OpenAI key**. Each run creates
pending source-backed candidates plus taxon-match, unmapped-code, trait-coverage,
and family-coverage audits as an artifact. It never commits downloaded data or
curates trait values automatically. Details are in
[`docs/bulk_trait_acquisition.md`](docs/bulk_trait_acquisition.md).

Low-confidence is not the same as wrong. Declared weak-source and taxonomic
inference candidates remain auditable in the broad track, while logical errors
such as inferring autonomous selfing from self-compatibility alone are quarantined.

The single scoreboard for how far acquisition has actually got is
`island-v2-acquisition-rate`, which unions every channel against the master
species set and reports the broad, source-backed coverage rate per trait. It is
the confirm-instrument for the rate-first workflow: land a source, rerun, read
the delta. Details are in
[`docs/acquisition_rate_report.md`](docs/acquisition_rate_report.md).

## Frozen v1

The original standalone R analysis is preserved without edits under
[`legacy/v1/`](legacy/v1/). It is retained for provenance and is not part of
the active v2 workflow.

## Current campaign rule

GBIF request catchments are used only to retrieve candidate records. Final
assignment is always against the original, exact island polygons.
