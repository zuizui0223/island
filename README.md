# Island floral syndrome — v2

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
→ GBIF download blocks and campaign ledger
→ archive collection and exact point-in-polygon assignment
→ flora and Bombus observation-process diagnostics
→ taxonomic / trait evidence tables
→ coverage and attrition audit
→ conditional island-level inference
```

The governing design documents are:

- [`docs/v2_pollination_regime_framework.md`](docs/v2_pollination_regime_framework.md) — scientific scope, island strata, preregistration, models, and continuation rules;
- [`docs/v2_channel_architecture.md`](docs/v2_channel_architecture.md) — compact channel-audit architecture;
- [`config/bombus_applicability.yml`](config/bombus_applicability.yml) — outcome-blind rule for the biological applicability of the Bombus hypothesis;
- [`config/bombus_observation_diagnostics.yml`](config/bombus_observation_diagnostics.yml) — effort-aware Bombus detection and non-detection policy.

## Active code

- `src/island_v2/` — v2 Python package
- `config/` — current data-acquisition and inference configuration
- `data/v2/` — external, staging, curated, and template data layers
- `docs/` — v2 architecture, data policy, and operational notes
- `.github/workflows/` — validation, GBIF submission/polling/collection

## Frozen v1

The original standalone R analysis is preserved without edits under
[`legacy/v1/`](legacy/v1/). It is retained for provenance and is not part of
the active v2 workflow.

## Current campaign rule

GBIF request catchments are used only to retrieve candidate records. Final
assignment is always against the original, exact island polygons.
