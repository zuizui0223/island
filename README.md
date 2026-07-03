# Island floral syndrome — v2

This repository is now organized around the reproducible v2 pipeline:

```text
frozen exact island universe
→ GBIF download blocks and campaign ledger
→ archive collection and exact point-in-polygon assignment
→ taxonomic / floral-trait evidence tables
→ island-level inference
```

## Active code

- `src/island_v2/` — v2 Python package
- `config/` — current data-acquisition and inference configuration
- `data/v2/` — external, staging, curated, and template data layers
- `docs/` — v2 architecture, data policy, and operational notes
- `.github/workflows/` — validation, GBIF submission/polling/collection

## Frozen v1

The original standalone R analysis is preserved without edits under
[`legacy/v1/`](legacy/v1/). It is retained for provenance and is not
part of the active v2 workflow.

## Current campaign rule

GBIF request catchments are used only to retrieve candidate records.
Final assignment is always against the original, exact island polygons.
