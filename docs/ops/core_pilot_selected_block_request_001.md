# Core-pilot selected GBIF block request 001

## Action requested

Run the manual workflow `collect_gbif_selected_block` with:

- `block_id`: `gbif_block_cell_w+120.0_s+30.0_003`
- `selection_rationale`: `Predeclared northern-temperate source-region coverage for Core-pilot intake; order override only, not based on traits, Bombus island occurrence or non-detection, reproductive outcomes, or model results.`

## Why this is safe

The block is already recorded in the frozen plant acquisition campaign as a `succeeded` GBIF block. The selected-block workflow uses the same frozen island universe, exact-island assignment rule, campaign ledger, and collection code as the scheduled collector. This request changes only collection order.

## Boundaries

This request must not be interpreted as accepting any of the following:

- taxonomy;
- island native or establishment status;
- source-region applicability;
- trait values;
- Bombus presence, absence, or non-detection;
- reproductive outcomes;
- analysis inclusion.

Any resulting raw exact-island plant labels remain non-curated inputs for downstream review queues only.

## Expected post-run checks

After the workflow completes, inspect `data/v2/staging/gbif/collected/collection_status.json` and `collection_manifest.csv`. If meaningful raw exact-island records are added, then rebuild or validate only non-curated downstream artifacts: taxon-intake review queues, Bombus applicability outputs if their reviewed source-region inputs exist, and Core-pilot nomination reports. Do not move any rows into curated evidence without human review.
