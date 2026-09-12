# PR183 — Database 2.0 alpha2 full taxonomy normalization

## Goal

Build a complete, auditable taxonomy crosswalk for all 115,328 provisional plant-name identities frozen in Database 2.0-alpha1, using GBIF v2 species matching against the current Catalogue of Life Extended Release (COL XR).

## Frozen source

- Database source: `2.0.0-alpha1`
- provisional taxa: 115,328
- source table: `data/v2/staging/gbif/collected/island_taxa.csv`
- source SHA-256: `c0586264d9a877b88c26e866264f2c77b20423cd93380a8eebe269582d53352e`
- alpha1 is immutable; this PR adds a separate crosswalk layer only.

## Qualified matcher

PR182 qualified the production GBIF v2 batch matcher against COL XR checklist
`7ddf754f-d193-4cc9-b351-99906754a03b`.

The 512-name deterministic pilot produced 500 EXACT matches and, after the synonym-aware accepted-target correction, 492/512 (96.1%) conservative automatic-resolution candidates. All automatic results remain candidates until the alpha2 layer itself is frozen.

## Full-run architecture

The source is split deterministically by its alphabetically stable row order:

- segments 0–10: 10,000 taxa each
- segment 11: 5,328 taxa
- total: 115,328 taxa
- maximum live segment concurrency: 2
- GBIF v2 POST batch size: 1,000 names

Each segment independently stores:

- compressed normalized crosswalk;
- COL XR metadata snapshot;
- raw request/response batches;
- segment manifest with source/checklist/metadata/crosswalk SHA-256 values;
- resolution summary.

The combine stage refuses to emit a full artifact unless all 12 segments are present, contiguous, non-overlapping, source-identical, checklist-identical, and metadata-identical.

## Conservative resolution policy

Automatic-resolution candidates require:

- `EXACT` match type;
- confidence >=95;
- no processing flags;
- no issues;
- Plantae classification;
- species-rank accepted target;
- either an accepted matched usage or an explicit `acceptedUsage` target from a matched synonym.

`VARIANT`, `CANONICAL`, `AMBIGUOUS`, `HIGHERRANK`, non-species targets, provisionally accepted usages, flagged matches, low-confidence matches, and unmatched names remain in the review queue.

## Merge gate

Merge only after the live full workflow produces:

1. 12 successful segment receipts;
2. exactly 115,328 unique submitted names in the combined crosswalk;
3. automatic + review-required + unmatched = 115,328;
4. one consistent COL XR metadata snapshot across all segments;
5. a full artifact containing `taxonomy_crosswalk.csv.gz`, `taxonomy_review_queue.csv.gz`, `gbif_colxr_metadata.json`, and `TAXONOMY_MANIFEST.json`;
6. no mutation of alpha1 candidate-flora or trait data.
