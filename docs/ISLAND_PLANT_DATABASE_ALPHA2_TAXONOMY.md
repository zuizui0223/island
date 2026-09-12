# Island Plant Database 2.0 alpha2 taxonomy normalization

## Purpose

Alpha1 freezes 115,328 provisional species-name identities derived from the completed
exact-island GBIF occurrence campaign. Alpha2 adds a **taxonomy crosswalk**; it does not
rewrite alpha1 in place.

The normalization target is GBIF's current Catalogue of Life Extended Release (COL XR):

- checklist key: `7ddf754f-d193-4cc9-b351-99906754a03b`
- GBIF v2 match endpoint: `https://api.gbif.org/v2/species/match`
- match metadata endpoint: `https://api.gbif.org/v2/species/match/metadata`
- COL XR registry record: `https://www.gbif.org/dataset/7ddf754f-d193-4cc9-b351-99906754a03b`

The legacy GBIF Backbone is not the normalization target because GBIF has discontinued
its updates and recommends migration to COL XR.

## Why a batch pilot first

The v2 matcher accepts up to 1,000 `NameUsageQuery` objects per POST request and returns
results in input order. This is operationally suitable for 115,328 taxa, but the matching
semantics must be qualified before the full run.

The first pilot uses an alphabetically stable, evenly spaced sample of 512 provisional
taxa. It records:

- exact query objects sent to GBIF;
- raw v2 responses;
- current matcher/index metadata;
- a normalized crosswalk;
- summary counts by match and resolution class.

No pilot result is allowed to mutate alpha1.

## Query policy

Each provisional alpha1 name is sent as:

- `scientificName`: the provisional species string;
- `taxonRank`: `SPECIES`;
- `kingdom`: `Plantae`;
- `strict`: `true`, to avoid silently accepting a higher-rank fallback;
- `genus` and `family`: retained occurrence-derived context when available;
- COL XR checklist key supplied to the endpoint.

The submitted genus/family are context, not independent accepted taxonomy.

## Promotion gate

A match becomes an **automatic-resolution candidate** only when all of the following hold:

1. v2 `matchType == EXACT`;
2. confidence is at least 95;
3. no v2 processing flags;
4. no v2 match issues;
5. the accepted target is rank `SPECIES`;
6. the accepted target status is `ACCEPTED`;
7. accepted classification is in kingdom `Plantae`.

An exact synonym may therefore point to an accepted species candidate, while preserving
both the matched synonym usage and accepted usage identifiers.

`VARIANT`, `CANONICAL`, `AMBIGUOUS`, `HIGHERRANK`, low-confidence, flagged, doubtful,
non-species, and otherwise non-qualified matches stay `review_required`. `NONE` stays
`unmatched`.

Even automatic-resolution candidates are not promoted into alpha1. They are candidates
for the next versioned taxonomy layer.

## Provenance and versioning

The pilot pins the alpha1 source summary SHA-256:

`c0586264d9a877b88c26e866264f2c77b20423cd93380a8eebe269582d53352e`

It also snapshots current GBIF match metadata, because COL XR is updated regularly and
an alpha2 release must identify the taxonomy/index actually used.

COL XR's GBIF registry record is licensed CC BY 4.0. The crosswalk preserves that source
licence and the checklist key. Database-wide redistribution decisions remain a separate
rights layer.

## Full-run gate

Do not run all 115,328 taxa until the pilot has established that:

- POST batch schema works on production;
- response cardinality is stable;
- exact/synonym/ambiguous/unmatched cases normalize correctly;
- metadata can be pinned;
- the automatic-resolution gate is conservative enough for observed data.

After qualification, the same code can run deterministic <=1,000-name shards with raw
response receipts and resumability.
