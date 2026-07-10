# v2 taxon intake: names, broad group, and island establishment before traits

## Purpose

The GBIF collector produces a broad Tracheophyta occurrence candidate universe.
It is not yet a verified island flora and it is not suitable as direct input to
floral-trait extraction. In particular, a collector taxon table can contain:

- genus-only or otherwise non-species names;
- synonyms or fuzzy backbone matches;
- ferns, lycophytes, and gymnosperms;
- cultivated or introduced plants recorded on an island; and
- genuinely native angiosperm candidates.

`island-v2-taxon-intake build` makes these distinctions explicit as **review
queues**. It never changes a curated table and never makes an island-native or
accepted-angiosperm decision automatically.

## Inputs

```text
island_taxa.csv
island_species_occurrences.csv
```

These are the cumulative exact-island products emitted by the GBIF collector.

## Outputs

```text
taxon_intake_review_queue.csv
island_establishment_review_queue.csv
taxon_scope_decisions_prefilled.csv
taxon_intake_summary.json
gbif_backbone_raw/*.json
```

`taxon_intake_review_queue.csv` has a row per raw collector taxon. The GBIF
Backbone fields are candidates only. `taxon_scope_review_status` remains
`pending` in every row.

`island_establishment_review_queue.csv` has a row per island x raw taxon. GBIF
occurrence records are not evidence of native establishment, so every row starts
as `island_establishment_status=unresolved` and `review_status=pending`.

`taxon_scope_decisions_prefilled.csv` is a convenience starting point for the
reviewed decision table consumed by `island-v2-taxon-scope build`; it must be
reviewed before use.

## Review order

```text
1. Resolve name rank and accepted taxon.
2. Confirm broad group: angiosperm / gymnosperm / non-seed vascular / unresolved.
3. Review island establishment using an island or regional flora source.
4. Accept only the relevant angiosperm taxa into the taxon-scope table.
5. Build accepted-angiosperm coverage.
6. Send only reviewed eligible taxa to trait evidence extraction.
```

The GBIF Backbone endpoint helps with steps 1-2. It does **not** substitute for
steps 3-4, and its response cannot establish island nativeness.

## Running the intake

Use the `build_taxon_intake` GitHub Actions workflow after collection has
produced a meaningful batch of raw taxa. It uploads a versioned artifact and
does not write to curated data. The CLI equivalent is:

```text
island-v2-taxon-intake build \
  --taxa-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --island-species-csv data/v2/staging/gbif/collected/island_species_occurrences.csv \
  --output-dir data/v2/staging/taxon_intake
```

For a reproducible offline rebuild using previously cached responses, add
`--offline`.
