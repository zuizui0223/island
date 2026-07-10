# v2 data layout

v2 begins with new, traceable data products and does not use v1 derived tables as analysis inputs.

```text
templates/
  taxa_template.csv                              accepted taxon batch for reviewable trait work
  accepted_evidence_template.csv                 one accepted human-reviewed trait record per row
  review_decisions_template.csv                  explicit accept/reject/adjudication decisions
  source_region_bombus_registry_template.csv    source-region applicability proposal/review table
  source_region_evidence_template.csv           Bombus-status evidence for a source region
  island_source_region_assignment_template.csv  frozen-island -> source-region assignment
  island_assignment_evidence_template.csv       geography/source-pool/dispersal evidence for that assignment

external/
  islands/                                       source island polygons and normalized island manifest
  gbif/                                          exact-polygon request manifests and downloaded archives

staging/
  gbif/                                          raw island-by-name candidates before taxonomy review
  pollinators/                                   exact-island Apidae/Bombus occurrence rows and diagnostics
  ...                                            raw source-discovery outputs and review queues

curated/                                         versioned, human-adjudicated trait and source-region evidence
```

## First v2 data product: GBIF island species-name candidates

1. Put `GlbIslands.gdb` / `BigIslands` or a documented replacement vector source in `external/islands/`.
2. Run `island-v2-gbif-flora prepare-islands` to create `islands_v2.gpkg` and `island_manifest.csv`.
3. Run `make-requests`, `submit`, `poll`, and `collect` in bounded batches.
4. Treat the output as occurrence-based candidate names, not a verified flora.
5. Normalize taxonomy, then audit establishment status and the observation process before trait evidence review.

The full policy and commands are in `docs/gbif_flora_v2.md`.

## Bombus applicability comes before pollinator outcomes

Before retrieving or interpreting Bombus records, populate and review the four
Phase 0.5 tables. Two distinct evidence tracks are mandatory:

- **source-region Bombus-status evidence** supports whether native Bombus occurs in the proposed source region;
- **island-assignment evidence** supports the island's geographic position, floristic source pool, and plausible dispersal connection.

These claims must never reuse the same evidence identifier merely because they
are both relevant to a later applicability decision. Agent-drafted rows remain
`unresolved` and cannot create a frozen `applicable` or
`structurally_not_applicable` classification until a human reviewer changes the
region, assignment, and their cited evidence to `accepted`.

```text
island-v2-bombus-applicability build
  --island-manifest-csv <frozen island_manifest.csv>
  --island-source-assignment-csv <reviewed assignments.csv>
  --source-region-registry-csv <reviewed source regions.csv>
  --source-region-evidence-csv <reviewed Bombus-status evidence.csv>
  --island-assignment-evidence-csv <reviewed geography/source-pool evidence.csv>
  --output-dir <applicability outputs>
```

The command produces a complete island-level registry. Missing assignments and
pending agent drafts are written explicitly as `unresolved`; they are not
omitted and cannot become Bombus-absence labels.

## Pollinator records are a separate raw data product

A separately declared Apidae GBIF campaign reuses the frozen islands and exact
point-in-polygon assignment. Its collector writes:

```text
island_pollinator_occurrences.csv
island_pollinator_record_coverage.csv
pollinator_collection_manifest.csv
pollinator_collection_status.json
```

These retain row-level taxon, source, date, coordinate, and record-basis data
for the effort-aware Bombus diagnostic. They are not a verified pollinator
service dataset and do not by themselves establish biological absence.

## Trait source discovery comes after Core-pilot nomination

`island-v2-trait-source-discovery` is a bounded, free source-discovery tool for
staged taxa from an island that has passed Core-pilot nomination. It may collect
bibliographic seeds, public-access receipts, and public-PDF candidate-page
locators. It never accepts a trait value, native/establishment status, Bombus
applicability, or analysis-inclusion decision. Those require the later review
and curation tables.
