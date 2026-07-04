# v2 data layout

v2 begins with new, traceable data products and does not use v1 derived tables as analysis inputs.

```text
templates/
  taxa_template.csv                              accepted taxon batch for web-search trait work
  accepted_evidence_template.csv                 one accepted human-reviewed trait record per row
  review_decisions_template.csv                  explicit accept/reject/adjudication decisions
  source_region_bombus_registry_template.csv    reviewed Bombus applicability at source-region level
  island_source_region_assignment_template.csv  reviewed frozen-island → source-region assignment
  source_region_evidence_template.csv           evidence supporting each source-region decision

external/
  islands/                                       source island polygons and normalized island manifest
  gbif/                                          exact-polygon request manifests and downloaded archives

staging/
  gbif/                                          raw island-by-name candidates before taxonomy review
  pollinators/                                   exact-island Apidae/Bombus occurrence rows and diagnostics
  ...                                            raw LLM candidate exports and review queues

curated/                                         versioned, human-adjudicated trait and source-region evidence
```

## First v2 data product: GBIF island species-name candidates

1. Put `GlbIslands.gdb` / `BigIslands` or a documented replacement vector source in `external/islands/`.
2. Run `island-v2-gbif-flora prepare-islands` to create `islands_v2.gpkg` and `island_manifest.csv`.
3. Run `make-requests`, `submit`, `poll`, and `collect` in bounded batches.
4. Treat the output as occurrence-based candidate names, not a verified flora.
5. Normalize taxonomy, then audit establishment status and the observation process before LLM trait collection.

The full policy and commands are in `docs/gbif_flora_v2.md`.

## Bombus applicability comes before pollinator outcomes

Before retrieving or interpreting Bombus records, populate and review the three
source-region templates. Then run:

```text
island-v2-bombus-applicability build
  --island-manifest-csv <frozen island_manifest.csv>
  --island-source-assignment-csv <reviewed assignments.csv>
  --source-region-registry-csv <reviewed source regions.csv>
  --source-region-evidence-csv <reviewed source evidence.csv>
  --output-dir <applicability outputs>
```

The command produces a complete island-level registry. Missing assignments are
written explicitly as `unresolved`; they are not omitted and cannot become a
Bombus-absence label.

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

## LLM trait workflow follows species-name review

After the candidate species list is normalized, run the LLM web-search workflow in batches. It returns source-cited candidates and never directly edits curated trait data. Genus/family distributions are generated only from reviewed evidence.
