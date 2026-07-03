# v2 GBIF island-flora rebuild

## What is being rebuilt

The v1 flora step used an island bounding box, a `scientificName` facet capped at 2,000 names, and treated request failure as an empty list. v2 does none of these.

For each standardized island polygon, v2 makes a GBIF Occurrence Download request with an exact `within` WKT geometry. GBIF's download service is asynchronous and returns a download key, status, link, and DOI. These identifiers are retained for each island-specific query.

## Why species names are acquired before trait coding

The first output is deliberately an **island-by-raw-name candidate table**, not a final flora and not yet the analysis matrix. This keeps three later decisions separate:

1. taxonomic normalization against a declared backbone;
2. evidence that a record represents an established member of the island flora rather than cultivation, misgeoreferencing, or transient occurrence; and
3. trait acquisition.

## Initial query policy

```text
exact island polygon, never a bounding box
angiosperm target clade under a declared taxonomy
hasCoordinate = true
hasGeospatialIssue = false
GBIF SPECIES_LIST download format
no first-pass year cut-off
no first-pass basis-of-record restriction
```

The absence of a year cut-off is intentional. A 2010–2025-only flora can turn uneven recent survey effort into false species absences. Record dates, data-source composition, basis of record, coordinate uncertainty, and cultivated/native evidence are retained for a later observation-process audit.

## Required input

The original `GlbIslands.gdb` / `BigIslands` source is not present in this repository. Place that source, or an explicitly documented replacement island polygon layer, at the path declared in `config/gbif_flora_v2.yml`.

The first command will produce a v2-normalized GeoPackage and manifest:

```bash
island-v2-gbif-flora prepare-islands \
  --input-path data/v2/external/islands/GlbIslands.gdb \
  --layer BigIslands \
  --output-dir data/v2/external/islands/prepared \
  --source-label global_islands_source \
  --min-area-km2 20
```

It repairs invalid geometry, splits multipart features into explicit polygon units, calculates equal-area km2, retains parent-feature identifiers, and assigns stable geometry-derived IDs.

## GBIF request lifecycle

```text
prepare-islands
    -> islands_v2.gpkg + island_manifest.csv
make-requests
    -> exact-WKT request_manifest.csv + taxon-resolution JSON
submit (bounded batches)
    -> GBIF download keys; failures explicitly visible
poll
    -> status, download link, DOI
collect
    -> island_species_raw.csv
```

GBIF credentials are required only at `submit` time and belong in `GBIF_USERNAME` and `GBIF_PASSWORD` environment variables or GitHub Action secrets. The account notification email may be supplied when requests are built.

## Minimum execution commands

```bash
island-v2-gbif-flora make-requests \
  --islands-gpkg data/v2/external/islands/prepared/islands_v2.gpkg \
  --output-csv data/v2/external/gbif/request_manifest.csv \
  --checklist-key 7ddf754f-d193-4cc9-b351-99906754a03b

export GBIF_USERNAME='your-gbif-username'
export GBIF_PASSWORD='your-gbif-password'
island-v2-gbif-flora submit \
  --request-manifest data/v2/external/gbif/request_manifest.csv \
  --max-requests 20

island-v2-gbif-flora poll \
  --request-manifest data/v2/external/gbif/request_manifest.csv

island-v2-gbif-flora collect \
  --request-manifest data/v2/external/gbif/request_manifest.csv \
  --download-dir data/v2/external/gbif/downloads \
  --output-species-csv data/v2/staging/gbif/island_species_raw.csv
```

## Interpretation boundary

A GBIF species-list download is occurrence evidence, not a verified island flora. Its role is a reproducible candidate set. The later data-process module will distinguish specimen support, other basis-of-record classes, spatial coverage, record-year structure, data-set concentration, cultivation/establishment status, and source-pool alternatives.
