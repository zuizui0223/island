# v2 GBIF island-flora rebuild

## What is being rebuilt

The v1 flora step used an island bounding box, a `scientificName` facet capped at 2,000 names, and treated request failure as an empty list. v2 does none of these.

For each standardized island polygon, v2 makes a GBIF Occurrence Download request with an exact `within` WKT geometry. GBIF's download service is asynchronous and returns a download key, status, link, and DOI. These identifiers are retained for each island-specific query.

## Public island source

v2 no longer depends on the locally stored legacy `GlbIslands.gdb`. It starts from **GSHHG 2.3.7**, the public global shoreline polygon dataset distributed by NOAA and the University of Hawai'i. GSHHG distributes ESRI shapefiles and encodes land-ocean boundaries as hierarchical closed polygons. The version, download URL, archive checksum, spatial resolution, area thresholds, and exact prepared geometries are written to the v2 manifest.

The first public-source command is:

```bash
island-v2-gshhg build \
  --config-path config/island_source_gshhg.yml \
  --output-dir data/v2/external/islands/gshhg
```

It downloads the pinned archive, extracts the high-resolution (`h`) L1 land-ocean layer, recombines split components that share a GSHHG polygon ID, calculates equal-area polygon size, and writes:

```text
data/v2/external/islands/gshhg/prepared/islands_v2.gpkg
data/v2/external/islands/gshhg/prepared/island_manifest.csv
data/v2/external/islands/gshhg/prepared/source_policy.json
```

The initial island universe retains L1 landmasses of at least 5 km2 (the primary operational threshold in `config/island_source_gshhg.yml`) and excludes continental-scale L1 landmasses above 7,000,000 km2. 1 km2 and 20 km2 are re-run as explicitly labelled sensitivity analyses. This threshold is explicitly versioned and will be assessed as a source-definition sensitivity analysis; it is not a hidden legacy classification.

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

The absence of a year cut-off is intentional. A 2010-2025-only flora can turn uneven recent survey effort into false species absences. Record dates, data-source composition, basis of record, coordinate uncertainty, and cultivated/native evidence are retained for a later observation-process audit.

## GBIF request lifecycle

```text
public GSHHG build
    -> islands_v2.gpkg + island_manifest.csv + source policy
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
  --islands-gpkg data/v2/external/islands/gshhg/prepared/islands_v2.gpkg \
  --output-csv data/v2/external/gbif/request_manifest.csv \
  --checklist-key 7ddf754f-d193-4cc9-b351-99906754a03b

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
