# v2 data layout

v2 begins with new, traceable data products and does not use v1 derived tables as analysis inputs.

```text
templates/
  taxa_template.csv                 accepted taxon batch for web-search trait work
  accepted_evidence_template.csv    one accepted human-reviewed trait record per row
  review_decisions_template.csv     explicit accept/reject/adjudication decisions

external/
  islands/                          source island polygons and normalized island manifest
  gbif/                             exact-polygon request manifests and downloaded archives

staging/
  gbif/                             raw island-by-name candidates before taxonomy review
  ...                               raw LLM candidate exports and review queues

curated/                            versioned, human-adjudicated trait evidence and decisions
```

## First v2 data product: GBIF island species-name candidates

1. Put `GlbIslands.gdb` / `BigIslands` or a documented replacement vector source in `external/islands/`.
2. Run `island-v2-gbif-flora prepare-islands` to create `islands_v2.gpkg` and `island_manifest.csv`.
3. Run `make-requests`, `submit`, `poll`, and `collect` in bounded batches.
4. Treat the output as occurrence-based candidate names, not a verified flora.
5. Normalize taxonomy, then audit establishment status and the observation process before LLM trait collection.

The full policy and commands are in `docs/gbif_flora_v2.md`.

## LLM trait workflow follows species-name review

After the candidate species list is normalized, run the LLM web-search workflow in batches. It returns source-cited candidates and never directly edits curated trait data. Genus/family distributions are generated only from reviewed evidence.
