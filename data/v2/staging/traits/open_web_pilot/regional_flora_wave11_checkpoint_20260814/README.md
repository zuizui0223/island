# Regional flora wave 11 checkpoint (2026-08-14)

This checkpoint changes the acquisition unit from repeated global batches to
complete, structured regional-flora surfaces.  It adds only reviewed
species-direct rows.  Genus inference is not stored here and is rebuilt later
by the shared `genus x trait_name` implementation.

## Acquired sources

- **Endemia New Caledonia**: all species-rank results from nine original
  structured flower-colour filters (including `Insignificant flower`, excluding
  `No flower`).  The acquisition covered 67 result pages, 629 colour
  assignments and 582 listed species.  Exact fixed-master name and family
  checks, existing-ledger deduplication and one-lineage-per-species treatment
  produced 364 new Medium rows.  Five exact-name records with a family conflict
  between Endemia and the fixed master were rejected.
- **Flora of Zimbabwe**: 251 family indexes (6,412 listings) were inventoried;
  493 relevant species pages were fetched and parsed with the existing
  measurement ontology.  Exact species treatments yielded 112 new Medium rows
  for 82 species.
- **New Zealand Plant Conservation Network (NZPCN)**: all 14 non-empty native
  flower-colour filters were enumerated (137 result pages, 2,564 colour
  assignments, 1,614 taxa).  After current-ledger deduplication, 283 original
  species pages were fetched.  All 283 passed exact page-name, family, Native
  biostatus and non-empty `Flower colours` field gates and are retained as High
  database evidence.

The committed wave contains 759 reviewed source rows for 729 species.  Every
row stores the original URL, an exact supporting excerpt or structured field,
retrieval timestamp, content hash, source lineage and cultivar status.  The
line-by-line audit has 759 accepted-correct rows out of 759 reviewed rows and
zero cultivar-contaminated rows.  This precision describes the reproducible
identity/value extraction audit; it is not a claim that a flora database has no
biological or taxonomic error.

## Guardrails

- fixed island master species and family must match exactly;
- only species-direct High or Medium evidence is emitted;
- multicolour states are retained as a state set;
- `pollen_vector_mode` and `reward_type` are not mapped into the structural
  axis;
- family inference and global fallback are absent;
- Validated Low is delegated to the common all-evidence rebuild and is joined
  by `genus x trait_name`, never `genus x axis`;
- coverage gains are not claimed from this checkpoint alone and must be read
  from the final workflow artifact.

## Rebuild

```bash
python -m island_v2.regional_flora_wave11_checkpoint \
  --source-snapshot-csv data/v2/staging/traits/open_web_pilot/regional_flora_wave11_checkpoint_20260814/regional_flora_wave11_source_rows_20260814.csv \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/high_yield_reproductive_wave10_checkpoint_20260814/combined_curated_evidence_20260814.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/high_yield_reproductive_wave10_checkpoint_20260814/combined_curated_manual_audit_20260814.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/regional_flora_wave11_checkpoint_20260814
```

`regional_flora_wave11_manifest_20260814.json` records hashes for the fixed
master, prior checkpoint, reviewed source snapshot and generated outputs.
