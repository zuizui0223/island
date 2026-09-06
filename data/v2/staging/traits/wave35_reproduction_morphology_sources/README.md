# Wave35 reproduction and morphology sources

This directory freezes the reviewed evidence used by the additive Wave35
overlay.  The fixed denominator is 106,295 accepted island species and three
axes (318,885 species-axis cells).  Evidence is resolved at
`accepted_species x trait_name`, source lineages are deduplicated, and quality
precedence is direct High, direct Medium, then validated genus Low.

## New source lanes

- **FloraWeb / BIOLFLOR**: the public 2026-06 scrape snapshot is fixed by URL
  and SHA-256 in `floraweb_trait_source_summary.json`.  Controlled BIOLFLOR
  reproductive values are High; flower-bound descriptive colour is Medium.
- **GloPL 2025**: Supplementary Data 7 of DOI
  `10.1038/s41467-025-61032-5` is fixed by workbook SHA-256.  Only autonomous
  selfing, self-incompatibility and floral symmetry enter the strict axes.
- **Ken Fern Useful Temperate/Tropical Plants and PFAF**: 6,779 public species
  pages were fetched with zero errors.  A positive statement for a wild,
  sexually reproducing species is mapped only to `self_incompatibility=SC` at
  Medium quality.  It is never mapped to autonomous selfing.  Provider
  redistributions share one Ken Fern species lineage.

`pollen_vector_mode`, `reward_type`, dichogamy, herkogamy, heterostyly and
other useful independent traits are retained in the independent ledgers but
are not substituted into the strict flower-structure or reproduction axes.

## Overlay inputs

- `wave35_resolved_direct_species_trait.csv.gz` contains conflict-resolved
  direct rows in which at least one Wave35 provider participated.
- `wave35_provider_touched_new_rule_audit.csv.gz` contains eligible
  `genus x axis x trait_name x inferred_state_set` rules absent from the frozen
  Wave34 frontier and touched by new direct evidence.
- `wave35_candidate_validated_low_species_trait.csv.gz` contains only
  still-unresolved Wave34 cells reached through those exact trait rules.
- `wave35_old_low_comparison.csv.gz` records the separate strict-rebuild audit.
  Invalidated old Low values are reported but are not silently subtracted from
  the frozen lossless secondary ledger.
- `wave35_source_manifest.json` fixes every source, formal GitHub artifact and
  file hash needed to audit the build.

Build the formal overlay after downloading the pinned Wave34 artifact:

```text
island-v2-wave35-trait-overlay build \
  --baseline-csv <wave34_secondary_species_axis_coverage.csv.gz> \
  --input-dir data/v2/staging/traits/wave35_reproduction_morphology_sources \
  --output-dir wave35-trait-coverage
```

The command rejects duplicate cells, unresolved direct conflicts, cross-trait
axis joins, support below three direct species, family inference, global
fallback, quality downgrades and any loss of an already-filled Wave34 cell.
