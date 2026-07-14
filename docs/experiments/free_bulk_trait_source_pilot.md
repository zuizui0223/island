# Free bulk trait source pilot

## Goal

Before spending retrieval budget on species-by-species search, exhaust open bulk trait
sources that can replace low-resolution genus/family/global cascade fills with direct
or curated species-level evidence.

The first live pilot is AusTraits, focusing on the `eFLOWER_Dun_2022` component.
The source metadata describes roughly 2,000 taxa, 36 floral traits, and more than
29,000 trait values derived from taxonomic literature.

## Live pipeline

The branch now contains a reproducible end-to-end experiment:

1. Query the same public Zenodo concept-record versions endpoint used by the official
   AusTraits R client.
2. Select the newest release by publication date.
3. Prefer a plain-text ZIP archive rather than an R-only serialized object.
4. Discover the long-format trait CSV from its column contract rather than from a
   hard-coded archive path.
5. Filter `eFLOWER_Dun_2022` and standardize only the transport columns required for
   auditing.
6. Preserve source values and source units exactly; do not infer categories or units.
7. Exact-match normalized scientific names to the current island species master.
8. Report the number of unique matched `species x mapped target trait` cells. This is
   the direct estimate of candidate fallback cells that a promoted source could replace.

The workflow is `.github/workflows/run-free-bulk-trait-pilot.yml`. Its artifact contains:

- `fetch_manifest.json`
- `coverage_report.json`
- the standardized `eFLOWER_Dun_2022.csv`

## Current mapping policy

The pilot intentionally maps only a small set of source traits while exposing the top
unmapped source traits for the next mapping pass:

- `flower_length` -> `flower_length_continuous_raw`
- `flower_diameter` -> `flower_diameter_continuous_raw`
- `flower_perianth_fusion` -> `flower_perianth_fusion_continuous_raw`
- `flower_structural_sex_type` -> `flower_structural_sex_type_raw`

The raw continuous names are deliberately unit-neutral. The standardized transport
file keeps `trait_unit`; unit normalization and any derived categorical trait must be
a separate reviewed transformation.

## Guardrails

- no fuzzy taxon matching;
- no silent synonym guessing;
- no continuous-to-category conversion during acquisition;
- no unit assignment when the source does not provide one;
- repeated observations do not inflate fallback-replacement potential: the promotion
  metric counts unique exact-master-matched `species x target_trait` cells;
- source rows remain traceable through `source_record_id`, dataset key, release DOI,
  release version, and fetch manifest.

## Promotion criteria

A source should move from `pilot` to the production bulk cascade only after the live
report establishes:

1. a non-trivial exact match to the island species master;
2. useful unique species-trait replacement coverage rather than merely many duplicate
   observations;
3. an explicit source codebook and reviewed mapping;
4. preserved provenance and units;
5. focused tests passing on the production adapter.

After promotion, the source should write the existing `candidate_long` contract under
`data/v2/staging/traits/bulk/**/trait_candidates.csv*`, then the current fill cascade can
consume it without changing the cascade algorithm. Direct evidence should replace
fallback values; fallback remains only for unresolved cells.

## Production promotion (2026-07-14)

The validated eFLOWER component is now promoted through
`.github/workflows/ingest-austraits-eflower-bulk.yml`. The production path verifies the
Zenodo-declared archive MD5, records a SHA-256, maps only directly corresponding APD
categories, writes review-pending source-backed candidates, and rebuilds the exhaustive
Island acquisition ledger.

The live v7.0.0 run against the 115,328-species master produced:

- 30,650 eFLOWER source rows and 21,114 master-matched rows;
- 2,589 source-backed candidate rows across 1,017 species;
- 938 species with `flower_primary_color` evidence;
- 594 species with `floral_symmetry` evidence;
- zero empty source URLs, evidence excerpts, or source record IDs;
- exactly 115,328 species rows and 1,499,264 species-trait rows in the resulting ledger.

Continuous flower measurements, `disymmetric` perianths, and flower-level structural-sex
states remain explicit unmapped audit rows. They are not coerced to flower-size classes,
the narrower symmetry ontology, or species-level `sex_system` values.

## Full-release target expansion

The eFLOWER-only production lane is complemented by
`.github/workflows/ingest-austraits-all-floral-bulk.yml`. This lane scans all 1,798,215
v7.0.0 records but retains only six inspected categorical source traits: flower colour,
perianth colour, perianth symmetry, whole-plant sex type, pollination syndrome, and
pollination system. The original AusTraits component dataset remains embedded in every
source record ID.

The full local production-equivalent join produced 31,647 candidates across 9,464 species:

- `sex_system`: 8,334 species;
- `pollination_functional_guild`: 3,953 species;
- `flower_primary_color`: 3,881 species;
- `pollen_vector_mode`: 3,577 species;
- `floral_symmetry`: 741 species.

All 65 observed `flower_colour` labels have explicit profile mappings. Unsupported broad
abiotic/self combinations, `disymmetric`, unknown sex type, rare pollination, and
unidentified specialised-biotic codes remain in the unmapped audit. No missing source
code is guessed.
