# Source-scale morphology acquisition (2026-08-09)

This wave replaces species-by-species Web search with complete official-flora
indexes. It preserves the fixed 106,295-species universe and writes only
species-direct evidence. Family inference, global fallback, and axis-only
extraction are not used.

## Sources and realized acquisition

| Source | Collection species profiles | Exact fixed-universe profiles | Profiles with descriptions | Novel species x trait | Novel species x axis |
|---|---:|---:|---:|---:|---:|
| Flora of Australia | 19,207 | 6,794 | 6,329 | 4,060 | 3,276 |
| PlantNET NSW FloraOnline | 7,505 | 4,236 | 4,052 | 758 | 699 |
| SANBI e-Flora of South Africa v1.42 | 21,586 taxon records | 2,811 | 2,811 | 2,548 | 1,947 |

PlantNET novelty is measured after excluding both the pre-wave formal ledger
and Flora of Australia. It therefore represents the second source's pure
candidate increment. The formal coverage workflow still recalculates direct
conflicts and every trait-specific genus rule before reporting net coverage.

All three sources were acquired without search API credentials or query cost. The
complete public collection index is read once, and each exact-match profile is
cached independently. A restart does not fetch a completed profile again.
SANBI is a versioned, CC BY 4.0 Darwin Core Archive rather than a profile API;
the acquisition pins the 15,561,518-byte archive at SHA-256
`157fe16ca6a5698df24fca4796aec2d69a8da20fd38a3bd36633eb0d5cdedfe4`.

## Evidence contract

The acquisition extracts these six traits concurrently from the original
species Description field:

- `flower_primary_color`
- `floral_form`
- `floral_symmetry`
- `flower_size_class`
- `tube_depth_class`
- `inflorescence_display`

Every selected row retains accepted species, exact family gate, source URL,
source citation, exact supporting excerpt, retrieval package, and source
lineage. Negation and measurement-owner gates prevent statements about bracts,
hairs, leaves, flower heads, hypanthia, capsules, spathes, opercula, filaments,
or comparison objects from being transferred to the flower.

Each source has a deterministic, trait-stratified 200-row manual audit.
Precision is `accepted_correct / all_reviewed`. Flora of Australia and PlantNET
are 200/200; SANBI is 198/200 (0.99 overall), all with zero cultivar
contamination. SANBI's flower-size stratum is 43/45 (0.9556); all other SANBI
strata are 1.0. The two rejected SANBI records are preserved in the review:
one component-organ measurement and one implausible source value. Flora of
Australia and SANBI meet the production sample minimum for all six traits.
PlantNET meets it for flower size, inflorescence display, and tube depth; its
13 form/symmetry/colour rows remain recorded but are not production-approved
in the separate PlantNET promotion.

## Reproduction

The raw HTTP caches are intentionally not committed. Given the prior formal
artifact files, the two resumable acquisitions are run as modules:

```powershell
py -3.13 -m analysis.acquire_flora_of_australia_traits_20260809 `
  --baseline-evidence-csv <prior-reviewed-source-lineages.csv.gz> `
  --universe-coverage-csv <prior-strict-species-axis-coverage.csv.gz> `
  --cache-dir <foa-cache> --output-dir <foa-output>

py -3.13 -m analysis.acquire_plantnet_nsw_traits_20260809 `
  --baseline-evidence-csv <prior-reviewed-source-lineages.csv.gz> `
  --baseline-evidence-csv <foa-output/flora_of_australia_evidence.csv.gz> `
  --universe-coverage-csv <prior-strict-species-axis-coverage.csv.gz> `
  --cache-dir <plantnet-cache> --output-dir <plantnet-output>

py -3.13 -m analysis.acquire_sanbi_eflora_traits_20260809 `
  --archive-path <flora_descriptions-v1.42.zip> `
  --dwca-dir <sanbi-dwca-cache> `
  --current-direct-ledger-csv <plantnet-direct-species-trait-ledger.csv.gz> `
  --strict-coverage-csv <plantnet-strict-species-axis-coverage.csv.gz> `
  --output-dir <sanbi-output>
```

Committed evidence, treatment inventories, deterministic reviews, summaries,
and SHA-256 manifests are under:

- `data/v2/staging/traits/direct_llm_pilot/20260809_flora_of_australia_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260809_plantnet_nsw_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260809_sanbi_eflora_source_acquisition/`

Formal promotion must run Flora of Australia, PlantNET, then SANBI against the
resulting artifact so later runs cannot rescan or double-count earlier evidence.
