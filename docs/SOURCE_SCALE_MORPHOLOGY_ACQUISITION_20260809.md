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
| Kew/WFO African floras (FTEA, FZ, FWTA) | 37,779 matched descriptions | 3,720 | 3,720 | 1,266 approved | 1,248 approved |
| WFO global six-flora wave | 22,981 matched descriptions | 8,631 | 8,631 | 4,742 approved | 3,788 approved |

PlantNET novelty is measured after excluding both the pre-wave formal ledger
and Flora of Australia. It therefore represents the second source's pure
candidate increment. The formal coverage workflow still recalculates direct
conflicts and every trait-specific genus rule before reporting net coverage.

All five source groups were acquired without search API credentials or query
cost. The complete public collection indexes are read once, and each
exact-match profile or pinned archive is processed deterministically. A restart
does not fetch a completed profile again.
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

The fourth source-scale wave reads three official Kew/WFO Darwin Core Archives
in one pass: Flora of Tropical East Africa, Flora Zambesiaca, and Flora of West
Tropical Africa. Names are reconciled through the pinned June 2026 WFO Plant
List using exact accepted names, exact fixed-universe synonyms, and their
family-consistent WFO synonym clusters. Of 1,491 novel candidate rows, an
independent deterministic holdout reviewed 200 source passages and accepted
197 (precision 0.985; cultivar contamination 0). The trait-level production
gate approves form (26/26), colour (98/99), size (13/13), and inflorescence
display (57/59). Symmetry (1 row) and tube depth (2 rows) remain stored but are
not scaled because each has fewer than ten reviewed records. After excluding
the three reviewed failures, 1,485 evidence rows remain, covering 1,266 novel
species x trait pairs and 1,248 species x axis pairs before formal conflict and
Validated-Low recalculation.

The fifth wave reads six additional pinned WFO-hosted archives in one local
pass: Brazilian Flora 2020, e-Flora Thailand, Flora of Panama, Flora
Neotropica, Memoirs of the New York Botanical Garden, and the Edinburgh
Rhododendron monograph. The prior formal direct ledger is used as an exclusion
set before extraction. The deterministic 200-row holdout accepted 187 rows
(0.935 overall; cultivar contamination 0). Approval is deliberately
trait-specific at the unchanged 0.95 threshold: form (30/30), symmetry
(10/10), flower size (30/30), and inflorescence display (40/40) scale to
production, while colour (64/70) and tube depth (13/20) remain preserved but
unpromoted. This yields 5,434 selected evidence rows, 4,742 species x trait,
and 3,788 species x axis candidates before conflict and Validated-Low
recalculation.

The sixth high-yield wave adds species treatments from the official Flora
Malesiana and Flora of the Guianas FlorML archives, four further WFO-hosted
Darwin Core archives, and a strictly re-gated subset of colour evidence that
the fifth wave preserved but did not promote. It excludes every
`accepted_species x trait_name` already present in formal run `31350993446`,
deduplicates equal provider-treatment lineages, and rejects colour statements
about fruits, indumentum, immature material, or dried material. Display counts
for sex-specific or partial inflorescence units are also excluded. The frozen
200-row review accepted 189 passages overall. The unchanged trait-specific
0.95 gate approves floral form (29/30), flower colour (59/60), and
inflorescence display (70/70), while flower size (31/40) remains stored but is
not promoted. The approved scopes contain 5,467 evidence rows, representing
5,231 species x trait and 4,731 species x axis candidates before formal
conflict resolution and all-evidence Validated-Low reconstruction.

The seventh wave extends that package with three more official WFO-hosted
floras: Flora of North America, Flora of China, and Flore d'Afrique Centrale.
Exact WFO accepted/synonym reconciliation found 14,819 fixed-universe species
treatments. A quote-completeness gate removes any display row that drops an
explicit second state (for example, capitula arranged in corymbs) and reruns
the deterministic provider x trait audit strata. The new 120-row review
accepted 116 passages (0.9667; cultivar contamination 0): form 20/20, colour
19/20, display 60/60, and size 17/20. Combined with the sixth-wave holdout,
the production gate is based on 320 reviews and approves form (49/50), colour
(78/80), and display (130/130); size remains blocked (48/60). The approved
combined package contains 9,337 evidence rows representing 9,049 species x
trait and 8,484 species x axis candidates before formal conflict resolution
and Validated-Low reconstruction. Plants of Nepal was inspected but excluded
because its scalable treatment fields did not contain target floral
descriptions.

A deterministic local replay against formal Web run `31350993446` resolved
8,945 new direct species x trait cells and 6,858 direct species x axis cells.
After rebuilding all trait-specific genus rules, strict coverage changed from
129,647 to 132,857 species x axis cells: 6,370 gross fills and 3,160 losses,
for a net gain of 3,210. The losses are retained scientific corrections, not
discarded acquisition failures: the new direct counterexamples invalidated
5,881 prior Low rows, added 3,462 new Low rows, and upgraded 2,090 Low cells to
direct evidence. Current-threshold eligible genus x trait rules increased from
2,698 to 2,837 (+139). Consequently, subsequent acquisition prioritizes
currently unresolved and all-axis-missing species; it does not suppress
counterexamples merely to maximize the headline coverage count.

## Reproduction

The raw HTTP caches are intentionally not committed. Given the prior formal
artifact files, the source-scale acquisitions are run as modules:

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

py -3.13 -m analysis.acquire_wfo_kew_africa_traits_20260810 `
  --backbone-zip <wfo-plant-list-2026-06.zip> `
  --archive-dir <three-kew-wfo-dwca-directory> `
  --strict-coverage-csv <sanbi-strict-species-axis-coverage.csv.gz> `
  --current-direct-ledger-csv <sanbi-direct-species-trait-ledger.csv.gz> `
  --output-dir <wfo-kew-africa-output>

py -3.13 -m analysis.acquire_wfo_global_flora_traits_20260810 `
  --backbone-zip <wfo-plant-list-2026-06.zip> `
  --archive-dir <six-wfo-provider-archives-directory> `
  --strict-coverage-csv <wfo-kew-africa-strict-species-axis-coverage.csv.gz> `
  --current-direct-ledger-csv <wfo-kew-africa-direct-species-trait-ledger.csv.gz> `
  --output-dir <wfo-global-six-output>

py -3.13 analysis/prepare_wfo_high_yield_source_package_20260810.py `
  --florml-evidence <malesiana-guianas-evidence.csv.gz> `
  --regional-evidence <four-additional-wfo-floras-evidence.csv.gz> `
  --prior-global-evidence <wfo-global-six-evidence.csv.gz> `
  --current-direct-ledger <run-31350993446-direct-ledger.csv.gz> `
  --audit-seed <frozen-high-yield-audit-seed.csv> `
  --output-dir <wfo-high-yield-output>

py -3.13 analysis/prepare_wfo_fna_foc_source_package_20260810.py `
  --extracted-evidence <fna-foc-africa-extracted-evidence.csv.gz> `
  --audit-seed <frozen-provider-trait-audit-seed.csv> `
  --manual-audit <reviewed-provider-trait-audit.csv> `
  --first-evidence <wfo-high-yield-output/evidence.csv.gz> `
  --first-audit <wfo-high-yield-output/manual-audit.csv> `
  --output-dir <wfo-combined-output>
```

Committed evidence, treatment inventories, deterministic reviews, summaries,
and SHA-256 manifests are under:

- `data/v2/staging/traits/direct_llm_pilot/20260809_flora_of_australia_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260809_plantnet_nsw_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260809_sanbi_eflora_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260810_wfo_kew_africa_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260810_wfo_global_six_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260810_wfo_high_yield_source_acquisition/`
- `data/v2/staging/traits/direct_llm_pilot/20260810_wfo_fna_foc_africa_source_acquisition/`

Formal promotion must run Flora of Australia, PlantNET, SANBI, then WFO/Kew
Africa, then the WFO global six-flora wave against the resulting artifact so
later runs cannot rescan or double-count earlier evidence.
The combined high-yield wave extends formal run `31350993446`; its reviewed
package is then passed through the same shared finalizer and trait-specific
genus-rule implementation rather than an independent Low implementation.
