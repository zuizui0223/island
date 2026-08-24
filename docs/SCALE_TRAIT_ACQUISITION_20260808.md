# Scale-first trait acquisition (2026-08-08)

This batch replaces per-species general-Web querying with source-scale
acquisition from structured floras.  It is incremental relative to reviewed
artifact run `31257281120`: a species x trait already represented by direct
High/Medium evidence is not emitted again.

## Sources and lineage policy

| Package | Original source | Strict traits | Evidence tier |
|---|---|---|---|
| Kew GrassBase | Royal Botanic Gardens, Kew species descriptions | `inflorescence_display` | High |
| NHM ADEPT: Flora of China | Flora of China, redistributed by NHM ADEPT | colour, form, inflorescence, tube depth | High |
| NHM ADEPT: Flora of North America | Flora of North America, redistributed by NHM ADEPT | colour, form, inflorescence, tube depth | High |
| NHM ADEPT: Flora of Pakistan | Flora of Pakistan, redistributed by NHM ADEPT | colour, form, symmetry, inflorescence, tube depth | High |
| Baseflor 2023.10 | Philippe Julve / Tela Botanica | colour, inflorescence; explicit `autogame` only | High except autogamy mapping Medium |
| Western Australian FloraBase | Western Australian Herbarium taxon profiles | flower colour | High |
| FloraWeb / BiolFlor snapshot | German BfN flora portal; fixed taxifydb snapshot | SI, mating system, explicit self-pollination/cleistogamy, floral form | Medium |
| Plants For A Future | Exact public species pages | explicit self-fertility, mapped only to SC | Medium |
| Useful Tropical / Temperate Plants | Ken Fern structured species properties | explicit wild-species self-fertility, mapped only to SC | Medium |
| Tree of Sex | Dryad published database and original `source_selfing` citations | SC, SI, variable | High |
| Alien Plants in Greece | Taxon-level Autogamy database records | discovery-only self-pollination | not admitted to strict axes |
| Direct reproductive studies | Original articles or exact quotations in official species assessments | autonomous selfing, SC | High |

NHM is a redistribution provider, not an independent original source.  When a
Flora of China or Flora of North America treatment is already present in the
ledger, the ADEPT row reuses that exact eFloras treatment lineage.  Otherwise
the lineage is the named original flora plus WFO taxon identifier.  Baseflor
duplicates are one lineage per species, so repeated taxonomic records cannot
be mistaken for independent agreement.

Plants For A Future and both Ken Fern databases are one Ken Fern lineage per
species. Tree of Sex rows use the underlying `source_selfing` citation as
lineage, not the Dryad or taxifydb redistribution. Purely dioecious Tree of
Sex records are excluded rather than translating sexual-system separation
into SI. PFAF statements mentioning apomixis or dioecy remain in the
acquisition inventory as category conflicts and are not evidence.

`pollen_vector_mode`, `reward_type`, and `sex_system` are not mapped to the
strict floral-structure or reproductive-assurance axes.  Baseflor
`entomogame`, `anemogame`, `hydrogame`, and `autogame` values remain outside the
strict three-axis ledger. They describe pollen transfer and do not by
themselves demonstrate unassisted fertilisation, fruit set, or seed set.
FloraWeb/BiolFlor `Selbstbestäubung` and Alien Plants in Greece `Autogamy
(self-pollination)` records follow the same fail-closed rule. They may be kept
for discovery or a future pollen-transfer ontology, but cannot populate
`autonomous_selfing_capacity`.

## Reproduction

The raw inputs are pinned by SHA-256 in
`scale_source_package_summary.json`.  The acquisition adapters cache every
download and are resumable.  Re-running the package builder creates an audit
template but never creates review decisions.

```bash
python analysis/acquire_kew_grassbase_traits_20260808.py \
  --strict strict_species_axis_coverage.csv.gz \
  --master data/v2/staging/gbif/collected/island_taxa.csv \
  --output grassbase

python analysis/acquire_florabase_colours_20260808.py \
  --strict strict_species_axis_coverage.csv.gz \
  --master data/v2/staging/gbif/collected/island_taxa.csv \
  --output florabase

python analysis/acquire_pfaf_reproductive_traits_20260808.py \
  --queue-csv pfaf_unresolved_reproduction.csv \
  --sitemap-xml pfaf-sitemap.xml \
  --robots-txt pfaf-robots.txt \
  --cache-dir pfaf-cache \
  --output-dir pfaf

python analysis/acquire_ken_fern_reproductive_traits_20260808.py \
  --queue-csv ken_fern_unresolved_reproduction.csv \
  --cache-dir ken-fern-cache \
  --output-dir ken-fern

python analysis/inventory_nhm_adept_traits_2026.py \
  --resource flora_of_china \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --output-dir nhm

python -m analysis.prepare_scale_trait_source_packages_20260808 \
  --grassbase-csv grassbase/fetch_inventory.csv.gz \
  --nhm-root nhm \
  --baseflor-xlsx baseflor-original-202310.xlsx \
  --florabase-csv florabase/fetch_inventory.csv.gz \
  --floraweb-csv floraweb_raw_2026-06-24.csv \
  --pfaf-csv pfaf/pfaf_self_fertile_candidates.csv.gz \
  --ken-fern-csv ken-fern/ken_fern_self_fertile_candidates.csv.gz \
  --tree-of-sex-csv tree_of_sex_taxifydb_2026-08.csv \
  --alien-plants-autogamy-html alien-plants-autogamy.html \
  --direct-literature-csv direct_reproductive_literature_manifest.csv \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --current-coverage-csv strict_species_axis_coverage.csv.gz \
  --current-lineages-csv reviewed_source_lineages.csv.gz \
  --output-dir scale-source-package
```

Repeat the NHM inventory command for `flora_of_north_america` and
`flora_of_pakistan`.

The committed audit is a deterministic stratified sample (`random_state =
20260808`) of 200 records from each large provider and every record from each
provider with at most 200 candidates, 1,622 total. Precision is always
calculated as `accepted_correct / all_reviewed`.  Production scaling still
passes through the common source-package gate: at least 200 reviewed records,
at least 10 per approved trait, precision at least 0.95, cultivar contamination
at most 0.02, exact species/synonym identity, complete quote and URL, and an
ontology-valid value.

The committed incremental package contains 13,227 evidence rows representing
12,622 previously absent species x trait pairs across 9,037 strict-universe
species.  The strict-universe input is mandatory and validated as exactly
106,295 species; two Baseflor gymnosperms present in the broader island master
were therefore excluded before audit and formalization.

The package also preserves 3,655 FloraWeb records for `pollen_vector_mode` and
`reward_type` in `scale_source_non_axis_candidates.csv.gz`.  These remain
trait-level data with a blank strict axis and cannot inflate the three-axis
coverage calculation.

## Formalization

Use `.github/workflows/open-web-review-promote.yml` with the committed
`scale_source_incremental_evidence.csv.gz` and
`scale_source_manual_audit.csv`.  The workflow appends the evidence to the
previous formal public-Web ledger, deduplicates source lineages, resolves
direct conflicts, rebuilds trait-specific Validated Low, and recalculates the
fixed 106,295 species x 3 axes exactly once.
