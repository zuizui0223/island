# Initial public trait-source audit

This is an inventory, not an approval to ingest any source. A source is approved
only after its downloadable table, trait definitions, value codebook, licence,
and name-matching behaviour have been checked against the v2 ontology.

## Current decision

EOL TraitBank is the first approved bulk-ingestion pilot because it has a stable
public all-traits archive and a dedicated conservative mapping file in
`config/eol_traitbank_terms.yml`. Approval is limited to pending candidate
generation: no EOL row is curated automatically, and unmapped terms remain in
the audit until reviewed.

In particular, BIEN is a valuable free validation and name-resolution resource,
but it is **not** the first floral-trait import source: its published high-volume
trait coverage is dominated by DBH, whole-plant height, growth form, seed mass,
wood density, and specific leaf area. Those do not establish coverage for flower
colour, floral architecture, pollen-vector mode, or reproductive assurance.

| Source | Scope | What is currently established | Decision | Required next action |
|---|---|---|---|---|
| EOL TraitBank all-traits archive | Global, all life | Zenodo record 13305577 provides `traits_all.zip`, an EOL graph export containing `pages.csv`, `traits.csv`, `metadata.csv`, `inferred.csv`, and `terms.csv`. | Approved as a free bulk-ingestion pilot for explicitly mapped predicate/value terms only. | Run `ingest_eol_traitbank_bulk` first with a smoke-test row cap, inspect `eol_traitbank_unmapped_terms.csv` and name-match audits, then expand `config/eol_traitbank_terms.yml` only from reviewed terms. |
| TRY Plant Trait Database | Global | The official portal reports free/open access, version 6 with 15,409,681 records across 305,594 taxa, and completion of the version 7 import in July 2025. Access is through the data portal rather than a fixed public bulk-download URL. | Screening only; potentially valuable if a requested output includes the required floral or reproductive fields. | Obtain the trait catalogue and a small requested export. Preserve release, terms, trait IDs, units, and source citation before writing any profile. |
| BIEN 4.2 | Global, strongest New World coverage | BIEN reports 25.9M trait observations for 54 standardized traits, but its described high-volume traits are DBH, height, growth form, seed mass, wood density, and SLA. Coverage is strongest for New World vascular plants. | Do not use as first floral-trait import. Retain for name validation, geographic checks, and possible non-floral covariates or sensitivity work. | Do not create a v2 floral profile unless a trait dictionary demonstrates direct ontology-mappable floral or pollen-vector fields. |
| GIFT | Global flora and trait archive | TRY reports trait information for 281,836 species and 109 traits, but the current reproducible download route, licence, and trait dictionary remain unverified here. | Screening only. | Verify current export route, licence, trait dictionary, and per-field coverage before selecting any value. |
| AusTraits | Australia | AusTraits is an open, CC-BY 4.0 harmonized database with nearly 500 traits across more than 30,000 Australian taxa, downloadable through Zenodo and an R interface. | Production regional enrichment for direct flower-colour and floral-symmetry evidence; not a global coverage source. | Continue ingesting checksum-verified releases and preserve all non-corresponding traits in the unmapped audit. |

## Decision rule

No source is allowed into `ingest_bulk_trait_source` until it has:

1. a stable download or reproducible export route;
2. a retained release citation and licence record;
3. a source-specific profile in `config/bulk_trait_sources.yml`;
4. a small-sample audit of exact-name, author-stripped, unmatched, and ambiguous matches; and
5. a documented mapping from source values to the v2 ontology.

The first operational test should use 1,000 island-master species and a limited
set of high-priority global traits. Reproductive-assurance and pollinator-guild
traits remain evidence-rich subset variables unless species-direct data coverage
is demonstrated.
