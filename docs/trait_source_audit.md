# Initial public trait-source audit

This is an inventory, not an approval to ingest any source. A source is approved
only after its downloadable table, trait definitions, value codebook, licence,
and name-matching behaviour have been checked against the v2 ontology.

## Current decision

EOL TraitBank, AusTraits, and the public GIFT API now have production acquisition
lanes. Approval is limited to pending candidate generation: no imported row is
curated automatically, and unmatched, excluded, or unmapped records remain in
audits until reviewed. GIFT is pinned to stable database version 3.2 and uses raw
records only (`derived=0`, `biasderiv=0`); restricted or absent public references
are not fetched.

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
| GIFT 3.2 | Global flora and trait archive | The official unauthenticated, versioned API exposes trait metadata, reference metadata, standardized species names, and reference-level raw trait records. The production lane selects direct flower colour, sexual system, and pollination-syndrome fields, filters to public unbiased reference/trait pairs, and archives raw responses plus primary citations. The GIFT disclaimer requires version and primary-source citation and warns about taxonomic harmonization and coverage gaps; this audit does not assert a separate permissive redistribution licence. | Approved for source-backed pending candidates from the public API only; no restricted API content and no derived values. | Re-run `Acquire GIFT direct floral traits`, inspect the exclusion/name-match audits, and cite GIFT 3.2 plus the retained `ref_long` source in downstream use. |
| AusTraits | Australia | AusTraits is an open, CC-BY 4.0 harmonized database with more than 500 released trait names across more than 30,000 Australian taxa, downloadable through Zenodo and an R interface. | Production regional enrichment for direct flower colour, floral symmetry, whole-plant sex system, pollen vector, and pollination functional guild; not a global coverage source. | Continue checksum-verified selected-trait ingests and preserve all non-corresponding values in the unmapped audit. |
| Meyer, Galloway & Eckert 2026 Dryad v4 | Global literature synthesis, 4,324 species rows | The CC0 output defines `sc`/`si` as largely self-compatible/self-incompatible. The pinned file hash is validated. Four rows derived only from a DIO (dioecy) code are excluded because dioecy does not demonstrate genetic SI. | Approved for pending species-direct `self_incompatibility` candidates after exact accepted-name or strict GBIF exact-synonym reconciliation. Unmatched/fuzzy/family-conflicting names remain audited. | Keep direct candidates separate from the low-confidence genus completion layer. Genus transfer requires at least five directly supported congeners with unanimous SC/SI and is emitted as `likely_SC`/`likely_SI`, never as direct evidence. |
| Razanajatovo et al. 2016 Supplementary Data 1 | Global experimental synthesis, 1,890 observations / 1,752 species | The CC-BY workbook reports continuous self-compatibility and autofertility indices with underlying reference keys. The workbook and extracted worksheet CSV are checksum pinned. | Approved for conservative pending candidates: consistently high compatibility supports SC, while low compatibility is never converted to SI; autofertility maps only to `autonomous_selfing_capacity`. | Retain intermediate/heterogeneous index sets as `mixed_or_variable`, preserve raw index values and reference keys, and keep exact/synonym name-match audits. |

## Decision rule

No source is allowed into `ingest_bulk_trait_source` until it has:

1. a stable download or reproducible export route;
2. a retained release citation and licence record;
3. a source-specific profile in `config/bulk_trait_sources.yml`;
4. a small-sample audit of exact-name, author-stripped, unmatched, and ambiguous matches; and
5. a documented mapping from source values to the v2 ontology.

The first operational test should use 1,000 island-master species and a limited
set of high-priority global traits. Reproductive-assurance traits remain an
evidence-rich subset unless direct coverage is demonstrated. Pollinator guild is
a deferred interaction-evidence layer, not an active floral-trait completion
target.
