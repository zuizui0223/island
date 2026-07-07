# Initial public trait-source audit

This is an inventory, not an approval to ingest any source. A source is approved
only after its downloadable table, trait definitions, value codebook, licence,
and name-matching behaviour have been checked against the v2 ontology.

| Source | Scope | What is currently established | Use status | Next step |
|---|---|---|---|---|
| TRY Plant Trait Database | Global | The official portal reports free/open access, v6 with 15,409,681 records across 305,594 taxa, and completion of the v7 import in July 2025. Data are requested through the portal rather than a stable direct bulk URL. | Screening only | Export the trait catalogue and a small requested sample. Map only documented values to the v2 ontology. |
| BIEN 4.2 | Global, strongest New World coverage | BIEN reports 25.9M trait observations and 54 standardized traits. It provides a data portal, R package, and validation services. | Screening only | Inspect the trait dictionary and export format; run name-overlap and codebook tests on 1,000 island-master species. |
| GIFT | Global flora and trait archive | TRY reports GIFT contains trait information for 281,836 species and 109 traits. | Screening only | Verify the current download route, licence, and trait dictionary before selecting any field. |
| AusTraits | Australia | AusTraits is an open, CC-BY 4.0 harmonized database with nearly 500 traits for more than 30,000 taxa. | Regional validation / enrichment only | Download a small sample and test its flower-trait fields and name match rate. It cannot be used as the global coverage source. |

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
