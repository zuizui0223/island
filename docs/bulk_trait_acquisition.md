# Bulk trait acquisition for the global island flora

## Purpose

The global species master is expected to grow beyond 100,000 accepted species.
It is therefore neither affordable nor scientifically sensible to make
species-by-species web retrieval the primary acquisition path.

The primary path is:

```text
public bulk trait export
-> source-specific codebook/profile
-> exact species-name matching to island_taxa.csv
-> pending source-backed evidence table
-> trait and family coverage audit
-> conservative / expanded sensitivity analyses
```

## What bulk import does and does not do

The importer accepts a downloaded long-format CSV or TSV. A YAML profile records:

- source columns containing scientific name, trait label, and trait value;
- any source-specific mapping to the v2 trait ontology;
- source reliability and source type; and
- the distinction between global-bulk priority traits and evidence-rich subset
  traits.

It does **not** infer trait values. It does not use a genus or family value for
an unmatched species. It does not convert a source code without a profile rule.
It does not automatically curate a candidate into accepted evidence.

Each source row becomes one `source_backed`, `species_direct`, `pending`
candidate only when all three conditions hold:

1. the source taxon matches the master accepted name exactly or differs only by
   an author citation;
2. the source trait label maps to a named v2 ontology trait; and
3. the source value maps to an allowed value for that trait.

Names containing `subsp.`, `var.`, or other infraspecific ranks are deliberately
not collapsed to a species-level master row. They remain in the name-match audit.

## Output layers

One run writes these separate products:

- `bulk_trait_candidates_*.csv`  -  pending, source-backed, species-direct rows;
- `bulk_trait_name_matches_*.csv`  -  every source-row taxonomic match decision;
- `bulk_trait_unmatched_taxa_*.csv`  -  records that require explicit taxonomic
  reconciliation rather than guessing;
- `bulk_trait_unmapped_values_*.csv`  -  source codes that require a reviewed
  profile mapping;
- `bulk_trait_coverage_by_trait_*.csv`  -  both species coverage and weighted
  island-species-pair coverage; and
- `bulk_trait_coverage_by_family_*.csv`  -  taxonomic bias diagnostics.

The manifest records the downloaded source checksum and all paths.

## Trait roles

`global_bulk_priority` traits are the traits for which direct broad coverage is
worth pursuing first: flower colour, symmetry, gross floral form, tube depth,
flower size, display, reward, pollen-vector mode, and sex system.

`evidence_rich_subset` traits are retained in the ontology but should not be
forced into a 100,000-species global table: pollination guild, SI/SC,
autonomous selfing, mating system, herkogamy, dichogamy, and cleistogamy.
They are analysed later as source-rich sensitivity panels, not imputed globally.

## Running an import

The GitHub Actions workflow `ingest_bulk_trait_source` is manual and needs no
OpenAI secret. It takes a direct dataset download URL, a stable source/version
name, a reviewed source profile, a dataset citation, and a stable dataset page
URL. It uploads only artifacts; it never commits downloaded data or candidates
into the repository.

Before running a large source, inspect a small sample and add a source-specific
profile to `config/bulk_trait_sources.yml`. A source profile is evidence
provenance, not mere plumbing: it is where trait-code semantics are recorded.

## EOL TraitBank all-traits lane

EOL's public "All trait data" Zenodo record provides `traits_all.zip`, an export
with `pages.csv`, `traits.csv`, `metadata.csv`, `inferred.csv`, and `terms.csv`.
That archive is not the canonical long CSV accepted by the generic importer, so
it has a dedicated preparation step:

```text
traits_all.zip
-> island-v2-bulk-traits prepare-eol
-> eol_traitbank_v2_long.csv
-> island-v2-bulk-traits ingest
```

`prepare-eol` streams `traits.csv` in chunks, resolves `page_id` through
`pages.csv`, resolves term labels through `terms.csv`, and exports only
predicate/value pairs declared in `config/eol_traitbank_terms.yml`. Unmapped EOL
terms are retained as counts in `eol_traitbank_unmapped_terms.csv`; they are the
review queue for expanding the mapping.

Use the manual workflow `ingest_eol_traitbank_bulk` for this lane. Start with a
small `max_trait_rows` smoke test, inspect the unmapped-term and name-match
audits, then rerun with `max_trait_rows = 0` for the full archive. `inferred.csv`
is deliberately not used for direct evidence candidates; only rows attached to
the EOL trait record itself enter the v2 pending-evidence path.

## GIFT 3.2 direct-record lane

The `Acquire GIFT direct floral traits` workflow uses the
[official versioned public API](https://biogeomacro.github.io/GIFT/articles/web_only/GIFT_API.html)
rather than an aggregated trait table. Its filtering follows the documented
[raw-trait contract](https://biogeomacro.github.io/GIFT/reference/GIFT_traits_raw.html)
and records the [GIFT data-quality disclaimer](https://biogeomacro.github.io/GIFT/articles/web_only/GIFT_Disclaimer.html):

```text
GIFT 3.2 metadata + references + species + traits_raw
-> deriv=0 / biasderiv=0 and public-reference checks
-> exact categorical ontology mappings + exclusion audit
-> version-pinned GIFT work-species names
-> generic canonical-long bulk importer
-> Island-master pending candidates
```

`island-v2-fetch-gift-traits` archives every returned raw row, the full species
and reference metadata used, query URLs, file checksums, and every exclusion
reason. The evidence excerpt is a compact serialization of the exact raw GIFT
fields, while `source_citation` retains the underlying bibliographic
`ref_long`. API rows leaking outside the requested reference/trait pair,
infraspecific records, unresolved names, restricted/missing references, and
ambiguous ontology mappings never become candidates.

The full GitHub production run
[`29338693288`](https://github.com/zuizui0223/island/actions/runs/29338693288)
archived 14,122 raw API rows and produced 19,243 Island-matched candidates
across 6,983 species. Required evidence fields were empty in zero candidates.
The resulting source coverage is 5,288 species for functional guild, 5,210 for
pollen-vector mode, 1,929 for sex system, and 1,613 for flower colour. Its
90-day artifact is `gift-direct-traits-29338693288` (artifact ID `8313108327`,
digest `sha256:ae660ca791bab115f21d8b2ce7c9fdf65446c66365b274b6e33a7af1ada9ccca`).

The three-source EOL + AusTraits + GIFT all-master build
[`29339087438`](https://github.com/zuizui0223/island/actions/runs/29339087438)
retains all 115,328 species and all 1,499,264 species-trait cells. It raises
source-backed any-trait coverage from 12,454 to 16,815 species and reduces the
eligible no-candidate count from 94,063 to 89,833. This remains incomplete
coverage, not a claim that the other species lack the traits.
