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
