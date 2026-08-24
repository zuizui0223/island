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

## Public-page and LLM gap-filling lanes

Bulk databases remain the first pass, but they are not the only configured
source class. Remaining gaps can be searched in public floras, illustrated
plant references, institutional plant sites, horticultural references,
Wikipedia, GBIF descriptions and exact-species scholarly metadata. The
[public-web shard campaign](public_web_trait_shard_campaign.md) records the
retrieved URL, bounded source excerpt, species identity and source-text hash.
A page contributes only the traits that its text states; descriptive colour or
shape is not converted into a pollination or mating-system claim.

Google and Brave are bounded URL-discovery options through their official APIs.
Search-result HTML is not scraped, and titles or snippets never count as trait
evidence. A discovered URL must be retrieved under the source site's access
policy and provide an exact species-matched excerpt before a source-backed row
can be created.

LLM processing has two separate acceptance paths, described in
[Reproducible LLM-assisted evidence extraction](v2_llm_evidence_extraction.md):

- quote-locked extraction can remain source-backed only when the frozen page,
  exact quote, species, ontology value and all packet hashes validate; and
- model-only completion has no source URL and is always stored as
  `evidence_type=inference`, `confidence=low`; direct `SI`/`SC` model outputs
  are weakened to `likely_SI`/`likely_SC` at export.

Raw model responses and generic candidate CSVs are not accepted as LLM inputs
to the all-master exporter. Source-backed rows outrank validated model-only
inference when both exist.

## Whole-master accounting contract

Acquisition coverage and output-row coverage are different quantities. Every
all-master build must emit exactly one nine-column row for each of the 115,328
unique master species. The current scope contains 106,295 applicable
family-classified angiosperms; all 9,033 non-angiosperm or family-unresolved
names remain present as explicit all-`unknown`, `confidence=low` rows. An
applicable all-`unknown` row means "not yet acquired", never biological absence.
The current requested-field counts and remaining gap count are reported with
the promoted public-web snapshot rather than described as complete coverage.

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

## Accounted eFloras production lane

The all-Flora workflow built complete inventories for Flora of North America,
Flora of China, and the Annotated Checklist of the Flowering Plants of Nepal,
then split every Island-matched treatment request across 12 shards. A page only
becomes evidence when `lblTaxonDesc` contains a confirmed species treatment
whose heading resolves back to the Island master. Genus pages, short or empty
records, and mismatching headings are rejected explicitly rather than reused as
species evidence.

The source run
[`29330032130`](https://github.com/zuizui0223/island/actions/runs/29330032130)
accounted for all 18,251 selected treatment requests: 14,734 confirmed treatment
rows and 3,517 rejected records, with zero pending requests and zero HTTP
errors. The confirmed records yielded 14,008 direct candidates for 8,668 Island
species: 7,834 flower-colour, 2,069 floral-form, 1,032 sex-system, 101
floral-symmetry, and 79 self-incompatibility species.

The original jobs returned a failure conclusion because their first contract
treated any rejected page as an incomplete run. The successful promotion run
[`29382678187`](https://github.com/zuizui0223/island/actions/runs/29382678187)
re-collected the immutable 12 shard artifacts without requesting eFloras again,
verified `confirmed + rejected = selected` for every shard, and added all 3,517
rejection records to `all_acquisition_failures.jsonl`. Its 90-day artifact is
`efloras-all-floras-combined` (artifact ID `8330314679`, digest
`sha256:9abe55a5c50c22c3bb272bf1509bf3b99c47e80fbd698660a4dd0382bfd8a896`).

The evidence-preserving re-extraction run
[`29409302362`](https://github.com/zuizui0223/island/actions/runs/29409302362)
reused those same 14,734 confirmed treatment texts and expanded only direct
floral-context terms (for example stellate, infundibuliform, hypocrateriform,
capitulum and asymmetric). It made no new taxonomic matches and no trait
inference. The promoted artifact contains 14,502 candidate rows: floral-form
coverage increased from 2,069 to 2,316 species and floral-symmetry coverage
from 101 to 141 species. Its artifact ID is `8340405765` and digest is
`sha256:adf0ed8529a22faf37e22873864469e29fc59dd5ece04fa4956926b82055a502`.

## USDA PLANTS direct Flower Color lane

The `Acquire USDA PLANTS direct floral traits` workflow snapshots the public
USDA PLANTS characteristic-result and characteristic-filter API responses. The
API is not versioned, so the pipeline archives the exact JSON fields, retrieval
time, query URLs, and SHA-256 checksums. It follows the database's
[official scope and citation guidance](https://plants.sc.egov.usda.gov/help).

Only the direct `Flower Color` values Blue, Brown, Green, Orange, Purple, Red,
White, and Yellow are mapped. Source IDs must also be explicitly classified as
Dicot or Monocot. Non-angiosperms, infraspecific names, hybrids, unknown values,
and unresolved identifiers remain in the exclusion audit; no value is inferred.
Obsolete IDs absent from the bulk result are resolved only through their exact
public `PlantProfile/{id}` response.

Production run
[`29382621509`](https://github.com/zuizui0223/island/actions/runs/29382621509)
snapshotted 2,117 source Flower Color IDs, prepared 1,868 conservative binomial
angiosperm rows, and retained 246 explicit exclusions. Exact Island-master
matching produced 1,451 flower-colour candidates across 1,451 species, with
zero empty required evidence fields. Its 90-day artifact is
`usda-plants-direct-traits-29382621509` (artifact ID `8330303160`, digest
`sha256:786f65fb8f10d786aa2ffb44a5f1f91522a9dc0f452b8331f42841eb2103bdf1`).

## Five-source all-master production ledger

The five-source production build
[`29382735518`](https://github.com/zuizui0223/island/actions/runs/29382735518)
materialized the exact EOL, AusTraits, GIFT, eFloras, and USDA artifacts above.
It retains all 115,328 master species and all 1,499,264 species-trait cells.
Every required evidence field is nonempty in every materialized source.

Broad direct-candidate coverage is 23,006 species (19.95%). The exhaustive
species ledger classifies 22,621 species as partially source-backed, 83,661 as
eligible but still without a configured-source candidate, 9,033 as outside the
taxonomic scope of floral traits, and 13 as having only non-source-backed
candidates. These are acquisition states, not biological absence claims.

Coverage by requested trait is 13,257 species for flower colour, 2,124 for
floral form, 855 for symmetry, 23 for flower-size class, 1 for tube depth,
8,674 for pollen-vector mode, 9,035 for pollination guild, 11,804 for sex
system, 37 for mating system, 178 for self-incompatibility, 308 for autonomous
selfing, 2 for cleistogamy, and 69 for dichogamy. The artifact
`all-master-trait-ledger-29382735518` has ID `8330348999` and digest
`sha256:00c3b7c0a12973bde27fe70f0387b48bcf410e3869402fd89a52e5c932f36c9a`.

## Artifact continuity

The monthly `Refresh durable trait source artifacts` controller dispatches
fresh EOL, AusTraits, eFloras, GIFT, and USDA snapshots before their Actions
artifacts expire. eFloras intermediate and combined artifacts use the same
90-day retention window as the other bulk lanes. An automatic all-master build
resolves the newest successful run that still contains the exact expected,
non-expired artifact name and records the resolved run IDs in its summary.
For a historical reconstruction, manually set `resolve_latest_sources=false`
and supply the pinned run IDs instead.

## Source-package-first reproductive acquisition

The Ferrer et al. global self-incompatibility workbook is acquired at workflow
runtime from its stable OUP article resolver and accepted only when its exact
SHA-256 is `8ef2f5ac99780ca19a15b847442272457634064577f8a12fbb9710f5521e5899`.
The copyrighted workbook is not committed. The artifact retains the download
hash, source citation, row identifier, exact structured supporting excerpt,
name-match audit, exclusions, conflicts, and derived candidates.

The importer deliberately does not treat the complete published SI column as
species-direct evidence. Reference-supported SC is retained; SI and partial SI
require a phenotype reaction plus a cited source. Self-sterility alone and the
authors' genus-transferred SI assignments are excluded. Exact accepted names
and strict, family-consistent GBIF synonyms are allowed. Multiple values that
resolve to one accepted species are quarantined rather than selected by row
order. Values agreeing with the Meyer compilation are assigned the Meyer
lineage so the two redistributed records cannot masquerade as independent
support.

This adapter emits no genus inference. Its High/Medium direct rows enter the
same all-master artifact as other bulk sources, and the integrated workflow
then rebuilds trait-specific Validated Low from all direct evidence. The
integrated workflow accepts an explicit all-master Run ID and artifact name, so
a formal coverage run can be pinned to the exact source-package build that it
consumed.

## Pladias Czech-flora source package

`island-v2-pladias-traits` downloads the hash-pinned official Pladias taxon
list, intersects exact species names with the island master, and fetches each
matched species page once. Page-level evidence records are cached and retained
with exact page identity, family breadcrumb, URL, retrieval time, content hash,
feature-specific quotation and Pladias feature citation. A rerun reuses every
completed page record instead of searching the species again.

The package emits only species-direct High candidates from the institutional
database. It preserves multistate flower colours, maps flower symmetry and
supported inflorescence terminology trait-by-trait, and keeps
self-incompatibility separate from mating system. Pollination syndrome is not
used as floral structure, and neither pollen-vector nor reward evidence is
projected onto the strict three axes. Apomixis and sterility are retained only
in the raw page record because the current ontology does not justify treating
them as reproductive assurance.

Every build also writes a fixed-seed, trait-stratified 200-species extraction
audit. It replays the exact-name identity gate, controlled-value mapping,
supporting quotation, content hash and cultivar exclusion, and reports
`accepted_correct / reviewed` plus cultivar contamination. This is an
extraction-contract audit, not an independent biological re-measurement of the
curated Pladias values; that distinction is recorded in the audit manifest.

Pladias requires database and original-source citation for scientific reuse;
the adapter records those citations and does not assert a Creative Commons
license. Because this package uses a large part of the database, its manifest
also records the requirement to contact the Pladias Governing Board before
publication. The all-master workflow materializes the package under
`data/v2/staging/traits/bulk/pladias_czech_flora/`, after which the common
integrator performs source-lineage deduplication, direct conflict handling and
the trait-specific all-evidence Validated Low rebuild.
