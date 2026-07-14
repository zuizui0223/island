# eFloras species-treatment acquisition

## Purpose

`scripts/fetch_efloras_traits.py` inventories taxa from the eFloras hierarchy,
matches species conservatively to the Island v2 master, downloads confirmed
species treatments, and emits source-backed trait candidates for review. It
does not treat rule-based extraction as curated evidence.

The implementation follows the public eFloras structure:

```text
eFloras home / flora landing page
  -> browse or volume entry
  -> family taxon ID
  -> genus taxon ID
  -> species taxon ID
  -> florataxon.aspx treatment / lblTaxonDesc
```

Flora IDs are discovered from the eFloras home page. A requested unlisted ID is
probed and recorded explicitly; IDs without a taxonomic browse structure are
not silently treated as successful empty sources.

## Matching and treatment validation

The Island master is applied after a Flora inventory is created. Match states
are kept distinct:

- `exact_accepted_name_match`
- `normalized_binomial_match` (author text stripped; never collapses an
  infraspecific name)
- `synonym_match` when an explicit synonym CSV is supplied
- `unmatched`
- `ambiguous`

Fuzzy matching is not used. Infraspecies and hybrids remain visible in the
inventory and are excluded from species-level acquisition.

A selected page is saved as a treatment only when `lblTaxonDesc` contains a
bold accepted heading of species rank, enough narrative text, and a botanical
description signal. Bibliography/distribution-only checklist leaves are logged
but are not counted as successful species treatments.

## Outputs

Each Flora (and shard, when enabled) writes:

- `discovered_taxa.csv`
- `master_match_report.csv`
- `species_treatments.csv`
- `trait_candidates.csv`
- `coverage_report.json`
- `errors.jsonl`
- `crawl_checkpoint.json`
- `treatment_checkpoint.json`
- optional raw `html/`

Candidate rows retain the accepted master name, original eFloras name, trait,
candidate value, Flora and taxon IDs, source URL, verbatim source excerpt,
record ID, evidence scope/status, extraction method, and name-match method.
The evidence status is always
`source_backed_rule_extracted_unreviewed` until human review.

The combine mode writes deduplicated cross-Flora treatment/candidate tables and
`combined_coverage.json`, including cross-Flora duplicate counts and coverage
by trait.

An empty result is never reported as success. `no_species_taxa_discovered`
means the published hierarchy ended above species rank (for example, a
family/genus-only checklist). This is distinct from `no_master_matches`, where
species names were enumerated but none matched the island master list.
`partial_unconfirmed_treatments` means at least one selected eFloras leaf was
rejected by the species-treatment contract; the report records its count as
`rejected_selected_treatments` instead of presenting partial coverage as fully
successful.

## Rate limiting and recovery

The default client uses randomized 1.5-3.5 second delays, a descriptive user
agent, a 120-second timeout, and up to six exponentially backed-off attempts.
The GitHub workflow uses 1-2 second delays and at most four Flora jobs in
parallel. It never launches a large single-host request burst.

Both browse state and treatment state are checkpointed. A resumed run reuses
completed treatment rows, retries transport-failed treatment requests, and
deduplicates output keys. Shard assignment is deterministic after conservative
master matching.

After three consecutive species requests each exhaust all configured retries,
the treatment circuit breaker stops that Flora run. Remaining species stay
pending in the checkpoint for a later resume, avoiding a large repeated burst
against a systematically failing source.

## Running

Discover the current published Flora set:

```bash
python scripts/fetch_efloras_traits.py \
  --discover-floras \
  --output-dir data/v2/staging/traits/bulk/efloras
```

Run a 20-species bounded acquisition for one Flora:

```bash
python scripts/fetch_efloras_traits.py \
  --flora-id 2 \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --max-species 20
```

Use the **Fetch eFloras species treatments** workflow for GitHub-hosted runs:

- `smoke`: at most 20 selected species per Flora
- `pilot`: at most 200 selected species per Flora
- `production`: no species limit

For production, first review the smoke/pilot `coverage_report.json` and
`errors.jsonl`, then dispatch `production`. Large production runs should use
`shard_index` and `shard_count`; run every index from zero through
`shard_count - 1` and combine their artifacts afterward.
