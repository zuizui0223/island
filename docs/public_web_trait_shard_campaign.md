# Resumable public-web trait shards

## Production contract

The public-web campaign follows one fixed path:

```text
106,295 family-classified angiosperm species
  -> deterministic scientific-name hash
  -> 128 shards (763-893 species per shard)
  -> species x provider checkpoint
  -> public-web/API search for remaining gaps
  -> frozen search evidence packet
  -> exact nine-column extraction
  -> enum and species-set validation
  -> shard checkpoint
  -> cache restore and resume
```

The v6 runner is `island-v2-web-trait-shards` and the workflow is
`run-public-web-trait-shards.yml`. No paid API or model key is used. Existing
global machine candidates are loaded first, so a source is queried only when a
target field remains missing.

The active acquisition fields are flower colour, flower shape/symmetry, mating
system and self-incompatibility. The legacy nine-column output still preserves
any pollination-guild evidence already found, but guild is no longer a field
whose absence triggers another provider lookup.

The enabled production providers are strict GBIF species descriptions,
Wikipedia, indexed flora, illustrated-plant and horticultural reference sites
(NCSU Plant Toolbox, NZPCN, PFAF and Useful Tropical Plants), and exact-species
Europe PMC title/abstract discovery. OpenAlex requires an API key and remains
disabled. World Flora Online remains disabled because its robots policy
disallows automated retrieval.

Search engines are discovery channels, not biological evidence. Google result
HTML is not scraped. An official Google Custom Search JSON or Brave Search API
adapter may write URL/title/snippet discovery records when credentials and a
hard query budget are supplied, but snippets cannot create a trait row. The
underlying page must still pass domain policy, species identity, retrieval and
quote validation before it can enter `source_texts.csv`.

These source classes are not treated as interchangeable. A field study or
review can support pollination or mating-system claims that a descriptive plant
page may not. Flora, encyclopedia and horticultural descriptions can support
the traits they state directly, but flower colour or form alone is never used
to infer a pollination guild. LLM processing is a separate downstream
extraction or low-confidence inference lane; it is not an additional public
source. See [Reproducible LLM-assisted evidence extraction](v2_llm_evidence_extraction.md).

The 106,295 denominator is produced by `config/angiosperm_scope.yml` from the
full 115,328-name master. All 9,033 non-angiosperm or family-unresolved names
remain in the final nine-column export as explicit `unknown`/`low` rows rather
than receiving floral values. An unexpected denominator or master-fingerprint
change still fails closed.

## Evidence packet

Every attempted batch writes an immutable packet containing:

- `species.csv`;
- `source_texts.csv` with retrieved text, URL, citation and source type;
- `lookup_errors.csv`;
- `v1_category_evidence.csv` with matched excerpts and rules;
- `v1_category_traits.csv` with exactly the nine requested columns;
- `v2_reported_candidates.csv` for the active v2 long-format pipeline;
- `packet_manifest.json` with file hashes, row counts and validation contract.

Public artifacts store bounded trait-oriented excerpts rather than whole plant
pages. The citation records the SHA-256 of the retrieved visible text and the
stored excerpt length; URLs and source attribution remain attached.

The nine-column table is retained because it matches the previous extraction
contract. The v2 candidate table is the analysis-facing output. Inferred values
remain explicitly `likely`/low-confidence and are not promoted to direct
pollinator evidence.

## Checkpoint semantics

Each shard has an independent `checkpoint.csv`:

- `pending`: not attempted;
- `running`: reserved by the current run;
- `retry`: a species-specific transient lookup failed;
- `completed`: search and validation completed, including valid all-`unknown`
  zero-hit rows;
- `exhausted`: lookup still failed after `max_attempts`.

An interrupted `running` row becomes `retry` on restore. A missing Wikipedia
page is a completed source miss, while a timeout/HTTP failure remains retryable.
Neither `unknown` nor `exhausted` is interpreted as biological absence.

v6 also has `provider_checkpoint.csv`, exactly one row per shard species and
configured provider. Provider states are `pending`, `running`, `hit`,
`no_hit`, `skipped_covered`, `retry`, `exhausted`, or migrated
`legacy_completed`; `skipped_covered` means seed evidence already supplied every
field targeted by that provider, whereas `no_hit` requires a completed lookup.
Adding or repairing a provider therefore requeues it without discarding
evidence already obtained from another source. A provider-wide sitemap/API
failure is fanned out to every affected species and cannot silently become a
completed zero-hit row.

## GitHub Actions persistence

The workflow runs at most four shards concurrently to bound traffic per public
site. It advances up to 500 species per shard by default and is scheduled every
six hours. Plant-site sitemaps are fetched once in the plan job, hash-frozen,
and shared by all 128 shards instead of being downloaded separately by every
job. A sitemap digest is retained as provenance but is not a provider
implementation version: an ordinary sitemap refresh therefore cannot reset all
106,295 species. If one sitemap is temporarily unavailable, the frozen partial
index proceeds and the affected provider/species rows remain retryable.

On `main`, a successful full 128-shard run automatically triggers exact
artifact promotion. Promotion dispatches the next full pass only when enabled
provider work remains and the completed pass advanced at least one species.
The successful promoted artifact in turn triggers the 115,328-row all-master
build. Pull-request validation uses a separate concurrency group and cannot
replace a pending production pass. Cache is only an accelerator; the next pass
can bootstrap from the newest complete 128-artifact run across branch scopes.

Wikipedia is queried directly by exact scientific-name title, with bounded
backoff and a per-species pause. GBIF accepts only exact rank=SPECIES matches or
an exact synonym whose accepted species identity is independently verified.
PFAF's ASP.NET form body is retained by a site-specific parser, and the former
silent 100-species cap for Useful Tropical Plants has been removed. Europe PMC
uses its REST API only, retains exact full-species title/abstract matches with
reproductive terms, and performs no trait inference.

The production cache uses contract v6. It migrates immutable v4/v5 artifacts
from run `29186029383`, preserves their evidence, reuses completed Wikipedia
accounting, and deliberately requeues GBIF and plant-site providers because
their identity/extraction contracts were repaired.

Actions cache stores only resume state and cumulative outputs:

- `checkpoint.csv`;
- `provider_checkpoint.csv`;
- `campaign_manifest.json` and `campaign_status.json`;
- `cumulative/trait_results.csv`;
- cumulative source/evidence/error/v2-candidate tables.

The current evidence packet is uploaded as a 90-day artifact. Packet directories
are deliberately excluded from the resume cache to prevent repeated source text
from inflating every cache generation.

## Operations

### Bounded official-search discovery

The standalone official-search CLI uses only official Google Custom Search JSON
and Brave Search APIs; it does not scrape Google, Brave, or other search-result
HTML. It neither provisions credentials nor bypasses provider quotas. Supply an
existing Google API key together with its Custom Search engine ID (`cx`), a
Brave Search API key, or both. A positive `--max-queries` is mandatory and is a
hard campaign-wide cap: every HTTP attempt, including a retry or Brave fallback,
consumes one unit.

Keep secrets in environment variables rather than command arguments:

```bash
export GOOGLE_CUSTOM_SEARCH_API_KEY="..."
export GOOGLE_CUSTOM_SEARCH_CX="..."
export BRAVE_SEARCH_API_KEY="..."  # optional fallback when Google is configured

island-v2-official-search-discovery \
  --species-csv data/species.csv.gz \
  --output-dir outputs/official-search-discovery \
  --max-queries 1000
```

The input must contain a `species` column; plain CSV and gzip-compressed CSV are
accepted. `discovery.csv` contains only species, provider, exact query, rank,
result URL, title, snippet, retrieval time and record hash. Search metadata and
snippets are discovery leads only and never biological evidence. A `result_url`
must subsequently be retrieved under the site's access policy, matched to the
exact species, and supported by a verified quotation from the source page before
it can become a trait candidate.

Run all shards manually:

```bash
gh workflow run run-public-web-trait-shards.yml \
  -f shard_start=0 \
  -f shard_end=127 \
  -f batch_size=500
```

Run a bounded diagnostic subset:

```bash
gh workflow run run-public-web-trait-shards.yml \
  -f shard_start=0 \
  -f shard_end=3 \
  -f batch_size=10
```

Normal reruns restore the newest cache for each shard and continue from the
first pending row. `retry_exhausted=true` is reserved for an intentional retry
after an external source or extraction bug has been fixed.

## First promoted public-web snapshot

Promotion run
[`29385927769`](https://github.com/zuizui0223/island/actions/runs/29385927769)
validated and combined all 128 immutable v4 artifacts from source run
`29186029383`. The accounted checkpoint contains 105,612 species: 6,670
completed, 98,842 pending and 100 retry. It contains results for 6,720 species,
1,383 species with at least one requested value, 3,560 evidence rows and 4,638
source-text rows. This is the initial snapshot already integrated into the
115,328-row all-master export; it is not claimed to be the completed campaign.
The promoted artifact is `public-web-traits-combined` (artifact ID
`8331432629`, digest
`sha256:c42941ccfa3b666a496da9fa458556a5be6dc07b6cd1583736746c6e433bd96e`).

Integration run
[`29386456144`](https://github.com/zuizui0223/island/actions/runs/29386456144)
verified exactly 115,328 unique species rows in the requested nine-column
table. Of the 106,295 applicable angiosperms, the historical snapshot had at
least one requested value for 18,764 species; 96,564 of all master rows still
had all five biological fields `unknown`. Requested-field coverage was 13,757
species for flower colour, 3,254 for flower shape, 8,612 for pollination guild,
21 for mating system and 193 for self-incompatibility. These are acquisition
coverage counts, not a review-complete claim. Future all-master builds accept
LLM rows only through the validated bundle contracts described in the LLM
document; legacy raw LLM JSONL is not an input.

The historical artifact is `all-master-trait-ledger-29386456144` (artifact ID
`8331626313`, digest
`sha256:d53e1212d0f8cbc614de7aab4ed78311250c7972dede0eed894c79595dfdcdde`).
Its nine-column table SHA-256 is
`51781e6b234cb4aef11e865abc1c5c0a49c32a74963162898d45f3681ab66863`.
