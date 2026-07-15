# Resumable public-web trait shards

## Production contract

The public-web campaign follows one fixed path:

```text
106,295 family-classified angiosperm species
  -> deterministic scientific-name hash
  -> 128 shards (763-893 species per shard)
  -> public-web search for remaining gaps
  -> frozen search evidence packet
  -> exact nine-column extraction
  -> enum and species-set validation
  -> shard checkpoint
  -> cache restore and resume
```

The runner is `island-v2-web-trait-shards` and the workflow is
`run-public-web-trait-shards.yml`. No paid API or model key is used. Existing
global machine candidates are loaded first, so a source is queried only when a
target field remains missing.

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

## GitHub Actions persistence

The workflow runs at most two shards concurrently to bound traffic per public
site. It advances 50 species per shard by default and is scheduled every six
hours. At that rate a first pass takes at most 18 scheduled advances for the
largest shard, plus retry passes.

Wikipedia is queried directly by exact scientific-name title, with bounded
backoff and a per-species pause. This avoids extra Wikidata calls. World Flora
Online is disabled in the production workflow while its TLS chain fails on the
GitHub-hosted runner. When explicitly enabled, that deterministic source failure
is retained in packet diagnostics but is not treated as a retryable species failure.
The production cache uses contract v5. It can migrate the immutable v4 shard
artifacts from run `29186029383`, preserving completed evidence, adding 711 newly
eligible angiosperms, and removing 28 now-inapplicable names. Pre-tuning pilot
caches remain intentionally ignored.

Actions cache stores only resume state and cumulative outputs:

- `checkpoint.csv`;
- `campaign_manifest.json` and `campaign_status.json`;
- `cumulative/trait_results.csv`;
- cumulative source/evidence/error/v2-candidate tables.

The current evidence packet is uploaded as a 90-day artifact. Packet directories
are deliberately excluded from the resume cache to prevent repeated source text
from inflating every cache generation.

## Operations

Run all shards manually:

```bash
gh workflow run run-public-web-trait-shards.yml \
  -f shard_start=0 \
  -f shard_end=127 \
  -f batch_size=50
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
