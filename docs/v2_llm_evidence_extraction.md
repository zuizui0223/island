# Reproducible LLM-assisted evidence extraction

## Purpose

The former v1 workflow obtained high coverage by giving species batches directly to ChatGPT. That workflow was useful for semantic extraction, but the input evidence, browsing path and exact model context were not frozen. Repeating the same species list could therefore produce a different undocumented evidence base.

v2 keeps the semantic advantage of an LLM but changes its role. The LLM is an extractor over a frozen source corpus, not a source of biological facts.

```text
bulk / public-source discovery
        ↓
source_texts.csv
        ↓
SHA-locked evidence packets
        ↓
LLM JSONL claims
        ↓
quote + ontology validation
        ↓
unreviewed v2 reported candidates
```

The packet builder and validators make no model call. Packet creation and
result validation are provider-neutral, so the same frozen packet can be
processed by ChatGPT, another hosted model or a local model. The exact model
identifier is recorded at ingest. An optional OpenAI adapter described below
also prepares requests and parses downloaded results offline; preparation alone
does not contact OpenAI or read an API key.

There are two deliberately separate LLM channels:

1. evidence-locked extraction from frozen public text, which can be source-backed
   only after species, ontology value, source chunk, exact quote and packet hashes
   all validate;
2. model-only completion, which is always non-source-backed
   `evidence_type=inference` and `confidence=low`.

Raw provider JSONL is never read by the all-species exporter. Model-only output
must first pass `island-v2-search-enabled-llm ingest`, including an exact
`job_id` and echoed `prompt_sha256`, and enter
`validated_trait_results.csv`. Evidence-locked output enters through
`v2_llm_reported_candidates.csv`, whose `source_chunk_id`, URL and quote are
retained in the final companion evidence ledger.

## Why this can recover more than regex extraction

The deterministic search lane already caches long flora, encyclopedia and institutional descriptions. Rule-based extraction can miss semantically equivalent wording even when the relevant statement is present. The LLM lane runs only after those texts are frozen and asks the model to interpret wording while citing the exact supporting passage.

This separates three operations that were mixed in the former manual workflow:

1. **discovery**: retrieve and cache source text;
2. **extraction**: map explicit source wording to a v2 ontology value;
3. **harmonization and validation**: verify the quote, species, requested trait and allowed value.

In the evidence-locked lane a model may improve recall in step 2, but it cannot
add a source-backed row without evidence from step 1. Model-only completion is
kept outside this path and remains explicitly non-source-backed.

## Build packets from cached source text

First obtain `source_texts.csv`, for example from `island-v2-v1-category-search`. Then build gap-only packets:

```bash
island-v2-llm-evidence prepare \
  --source-csv /tmp/public_flora_scout/source_texts.csv \
  --species-csv data/v2/staging/traits/validation_pilot/validation_pilot_species.csv \
  --candidate-csv /tmp/public_flora_scout/v2_reported_candidates.csv \
  --output-dir /tmp/v2_llm_evidence_packets \
  --batch-size 10
```

Existing candidate CSVs can be repeated. Traits already covered for a species are removed from that species' extraction task.

Each packet contains:

- `packet_input.json`: species, unresolved target traits, allowed ontology values and frozen source chunks;
- `prompt.txt`: the versioned extraction instructions;
- `result_template.jsonl`: an empty JSONL target;
- `packet_manifest.json`: SHA-256 hashes for the prompt, packet input, ontology and source CSV.

Only `species_direct` source text enters a packet. Indirect genus-level or background text is excluded.

## LLM output contract

The model returns one JSON object per supported claim:

```json
{"accepted_species":"Campanula example","trait_name":"floral_form","standardized_value":"bell_campanulate","source_chunk_id":"src_...","evidence_quote":"corolla broadly campanulate","confidence":"high"}
```

Unsupported traits are omitted. The model must not emit `unresolved`; no output means unresolved.

The default extraction target covers flower colour, form, symmetry, tube depth,
size, inflorescence display and reward; sex system, compatibility, selfing,
mating system, herkogamy, dichogamy and cleistogamy; and coarse pollen-vector
mode. Pollinator functional guild is not part of the exhaustive target.

The model must not:

- browse or use prior knowledge;
- infer effective pollinators from flower colour or form;
- infer biological absence from missing text;
- collapse sex system, self-incompatibility, autonomous selfing, mating system,
  herkogamy, dichogamy or cleistogamy into one state;
- infer selfing, outcrossing, compatibility or pollen vector from
  hermaphroditism or monoecy;
- convert self-compatibility into autonomous selfing.

## Validate a completed packet

```bash
island-v2-llm-evidence ingest \
  --packet-dir /tmp/v2_llm_evidence_packets/batch_00001 \
  --result-jsonl result.jsonl \
  --output-dir /tmp/v2_llm_ingest/batch_00001 \
  --model-provider openai \
  --model-id '<exact model identifier>'
```

A claim is accepted only when:

- the species occurs in the packet;
- the trait was an unresolved target for that species;
- the standardized value is in the packet's frozen ontology vocabulary;
- the cited source chunk belongs to that species and is `species_direct`;
- the evidence quote occurs in the frozen source chunk after whitespace normalization;
- the packet and prompt hashes still match their manifest.

The output `v2_llm_reported_candidates.csv` remains `review_status=unreviewed_llm_extraction`. It is source-backed machine extraction, not curated truth.

Invalid claims are written to `llm_rejected_claims.csv`. Strict mode fails the
command when any invalid claim is present, removes any accepted-candidate CSV
from that output directory and writes a failed ingest manifest.
`--allow-rejected` keeps valid claims while retaining the rejection audit.

The all-master exporter does not accept this candidate file by itself. A
publishable evidence-locked bundle contains the validated
`v2_llm_reported_candidates.csv`, `llm_ingest_manifest.json`, `packet_input.json`,
`packet_manifest.json`, `prompt.txt` and the exact `result.jsonl`. The exporter
recomputes all file and extraction hashes, checks the source URL, species-direct
scope and exact quote again, and accepts the directory only through
`--llm-extracted-bundle`. This prevents a candidate CSV from being detached
from the frozen evidence that justified it.

## Reproducibility record

Every validated row carries:

- source URL and source type;
- exact evidence quote;
- source chunk ID and source-text SHA-256;
- prompt version and prompt SHA-256;
- packet-input SHA-256;
- ontology SHA-256;
- model provider and exact model ID;
- deterministic extraction-run ID.

Therefore a later model or prompt can be evaluated against the same frozen corpus without repeating web discovery. Coverage gains can be attributed separately to better discovery versus better semantic extraction.

## All-species job preparation

`prepare-search-enabled-llm-jobs.yml` downloads one immutable all-master
artifact, keeps only the 106,295 applicable angiosperms, removes species whose
five requested fields are already complete, and writes deterministic 128-shard
jobs only for the remaining unknown fields. It makes no model call and uses no
API key. A response is rejected if it fills a field that was not a target of its
frozen job.

Each job binds the species, current taxonomy, target fields, prompt version and
prompt-policy hashes. Provider output must echo the exact `job_id` and
`prompt_sha256`; the optional `job_sha256`, when supplied, must also match. The
ingest command validates the requested nine-column CSV row, rejects values in
non-target fields, and records a canonical result hash and completed job-history
row. Stale rows are pruned when the master, applicability table, current trait
table or prompt contract changes.

## Offline OpenAI Batch adapter

`island-v2-openai-batch-llm` is an optional, fixed-model adapter for the
model-only lane. It uses the `gpt-5.4-nano-2026-03-17` snapshot and prepares
structured `/v1/responses` requests with `store=false`, no tools and no web
search. Its instructions override the search-oriented campaign wording: the
model must not claim that it browsed or verified a source, and every accepted
value remains `inference`/`low` downstream.

Preparing Batch files is offline:

```bash
island-v2-openai-batch-llm prepare \
  --jobs-jsonl /tmp/search_enabled_llm_campaign/jobs_shard_0000.jsonl \
  --output-dir /tmp/openai_batch \
  --max-usd 10 \
  --max-output-tokens 256
```

Repeat `--jobs-jsonl` for additional shards. The command independently
revalidates every campaign job and writes deterministic
`batch_requests_*.jsonl` files plus `request_manifest.json`. It follows the
official Batch limits of at most 50,000 requests and less than 200 MB per input
file. See the [OpenAI Batch API guide](https://developers.openai.com/api/docs/guides/batch).

`--max-usd` is required. Preparation counts every UTF-8 byte in each canonical
request body as one input token, adds the configured maximum output tokens, and
refuses to write requests when that conservative worst-case estimate exceeds
the supplied ceiling. The manifest records the estimate and the pinned model's
Batch rates. This is a fail-closed preparation guard under the recorded prices,
not an OpenAI account billing limit; pricing and account limits must still be
checked before a separately authorized upload. The current rates and model
snapshot are linked from the
[official GPT-5.4 nano model page](https://developers.openai.com/api/docs/models/gpt-5.4-nano).

After an authorized Batch has completed and its raw output has been downloaded,
parsing is also offline and strict by default:

```bash
island-v2-openai-batch-llm parse-output \
  --request-manifest /tmp/openai_batch/request_manifest.json \
  --batch-output-jsonl /tmp/openai_batch/raw_output.jsonl \
  --output-dir /tmp/openai_parsed
```

Missing, duplicate, refused, incomplete, failed or wrong-model responses are
recorded in `parse_audit.csv`. In strict mode any audit error suppresses
`provider_results.jsonl`; a clean file contains only
`job_id`, `prompt_sha256` and the model's `result` string.

## Model-only publication boundary

`ingest-search-enabled-llm-results.yml` binds one immutable prepared-job
artifact to one or more successful provider-result artifacts and runs
`island-v2-search-enabled-llm ingest`. It publishes `validated-llm-traits` only
after job, prompt and canonical result hashes pass. The all-master exporter then
requires these three files together:

- `validated_trait_results.csv`;
- `job_history.csv`;
- `validated_llm_manifest.json`.

Their hashes, row counts, current master taxonomy, target fields, prompt/job
hashes, canonical result hashes and completed history are all revalidated
through `--llm-inference-bundle`. Even a valid bundle cannot create
source-backed evidence: it is forced to `evidence_type=inference` and
`confidence=low`, has no source URL, and direct model-only `SI`/`SC` values are
weakened to `likely_SI`/`likely_SC`. Raw provider JSONL is never accepted by the
all-master exporter.
