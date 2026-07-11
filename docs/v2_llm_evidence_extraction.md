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

No model call is made by the repository. Packet creation and result validation are provider-neutral, so the same frozen packet can be processed by ChatGPT, another hosted model or a local model. The exact model identifier is recorded at ingest.

## Why this can recover more than regex extraction

The deterministic search lane already caches long flora, encyclopedia and institutional descriptions. Rule-based extraction can miss semantically equivalent wording even when the relevant statement is present. The LLM lane runs only after those texts are frozen and asks the model to interpret wording while citing the exact supporting passage.

This separates three operations that were mixed in the former manual workflow:

1. **discovery**: retrieve and cache source text;
2. **extraction**: map explicit source wording to a v2 ontology value;
3. **harmonization and validation**: verify the quote, species, requested trait and allowed value.

A model may improve recall in step 2, but it cannot add a row without evidence from step 1.

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

The model must not:

- browse or use prior knowledge;
- infer effective pollinators from flower colour or form;
- infer biological absence from missing text;
- collapse self-incompatibility, autonomous selfing, mating system or cleistogamy into one state;
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

Invalid claims are written to `llm_rejected_claims.csv`. Strict mode fails the command when any invalid claim is present; `--allow-rejected` keeps valid claims while retaining the rejection audit.

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
