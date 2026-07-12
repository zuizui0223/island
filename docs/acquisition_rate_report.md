# Broad trait acquisition-rate report

## Why this exists

v2 runs several trait-acquisition channels in parallel — human-reviewed
curation, the machine wave campaign, bulk source joins, and search-enabled LLM
extraction — and each writes its own staging output. Until now no single view
answered the operational question that decides whether the pipeline is working:

> Of the whole master species universe, what fraction now has a source-backed
> candidate for each trait, and which channel supplied it?

`island-v2-acquisition-rate` unions every channel against the frozen master
species set (`data/v2/staging/gbif/collected/island_taxa.csv`) and reports the
**broad acquisition rate**: the fraction of the master with a source-backed
candidate. This is the confirm-instrument for a rate-first workflow — run it,
land a new source, run it again, and read the delta.

## What "broad" means

- **Broad-acquired** = any channel marked `counts_as_broad: true` supplied a
  non-empty value for a species-trait. It is candidate *presence*, not an
  accepted trait value and not a biological absence for the missing remainder.
- Inference/proxy channels (currently the search-enabled LLM track) are still
  counted and shown by track, but never fold into the broad rate.
- A species-trait with only empty/non-committal values (`unknown`, `na`, …)
  does not count as acquired.

The scoreboard is deliberately upstream of adjudication. Quality is enforced
later, in the curated review track; this report measures reach.

## Run it

```bash
island-v2-acquisition-rate run --output-dir data/v2/staging/traits/acquisition_rate
```

Outputs:

- `acquisition_rate_summary.json` — headline colour, form, core-minimum
  (colour **and** form), any-trait broad rates, plus off-master drops per source.
- `acquisition_rate_by_trait.csv` — per-trait broad vs. any-channel species
  coverage and rate.
- `acquisition_rate_by_track.csv` — per-channel species reach and species-trait
  pair counts.

## Configuration

`config/acquisition_rate.yml` declares the master denominator, the core trait
sets, the reported trait order, and the channels. Each channel names an adapter
and whether it counts toward the broad rate:

| adapter | reads | notes |
|---|---|---|
| `candidate_index` | one deduplicated machine candidate index `.csv.gz` | the wave campaign |
| `candidate_long` | generic long trait CSV/gz by glob | bulk joins, new sources |
| `curated_evidence` | accepted human-reviewed evidence by glob | keeps only `review_status=accepted` |
| `search_enabled_llm` | LLM result `batch_*.jsonl` by glob | inference track, not broad |

A channel whose files do not exist yet contributes nothing and is safe to leave
declared. A newly landed bulk source (EOL TraitBank, Wikidata, other free bulk
downloads) appears in the rate as soon as it writes to its staging path in the
`candidate_long` layout — no code change needed.

## Reading the baseline

The current checked-in snapshot shows how far the pipeline still has to go and
is the reason the rate-first workflow exists: the machine wave campaign is the
only channel with real reach so far, and even flower colour — the single most
important target — sits far below 1% of the master. Landing bulk sources is the
lever that moves these numbers, and this report is how each landing is judged.
