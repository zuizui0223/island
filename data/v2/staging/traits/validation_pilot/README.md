# Trait-evidence validation pilot (Stage 6)

A **stratified gold-standard validation pilot** for the trait-evidence
review-packet pipeline. It is a candidate species list for **human review**, not
an analysis input. No trait value, native/establishment status, Bombus
applicability, or analysis inclusion is decided here.

## What this is

- `validation_pilot_species.csv` — 50 species, one from each of 50 major
  angiosperm families, selected by the Stage-6 stratified selector
  (`island-v2-trait-pilot-eval build`).
- `validation_pilot_summary.json` — selection provenance (pool size, strata,
  per-stratum counts).
- `reported_proxy_candidates/` — unreviewed reported-source and rule-based
  proxy candidates from GitHub Actions run `29060068568` on commit
  `841354972a2bae80b407d593228715c2a56b1702`.
- `globi_interaction_evidence/` — unreviewed explicit GloBI interaction records
  from GitHub Actions run `29017098053` on commit
  `8f586d420178cec07c049641a7f2f1168a0f5699`.

Rebuild with:

```
python scripts/build_validation_pilot.py \
    --taxa data/v2/staging/gbif/collected/island_taxa.csv \
    --output-dir data/v2/staging/traits/validation_pilot
```

## How it was selected (and what it deliberately avoids)

- Pool: the collected `island_taxa.csv` (115,328 species) restricted to
  angiosperms — ferns, lycophytes, and gymnosperms are excluded because the study
  is about floral syndromes (105,795 species / 479 families).
- Selection spreads **evenly across families** (round-robin, largest families
  first, deterministic). It is **not** driven by GBIF record count, which would
  bias the pilot toward human-associated / cultivated / easily observed species
  and hide exactly the off-target-literature and wrong-transfer failure modes the
  pilot is meant to measure. In the current pilot, species with a single GBIF
  record sit alongside species with thousands.

## Fail-closed status of the stratification axes

Only `family` is populated from data. The trait-derived axes —
`native_status`, `growth_form`, `flower_conspicuousness` — are seeded `unknown`
on purpose: those labels are exactly what this pilot's human review supplies, so
they are never auto-decided. Two ways to enrich them before review:

1. `native_status` can be **derived from GBIF `establishmentMeans`** already in
   the collected occurrences (data, not a guess).
2. `growth_form` / `flower_conspicuousness` remain review targets, or can be
   pulled from a free source (e.g. GBIF/Wikidata) as candidates, still requiring
   review.

Once these axes carry values, rerunning the builder stratifies across all four
axes (native/introduced × woody/herb/vine/epiphyte × family × conspicuousness).

## Generated unreviewed candidate artifacts

The checked-in candidate artifacts are staging products for review only. They do
not decide trait values, pollinator effectiveness, native status, Bombus
applicability, or analysis inclusion.

- `reported_proxy_candidates/web_reported_scout/web_reported_candidates.csv`:
  27 Wikipedia reported candidates across 15 species; lookup errors: 0.
- `reported_proxy_candidates/m0_descriptive_scout/m0_descriptive_candidates.csv`:
  1 GBIF species-description candidate; lookup errors: 0.
- `reported_proxy_candidates/trait_proxies/trait_proxy_candidates.csv`:
  17 rule-based `floral_syndrome_proxy` candidates across 17 species.
- `globi_interaction_evidence/globi_interaction_evidence.csv`:
  21 explicit GloBI pollination interaction claims across 5 queried taxa; query
  errors: 0.

## Next step: source discovery

Literature source discovery (Stage 1, `island-v2-trait-source-discovery`) runs
against the public OpenAlex/Crossref/Unpaywall APIs. Those hosts are **blocked by
egress policy in the Claude Code session environment**, so discovery is run from
the GitHub Actions workflow (`discover-core-pilot-trait-sources.yml`), whose
runners have network access. The pipeline is fail-closed: an unreachable source
yields zero leads with a logged error, never a fabricated lead.
