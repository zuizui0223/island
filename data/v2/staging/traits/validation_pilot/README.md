# Trait-evidence validation pilot (Stage 6)

A **stratified gold-standard validation pilot** for the trait-evidence
review-packet pipeline. It is a candidate species list for **human review**, not
an analysis input. No trait value, native/establishment status, Bombus
applicability, or analysis inclusion is decided here.

## What this is

- `validation_pilot_species.csv`  -  50 species, one from each of 50 major
  angiosperm families, selected by the Stage-6 stratified selector
  (`island-v2-trait-pilot-eval build`).
- `validation_pilot_summary.json`  -  selection provenance (pool size, strata,
  per-stratum counts).
- `reported_proxy_candidates/`  -  unreviewed reported-source and rule-based
  proxy candidates from GitHub Actions run `29060068568` on commit
  `841354972a2bae80b407d593228715c2a56b1702`.
- `candidate_review_queue/`  -  prioritized adjudication queue generated from the
  checked-in reported/proxy/descriptive candidates. Decision columns are blank;
  no candidate is curated by queue generation.
- `globi_interaction_evidence/`  -  unreviewed explicit GloBI interaction records
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
  angiosperms  -  ferns, lycophytes, and gymnosperms are excluded because the study
  is about floral syndromes (105,795 species / 479 families).
- Selection spreads **evenly across families** (round-robin, largest families
  first, deterministic). It is **not** driven by GBIF record count, which would
  bias the pilot toward human-associated / cultivated / easily observed species
  and hide exactly the off-target-literature and wrong-transfer failure modes the
  pilot is meant to measure. In the current pilot, species with a single GBIF
  record sit alongside species with thousands.

## Fail-closed status of the stratification axes

Only `family` is populated from data. The trait-derived axes -
`native_status`, `growth_form`, `flower_conspicuousness` - are seeded `unknown`
on purpose: those labels are exactly what this pilot's human review supplies, so
they are never auto-decided. Two ways to enrich them before review:

1. `native_status` can be **derived from GBIF `establishmentMeans`** already in
   the collected occurrences (data, not a guess).
2. `growth_form` / `flower_conspicuousness` remain review targets, or can be
   pulled from a free source (e.g. GBIF/Wikidata) as candidates, still requiring
   review.

Once these axes carry values, rerunning the builder stratifies across all four
axes (native/introduced x woody/herb/vine/epiphyte x family x conspicuousness).

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
- `candidate_review_queue/trait_candidate_review_queue.csv`:
  45 review rows across 17 species; 23 direct reported rows first, 4 indirect
  source-check rows second, 17 proxy rows third, and 1 descriptive source-check
  row last. Allowed adjudication values are `accepted`, `rejected`, and
  `needs_source_check`.
- `globi_interaction_evidence/globi_interaction_evidence.csv`:
  21 explicit GloBI pollination interaction claims across 5 queried taxa; query
  errors: 0.

## Next step: adjudicate candidates

Work through `candidate_review_queue/trait_candidate_review_queue.csv` in
`review_priority` order:

1. Fill `adjudication_decision` for `P1_species_direct_reported` rows first.
2. Use `needs_source_check` for rows that require source-page inspection before
   acceptance or rejection.
3. Keep `P2_species_indirect_source_check` strict: accept only when the source
   truly applies to the accepted species, not just a genus or related taxon.
4. Treat `P3_proxy_candidate` rows as proxies only. They are not explicit
   reported values.
5. Promote only rows with `adjudication_decision=accepted` and nonblank
   `final_value` into curated outputs.

After editing the queue, validate the adjudication sheet before curation:

```
python -m island_v2.trait_candidate_review_queue validate \
    --review-queue-csv data/v2/staging/traits/validation_pilot/candidate_review_queue/trait_candidate_review_queue.csv \
    --output-dir data/v2/staging/traits/validation_pilot/candidate_review_queue
```

The validator fails on unknown decision values, accepted rows missing
`final_value` / reviewer provenance, accepted ontology/proxy rows whose
`final_value` is outside the controlled vocabulary, or rejected / source-check
rows that still carry a `final_value`. It also writes
`accepted_trait_candidates_for_curation.csv`, preserving source type, evidence
scope, raw excerpt, and proxy basis fields for accepted rows.

### Optional one-row review app

For a less error-prone review loop, use the Streamlit app:

```
pip install -e .
streamlit run src/island_v2/trait_review_app.py
```

The app opens the next undecided row in priority order, shows the source excerpt
and link, suggests controlled `final_value` choices, writes the decision back to
`trait_candidate_review_queue.csv`, and runs the same validator before saving.
It still records only human review decisions; it does not curate or infer trait
values by itself.

### Optional multimodal candidate queue

Image-derived visible traits and LLM-extracted ecology statements should stay in
separate lanes before review:

- image rows are only for visible M0 floral traits such as colour, form,
  symmetry, tube depth, size class, or inflorescence display;
- reproductive system, selfing capacity, pollen vector, and visitor evidence
  must come from source text excerpts, not from image inference;
- image files are not downloaded by the queue builder. Store an `image_url` and
  optional model metadata, then let the review app display one image at a time.

Create input templates:

```
python -m island_v2.trait_multimodal_candidates templates \
    --output-dir data/v2/staging/traits/validation_pilot/multimodal_candidate_inputs
```

After an image model, manual image pass, or LLM text-extraction pass fills those
templates, build a review queue:

```
python -m island_v2.trait_multimodal_candidates build \
    --visual-candidates-csv data/v2/staging/traits/validation_pilot/multimodal_candidate_inputs/visual_trait_candidates_template.csv \
    --reported-ecology-candidates-csv data/v2/staging/traits/validation_pilot/multimodal_candidate_inputs/reported_ecology_candidates_template.csv \
    --output-dir data/v2/staging/traits/validation_pilot/multimodal_candidate_review_queue
```

The resulting queue uses the same `accepted` / `rejected` /
`needs_source_check` decisions and the same Streamlit review app. Visual rows
are grouped as `P1_visual_trait_candidate`; LLM text rows are grouped as
`P1_reported_ecology_llm_source_check`.

### No-human-review machine evaluation

When manual adjudication is not part of the workflow, keep the outputs as
machine layers instead of writing them to curated accepted traits:

```
python -m island_v2.trait_machine_method_eval evaluate \
    --species-csv data/v2/staging/traits/validation_pilot/validation_pilot_species.csv \
    --review-queue-csv data/v2/staging/traits/validation_pilot/machine_method_eval/source_refresh/candidate_review_queue/trait_candidate_review_queue.csv \
    --globi-csv data/v2/staging/traits/validation_pilot/machine_method_eval/source_refresh/globi_interaction_evidence/globi_interaction_evidence.csv \
    --output-dir data/v2/staging/traits/validation_pilot/machine_method_eval/machine_outputs
```

This writes:

- `machine_trait_selected.csv`: main no-review trait rows selected from
  species-direct explicit text;
- `machine_trait_sensitivity.csv`: indirect, descriptive, visual, or proxy rows
  kept out of the main no-review layer;
- `machine_interaction_claims.csv`: source-backed interaction claims, not
  effectiveness or guild decisions;
- `method_comparison.md/json/csv`: method counts and recommended roles.

The current validation-pilot run recommends `web_reported_species_direct` as the
main no-review trait lane, `globi_interaction_claims` as the interaction-claim
lane, proxy rows as sensitivity only, and LLM excerpt extraction as the next lane
to build for reproductive mode, selfing capacity, and visitor statements.

## Source discovery

Literature source discovery (Stage 1, `island-v2-trait-source-discovery`) runs
against the public OpenAlex/Crossref/Unpaywall APIs. Those hosts are **blocked by
egress policy in the Claude Code session environment**, so discovery is run from
the GitHub Actions workflow (`discover-core-pilot-trait-sources.yml`), whose
runners have network access. The pipeline is fail-closed: an unreachable source
yields zero leads with a logged error, never a fabricated lead.
