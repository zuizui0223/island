# Trait evidence & source-discovery design (v2)

The trait workflow is a **source-lead / literature-scout and evidence-review
preprocessor**, not a trait extractor. It never auto-decides a trait value,
island establishment, Bombus applicability, or analysis inclusion. This document
records the target design and its staged status.

## Retained fail-closed principles (do not regress)

- A source lead is **not** trait evidence.
- Blank / `unresolved` is **never** biological absence.
- No analysis release until **accepted taxon scope AND accepted island
  establishment**.
- Bombus applicability is **never** decided on the trait-extraction side.
- Model / retrieval output stays a **candidate database**; human review required.

## Stage 1 — three-lane scouting  (by source type, feeding different trait layers)

Retrieval is split into three lanes because a single source type cannot serve every
trait layer, and obscure island endemics with little primary literature still have
descriptive coverage:

- **Stage 1A — literature scout** (`trait_source_discovery.py`; OpenAlex / Crossref /
  Unpaywall) → **M1/M2 species-direct evidence**.  ✅ implemented (v2.1)
- **Stage 1B — descriptive / flora scout** (`trait_descriptive_scout.py`,
  `island-v2-trait-descriptive-scout scout`; GBIF species descriptions first, then
  POWO/WCVP → flora & horticulture DBs → Wikidata/Wikipedia) → **M0 floral-phenotype
  candidate extraction**.  ✅ implemented (GBIF lane); further sources incremental.
  It retrieves free-text descriptions and emits an M0 candidate **only where a
  controlled-vocabulary term appears verbatim** (`config/m0_descriptive_keywords.yml`),
  each with a source excerpt. Output is long-format
  (`trait_name` / `provisional_candidate_value`, `review_status=unreviewed`), never a
  finalized wide column; multiple colours yield multiple candidates (nothing collapsed
  to a single primary); a blank is never biological absence.
- **Stage 1C — interaction evidence** (`interaction_evidence.py`; GloBI) → **M2
  subset** (explicit flower-visit / pollination claims).  ✅ implemented

### Stage 1A v2.1 — trait-group retrieval + taxon relevance tiers

Retrieval is **species × trait-group** (A floral morphology/colour, B pollination/
pollen vector, C reproductive assurance) with the binomial quoted as an **exact
phrase**, so retrieval biases toward on-species work rather than anything sharing the
genus or the trait keywords. Each lead carries `query_trait_group`, `query_template`,
`taxon_relevance_score`, `taxon_match_kind`, `taxon_relevance_tier`,
`title_relevance_score`, `abstract_relevance_score`, `trait_keywords_matched`,
`likely_evidence_type`, and `provisional_source_reliability_hint`.

**Taxon relevance tiers.** The first validation run showed trait relevance is high
but does not mean a lead is about the *target species* (767 leads, 642 trait-relevant,
but only ~24 named the species). `score_taxon_relevance` checks whether the queried
binomial (or its genus/epithet) appears in the title/abstract, and the score maps to a
tier: **S** (species named with confidence, score ≥ 3), **A** (genus-only / flora
rescue, score 1–2), **B** (species never named, background). Leads rank by taxon
relevance first, then trait relevance, then source grade. The report includes
`n_taxon_matched_species` / `n_zero_taxon_species` (how many pilot species got any
on-species lead at all).

### Tier-routed compression (`trait_lead_packet.py`, `island-v2-trait-lead-packet compress`)

Leads are routed by tier, not merged: **Tier S** → the bounded review packet
(round-robin per species/trait-group, `config/lead_review_packet.yml`, default
150–300); **Tier A** → a quarantined genus/flora rescue CSV; **Tier B** → a retained
background-literature CSV that is never reviewed. **Only the Tier-S packet** feeds a
human `irrelevant_literature_rate` evaluation (Stage 6).

## Stage 2 — M0 / M1 / M2 trait layers  ✅ implemented

`config/trait_layers.yml` + `src/island_v2/trait_layers.py`
(`island-v2-trait-layers annotate`). Annotates a reviewed trait-candidate table
with `trait_layer` (M0/M1/M2) and `analysis_track` (`main_direct` /
`expanded_sensitivity` / `excluded`). It governs analysis eligibility only and
never decides a trait value. M1 (and M2 guild) main-analysis values are
restricted to species-direct evidence; any hierarchical (genus/family) inference
is always routed to `expanded_sensitivity`, never a primary result. The layer
membership is exactly:

- **M0 floral phenotype**: `flower_primary_color`, `floral_symmetry`,
  `floral_form`, `tube_depth_class`, `flower_size_class`, `inflorescence_display`.
  May be candidated broadly from floras, descriptions, herbarium and
  botanic-garden pages.
- **M1 reproductive assurance / pollen vector**: `pollen_vector_mode`,
  `self_incompatibility`, `autonomous_selfing_capacity`, `mating_system`,
  `herkogamy`, `dichogamy`, `cleistogamy`. **Main-analysis values restricted to
  species-direct evidence.**
- **M2 pollinator functional evidence**: pollinator guild, evidence type,
  effectiveness evidence.

Genus/family inference stays allowed but **isolated as `hierarchical_inference`
for an expanded sensitivity track only** — never in the main analysis — for:
`self_incompatibility`, `autonomous_selfing_capacity`, `mating_system`,
`herkogamy`, `dichogamy`, `pollination_functional_guild`.

## Stage 3 — pollinator functional evidence (long format)  ✅ implemented

`config/pollinator_evidence_schema.yml` + `src/island_v2/pollinator_evidence.py`
(`island-v2-pollinator-evidence validate|index`). Validates the long-format
table and derives a functional-replacement index from **accepted** rows only,
without mutating the raw table. A collapsed `mixed_or_generalist` column is
rejected. A functional-replacement candidate has 2+ guilds each confirmed as an
effective pollinator — never a single generalist label.

The single categorical `mixed_or_generalist` is too coarse to test alternative
pollination-channel hypotheses. The **long-format** evidence table is not
collapsed at the raw stage:

`accepted_species, pollinator_guild, evidence_type {floral_visitor,
pollen_contact, pollen_deposition, fruit_or_seed_set_contribution,
exclusion_experiment}, effectiveness_class, source_citation, source_url,
study_location, season, wild_or_cultivated, review_status`.

A functional-replacement index may be derived later, but only from this raw
long-format table.

## Stage 4 — measurement-first size / tube depth  ✅ implemented

`config/measurement_classification.yml` + `src/island_v2/trait_measurements.py`
(`island-v2-trait-measurements annotate`). Stores the measurement first and
derives the class by a **versioned** rule (`classification_rule_version`), so the
class is reproducible and reviewer-independent; the raw fields are never
overwritten. A range (e.g. `10-20`) uses its midpoint; an unknown unit or an
unparseable value yields a blank (unresolved) class, never a guess.

`raw_measurement, raw_unit, measurement_structure, measurement_source_text,
derived_class, classification_rule_version`.

## Stage 5 — structured flower colour  ✅ implemented

`config/colour_schema.yml` + `src/island_v2/trait_colour.py`
(`island-v2-trait-colour validate|annotate`). Captures the structured fields,
validates the controlled vocabulary, and appends
`species_level_primary_colour_admissible`. Cultivar, geographic variation, and
anthesis colour change are **retained but flagged non-admissible**, never merged
into a species-level colour.

`raw_colour_description, primary_colour, secondary_colour,
within_species_variation, geographic_scope, wild_or_cultivated,
colour_evidence_type`.

## Stage 6 — stratified gold-standard validation pilot  ✅ implemented

`config/pilot_stratification.yml` + `src/island_v2/trait_pilot_eval.py`
(`island-v2-trait-pilot-eval build|evaluate`). `build` selects a 30–50 species
pilot by spreading selections **evenly across strata** (round-robin, largest
strata first, deterministic) — never by GBIF record count, which biases toward
human-associated / cultivated / easily observed species. Strata combine: native
vs introduced; woody / herb / vine / epiphyte; major angiosperm families;
conspicuous / inconspicuous / wind-pollination-candidate flowers.

`evaluate` compares human review vs automatic candidates and computes:
`species_direct_evidence_proportion`, `irrelevant_literature_rate`,
`source_to_trait_relevance_rate`, `trait_specific_unresolved_rate`,
`wrong_hierarchical_transfer_rate`. **`irrelevant_literature_rate` (how much
off-target literature is pulled in) is a primary metric — not just agreement.**
Each rate counts only explicit judgements; a blank is never a judgement and a
zero-denominator rate is reported as `null` (undefined), never fabricated.
