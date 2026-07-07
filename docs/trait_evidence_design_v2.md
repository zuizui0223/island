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

## Stage 1 — trait-group-targeted retrieval  ✅ implemented

`src/island_v2/trait_source_discovery.py`. Retrieval is **species × trait-group**
(not species-only), so high-quality but trait-irrelevant papers
(micropropagation, tissue culture, nutrient/bromelain content, nanoparticles,
general genomics/pathology) are demoted rather than surfaced.

Three query groups: A floral morphology/colour, B pollination/pollen vector,
C reproductive assurance. Each lead now carries `query_trait_group`,
`query_template`, `title_relevance_score`, `abstract_relevance_score`,
`trait_keywords_matched`, `likely_evidence_type`, plus the existing
`provisional_source_reliability_hint`. Leads are ranked within (taxon, group) by
relevance then source grade; a `title_relevance_score == 0` lead can never head a
taxon's trait review. Source grade (A/B/C/D) is a **source-quality** hint;
relevance is a separate **trait-relevance** axis — both are reported.

## Stage 2 — M0 / M1 / M2 trait layers  ⏳ next

Separate the trait ontology into three explicit layers so acquisition and
analysis-eligibility rules differ per layer:

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

## Stage 3 — pollinator functional evidence (long format)  ⏳ next

The single categorical `mixed_or_generalist` is too coarse to test alternative
pollination-channel hypotheses. Keep a **long-format** evidence table; do not
collapse to `mixed_or_generalist` at the raw stage:

`accepted_species, pollinator_guild, evidence_type {floral_visitor,
pollen_contact, pollen_deposition, fruit_or_seed_set_contribution,
exclusion_experiment}, effectiveness_class, source_citation, source_url,
study_location, season, wild_or_cultivated, review_status`.

A functional-replacement index may be derived later, but only from this raw
long-format table.

## Stage 4 — measurement-first size / tube depth  ⏳ next

Do not store only `small/medium/large` or `shallow/deep` classes (reviewer- and
taxon-dependent). Store the measurement first and derive the class by rule:

`raw_measurement, raw_unit, measurement_structure, measurement_source_text,
derived_class, classification_rule_version`.

## Stage 5 — structured flower colour  ⏳ next

`flower_primary_color` alone is insufficient. Capture:
`raw_colour_description, primary_colour, secondary_colour,
within_species_variation, geographic_scope, wild_or_cultivated,
colour_evidence_type`. Cultivar, geographic variation, and anthesis colour change
must not be merged into a species-level colour.

## Stage 6 — stratified gold-standard validation pilot  ⏳ next

Do **not** select the first pilot by GBIF record count (biases toward
human-associated / cultivated / easily observed species). Stratify a 30–50
species pilot by at least: native vs introduced candidate; woody / herb / vine /
epiphyte; major angiosperm families; conspicuous / inconspicuous /
wind-pollination-candidate flowers.

Compare human review vs automatic candidates and evaluate:
`species_direct_evidence_proportion`, `irrelevant_literature_rate`,
`source_to_trait_relevance_rate`, `trait_specific_unresolved_rate`,
`wrong_hierarchical_transfer_rate`. **`irrelevant_literature_rate` (how much
off-target literature is pulled in) is a primary metric — not just agreement.**
