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

## Evidence policy  -  reported vs proxy (cross-layer, M0-M4)

Trait evidence is captured in **two strictly separated classes**, kept in
**separate columns**, never merged. A `candidate_class` field records which:

### 1. `reported`  -  an explicit statement exists

For any trait  -  including reproductive system, self-incompatibility /
self-compatibility (SI/SC), and pollination syndrome  -  sources are **not limited to
the primary literature**. Any source with an **explicit statement** is used and kept:
peer-reviewed literature, Floras, taxonomic/biodiversity databases, Wikipedia,
Wikidata, botanic-garden pages, horticulture sites, and other web descriptions. A
web-derived explicit statement is **not discarded**; it is recorded with its source
type so it can be weighted at review time.

Every `reported` candidate retains, at minimum:
`candidate_class=reported`, `source_type`, `source_url`, `raw_description` (the
verbatim excerpt), `evidence_scope` (species_direct / species_indirect / genus / ...),
`wild_or_cultivated`, `review_status`. These are recorded as a
**`reported_trait_candidate`**.

### 2. `proxy`  -  no explicit statement; inferred to fill a gap

When no explicit statement exists, a value **may** be suggested from floral
phenotype, taxonomic tendency, life history, distribution, or insularity  -  but only
as a **proxy**, never a decided value. Proxy traits carry a distinct name
(`*_proxy`), e.g. `floral_syndrome_proxy`, `reproductive_assurance_proxy`,
`compatibility_system_proxy`. A proxy is a phenotype/life-history-based candidate for
gap-filling; it is **not real data**.

### Distinctions that must never be conflated

| reported (explicit) | proxy (inferred) |
| --- | --- |
| `self_compatible` reported | `likely_self_compatible_proxy` |
| `Bombus` visitor reported | `large_bee_or_Bombus_like_floral_phenotype_proxy` |
| `autonomous_selfing` reported | `reproductive_assurance_like_proxy` |

### Design stance

Keep the reliability model **simple and coverage-first**  -  do not over-engineer
confidence scoring. The one non-negotiable is provenance: always retain the
**raw text** and **source type** so a reviewer can later assess each candidate.
In **every trait evidence layer and every analysis model (M0-M4)**, `reported`
evidence and `proxy` are **separate columns**; proxies never substitute for reported
values and never enter a main analysis in place of species-direct reported evidence.

> **Two different "M" namespaces  -  do not conflate:**
> - **Trait evidence layers** (data organization, Stage 2): *floral phenotype*,
>   *reproductive assurance / pollen vector*, *pollinator functional evidence*.
>   (`config/trait_layers.yml` currently keys these M0/M1/M2; treat them as evidence
>   *layers*, not the models below.)
> - **Analysis model ladder M0-M4** (the hypothesis sequence, Stage 7): baseline ->
>   Bombus association -> mediation -> joint -> compensation. Defined in
>   "Analysis model ladder" below and `config/analysis_models.yml`.
>
> The reported-vs-proxy policy applies to both namespaces.

## Stage 1  -  three-lane scouting  (by source type, feeding different trait layers)

Retrieval is split into three lanes because a single source type cannot serve every
trait layer, and obscure island endemics with little primary literature still have
descriptive coverage:

- **Stage 1A  -  literature scout** (`trait_source_discovery.py`; OpenAlex / Crossref /
  Unpaywall) -> **M1/M2 species-direct evidence**.  implemented (v2.1)
- **Stage 1B  -  descriptive / flora scout** (`trait_descriptive_scout.py`,
  `island-v2-trait-descriptive-scout scout`; GBIF species descriptions first, then
  POWO/WCVP -> flora & horticulture DBs -> Wikidata/Wikipedia) -> **M0 floral-phenotype
  candidate extraction**.  implemented (GBIF lane); further sources incremental.
  It retrieves free-text descriptions and emits an M0 candidate **only where a
  controlled-vocabulary term appears verbatim** (`config/m0_descriptive_keywords.yml`),
  each with a source excerpt. Output is long-format
  (`trait_name` / `provisional_candidate_value`, `review_status=unreviewed`), never a
  finalized wide column; multiple colours yield multiple candidates (nothing collapsed
  to a single primary); a blank is never biological absence.
  - **1B web-reported sub-lane** (`trait_web_reported_scout.py`,
    `island-v2-trait-web-reported-scout scout`; Wikipedia + Wikidata): the first pilot
    showed GBIF descriptions are nearly empty for obscure endemics (1 M0 candidate /
    50 species), so this sub-lane resolves each species to its Wikidata item and English
    Wikipedia article and emits **`reported`** M0/M1 candidates wherever a controlled
    term appears verbatim (`config/m0_descriptive_keywords.yml` +
    `config/m1_reported_keywords.yml`). Web/DB explicit statements are kept with full
    provenance (`candidate_class=reported`, `source_type`, `source_url`,
    `raw_description`, `evidence_scope`, `wild_or_cultivated`, `review_status`).  implemented
- **Stage 1C  -  interaction evidence** (`interaction_evidence.py`; GloBI) -> **M2
  subset** (explicit flower-visit / pollination claims).  implemented (workflow pending)

- **Stage 1D  -  optional multimodal staging**
  (`trait_multimodal_candidates.py`,
  `island-v2-trait-multimodal-candidates templates|build`) -> **review queue only**.
  This lane keeps image and LLM evidence separate. Image inputs may propose only
  visible M0 floral traits (`flower_primary_color`, `floral_form`,
  `floral_symmetry`, `tube_depth_class`, `flower_size_class`,
  `inflorescence_display`) and store `image_url` rather than downloading image files.
  LLM/text inputs may propose M1/M2 ecology traits only when a source excerpt is
  preserved. Neither lane curates a value; both are converted into the same
  human-adjudicated queue as the text scouts.

- **Stage 1E  -  no-human-review machine method evaluation**
  (`trait_machine_method_eval.py`,
  `island-v2-trait-machine-method-eval evaluate`) -> **machine layer outputs**.
  This path is for runs where manual review is intentionally skipped. It maps
  controlled-vocabulary values, selects only species-direct explicit-text rows
  for the main no-review trait layer, keeps indirect/descriptive/proxy rows in
  sensitivity layers, and preserves GloBI records as interaction claims rather
  than effectiveness or guild decisions.

- **Stage 1F  -  no-human-review interaction-guild machine layer**
  (`interaction_guild_machine.py`,
  `island-v2-interaction-guild-machine build`) maps source-backed GloBI partner
  taxon names to configured pollinator guild candidates. Unmapped partners are
  written to holdouts. The schema-compatible long-format output remains
  `review_status=unreviewed`; a separate machine index records multi-guild
  functional-replacement signals for no-review sensitivity analyses.

### Stage 1A v2.1  -  trait-group retrieval + taxon relevance tiers

Retrieval is **species x trait-group** (A floral morphology/colour, B pollination/
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
tier: **S** (species named with confidence, score >= 3), **A** (genus-only / flora
rescue, score 1-2), **B** (species never named, background). Leads rank by taxon
relevance first, then trait relevance, then source grade. The report includes
`n_taxon_matched_species` / `n_zero_taxon_species` (how many pilot species got any
on-species lead at all).

### Tier-routed compression (`trait_lead_packet.py`, `island-v2-trait-lead-packet compress`)

Leads are routed by tier, not merged: **Tier S** -> the bounded review packet
(round-robin per species/trait-group, `config/lead_review_packet.yml`, default
150-300); **Tier A** -> a quarantined genus/flora rescue CSV; **Tier B** -> a retained
background-literature CSV that is never reviewed. **Only the Tier-S packet** feeds a
human `irrelevant_literature_rate` evaluation (Stage 6).

## Stage 2  -  M0 / M1 / M2 trait layers  implemented

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
for an expanded sensitivity track only**  -  never in the main analysis  -  for:
`self_incompatibility`, `autonomous_selfing_capacity`, `mating_system`,
`herkogamy`, `dichogamy`, `pollination_functional_guild`.

## Stage 3  -  pollinator functional evidence (long format)  implemented

`config/pollinator_evidence_schema.yml` + `src/island_v2/pollinator_evidence.py`
(`island-v2-pollinator-evidence validate|index`). Validates the long-format
table and derives a functional-replacement index from **accepted** rows only,
without mutating the raw table. A collapsed `mixed_or_generalist` column is
rejected. A functional-replacement candidate has 2+ guilds each confirmed as an
effective pollinator  -  never a single generalist label.

The single categorical `mixed_or_generalist` is too coarse to test alternative
pollination-channel hypotheses. The **long-format** evidence table is not
collapsed at the raw stage:

`accepted_species, pollinator_guild, evidence_type {floral_visitor,
pollen_contact, pollen_deposition, fruit_or_seed_set_contribution,
exclusion_experiment}, effectiveness_class, source_citation, source_url,
study_location, season, wild_or_cultivated, review_status`.

A functional-replacement index may be derived later, but only from this raw
long-format table.

## Stage 4  -  measurement-first size / tube depth  implemented

`config/measurement_classification.yml` + `src/island_v2/trait_measurements.py`
(`island-v2-trait-measurements annotate`). Stores the measurement first and
derives the class by a **versioned** rule (`classification_rule_version`), so the
class is reproducible and reviewer-independent; the raw fields are never
overwritten. A range (e.g. `10-20`) uses its midpoint; an unknown unit or an
unparseable value yields a blank (unresolved) class, never a guess.

`raw_measurement, raw_unit, measurement_structure, measurement_source_text,
derived_class, classification_rule_version`.

## Stage 5  -  structured flower colour  implemented

`config/colour_schema.yml` + `src/island_v2/trait_colour.py`
(`island-v2-trait-colour validate|annotate`). Captures the structured fields,
validates the controlled vocabulary, and appends
`species_level_primary_colour_admissible`. Cultivar, geographic variation, and
anthesis colour change are **retained but flagged non-admissible**, never merged
into a species-level colour.

`raw_colour_description, primary_colour, secondary_colour,
within_species_variation, geographic_scope, wild_or_cultivated,
colour_evidence_type`.

## Stage 6  -  stratified gold-standard validation pilot  implemented

`config/pilot_stratification.yml` + `src/island_v2/trait_pilot_eval.py`
(`island-v2-trait-pilot-eval build|evaluate`). `build` selects a 30-50 species
pilot by spreading selections **evenly across strata** (round-robin, largest
strata first, deterministic)  -  never by GBIF record count, which biases toward
human-associated / cultivated / easily observed species. Strata combine: native
vs introduced; woody / herb / vine / epiphyte; major angiosperm families;
conspicuous / inconspicuous / wind-pollination-candidate flowers.

`evaluate` compares human review vs automatic candidates and computes:
`species_direct_evidence_proportion`, `irrelevant_literature_rate`,
`source_to_trait_relevance_rate`, `trait_specific_unresolved_rate`,
`wrong_hierarchical_transfer_rate`. **`irrelevant_literature_rate` (how much
off-target literature is pulled in) is a primary metric  -  not just agreement.**
Each rate counts only explicit judgements; a blank is never a judgement and a
zero-denominator rate is reported as `null` (undefined), never fabricated.

## Analysis model ladder (M0-M4)  (Stage 7  -  analysis stage, gated)

The confirmatory/exploratory **model ladder**. These are analysis models, distinct
from the trait evidence layers above. They are **gated**: no model is fitted or
released until accepted taxon scope AND island establishment AND the region-scoped
Bombus applicability freeze hold (see fail-closed principles). The outcome throughout
is a **floral-phenotype-type proxy** (e.g. `large_bee_or_Bombus_like_floral_phenotype_proxy`),
never an assertion of the actual pollinator. `reported` evidence and `proxy` stay in
separate columns in every model.

- **M0  -  baseline model.** Covariates only: geography, source pool, lineage
  composition, island area, isolation, climate, establishment, observation process.
  Question: can the baseline alone explain the composition of floral-phenotype,
  reproductive-assurance, and pollination-channel **proxies**?

- **M1  -  Bombus-channel association.** Add **Bombus channel state** to the M0
  covariates and test whether it explains the large-bee/Bombus-like floral phenotype
  proxy outcome. This is a comparison of *Bombus-associated floral phenotype*, **not**
  confirmation of the real pollinator.

- **M2  -  reproductive-assurance-mediated route.** Test whether Bombus channel state
  relates to the floral-phenotype outcome **via** reproductive assurance /
  compatibility / selfing  -  using either `reported` evidence or `proxy` for the
  mediator. (Mediation path: Bombus channel -> reproductive assurance -> phenotype.)

- **M3  -  joint direct-plus-mediated route.** A single model carrying **both** paths
  simultaneously:
  - Bombus channel -> floral phenotype proxy (direct), and
  - Bombus channel -> reproductive assurance (reported/proxy) -> floral phenotype proxy
    (mediated).

- **M4  -  alternative functional replacement / compensation route.** On islands where
  the Bombus channel is weak, absent, or not applicable, test whether **other channel
  types compensate**: `open_or_generalist_insect_like`, `small_bee_or_fly_like`,
  `butterfly_or_moth_like`, `bird_or_bat_like`, `wind_like`. The main term is the
  **Bombus channel state x alternative_functional_replacement_proxy interaction** -
  i.e. whether the Bombus-like phenotype reduction on Bombus-absent islands **weakens,
  disappears, or reverses** because another channel is present.

Guardrails (all models): keep `reported` and `proxy` in separate columns; a proxy is a
candidate phenotype-type outcome built from web/DB/literature statements or trait
information, **not real data**; never label a taxon "Bombus-pollinated"  -  use
`large_bee_or_Bombus_like_floral_phenotype_proxy`.
