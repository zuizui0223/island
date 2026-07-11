# Category-first trait collection for the v1 analysis

## Principle

Trait collection and statistical coding are separate steps.

The LLM collects the following nine columns without reducing them to binary
values at acquisition time:

```text
species,flower_color,flower_shape,pollination_guild,pollination_notes,mating_system,self_incompatibility,evidence_type,confidence
```

The v1 analysis then derives coarse variables from the validated rows.

## v1 derived categories

### Flower colour

- `plain`: white, cream/ivory, green, brown, dull or inconspicuous descriptions
- `conspicuous`: blue/purple, red/pink or yellow/orange descriptions
- `unknown`: no recognized colour description

Where a description contains both plain and conspicuous colours, the row is
classified as conspicuous because a conspicuous signal is present.

### Floral form

- `generalized`: actinomorphic/radial, open/bowl/dish/cup, brush/puff,
  composite head/capitulum, or reduced wind-type descriptions
- `specialized`: zygomorphic/bilateral/bilabiate, tubular, salverform,
  funnel/trumpet, urceolate, papilionaceous or spurred descriptions
- `unknown`: no recognized form description

When both generalized and specialized terms occur, specialized takes precedence
because the restrictive tube or bilateral architecture is biologically relevant
to pollinator access.

### Self-compatibility

- `self_compatible`: `SC`, `likely_SC`
- `self_incompatible`: `SI`, `likely_SI`
- `unknown`: `unknown`

### Pollination filter

- animal-pollinated: bees, bumblebees, flies, butterflies, moths or birds
- non-animal: wind or self
- mixed: retained separately
- unknown: retained as missing evidence

### Mating-system summary

- outcrossing: `obligate_outcrossing`
- mixed: `mixed_mating`
- selfing: `mainly_selfing`, `obligate_selfing`
- unknown: `unknown`

## Evidence layers

Primary analysis eligibility is assigned separately for colour, form and
self-compatibility. A row is primary evidence when:

1. evidence includes `field_study`, `review` or `flora`; and
2. confidence is `high` or `medium`; and
3. the relevant trait is not unknown.

`horticulture`, inference-only and low-confidence rows are retained for secondary
or sensitivity analyses. They are not silently discarded or promoted.

## Prompt packets

Generate deterministic prompt packets from the global island plant master:

```bash
island-v2-v1-category-traits prepare \
  --species-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --output-dir /tmp/v1_category_packets \
  --species-column accepted_species \
  --batch-size 25 \
  --max-species 1000
```

Each packet contains:

- `prompt.txt`
- `expected_species.csv`
- `result_template.csv`

The manual GitHub Actions workflow `build_v1_category_prompt_packets` creates the
same packet structure as a downloadable artifact.

## Free automatic search baseline

This is a legacy-comparison bridge, not the active v2 analysis contract. The
same run also writes `v2_reported_candidates.csv`; active acquisition consumes
that long-format table and treats the nine-column v1 files as diagnostics.

The legacy nine-column table can also be populated without an LLM or paid API.
Use bulk candidates first, then query public text only for the remaining gaps:

```bash
island-v2-v1-category-search search \
  --species-csv data/v2/staging/traits/validation_pilot/validation_pilot_species.csv \
  --candidate-csv data/v2/staging/traits/eol_traitbank/trait_candidates_campaign.csv.gz \
  --output-dir /tmp/v1_category_search \
  --max-taxa 25
```

This bounded command searches GBIF species descriptions, Wikimedia text,
World Flora Online, OpenAlex title/abstract metadata, NCSU Plant Toolbox, New
Zealand Plant Conservation Network, PFAF and Useful Tropical Plants. It writes the exact v1 table alongside
`source_texts.csv` and `v1_category_evidence.csv`, which retain every accepted
match, source URL and excerpt. Only species-direct explicit statements are
collapsed. Missing evidence remains `unknown`; the command does not call an LLM
and does not treat a failed search as biological absence.

Retrieval is a gap-filling cascade. Bulk candidates and flora/monograph sources
are screened first, followed by literature, biodiversity databases and
encyclopedias, then horticulture pages. A lower tier is queried only for traits
that remain missing, and it cannot replace a value already found at a higher
tier. The long evidence table keeps all candidates with `source_tier` values
`A_flora`, `A_literature`, `B_database`, `B_encyclopedia`, `C_horticulture` or
`D_inference`, so downstream analyses can choose strict or coverage-first data.

Cached text can be re-extracted after rule improvements without repeating any
network requests:

```bash
island-v2-v1-category-search extract-text \
  --source-csv /tmp/v1_category_search/source_texts.csv \
  --species-csv data/v2/staging/traits/validation_pilot/validation_pilot_species.csv \
  --output-dir /tmp/v1_category_search_reextract
```

## Direct evidence and likely proxies

The legacy nine-column file remains the direct-evidence contract. Descriptions
such as `pollinated by birds` can populate `pollination_guild`, but weaker
statements such as `attracts hummingbirds` cannot. All proxy results are written
separately and every inferred value starts with `likely`:

- `v1_likely_traits.csv`: one inference per row with basis fields, values, URLs,
  rule ID and confidence;
- `v1_likely_traits_wide.csv`: one species per row with the likely values;
- `v1_priority_traits.csv`: direct-first analysis table with a `direct`, `likely`
  or `unknown` mode for pollination, selfing and self-incompatibility.

The bounded proxy rules are intentionally conservative:

- red, pink or orange plus a tubular/funnel/trumpet flower gives
  `likely_birds_or_butterflies`;
- blue, purple or yellow plus a restrictive bee-associated form gives
  `likely_bees_or_bumblebees`;
- an explicit visitor/attraction statement gives `likely_<visitor>` when no
  direct pollination statement exists;
- autonomous selfing, autogamy, autofertility or cleistogamy gives
  `likely_selfing_capable` and, only when compatibility is otherwise unknown,
  `likely_SC`.

Self-compatibility alone is not converted to mainly selfing, and a floral
syndrome is not treated as proof of an effective pollinator. Island-level
Bombus presence or absence is never an input to these rules. It remains an
independent downstream explanatory variable, avoiding circular support for a
bird/butterfly replacement hypothesis.

For roughly 100,000 species, run bulk sources and scientific-name joins first,
then extract cached flora/encyclopedia text, and only query websites for the
remaining flower-colour and flower-shape gaps. The deterministic likely pass is
local and is run after direct extraction; it does not add network requests.

## Validate and recode a returned CSV

```bash
island-v2-v1-category-traits ingest \
  --result-csv result.csv \
  --expected-species-csv expected_species.csv \
  --output-csv validated_v1_traits.csv
```

The command rejects column-order changes, invalid category labels, duplicate or
missing species, blank evidence types and non-low confidence on fully unknown
rows. It writes the raw nine columns plus the derived v1 categories and a
coverage manifest.
