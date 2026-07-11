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
