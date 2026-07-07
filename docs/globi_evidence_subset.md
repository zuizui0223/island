# GloBI interaction evidence subset

## Purpose

`island-v2-interaction-evidence` retrieves a bounded list of explicit GloBI
interaction claims with one of four relations:

- `pollinates`
- `pollinated by`
- `visits flowers of`
- `flowers visited by`

The output is a long-format, source-backed **M2 evidence subset** for evaluating
alternative pollinator functional replacement. Each record preserves source and
target names, the relation, study and dataset citations when supplied by GloBI,
and available locality/date fields.

## What it does not infer

The collector never:

- marks a species as non-animal-pollinated because GloBI has no record;
- treats a flower-visit record as evidence of pollen deposition, seed set, or
  effective pollination;
- assigns a functional guild from morphology, colour, or taxonomic family;
- changes island establishment, Bombus applicability, or analytical inclusion.

`effectiveness_class` is therefore always `not_assessed` at retrieval. A later
review step may distinguish visitor evidence, pollen-contact evidence,
pollen-deposition evidence, fruit/seed-set evidence, and exclusion experiments.

## Bounded use

GloBI's own guidance warns that large custom taxon-list searches can be unstable.
The collector limits one run to 100 predeclared taxa and queries both source and
target orientations. It is intended for selected evidence-rich panels, not for
filling the whole global flora.

## Usage

```bash
island-v2-interaction-evidence collect \
  --taxa-csv path/to/predeclared_taxa.csv \
  --output-dir data/v2/staging/globi_pilot \
  --max-taxa 100
```

The input must contain `accepted_species`. The output remains unreviewed until
its taxonomic match, geographic relevance, relation direction, and study
provenance are checked.
