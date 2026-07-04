# v2 taxonomic scope: raw vascular acquisition, primary angiosperm analysis

## Decision

The primary v2 analysis universe is the **reviewed angiosperm flora of each
island**, not all vascular plants.

The present GBIF campaign queries `Tracheophyta` deliberately as a broad raw
occurrence universe. Those records are useful for spatial assignment, source
structure, historical coverage, and broad source-pool context. They are not a
valid shared denominator for floral phenotype, pollen-vector mode, or
reproductive-assurance composition.

## Why the denominator must be angiosperms

The primary v2 responses concern floral signal, floral architecture,
pollen-vector mode, and reproductive assurance. Non-seed vascular plants do
not have flowers or an angiosperm-compatible pollen-vector system. Gymnosperms
also do not supply a comparable floral-trait observation unit. Including either
in a "whole flora" denominator would mechanically change the apparent wind/
biotic composition between islands for reasons outside the stated hypothesis.

Use these terms precisely:

```text
raw vascular-plant occurrence universe
whole angiosperm flora              # primary composition denominator
animal-pollinated angiosperm subset # conditional Bombus-domain subset
```

Avoid calling the primary denominator simply "whole flora."

## Operational sequence

```text
GBIF Tracheophyta raw occurrences
  -> exact original-island assignment
  -> raw observation-process audit
  -> species-level taxonomic-scope review
  -> accepted-angiosperm island species universe
  -> angiosperm coverage gate
  -> pollen-vector / floral / RA trait evidence and Phase 1 attrition
```

`island-v2-taxon-scope build` consumes the collector's
`island_species_occurrences.csv` and a reviewed scope decision table. It writes:

```text
island_species_taxonomic_scope.csv
island_angiosperm_species.csv
island_angiosperm_coverage.csv
taxon_scope_review_queue.csv
taxon_scope_summary.json
```

For confirmatory floral analyses, run the observation diagnostics with the
reviewed angiosperm coverage table:

```text
island-v2-diagnostics run \
  --effort-csv <island_observation_effort.csv> \
  --angiosperm-coverage-csv <island_angiosperm_coverage.csv> \
  --output-dir <diagnostics directory>
```

Without `--angiosperm-coverage-csv`, diagnostics are explicitly labelled
`raw_tracheophyta_screening_only`; they can describe collection effort but must
not determine the confirmatory floral-analysis island set.

## What remains informative outside the primary denominator

Raw vascular richness and the proportions of reviewed angiosperms,
gymnosperms, non-seed vascular taxa, and unresolved taxa can be retained as
source-pool/lineage-composition descriptions or M0 contextual diagnostics. They
are not direct floral outcomes.
