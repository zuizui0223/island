# Free bulk trait source pilot

This experiment tests whether curated, openly accessible bulk plant-trait datasets can replace low-resolution genus/family/global fallback values before species-by-species web retrieval.

## Pilot source

AusTraits / eFLOWER_Dun_2022 is the first source inspected because its metadata describes roughly 2,000 taxa, 36 floral traits, and more than 29,000 trait values scored from taxonomic literature.

## Decision rule

A bulk source is promoted into the production cascade only when:

1. its licence and citation are recorded;
2. source trait definitions are reviewed explicitly;
3. taxon names can be reconciled to the island master without silent fuzzy matching;
4. source values are preserved in a lossless raw layer before any v2 categorical derivation;
5. direct/curated evidence replaces fallback values but never the reverse;
6. unmatched and unmapped records are reported, not guessed.

## Pilot metrics

The adapter reports:

- input records;
- unique source taxa;
- exact master-name matches;
- unmatched taxa;
- mapped and unmapped source traits;
- records eligible to replace each existing fill tier.

## Initial mapping target

The pilot intentionally starts with source traits that can be retained without lossy biological interpretation:

- flower_length -> flower_length_mm
- flower_diameter -> flower_diameter_mm
- flower_perianth_fusion -> flower_perianth_fusion_fraction
- flower_structural_sex_type -> raw structural-sex evidence

The first three are stored as continuous raw traits. Existing categorical traits such as `flower_size_class` or `floral_form` should be derived later using preregistered thresholds or explicit rules, rather than overwriting the continuous measurements.

## Scope

This PR is an experiment, not a claim that AusTraits alone provides global coverage. The intended workflow is:

open bulk sources -> source-specific adapters -> taxon reconciliation -> evidence ledger -> fallback replacement audit -> only then species-by-species retrieval for remaining high-value gaps.
