# Wave57 reproductive-assurance acquisition

## Verified Wave56 result

Run 34029906166 succeeded operationally. Its final artifact contains 94 raw records for 49 taxon labels, not 94 adopted reproductive records: 93 are off the reproductive axis; the sole reproductive row is genus-only `Cymbalaria`, incorrectly tagged species-direct/High from an urban-community generalization. It is rejected. None is promoted. Of the previous 300 target rows, 21 fail the conservative species-name/axis/public-universe gate. The frozen analysis denominator is NOT reduced by this discovery.

Wave56's four manually reviewed CSVs (17 rows / 13 species) are separate from that machine output and are inherited unchanged on this branch. Their current-canonical collision audit remains incomplete.

## New source-reviewed packet

Four trait rows describe two public-Wave52-unresolved species: `Vandellia micrantha` and `Vandellia setulosa`, each with SC and autonomous self-pollination observations. These are source-reviewed candidates, not four new coverage cells. The source is Higashiyama et al. (2006), DOI 10.1104/pp.106.083832, Materials and Methods. The original names, scope and shared source lineage are retained. SC is not used to infer autonomy: the source has separate explicit statements for each.

The current POWO name for the frozen-master `Vandellia micrantha` is `Torenia micrantha`; both it and `Lindernia micrantha` are documented homotypic synonyms. Do not count aliases as multiple cells or use this record as a Vandellia genus-rule vote. All four rows have promotion_allowed=false and genus_rule_training_allowed=false. Observations are growth-chamber plant-material observations, not estimated natural-population selfing rates.

## Executable next wave

`scripts/acquire_wave57_reproductive.py` audits the previous artifact and performs 40 genus-centred reproductive literature searches across Europe PMC and Crossref, scanning all exact species names in returned titles/abstracts rather than just the first alphabetic target species. Provider failures, truncation and zero-yield queries are recorded separately. Three consecutive failures stop that provider. Abstract co-occurrence remains an unreviewed lead: no normalized value, High grade, genus-rule vote, or automatic promotion is generated. Search retrieval is bounded and is not exhaustive. Raw metadata cache is transient and excluded from uploaded artifacts.

The push-triggered workflow `.github/workflows/acquire-wave57-reproductive-literature.yml` uses immutable public Wave52/Wave56 run artifacts. Seven regression tests cover the genus-only false-positive, off-axis output, synonym mapping and fail-closed grading.

## Canonical boundary

Reproduction remains **48,527 / 106,295**; total Chapter 1 remains **222,759 / 318,885 (69.8556%)**. These are the documented private-plus-public checkpoints, not recalculated by Wave57. Public-Wave52 presence checks cannot establish net gain against private TRY plus Waves53-55. The remaining promotion gate is a synonym-aware, species-axis and trait-conflict collision audit against that current canonical ledger. No private TRY rows are stored here.
