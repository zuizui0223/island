# Rodger et al. 2021 pollinator-contribution snapshot

This directory records the exact structured source used for the
`autonomous_selfing_capacity` checkpoint.

- Article: Rodger et al. (2021), *Widespread vulnerability of flowering plant
  seed production to pollinator declines*, DOI `10.1126/sciadv.abd3524`.
- Dataset: Figshare `10.6084/m9.figshare.14607882.v1`, file ID `28041726`,
  `S3_pc_publish.csv`, CC BY 4.0.
- Retrieved: 2026-08-11.
- Original file SHA-256:
  `38ffb9b6438e8e0f683c0f1504ff500c59603b641dcf6936e4e29f232bc66e6a`.

`rodger_2021_pollinator_contribution_snapshot.csv.gz` is a deterministic
UTF-8, gzip-mtime-zero projection of the 1,528 original rows and the columns
needed to reproduce the checkpoint.
`rodger_2021_two_backbone_synonyms_20260811.csv` contains only 48 submitted
names for which exact GBIF API and WFO Plant List 2026-06 resolutions agreed
on the same accepted species and family. The checkpoint does not accept
one-backbone, fuzzy, or family-conflicting names.

The pollinator-exclusion fields (`auto.*`) measure fruit or seed production
without pollinator access. Positive measured output is evidence of autonomous
capacity; measured zero is evidence of absence. GloPL hand-outcrossed bagged
flowers, self-compatibility, and mating-system labels are never substituted.
