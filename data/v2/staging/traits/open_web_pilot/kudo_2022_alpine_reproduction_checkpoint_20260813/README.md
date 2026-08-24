# Kudo 2022 alpine reproductive evidence checkpoint

This checkpoint audits all 46 rows in Table 1 of Kudo (2022), *Outcrossing
syndrome in alpine plants: Implications for flowering phenology and pollination
success*, Ecological Research 37:288–300,
[doi:10.1111/1440-1703.12314](https://doi.org/10.1111/1440-1703.12314).

The public article reports controlled self-pollination and/or pollinator-
exclusion experiments for 41 taxa. The rendered Table 1 was retrieved on
2026-08-13 JST and its normalized inner text has SHA-256
`24737a195d4aca322275119f53b949a9adef070edef7c64a176afffe0e0040fe`.

Of 46 table rows, 23 are accepted as High species-direct evidence. Five rows
contain no trait statement, eight are below species rank, seven lack a strict
target-master match, two have a source/master family conflict, and one contains
an unresolved source-name typo. No fuzzy match or rank roll-up is used.
`Loiseleuria procumbens` is the sole synonym promotion, because the committed
GBIF/WFO snapshot independently agrees on `Kalmia procumbens`.

Mappings remain trait-specific:

- `Self-incompatibility` -> `self_incompatibility=SI`
- `Autonomous selfer` -> `autonomous_selfing_capacity=autonomous`
- `Mixed mating` -> `mating_system=mixed_mating`
- `Obligate outcrosser` -> `mating_system=predominantly_outcrossing`

These records do not convert compatibility into autonomy or mating system.
They also do not emit genus inference. The shared all-evidence rebuild performs
lineage deduplication, direct-conflict handling, genus x trait validation, and
strict species x axis coverage.
