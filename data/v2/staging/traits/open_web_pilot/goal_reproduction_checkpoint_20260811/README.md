# Goal reproduction checkpoint (2026-08-11)

This checkpoint extends the individually reviewed curated ledger used by
`open-web-review-promote.yml`. It contains 78 new species-direct records on
top of the 170-row NParks/high-leverage checkpoint:

- 49 additional NParks Flora & Fauna Web morphology records from complete
  near-rule genus inventories;
- two primary-source records for *Silene colorata*;
- 23 exact-name population mating-system records from Whitehead et al. (2018),
  DOI `10.3389/fevo.2018.00038`;
- one published quantitative outcrossing-rate report for *Magnolia
  salicifolia*;
- three primary-study mating-system statements for *Ficus arpazusa*, *Ficus
  petiolaris*, and *Dioscorea alata* (DOIs
  `10.1007/s10592-008-9776-x`, `10.3732/ajb.1100472`, and
  `10.3390/plants10071412`).

The Magnolia record is species-direct Medium evidence; the three added
primary-study records are species-direct High evidence. None is a genus
inference itself. The *Dioscorea alata* record uses an explicit statement that
cross-pollination is essential, without substituting its sex-system field for
mating system. A FloraWeb/BiolFlor pollen-vector record for *Hieracium
umbellatum* was reviewed but rejected from autonomous-selfing evidence because
the source explicitly treats pollination independently from fertilization; the
same cross-trait problem was also excluded for the two existing Hieracium rows.
The shared all-evidence rebuild may form Validated Low only from
`genus x axis x trait_name` rules with at least three direct species and the
configured dominance and masked-validation gates. Family inference, global
fallback, n=2 adoption, and substitution among reproductive traits are not
used.

`whitehead_2018_supplementary_data.zip` is the downloaded source supplement.
Every accepted record retains its exact source excerpt, stable URL, retrieval
hash, source lineage, identity decision, and line-by-line audit. File hashes,
row counts, mappings, and policy flags are pinned in
`goal_reproduction_checkpoint_manifest_20260811.json`.
