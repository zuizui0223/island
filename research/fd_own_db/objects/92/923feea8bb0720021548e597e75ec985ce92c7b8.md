# Completed GIFT control-trait acquisition — 25 November 2026

Protected [run 45530336990](https://github.com/zuizui0223/island/actions/runs/45530336990) SUCCESS for both jobs. Exact source logs reproduced in [control_trait_coverage.json](control_trait_coverage.json); artifact digests are GitHub metadata, not independently redownloaded archives. Full original 4,505 botanical island / 8,264 geography frame. 4,504 geography-linked botanical floras, 1,558 tropical. Public raw nonderived records with matched/resolved name, unrestricted cited source, no uncertain/infraspecific scope, single literal state; normalized exact accepted name match; contradictory species states removed.

| Trait | Public cited references read | Raw rows | Candidate rows | Unambiguous source names | Original flora matched |
|---|---:|---:|---:|---:|---:|
| Fruit dryness 3.18.1 | 3/3 | 2,751 | 2,689 | 2,530 | 2,183 |
| Growth form 1.2.1 | 167/167 | 240,714 | 216,251 | 123,825 | 48,686 |

Fruit 1190 fleshy, 993 dry names; growth 26,688 herb, 13,005 tree, 8,911 shrub, 82 other. Excluded contradictory states: 6 fruit, 5016 growth names. Exact name normalization is not an independent taxonomic review. GIFT metadata public count is different from fetched raw rows; metadata reference file hashes differed between the two jobs. Preserve both original hashes, do not claim identical reference snapshot across traits.

| Scope | Fruit ≥5 / ≥30 | Growth ≥5 / ≥30 | Fruit ≥5 and ≥50% roster coverage | Growth ≥5 and ≥50% roster coverage |
|---|---|---|---:|---:|
| Global all-recorded | 2582 /1533 | 3407 /2418 | 43 | 3097 |
| Tropical all-recorded | 710 /270 | 1035 /617 | 1 | 785 |
| Tropical recorded-native | 102 /43 | 132 /98 | 0 | 86 |

This is a **candidate-trait coverage audit**, not a finalized evidence-grade ledger or independently verified whole-flora representation. Tropical fruit coverage is especially sparse relative to the total flora; only 1 of the 710 min5 islands reaches 50%. The direct fruit endpoint may meet nominal sample size if adequate animal maps arrive, but not an unqualified whole-flora inference. No E1/E2, no association with isolation or animal diets, no change to frozen primary gate. Need protected repeat with source-reference identity, fixed grading rule, exact matched-species roster, genus-block variation and observation-mask diagnostics. Previous job released aggregate coverage only; no reusable protected species table was archived.

## Independent bat preparation

The official [Marsh et al. 2022 MDD v1.2 archive](https://zenodo.org/records/6644198) lists MDD_Chiroptera.zip (2.9 GB, MD5 3ecf49e018f4ab9654e22258776076fc) and mapped species metadata mdd_spList_wFamilieswOrders_mapped_6362species.csv (MD5 4a800c367f7f8d767779d2f68df40f93). Expert native range maps; not equivalent to IUCN code-filtered ranges. The published [mdd package](https://alrobles.github.io/mdd/) simplifies polygons at 1000 m tolerance; do not substitute simplified maps for 1 km² island primary comparisons. EltonTraits mammals [original Figshare file 5631084](https://figshare.com/ndownloader/files/5631084) uses older MSW3 taxonomy and contains taxonomically interpolated diet attributes. Source read and exact-name overlap census are separate from range-map analysis. Prior independent island mammal work may have a relevant archive, but no verified transferable bat incidence located; no held-out data decoded.

A new workflow checks the official mapped-species list and original mammal diet file before any range acquisition or trait effects. Report exact-name match fractions and candidate diet classes, missing taxonomy and diet provenance fields. Unmatched species never become diet zero. A global list match rate is not the actual incidence-conditioned 5% matching gate. No bat polygons or animal × plant association in this stage.
