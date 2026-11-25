# Functional floral form and GloBI visitor channels: finer ecological qualification

**11 November 2026 · fd Draft PR #14 · original-source-pinned posthoc analyses.** No original botanical genus, species, island, pollinator taxon or GloBI provider/reference identity is released. The endpoint is **global GloBI documented flower-interaction reference shares**, NOT true local pollinator abundance, effective pollination, or a named pollination syndrome.

## Provenance and cohort

We reused the frozen source-scored Direct High/Medium, single-state flower colour, form and symmetry evidence and GloBI v0.9 source-genus breadth. The original private progressive `island` ZIP SHA256 is `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`; the original 2026-06-04 GloBI source-genus breadth ZIP is `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`. An unambiguous frozen GBIF family crosswalk provides **484 matched genera in 88 families**, representing **1,590 direct floral morphology species**. Only the percent open-radial floral form per genus is added to the three original GloBI source-documentation controls. Each of the 88 families is omitted in full from the ridge training data for its test genera (alpha=100, train-only normalization, five-share simplex projection).

## 1. Bird documentation is the clearest family-equal individual response

MSE gain = source-reference/plant-coverage-only prediction error minus prediction error after adding the proportion of open-radial floral forms. Positive values are better prediction. Conditional 4,999 family-equal bootstrap intervals **give each family one vote** rather than overweighting large plant families.

| GloBI documented visitor category | Genus-equal MSE gain | Family-equal MSE gain | Family-equal 95% interval | Positive families |
| --- | ---: | ---: | ---: | ---: |
| Bombus | +0.000272 | −0.000265 | [−0.001749,+0.001172] | 42/88 |
| Other bees | +0.001348 | +0.001717 | [−0.000947,+0.004345] | 52/88 |
| Lepidoptera (moths AND butterflies, unseparated) | +0.001054 | +0.000164 | [−0.002718,+0.002920] | 47/88 |
| **Flower-visiting birds** | **+0.000726** | **+0.002091** | **[+0.000927,+0.003533]** | 53/88 |
| Diptera | +0.002334 | +0.001741 | [−0.000362,+0.004009] | 55/88 |
| **Five-channel equal-weight average** | **+0.001147** | **+0.001090** | **[+0.000122,+0.002077]** | 54/88 |

Even if the five compositional fractions are analyzed individually, they are not independent outcomes. The family-equal bird share is the only individually positive interval in this exploratory study. The 5-channel average avoids cherry-picking a favorable guild. Excluding the two largest plant families from *test scoring only* (377 remaining genera) still gives a positive pooled MSE gain **+0.001571**; this is not an independent study.

## 2. Intra-family morphological negative control

Shuffle the **open-form feature only within each family** while retaining its genus membership, 3 original source-effort predictors, all five observed GloBI channel proportions and the same whole-family leaveout folds. Out of 88 families, **32 / 308 genera** actually have two or more different open-form values permitting a nontrivial permutation. There are 999 fixed-seed shuffles.

For the five-channel average, the observed genus-equal gain **+0.001147** exceeds the shuffled 95% range **[+0.000201,+0.000751]**, posthoc upper tail **0.001**. The documented *bird* proportion: observed +0.000726 versus null 95% [+0.000057,+0.000546], upper tail 0.001. Diptera observed +0.002334 versus null [+0.000395,+0.001561], upper tail 0.001. Other bees have within-family null tail 0.075 and a family-equal interval crossing zero.

This is NOT a phylogenetically disjoint, subfamily-disjoint, GloBI reference-provider-independent, preregistered or effective-pollination test. Limited within-family trait variation strongly limits what the permutation can identify. The six guild/pooled aggregates and all 999 null values were recomputed using a separately written function and agreed to machine precision for observed genus/family gain and permutation summaries. Public tests use synthetic plant profiles only.

## 3. Critical sensitivity: does explicit pollination evidence qualify the original global GloBI genus?

The original GloBI `n_pollination_rows` and `n_flower_visit_rows` are **per genus**, not per individual guild. Of 484 selected genera **151 have zero explicitly classified pollination rows** despite being included with >=3 total independent flower-interaction references. The dataset has **121,397** broadly classified pollination rows versus **170,565** generic flower-visit rows in these genera. These counts are source-documentation types, not validated pollen transport.

To avoid the shortcut of treating visits as proven pollination, keep only genera having a minimum number of separately classified `pollination` interaction rows (the ORIGINAL five-share outcome still mixes `pollination` and `flower_visit` across each channel). Rerun the same leave-one-family-out model and a conditional family-equal bootstrap:

| Explicit `pollination` rows per genus required | Tested genera / families | Family-equal bird-share MSE gain | Conditional family 95% interval | Five-channel average gain |
| --- | --- | ---: | --- | ---: |
| >=1 | 333 / 72 | +0.002095 | [+0.000689,+0.003721] | +0.001343 |
| >=3 | **294 / 68** | **+0.002236** | **[+0.000825,+0.004049]** | **+0.001326** |
| >=10 | 248 / 62 | +0.001077 | [+0.000304,+0.002035] | +0.001049 |
| >=20 | 210 / 57 | +0.000909 | [+0.000164,+0.001836] | +0.000689 |
| >=50 | 173 / 50 | +0.000266 | [−0.000020,+0.000640] | +0.000676 |

The bird documented-share directional improvement remains positive through >=20, but the >=50 subgroup no longer has a bootstrap interval excluding zero. Across thresholds, the **sampling set and family identity change**; they are nested, post-hoc sensitivity checks, not five independent effects. Some recorded pollination rows may not be from the target bird channel at all; this filtering does **not** change the reference shares into validated guild-specific effective pollination. The Diptera family-equal intervals still cross zero at all shown thresholds. GloBI's original data integration documentation emphasizes geographical/taxonomic source bias and independent attribution to original datasets.

## 4. Source-provider blocking persists

The frozen 2,228-genus GloBI aggregate does **not** contain per-plant-genus × channel × reference original study/provider IDs; the separate pollinator-species `reference_examples` catalog is not joinable back to the scored plant-genus outcomes. The original GloBI v0.9 `interactions.tsv.gz` (Zenodo `10.5281/zenodo.20546682`, SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`) is necessary for a truly provider-heldout estimate. **No provider-heldout positive empirical result is claimed.**

## Decision

**The noncausal finding currently worth retaining:** in globally documented GloBI five-channel flower-interaction reference shares, simple directly recorded floral openness exhibits a small, heterogeneous predictive signal across unseen botanical families, with the clearest family-equal evidence for the documented bird channel. Its magnitude and interpretation depend on a selected subset of taxa and reference classes.

**Not warranted:** that islands with open flowers have more effective bird pollination, that birds or Diptera shaped flower evolution, that moth/butterfly visitation has shifted, or that a universal pollination syndrome exists. Confirmatory next tests need per-guild, per-provider source-original interactions and independent local flower-visitor / pollen-transfer observations.

## Reproducibility

- `src/fd/globi_family_guild_audit.py` / `tests/test_globi_family_guild_audit.py`: fail-closed genus-family IDs, one-family-one-vote estimates, whole-family train exclusion, matched 999-draw within-family morphological label null, aggregate-only outputs.
- `per_guild_88_family_within_family_morphology_null.csv`: six anonymous original-data summary rows.
- `pollination_evidence_family_coverage_sensitivity.csv`: eighteen anonymous rows summarizing bird, Diptera and five-channel pooled thresholds.
- GloBI and `island` original taxon/provider rows were never committed to the public fd branch. All empirical results are exploratory; public CI validates synthetic rules and no raw genus names.


## 5. Independent pollinator-taxon catalog evidence-class audit

The independently frozen original `channel_taxon_catalog_confirmatory.csv` in the GloBI qualified-interactions artifact **10291643834** (separate from genus→guild modeled input) has one row per documented visitor taxon per channel, with raw `n_pollination_rows` and `n_flower_visit_rows` counts. Summing those source-v0.9 catalog rows yields:

| Documented visitor channel | Confirmatory visitor taxa | `pollination` labeled interaction rows | `flower_visit` labeled rows | Fraction explicitly labeled `pollination` |
| --- | ---: | ---: | ---: | ---: |
| Bombus | 146 | 41,495 | 135,466 | **23.45%** |
| Other bees | 2,206 | 143,049 | 115,980 | **55.23%** |
| Lepidoptera | 2,360 | 8,515 | 167,712 | **4.83%** |
| Flower-visiting birds | 370 | 2,670 | 8,651 | **23.58%** |
| Diptera | 1,846 | 34,753 | 46,800 | **42.61%** |

The source contains **visitor-taxon-level counts**, not the exact same 484 botanical genus × GloBI channel training records. Therefore these global fractions **cannot** be used to back-correct the five-channel genus MSE or infer local island pollination effectiveness. Even a record labeled `pollination` is not equivalent to independently documented successful pollen transfer/seed set. For Lepidoptera in particular, the literature reference composition heavily reflects plain flower-visit records, and the original GloBI group does not distinguish daytime butterflies from nocturnal moths.

Released grouped original-data source: [`globi_qualified_pollinator_catalog_evidence_classes.csv`](globi_qualified_pollinator_catalog_evidence_classes.csv). Original pollinator species names, precise reference citations and provider IDs were NOT exported to fd.
