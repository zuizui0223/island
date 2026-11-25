# 2026-11-16 — Independent NHM numeric corolla-tube length: original GIFT source checks and genus-heldout deep-tube prediction

**Status: completed original-data analyses, exploratory / SOURCE-PROVENANCE-limited, NOT a pollination-shift or adaptive field-phenotype result.** Original source-hash verified [NHM numeric-length opportunity run #43874963541](https://github.com/zuizui0223/island/actions/runs/43874963541) SUCCESS, [anonymous artifact #14228584961](https://github.com/zuizui0223/island/actions/runs/43874963541/artifacts/14228584961). Follow-up [source-only genus-heldout numeric-to-Direct-depth predictor run #43883317272](https://github.com/zuizui0223/island/actions/runs/43883317272) SUCCESS, 14 combined synthetic tests and source-original GIFT/NHM/GSHHG checks, [anonymous artifact #14232618500](https://github.com/zuizui0223/island/actions/runs/43883317272/artifacts/14232618500). The published source and protected island botanical names, coordinates and exact GIFT candidate region IDs are NOT released.

## New source evidence: actual centimetre-length data, not another imputed shape category

The immutable [Natural History Museum plant-traits repository at commit 4843c22b](https://github.com/NaturalHistoryMuseum/plant-traits/tree/4843c22bb2fc91018d9828dcddebf067c83b1097) contains OCR-extracted `corolla tube min. length [cm]` / `corolla tube max. length [cm]` in its original five constituent descriptive sources. Earlier fd NHM scripts used flower colour, form and symmetry only and overlooked this numerical axis. We reconstructed all FIVE original provider CSVs independently, each checked against its archived Git blob. The NHM `combined.csv` merely equals the union of original exact species (11,262) and was NOT counted as a sixth independent corpus. Its row-level botanical study provenance is missing. We parse only explicit positive finite centimetre measurements, accept consistent min≤max, collapse multiple rows/corpora per species by median, and separately report source disagreement (a 1.5-fold between-corpus threshold is only a sensitivity, not study-level independence).

**Actual independent-release NHM numeric tube measurements**: 11,262 accepted exact species; 469 have >=2 different originally named corpora with consistent lengths (197 multi-source disagreements >1.5x). GIFT candidate native mainland original accepted plant rosters include 9,492 of the source-only NHM numeric species; 3,410 are also directly described in the original `fd` ordinal `tube_depth_class` evidence. Original 4,505 GIFT island floras remain the study population.

| Original categorical Direct depth | EXACT plant names with external tube cm | NHM per-species length median | NHM IQR cm |
|---|---:|---:|---|
| shallow | **1,185** | 0.35 cm | 0.25–0.445 |
| intermediate | **2,052** | 1.00 cm | 0.70–1.35 |
| deep | **807** | 3.50 cm | 2.54–5.00 |

These are positive **cross-source agreement in rank of functional corolla-tube length** among species, but the ordinal labels may derive from overlapping literature, NHM numerical values are OCR-extracted historical descriptions rather than matched island/mainland field populations, and no within-species floral evolution is measured. An ordinal `deep` is not universally one fixed centimeter cutoff; we never threshold NHM lengths into fictive independent `deep` labels.

**Original source-region congener opportunity:** When the focal island's species names (but NOT every species found on ANY island globally) are excluded from geographically closest five and climate-nearest ten native GIFT source rosters, original Direct ordinal depth with >=5 target carriers and >=2 matching genuine botanical genera supported by >=10 source numeric species in BOTH modes is available on:
- tropical: all-recorded 436 islands (78 geographical blocks), known native 75 (30 blocks), known native nonendemic 74 (29 blocks);
- northern midlatitude: all-recorded 1,428 (79 blocks), native 189 (33 blocks);
- southern extratropical: all-recorded 34 (17 blocks), native 3 (insufficient).
Under >=2 original independent corpus agreement, tropical all-recorded 211, known native 37. This is an OPPORTUNITY result, not prediction or causal evidence.

## Critical failed strict calibration and justified genus-heldout substitute

A more stringent calibration from **mainland GIFT species recorded on NONE of the 4,505 islands** and having BOTH NHM numerical tube length and original Direct ordinal depth gave **0 species**. The initial [Actions #43882350699](https://github.com/zuizui0223/island/actions/runs/43882350699) correctly FAILED CLOSED at `independent external source calibration missing`. A separate intermediate run calculated results but failed the stale CI assertion, not the data calculation. We do NOT call the substitute an all-island-name-independent study.

Instead, for EACH original focal botanical genus `g`, train a logistic calibration `logit Pr(Direct deep) = alpha_{−g} + beta_{−g} log(1+NHM corolla tube length in cm)` from the **3,410** mainland-candidate accepted species where both independent-release NHM length and original Direct depth are recorded, excluding ALL Direct ordinal records of the WHOLE focal genus from that calibration. There are 689 deep and 2,721 not-deep in this original candidate calibration universe. This avoids using the predicted genus's own ordinal target phenotype labels but includes other genera that can occur on other island floras. Original source NHM and Direct labels might share botanical literature; neither the source-specific calibration nor NHM constituent corpora are independent experimental population samples.

Each focal island then predicts the original directly documented deep or not-deep status of the **SAME target accepted species** under two models:
1. **GIFT-source congener length**: apply genus-heldout calibration to different accepted source congener species' NHM cm, and use genus-specific mean predicted deep probability for every focal island plant in that genus;
2. **SAME source donors pooled**: use exactly those different accepted native mainland numeric species, pooled without the focal island's own genus composition.
Positive score = pooled Bernoulli Brier MSE minus congener-genus Brier MSE, both on identical targeted original island species. 999 original 10-degree geographic-block resamples and first-source-group resamples; the original mainland geo5/climate10 modes are NOT independent biological replications. A second support gate requires the documented proportion of original island Direct depth and original source NHM numeric length in the SAME matched genera to differ by <=10 percentage points in BOTH source modes. Status-known native and two-corpus agreement form selected nested subsets.

## New numerical source-only heldout-genus prediction

| Region and original botanical status | NHM numeric donor support | Original qualified island floras | Geo5 Brier gain, 10-degree-block equal [95%] | Climate10 Brier gain [95%] | Same-genus Direct documentation balanced in both modes |
|---|---|---:|---|---|---:|
| **Tropical all-recorded** | >=1 published numeric source | **425** | **+0.10803 [+0.09526,+0.12027]** | **+0.12198 [+0.11130,+0.13242]** | **8** |
| Tropical all-recorded | >=2 consistent named NHM corpora | **116** | **+0.09407 [+0.06271,+0.12218]** | **+0.09145 [+0.05381,+0.12341]** | **0** |
| Tropical known-native | any numeric source | **68** | +0.04797 [−0.00578,+0.09544] | +0.08722 (conditional) | **3** |
| Northern midlatitude all-recorded | any numeric source | **1,422** | +0.01305 [−0.00292,+0.02869] | +0.01823 [+0.00466,+0.03236] | **134** |
| Northern high latitude all-recorded | any numeric source | **136** | **−0.01517 [−0.02447,−0.00603]** | −0.02241 [−0.03200,−0.01320] | **8** |

Southern extratropical all-recorded has only 22 islands with source-only congener prediction in BOTH modes, 8 geographic blocks, below the fixed eligibility gate; do not announce any negative or positive effect there. The 15/30 source-opportunity panels from the first audit and **11/30** qualified 30-island/10-block actual prediction panels are **not the same denominator**: latter restricts to the identical paired original island ordinal-known plants and source genus evidence in both modes.

**Defensible positive new result:** In tropical sampled floras, source-only centimeter measurements from different accepted mainland congener plants, mapped to depth categories with the focal botanical genus completely withheld from calibration, are considerably more predictive of original island ordinal deep-tube membership than pooling those same source numeric specimens without genus identity. This independently supports a functional-length axis and genus-associated botanical assemblage composition, not a unique bird/butterfly/moth floral syndrome. The score improvement is a selected documented cohort of 425 original tropical islands with strong source and target documentation asymmetry: **ONLY 8/425** meet original <=10pp annotation balance, and known-native has only **3/68**. Thus it is **NOT robustly identified** as a general deep-tube ecological effect among native tropical island plants. Never use a positive cluster CI as a causal guild-effect confidence interval.

## What to do next

Priority 1: quantitative **original source numeric length by botanical FAMILY and plant-provider heldout** on the exact original paired cohorts to assess study/lineage leakage; county-level/source-list potential opportunity not equal actual seed arrival. Priority 2: seek truly comparable independent island/mainland source population COROLLA LENGTH and actual visitor measurements (birds, butterflies separately from hawkmoths) linked to accepted species; NHM is not local phenotype. Priority 3: if useful, compare actual real symmetry+length within a FULLY independent trait provider, but earlier NHM symmetry-only original cross-corpus support failed; do not fill it from original Direct species ordinal/depth to claim independent two-trait validation.

Neither this output nor current source-only congener flower-pair forecasts identify historical pollinator turnover or island floral adaptation. Maintain [fd PR #16](https://github.com/zuizui0223/fd/pull/16) DRAFT.
