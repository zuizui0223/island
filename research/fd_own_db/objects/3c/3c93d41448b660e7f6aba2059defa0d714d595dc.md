# Falsification contract — GIFT 4,505 island floral category prediction without focal genus oracle

**Frozen 19 November 2026 BEFORE original-data results from the new model.** This test probes how much of earlier successful *taxon-conditional* interregional prediction depends on knowing the target island's true botanical genus repertoire. It is designed as a negative-control extension to the 7/13 geographically robust, genus-conditioned botanical source result, not an adaptive pollination claim.

## New no-oracle predictor

Source botanical parent: every original **4,505** GIFT accepted island floras, each directly High/Medium species-documented floral endpoint, not historical triple-trait complete-case 928 islands. Original 10 single traits and 3 TRUE same-accepted-species categorical two-trait pairs; all-recorded and strict known-native original roster separately; all FOUR GSHHG macroregions and three training-deletion modes.

Predict a full categorical distribution of flower states for each target island WITHOUT looking up its real genus roster, genus frequencies, target flower category frequencies or species identities. The no-oracle predictor learns:

1. **GIFT target-region OTHER-ISLAND genus frequencies** from the number of region-exclusive accepted botanical species with direct state evidence per genus, excluding the focal entire accepted island flora names and all those present on ANY original island in the focal geographic 10-degree block, same first geo_k5 candidate mainland source entity or BOTH. This is an original directly documented **source species-pool genus weight**, not an estimated target island genus frequency or proof of colonization.
2. **GIFT OTHER-REGION species floral states** in the same genus, with outcome-blind hash-selected equal local/external donor counts up to 20, genus shrinkage alpha=5 onto the exact external regional categorical prior. The actual source names are disjoint by broad region and entire focal source blacklist; none is an independent ancestral mainland specimen.
3. Sum over weighted source-region genus probabilities to obtain ONE constant no-oracle floral category probability distribution per original target-region/geographic training group. A source genus without external trait-bearing congeners uses the outside-region genus-BLIND category prior; no invented category states or source taxa are used.

Source-training names in both regional genus weights and outside-region Direct category models obey the **same** strict geographic/source group blacklist. Original focal genus names are NEVER used to calculate any no-oracle source genus weight. The original region identity and spatial group needed for provenance holdout are known inputs, so this is NOT a completely geography-free prospective island SDM.

## Fair comparators: same original focal Direct plants, no posthoc cohort changing

Score **every** actual focal high/medium Direct accepted plant (minimum five per original target island), including genera without external donor data, by proper multiclass per-plant Brier/log loss. Before/after and reference models use EXACT same focal plant specimen roster.

Three source-trained predictions:
- **External GENUS-BLIND:** the pooled category prior of OTHER-REGION source plants after the same name/geography/source exclusions.
- **Regional local GENUS-BLIND:** the category prior of distant other GIFT island plants in the target macroregion with the same source blacklist. This is a strong empirical geography-only source-frequency benchmark, which may outperform the genus model.
- **No-oracle regional genus mix × external genus trait profiles** described above.

A FOURTH, explicitly labelled **oracle-known-focal-genus** comparator directly uses each focal plant's true genus identity and the same externally trained genus categorical model. It quantifies dependence on taxonomic foreknowledge; it is not silently included in no-oracle inference.

Freeze original source overlap gate before source analysis: at least 20 original Direct eligible distinct accepted species each in regional and other-region donor pools, two shared botanical genera and >=10% of original local source Direct species having an external-region trained congener. Original respondent sources failing this are **not identifiable** (do not infer absence of floristic or pollinator effects).

Require >=30 focal botanical island floras and >=10 original physical 10-degree blocks and >=10 distinct geo_k5 first candidate mainland-source groups for original source-block grouped 499-draw exploratory confidence intervals. A candidate source gain is robust only if Brier AND log-loss 95% intervals exclude zero in the same direction in BOTH groupings. Record **negative** and source-insufficient endpoints, not only positives. For a cross-region generality result demand the same positive sign in northern midlatitudes, tropics and southern extratropics under ALL THREE frozen training group exclusions. A separate, stronger claim must also beat the geographically regional blind prior. Original known-native subcohort cannot be treated as independent from all-recorded, and independent source regions may share lineage/trait coding.

## Scientific ceiling

If the no-oracle genus-mixture predictor fails against one or both source-blind baselines, this does NOT falsify the earlier confirmed **7/13 genus-conditioned** result. It demonstrates dependence on knowing the focal genus composition and/or regional phylogenetic assortment and limits claims of prospective **functional floral repertoire assembly** prediction.

If it succeeds, the result is STILL **source-pool-level genus-weight transfer based on different plant species recorded on other original GIFT ISLAND floras**. It does not identify ancestral mainland phenotypes (independent strict GIFT candidate mainland holdout found zero source novel Direct colour×form plant names outside all 4,505 island floras), actual within-lineage evolution, oceanic versus continental geological causal selection, animal guild disappearance, visitor preference, floral local phenotype shift or effective pollen/seed fitness.

**Pre-result original source protocol:** [oracle-free original model](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/fd_global_original_oracle_free_genus_mix.py), [synthetic negative/positive tests](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/tests/test_fd_global_original_oracle_free_genus_mix.py), [source-locked execution](https://github.com/zuizui0223/island/actions/runs/44474030421). Until original source CI success and signed anonymous result verification, no new numerical outcome belongs on main.
