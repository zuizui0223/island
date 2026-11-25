# Original GIFT 4,505 island genus composition forecast: geography, climate, area, isolation and candidate source opportunity

**Frozen pre-result protocol, 19 November 2026.** The global 4,505 GIFT island floral analysis showed 7/13 conditional traits transferable between different species of the same genus **when the actual target island genus roster is supplied**. A simple no-target-genus model has only 1/13 versus external blind, zero versus the in-region floral trait frequency baseline. The next source-original independent ecological check explicitly forecasts island floristic GENUS composition from geographic features, before reprojecting floral traits.

The outcome is **whole-flora accepted botanical species frequency by genus** on each GIFT island, all 4,505 floras retained as parent. That is a checklisted genus composition, not an observed colonization probability, and checklist omission is not proof of genus absence. Strict native accepted GIFT floras are a nested later sensitivity.

## Data and no-genus-oracle prediction protocol

Original GSHHG geography: 10-degree blocks, lat/long on the sphere, four original climate PCs, log(1+area km²), log(1+distance to continent km), source first geo_k5 GIFT mainland candidate entity and full top-five ranked candidate entities. Normalization of geographic features uses their original geography-only median/IQR, not plant genus or morphology labels.

For each original target island, exclude ALL other GIFT island donor communities in the **SAME 10-degree physical block** and ALL those with the **SAME first GIFT candidate mainland source entity** from TRAINING. The target's true genus or plant identities are used ONLY for evaluation. Require at least 30 donor floras from at least 10 blocks and 10 source groups. Focal must have at least 10 accepted botanical species, 5 recognized genera, and complete original geography; missing/uncertain source availability is accounted rather than interpreted as biological nonarrival.

The five frozen genus-frequency predictors all use the SAME admitted source islands and identical plant/genus category outcome:
- **Region-equal baseline:** each admitted GIFT donor island contributes its normalized whole-flora genus species composition equally.
- **Geographic nearest 40:** spatial WGS84 spherical nearest neighboring source floras, rank-weighted and shrunk 20% toward the common region mean.
- **Climate nearest 40:** original four climate PC Euclidean standardized similarity, same neighbor count, weighting and shrinkage.
- **Environment nearest 40:** four climate PCs, log(area), log(distance to continent) with the same frozen regularization.
- **Environment + candidate-source similarity:** environment metric plus 0.75 × (1−Jaccard similarity between original five GIFT candidate mainland source entities). Donors sharing the FIRST-ranked source are STILL forbidden; this is opportunity similarity, not observed dispersal.

Unknown target genera absent from all original eligible source floras are grouped into one COMMON unknown-genus category with the same fixed 0.02 probability across models; no focal genus gets an individualized source probability or trained model coefficient. All scored target source accepted plant taxa share identical original test specimens and model donor sets.

## Strict numerical confidence and scientific ceiling

Evaluate proper multiclass categorical Brier and log loss of actual whole-flora accepted species genus distribution, versus the region-equal baseline. Retain negatives and source-ineligible original islands. Calculate separately original 10-degree physical-block-equal and first GIFT source-group-equal 399 bootstrap 95% exploratory intervals; require both proper metrics and both group signs to agree for a supported region. Report all four original GSHHG regimes and original all-recorded plus nested literal-known-native scope. Region/source pooled groups are not independent evolutionary events, intervals not multiplicity-corrected.

A reproducible environmental advantage over region-equal genus frequencies would support SOURCE/ENVIRONMENT-related floristic genus sorting, not independent mainland source ancestry. The geospatial source-neighbor benchmark must also be shown; any claim of additional source Jaccard importance requires direct PAIRED source-model differences and group intervals, not just two source models independently beating the baseline.

**Not measured:** actual genus arrival, local genus dispersal, source-congeneric phylogenetic ancestry, actual mainland flower traits, native-specific evolutionary adaptation, insect guild replacement, pollination effectiveness, flower-specific seed/fruit fitness. Published original top-five candidate mainland regions are source POTENTIAL, not actual observed dispersing seeds. The first-stage taxonomic genus prediction must be separately connected to the original species-Direct floral traits before claiming improved island functional diversity predictions.

**Status:** original SHA-locked GIFT source Actions [#44480529633](https://github.com/zuizui0223/island/actions/runs/44480529633) running; synthetic source coverage and genus loss tests enforce no-target-genus selection, source geographic holdout and negative score retention. Never merge an unverified numeric biological conclusion to main.
