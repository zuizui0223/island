# Outcome-blind full-flora genus-disjoint matched floral functional beta — original source run 18 Nov 2026

**SOURCE-VERIFIED ORIGINAL RESULTS; observational and exploratory, not independent evolutionary replication.** Primary strict full-roster genus-disjoint category overlap is in [first-pass report](README.md). This follow-up **does not** decrease the global botanical source population from 4,505 GIFT floras or invent Direct phenotypes. Code [genus_disjoint_matched_beta.py](../../src/fd/genus_disjoint_matched_beta.py) and [synthetic tests](../../tests/test_genus_disjoint_matched_beta.py), protected original [island runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_genus_disjoint_matched.py).

## Original-source execution and provenance

- SHA256-locked original high/medium Direct accepted-species ledger and 4,505 GIFT island×plant accepted floras; immutable original published GSHHG and Weigelt geography crosswalk. **All-recorded 2,324** genuine colour×floral-form Direct-n≥5 floras; separately native-known **327**, both exact originals.
- Initial matched source [island Actions #44293839934](https://github.com/zuizui0223/island/actions/runs/44293839934) **SUCCESS**, synthetic integrity regressions **16 passed**, anonymous artifact **14421444251**, signed archive SHA256 \`9d6f79013f455c105e580a4945edc9c4957f068c3ba2c21af4292d28c60e81fa\`.
- Stronger 3-block-disjoint source run [#44295181899](https://github.com/zuizui0223/island/actions/runs/44295181899) **SUCCESS**, expanded original synthetic integrity **17 passed**, anonymous artifact **14422213576**, signed archive SHA256 \`77af3dfa62fc617c0161fbf50f6e35b1c33ad5124510f586abd35b1afc77583e\`.
- Anonymous result JSONs contain original regional counts and numerical summaries, **not** island, accepted plant, botanical genus, exact locality or mainland source identifiers. Both runs reproduced 4,505/2,324/327 denominator gates. The duplicate jobs use the SAME source and are **not** independent empirical samples.

## Frozen estimand and support design

For each *same focal botanical island*, draw one treatment candidate **G0: no accepted botanical genera shared across either entire selected floristic roster** versus one control **G1: zero shared full-roster accepted species but >=1 shared accepted genus**. All three island rosters have certifiable accepted-genus names; same original broad region; focal–donor pair must lie in distinct physical 10-degree blocks; donorship and covariates selected **without inspecting flower category outcomes**.

The outcome is exact species-frequency total-variation floral category overlap \`O = Σ min(p_i, p_j)\`. Within-focal contrast is \`O(focal, G0) − O(focal, G1)\`. Donor-to-donor calipers are focal–donor great-circle distance (within 2× primary), Direct count (2×), full-roster species count (2×), genus richness (2×), Direct coverage fraction (2×), area (within 4×); 4× broader caliper source mode is a sensitivity, not a trial selected by outcomes. Additional modes demand treatment versus control donor lie in distinct blocks or climate PC1 difference <=1.5 (when observed). One matched triplet per focal island, with unmatched focals excluded *as missing support*, never assigned zero effect. This is **not** a randomized experiment or direct observed GIFT source-regional phylogenetic comparison.

The source-independent first GIFT mainland donor-region ID is unavailable in the independently published geographic crosswalk used here, so **no independent first-GIFT-source reweighting/bootstrap** is claimed.

## Verified primary 2× original source results (all-recorded, different focal–donor 10-degree blocks)

| Geography | Matched focal islands | Independent focal 10° blocks | Observed G0 overlap | Matched G1 overlap | Focal-equal G0−G1 | 10° focal-block-equal contrast [exploratory 95% focal-block bootstrap] |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| All eligible within own region | **1,112** | **139** | **0.215415** | **0.233675** | **−0.018260** | **−0.006463 [−0.021426,+0.008933]** |
| Northern midlatitude | **942** | **75** | 0.208659 | 0.228364 | −0.019705 | −0.003137 [−0.020881,+0.015756] |
| Tropical | **128** | **52** | 0.273872 | 0.279019 | −0.005147 | −0.007543 [−0.042162,+0.027060] |
| Southern extratropical | **42** | **19** | 0.188781 | 0.214594 | −0.025813 | −0.015544 [−0.065861,+0.039978] |

*The two donor groups are observationally taxonomically different. Even after matching, their source pool, ancestry, archipelago, third-party observer and effective guild histories are NOT demonstrated exchangeable.* All matched focal–donor block networks per geography contain **ONE connected component**. Therefore this 10° focal-block bootstrap ignores extensive donor re-use and is **not a conservative independent-replication confidence interval**.

A different 2× test additionally requires G0 and G1 DONORS to be in distinct original geographical blocks: **global 1,109 focals, 138 focal blocks**, G0−G1 focal-equal **−0.024525**, focal-block-equal **−0.018519 [−0.037257,−0.002270]**. That apparent negative interval **must not be elevated to ecological selection**: the donor network remains ONE connected block component, matched structure changes and the stricter disjoint three-block reselection below is not consistently negative.

## New outcome-blind disjoint three-block hypergraph stress

To eliminate reused physical blocks **within each reselected subset**, each selected source triple must use THREE distinct block IDs and no original 10° block may appear in ANY other selected triple. Random selection priorities are set *before floral responses*; 199 seeded **numerical** greedy re-selections, each choosing one block-disjoint subset, are **not 199 biological experiments**. The 5%–95% spread is a **reselection sensitivity range, not a confidence interval or permutation p-value**.

| Source/matched geography and primary 2× match mode | Three-block-eligible candidate triplets | Physically disjoint triplets per reselection (min / median / max) | Median observed G0−G1 overlap | 5%–95% outcome-blind reselection spread |
| --- | ---: | --- | ---: | --- |
| Global same-region matches | 810 | **19 / 22 / 26** | **−0.01283** | **[−0.06299,+0.04528]** |
| Northern midlatitude | 665 | **6 / 9 / 11** | **−0.02995** | [−0.10078,+0.03851] |
| Tropical | 109 | **8 / 11 / 13** | **−0.00076** | [−0.05844,+0.07493] |
| Southern extratropical | 36 | **3 / 3 / 5** | **−0.03889** | [−0.19524,+0.08993] |
| Global, both donors in independent blocks by design | 1,109 | **19 / 24 / 27** | −0.02174 | [−0.06581,+0.03025] |
| Tropical, both donors in independent blocks | 128 | **8 / 11 / 13** | −0.01123 | [−0.08356,+0.06405] |
| Global, climate PC1 donor caliper | 785 | **19 / 23 / 26** | −0.00223 | [−0.05236,+0.04347] |
| Tropical, climate PC1 donor caliper | 105 | **8 / 11 / 13** | +0.00256 | [−0.07448,+0.06682] |

The truly block-disjoint matching opportunities drop from 1,112 focal comparisons to only 19–26 retained triplets per global draw, 8–13 tropical: the region-level **geographically independent support is thin**, even if original pair counts seem massive. Full selection requires matching true source-region/archipelago histories before calling them independent evolutionary realizations.

## Complete source-mask/caliper/native sensitivity, verified follow-up

A new source-protected, same-code re-execution [island #44296122505](https://github.com/zuizui0223/island/actions/runs/44296122505) **SUCCESS** (17 synthetic tests; independent CI artifact **14423576816**, SHA256 \`7ab06ddfee363ed754ae1f286b1950dbd6168b604ec29edf90929f8e44227951\`), printed **all 36** original matched numerical panels meeting >=10 focal islands and >=5 blocks, including frozen 4× broad calipers and independently status-known native cohorts. Original 4,505/2,324/327 botanical population and trait-mask checks remain unchanged. The 3 successful source re-runs are numerically identical when the selection mode is identical, not separate biological replications.

**Broader 4× primary donor caliper** (same within-focal outcome-blind design) changes the all-recorded global matched focal count from 1,112 to **1,123**, primary focal-equal genus-disjoint minus genus-sharing overlap from **−0.018260** to **−0.015773**, and block-equal value from **−0.006463** to **−0.000246** [exploratory focal-block 95% **−0.017542,+0.016936**]. The 3-block-disjoint global median remains **21 triplets**, effect median **−0.012011**, with 5%–95% outcome-blind reselection spread **[−0.060282,+0.038178]**. These are **different donor matches of the SAME data**, not replication.

**Tropical caliper sensitivity:** primary 2× **128** focal islands, focal-equal difference **−0.005147**; primary 4× **139** focal islands, focal-equal difference **+0.007939**. Yet 4× geographic focal-block-equal difference is still near zero **−0.002098** with exploratory interval **[−0.037998,+0.031983]**; 4× three-block-disjoint median 10 independent triples, median effect **+0.001620**, outcome-blind 5%–95% **[−0.067005,+0.069752]**. **The raw focal-equal tropical sign reverses when loosening the prespecified covariate caliper**, contradicting any stable directional tropical G0–G1 ecological effect.

**Strict original known-native** is much narrower even though it uses the independent native-reported GIFT status rather than declaring unresolved plants native. Globally all eligible native-known Direct **327** islands remain the source parent, but only **71 focal islands / 19 focal blocks** support a 2× G0 vs G1 matched comparison; all source-qualified focal matches in this output come from **northern midlatitude**. Native 2× primary focal-equal contrast **−0.024396**, focal-block-equal **−0.055031** [exploratory **−0.109765,+0.001263**]. Physically nonoverlapping three-block matchings retain a median only **3 independent triplets**, with outcome-blind selection 5%–95% effect span **[−0.25,+0.1465]**. There is **no tropical strict-native matched effect to report**. Thus first-pass tropical all-recorded category overlaps must not be transformed into native flora adaptive-guild conclusions.

The strict donor-block-only 2× focal-block bootstrap occasionally excludes zero, but its connected-source network and the very small geographically disjoint subsets fail to establish independent causal geology or pollinator effects. Do not select that one conditional interval over the full preregistered caliper/independence sensitivity family.

## Scientific interpretation and main manuscript treatment

1. **Positive direct recorded botanical phenomenon retained:** original genuinely different floristic genera can nevertheless occupy some of the same recorded whole-flower colour×form categories; first-pass unadjusted strict genus-disjoint global overlap **0.214879** across 155,569 original source island pairs.
2. **Falsification of “strong above-genus-null cross-lineage convergence” remains:** earlier global fixed-coverage genus-preserving source-shuffle distance expectation was 0.256538 vs observed 0.280641 among those 155,569 pairs (greater observed divergence, NOT extra similarity). After within-region km/effort/genus-richness/area matching, there is no robust large additional positive G0-over-G1 functional similarity. The tropical response is small and variable, especially on independent block subsets.
3. **Partial identification, not mechanistic rejection:** the 2× focal-block intervals can cross zero even with thousands of correlated focal matches; geographical three-block reselection gives extremely few independent units and broad sensitivity spreads crossing zero. These cannot reject subtle cross-lineage adaptive convergence with known power, nor demonstrate mainland source colonization or pollinator guild loss.
4. No original tropical island-local butterfly/Bombus/bird/hawkmoth effective pollen deposition, true plant population-level island vs mainland phenotype, source dispersal trial, phylogenetic independent sister-lineage comparator or GIFT-mainland-source-group bootstrap is included.
5. **Main story remains global floral functional repertoire recurrence + regionally heterogeneous source-relative paired morphology + strong genus-composition contribution.** This matched result is a robustness/identification limit, not a reason to discard the large global island survey or to return to the invalid 17D KDE volume.

**Paper status remains DRAFT** [PR #18](https://github.com/zuizui0223/fd/pull/18), along with broader [PR #16](https://github.com/zuizui0223/fd/pull/16). Next decisive target is **actual independent botanical-source and regional guild observations**, not unadjusted extra exploratory comparisons.
