# 13 November 2026: pollinator opportunity versus **island floral functional diversity**

**Original-data macro-ecological analysis completed; not an estimate of the causal impact of pollinator decline.** Work branches: [fd PR #15](https://github.com/zuizui0223/fd/pull/15), [island original-source PR #325](https://github.com/zuizui0223/island/pull/325). Source-backed GitHub Actions [#43339791418](https://github.com/zuizui0223/island/actions/runs/43339791418) **SUCCESS** and only anonymous [aggregate artifact #13998858022](https://github.com/zuizui0223/island/actions/runs/43339791418/artifacts/13998858022) was uploaded. Earlier validated runs #43338885931 and #43339432706 separately completed. Synthetic tests: 4/4 passed. The original species, island, geography, GIFT source entities and raw GBIF channel observations were never exported.

## Scientific question and answer

**Does a source-relative contraction/reorganization of recorded floral functional-proxy diversity occur on islands, and can it be attributed to local pollinator channel shortage or replacement?**

- **YES, descriptively, in a selected observational sample:** source-relative floral diversity and structural morphology differ geographically. In the **743 all-recorded northern-midlatitude islands**, both marginal and *joint* recorded floral diversity are lower than their corresponding five-nearest GIFT candidate native-mainland source plant trait pool under spatial-block and independent source-signature sensitivity (but these are only two overlapping, conditional spatial groupings, not source-independent plant evolutionary histories). Open-floral-form and actinomorphy prevalence are higher.
- **NO, guild disruption effects are NOT IDENTIFIED:** none of the **20** five-guild × four-region comparisons passes outcome-blind support/positivity requirements. The source-positive/adequately-nondetected guild state is NOT historical decline, a partner guild extinction or pollination-service loss. No per-plant pollen or fitness effect is measured.
- **NO, a single globally conserved contraction pattern is NOT established:** tropical broad-flora contraction lacks block/source supported uncertainty, southern-extratropical sampled floras show a different direction, and the strict native-qualified island subset substantially differs from the broad 743-island population.
- **NO, an actual pollinator shift is NOT identified:** the 21 source-positive Bombus adequately nondetected / bird-detected candidate islands show no source-relative floral FD or coordination increase robust to 14-spatial-block resampling.

The main finding is currently **source-relative observational floral assemblage reorganization under nonrandom ascertainment, with insufficient ecological partner-loss identification**, not proven pollinator-induced floral simplification.

## Locked outcome construction

Original source-verified Direct high/medium *species* evidence only: 2,234 accepted plant species with a single recorded flower colour, floral form, and symmetry state. 928 island floras have >=30 such species; 926 also have >=30 candidate source-pool species. Direct-only native records support 87 islands with >=30 complete species (67 northern midlatitude and 17 tropical, 3 other islands suppressed by minimum region/block support). Each island/source has 30 randomly sampled species without replacement, **30 rarefactions** at the verified final runner.

- `Q_marginal=(1/3)*sum_trait(1-sum_state(p^2))`, an **additive**, trait-pairing-blind categorical Rao Q.
- `Q_joint=1-sum_{observed colour×form×symmetry profile}(p^2)`, and `N2_joint=1/sum(profile p^2)` — respond to actual joint state prevalence (never expand ambiguous multi-state entries into hypothetical cross-products).
- `I(colour;form×symmetry) - exact fixed-margin permutation expectation` (the expected mutual information computed exactly from each cell's hypergeometric count, not noisy shuffled labels). Holding the form×symmetry pair intact avoids breaking ontological structure of floral architecture.
- Open-radial and actinomorphic species proportions as visitor-access *proxies*, not named pollination syndromes.
- Source-relative `delta=island - candidate GIFT pool` on the same observed Direct profiles and equal rarefaction size. The five geographically nearest GIFT flora regions are **plausible sources, not proven ancestral population or colonization pools**.

## Regional all-recorded floral outcomes: controlled observational source differences

**The following are block-equal mean island-minus-GIFT-source differences**; interval = 1,999 whole-spatial-block conditional bootstrap draws. A second separate 1,999 source-signature-equal draw on 66 nearest-mainland-source signatures is shown for the main North-midlat group. No genus/phylogeny/environment/availability-conditioned evolutionary selection inference is made.

| Region, broad recorded flora | Physical islands / spatial blocks | Marginal Q difference [block interval] | Joint Q difference [block interval] | Open-radial share difference [block interval] |
| --- | ---: | ---: | ---: | ---: |
| **Northern midlatitude** | **743 / 53** | **−0.02571 [−0.03707,−0.01508]** | **−0.01420 [−0.02189,−0.00687]** | **+0.04582 [+0.02357,+0.06756]** |
| Northern high latitude | 53 / 7 | −0.08249 [−0.11278,−0.05280] | −0.04614 [−0.06795,−0.02327] | +0.13718 [+0.08891,+0.18896] |
| Tropical | 106 / 42 | −0.00417 [−0.01175,+0.00335] | +0.00268 [−0.00265,+0.00859] | +0.03920 [+0.00415,+0.07743] |
| Southern extratropical | 26 / 13 | +0.04309 [+0.02395,+0.06310] | +0.03028 [+0.01486,+0.04510] | +0.17101 [+0.13144,+0.21069] |

The northern high cohort only includes **3 distinct nearest GIFT source signatures**, so block CI cannot be treated as independent source robustness; similarly the southern 26-island cohort is small. The tropical positive open-flower association does not survive the source-signature-equal bootstrap.

**66-nearest-source-signature-equal sensitivity (northern midlat):**
- Marginal Q mean **−0.01746 [−0.02781,−0.00771]**.
- Joint Q **−0.00985 [−0.01744,−0.00263]**.
- Open-radial share **+0.02904 [+0.00558,+0.05200]**.
- Actinomorphic share **+0.03965 [+0.02464,+0.05555]**.
- Colour×architecture **excess MI +0.00501 [−0.01103,+0.02263]** — **not source-signature robustly directional**, even though its island-equal median difference is negative. Thus *all dimensions do not collapse to one scalar syndrome or change uniformly*.

This is compatible with source-relative **floral assemblage filtering** but can reflect floristic turnover, taxa/provenance, mainland pool composition, unequal detection, environment and shared source species. It is NOT evidence that contemporary pollinators selected open flowers or that flower colour evolved locally.

## Falsification of universal or native-specific claims: identical-island comparison

Only **67 northern-midlatitude** islands have >=30 *strictly source-native* Direct triple-complete species. They are a highly selected subset of the 743 all-recorded islands. On the **SAME 67 ISLANDS**, the marginal Q source-relative median is **+0.01304 all-recorded** vs **+0.01546 strictly native** (not the −0.03877 median of the entire 743-island population). The within-island native-minus-all difference has a **block-equal mean +0.00136 [−0.00553,+0.00786]**.

For the same 67, joint Q difference native-minus-all **−0.00314 [−0.01237,+0.00663]**, open-radial share **+0.01627 [−0.00490,+0.03970]**; these intervals include zero. Consequently, the missing 743-island pattern in native-only analysis is **primarily a cohort-selection/coverage problem**, NOT demonstrated reversal caused by non-native botanical states. The same-island contrast is a posthoc record-origin sensitivity, not an introduced-plant experiment.

Within **17 tropical native-qualified islands**, native-minus-all marginal Q **−0.01094 [−0.01980,−0.00135]** across 10 blocks, but the cohort is small and status-selective; do not generalize a global native effect.

## Qualified pollinator channel opportunity: empirical support gate

| Qualified GBIF channel | Northern midlatitude (retained / adequately nondetected) | Tropical (retained / adequately nondetected) | Eligible guild×regime contrasts across 4 regimes |
| --- | ---: | ---: | ---: |
| **Bombus** | **480 / 5** | **2 / 16** | 0 |
| Other bees | 466 / 6 | 56 / 0 | 0 |
| Lepidoptera (moth+butterfly pooled) | 646 / 0 | 60 / 0 | 0 |
| Flower-visiting birds | 654 / 0 | 58 / 1 | 0 |
| Diptera | 562 / 1 | 53 / 2 | 0 |

Prospective minimum *associational* positivity: >=10 retained and >=10 adequately nondetected islands, each in >=10 source-geographic spatial blocks, and >=3 nearest GIFT source signatures supporting BOTH states in the same regime. **0/20** pass. Cross-regime extrapolation, zero-imputing unresolved channels and treating record counts as abundance are prohibited. This is a quantitative **causal-identification failure**, not a null ecological effect.

## Hypothetical Bombus→bird replacement opportunity: 21 candidate islands

Pre-frozen source-positive Bombus adequately nondetected AND island bird-detected cohort: **21 islands / 14 physical spatial blocks**. Equal-30 source-relative results, block-equal differences and 95% conditional intervals:

| Floral response | Block mean [95% conditional interval] |
| --- | ---: |
| Marginal Rao Q | **−0.00035 [−0.01644,+0.01542]** |
| Joint Rao Q | **+0.00689 [−0.00176,+0.01651]** |
| Joint Hill N2 | **+0.4766 [−0.5626,+1.5050]** |
| Colour×architecture excess MI | **−0.01729 [−0.05171,+0.01683]** |
| Open-radial fraction | **+0.01150 [−0.04897,+0.07133]** |
| Actinomorphy fraction | **+0.02672 [−0.01795,+0.07455]** |

**All intervals cross zero.** This neither confirms floral functional compensation nor says pollinators have no ecological effect: GIFT is not actual ancestor, bird detections need not be floral visitation, and the guild contrast lacks counterfactual overlap.

## Next empirical boundary

A real pollinator-*decline* study requires **dated**, within-island sampling effort and repeat guild detection/abundance with adequate nondetections; a *shift* study requires the same plant in island and potential mainland source **actual visitors** and time or comparable independent sites. Original GloBI v0.9 lacks event dates and 928-island exact-species GloBI coverage falls from 109 to 28 islands after the largest coarse archive-source omission; actual effective pollination remains unknown. Prior GloPL all-site isolation–pollen-limitation slope is a separate, globally supported but geographically/confounder-limited response, not an estimate of the island floral FD mechanism.

Species-fixed flower trait states only permit **floristic selection/assemblage reorganization**, not island-level within-species phenotype evolution. Build the future encounter histories and revisit H1–H3 with source-, taxon-, and observation-controlled repeated islands rather than promoting the currently unsupported pollinator shortage slope.

## Reproducibility

- [fd analysis source](../../src/fd/pollinator_functional_diversity.py) and [privacy-safe original-input runner](../../scripts/run_pollinator_functional_diversity.py).
- [synthetic null/guard tests](../../tests/test_pollinator_functional_diversity.py), including exhaustive fixed-margin shuffle expectation and equal-marginal/different-joint data.
- [original source workflow](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/.github/workflows/run-fd-functional-guild-impact.yml), successful [#43339791418](https://github.com/zuizui0223/island/actions/runs/43339791418), original aggregate [artifact #13998858022](https://github.com/zuizui0223/island/actions/runs/43339791418/artifacts/13998858022).

**PRs remain Draft. Scientific claims above are observational and conditional on the recorded taxon/flora evidence.**
