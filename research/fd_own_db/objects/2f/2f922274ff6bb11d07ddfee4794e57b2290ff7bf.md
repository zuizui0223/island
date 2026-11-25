# Predicting floral repertoires with taxon-held-out mainland and genus references

**8 November 2026.** Independent [fd](https://github.com/zuizui0223/fd) analysis on the source-audited [island PR #321](https://github.com/zuizui0223/island/pull/321) floral trait universe and geographically paired [continental/oceanic island study](../20261108_matched_flower_repertoire/README.md).

**Research claim boundary:** The empirical result is *conditional prediction of flower-colour × form state-set profiles given known genus membership*. It is NOT prediction of colonization, causal geological selection, within-lineage evolution, or equality of realized pollination function.

## 1. Why this predictive step is new and useful

The matched-island study found that exact recorded floral profile sets differ across geological island classes (All/Direct exact L2 under label swaps BY q≈0.00042), but its original categorical trait-similarity kernel did not reject. A stricter globally coherent **genus × regional-availability profile permutation** [did not support an excess global geological contrast](../20261108_genus_source_control/README.md), although simple subtraction of the GIFT nearby mainland category spectrum retained a descriptive exact-profile residual.

These three results address distinct estimands. Instead of treating a null non-rejection as mechanistic proof, we now ask: **Do genus and possible mainland source information predict previously unseen taxon-name floral profiles? How much skill is added by each information layer?** This requires a genuine name-level train/test split, rather than adding further post-hoc category p-values.

## 2. Evaluation population and leakage audit

Use the SAME 37 matched continental–oceanic island pairs (23 northern-midlatitude, 14 tropical; **74 unique islands**) from the previous fd research line. Four source models are evaluated in nested All/Direct evidence scopes.

| Quantity | All | Direct |
|---|---:|---:|
| Tested island × complete taxon-name occurrences | 13,141 | 11,664 |
| Distinct taxon names held out from EVERY GIFT donor pool | 3,316 | 2,605 |
| Remaining distinct complete-floral GIFT mainland reference taxa | 6,321 | 4,172 |
| Globally observed floral profile reference types | 496 | 496 (467 present in Direct) |
| Donor genera in GIFT reference | 1,517 | 1,466 |

The original native-mainland dataset has 748 geographic entities and 1,666,186 source-region × taxon rows. All taxon names appearing on any tested island are removed **globally from the training reference within the same evidence scope**, even when the same taxon is listed in multiple mainland entities. Test records never enter genus donor counts or the global profile prior. This is an out-of-*taxon-name* prediction check, not an evolutionary or phylogenetically independent external validation; donor species are still close relatives.

Complete island species-name support is fixed, so models condition on **which genera were actually recorded on islands**, not on genera predicted to arrive.

## 3. Four fixed predictor models

A tiny fixed pseudo-count of 0.02 per global profile avoids zero probability for unseen categories. Main hierarchical shrinkage is 8 pseudo donor taxa. Recorded multistate colour/form descriptions remain **single set-valued profiles**, never expanded to unobserved colour × form combinations.

| Model | Mainland reference for each held-out island taxon |
|---|---|
| Global mainland | GIFT-native flora pooled globally, without genus |
| Nearest mainland | Profiles in the island's five nearest GIFT native-mainland source entities, without genus |
| Global genus | Profiles of **other GIFT mainland species in the same genus**, shrunk to global prior |
| Nearest-mainland genus | Same-genus profiles among the five nearest mainland entities, shrunk to global-genus model if sparse |

The main evaluation uses two complementary losses on the *actual* held-out colour/form profile: **categorical proper log loss** (smaller means more exact-profile predictive information) and **expected original set-Jaccard trait distance** (smaller means more morphologically similar predictions). We also compare approximate predicted island unique-profile incidence with observed incidence (L2 and categorical-kernel errors). The latter uses normalized expected occupancy and is explicitly an approximation.

## 4. Model performance: actual held-out data

| Scope | Model | Held-out log loss ↓ | Original colour/form distance loss ↓ |
|---|---|---:|---:|
| All | Global mainland | 5.3563 | 0.7698 |
| All | Nearest mainland | 6.2872 | 0.7658 |
| All | **Global genus** | **4.6474** | 0.6256 |
| All | Nearest-mainland genus | 4.6618 | **0.6115** |
| Direct | Global mainland | 5.4687 | 0.7792 |
| Direct | Nearest mainland | 6.5307 | 0.7710 |
| Direct | **Global genus** | **5.1305** | 0.6950 |
| Direct | Nearest-mainland genus | 5.1498 | **0.6843** |

### Improvement over global mainland profile frequencies

37 matched pairs are resampled 4,999 times to form *pair-conditional* percentile intervals (pairs can share species; they are not phylogenetic independent replicates).

| Skill improvement (positive = lower error) | All | Direct |
|---|---|---|
| **Global genus: log-loss improvement** | **+0.709** [0.644, 0.774] | **+0.338** [0.309, 0.369] |
| **Global genus: trait-distance improvement** | **+0.144** [0.138, 0.151] | **+0.084** [0.080, 0.089] |
| Nearest genus vs global genus: log-loss improvement | **−0.0145** [−0.0212, −0.0077] | **−0.0194** [−0.0230, −0.0156] |
| Nearest genus vs global genus: trait-distance improvement | **+0.0140** [0.0123, 0.0158] | **+0.0107** [0.0089, 0.0125] |

This is a quantitative, independently held-out **genus predictive-skill result**. Geographic mainland proximity on top of genus membership is much less reliably informative about **the exact** floral profile, even when it makes predicted profiles somewhat closer in categorical trait space. The identity-level match can worsen while trait-similarity prediction improves; these are not equivalent endpoints.

### Geographic source coverage is limiting

| Reference genus available | All continental/oceanic | Direct continental/oceanic |
|---|---|---|
| Across all GIFT native-mainland reference species | 73.7% / 72.7% | 66.3% / 65.0% |
| In five nearest mainland entities | **38.6% / 36.6%** | **29.2% / 27.4%** |

Thus much of the nearest-mainland genus model **falls back to global genus**. This is an observational source-coverage limitation, not proof that a genus is actually absent from the historical mainland source region. It also explains why the earlier ~92% overlap of broad flower-profile types does **not** imply adequate lineage-conditioned donor coverage.

## 5. Prespecified family of predictive residual tests

The original geological contrast is compared with **observed minus model-predicted normalized profile incidence** for each island; within-pair label-sign swaps are exploratory and assume conditional exchangeability of geological classes after geographic matching. The stated model is nearest-mainland genus, four tests (All/Direct × exact profile L2/category similarity), 9,999 swaps.

| Residual geological contrast | All raw p | Direct raw p | Inference |
|---|---:|---:|---|
| Exact recorded profile L2 | **0.0001** | **0.0001** | BH q=0.0002 and BY q≈0.00042 (four-test family) |
| Original category similarity kernel | 0.1109 | 0.1302 | No adjusted support |

The **exact-category residual persists** in this predictive model, but the original taxon/genus permutation null's pooled global non-rejection remains valid. These are distinct references: a source-trained, incompletely covered **predictive mean** need not reproduce the null distribution from coherent shuffling among the observed island taxa. Do not call this result source/genus-independent evolutionary differentiation or geological causality. Prediction error/coverage bias could generate a systematic residual.

## 6. Native-only and source/smoothing sensitivities

A fixed-status sensitivity keeps **28** original matched pairs (20 north-midlatitude, 8 tropical), with at least 20 explicitly native complete names on every endpoint under All/Direct. Retraining the models after holding out those names gives global-genus log-loss improvement of **+0.875** All and **+0.507** Direct; trait-distance improvements **+0.164** and **+0.105**. However this smaller, status-selected population has strongly unequal observed native name counts in continental versus oceanic islands and cannot be generalized to all island floras.

For the full 37 pairs, changing the source from 5 to 10 and 20 nearest GIFT entities **increases the incremental trait-distance benefit** of nearest-mainland genus profiles (All from +0.0140 to +0.0174 and +0.0197; Direct +0.0107 to +0.0136 and +0.0154), but the **exact-profile log-loss penalty persists or increases** under shrinkage=8. Crucially, changing shrinkage from 8 to 20 **reverses the small log-loss penalty into a small benefit** (+0.0150 All, +0.0063 Direct), while shrinkage=2 worsens it. Therefore **claim only robust genus-global predictive gains**, not a universally beneficial or harmful mainland-distance effect.

All 5/10/20-source, 2/8/20-pseudo-count and native-only sensitivity values are in [model_sensitivity.csv](model_sensitivity.csv). These are **post-hoc correlated sensitivity checks**, not separately corrected confirmations.

## 7. What this resolves and what remains

**Established as a conditional prediction result:** Floral colour/form profiles exhibit strong **genus-associated predictability** across these 74 island assemblages, on taxon names withheld from the GIFT training records. The added geographic source layer can improve trait-distance approximation, but its exact-profile value is highly smoothing-dependent and constrained by missing donors.

**Not established:** Which genus actually arrives or goes extinct; historical source ancestry; plant pollinator-channel effects; geological-origin selective pressure; phylogenetically independent floral trait convergence; equivalence of reproductive function; that genus composition uniquely *causes* the observed oceanic–continental difference. The older globally genus-conditioned null did not show clear additional origin separation.

**Most informative next test:** A full plant-arrival/establishment model with independently sourced dispersal opportunity and archipelago/phylogeny structure, coupled to separately measured native status and visitation/reproduction. The present data identify the **predictive constraint of known lineage composition** but not its colonization mechanism.

## 8. Reproducibility and rights

- Executable, private-input-only module: [`src/fd/heldout_mainland_prediction.py`](../../src/fd/heldout_mainland_prediction.py).
- Strict [methods and held-out-scope contract](../../docs/HELDOUT_MAINLAND_PREDICTION.md).
- Source statistics: [main model scores](main_model_scores.csv), [paired improvements](paired_score_improvements.csv), [GIFT genus coverage](source_genus_coverage.csv), [origin residual diagnostics](origin_predictive_residuals.csv), and [post-hoc sensitivity](model_sensitivity.csv).
- Source inputs: original `island` PR321 frozen trait/archive, GIFT 3.2 native mainland lists and precomputed geography-only source assignments, and independent Roeble/Zizka geological-origin crosswalk. **Raw identity-level source records were processed locally and never published or passed to public Actions.**
- Taxon-native candidate record identities, original island IDs, source entity IDs and original provider data remain private according to the [island redistribution-rights audit](https://github.com/zuizui0223/island/blob/3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5/docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md). Public GitHub CI validates the implementation only with synthetic fixtures.
