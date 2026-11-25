# 2026-11-16 — ORIGINAL GIFT 4,505-island source-only congener heldout prediction of genuine two-trait floral architecture

**Completed source-protected original data run:** [island Actions #43868038587](https://github.com/zuizui0223/island/actions/runs/43868038587) SUCCESS, 16 independent synthetic tests passed and original accepted GIFT island botanical flora / source region / high-medium Direct trait SHA validation completed. Anonymous original [artifact #14226687228](https://github.com/zuizui0223/island/actions/runs/43868038587/artifacts/14226687228), exact extracted source JSON SHA256 `1bd107e79e937b69a57684ec0b455ef77677a94b8a7aa864315951834cd5643c` (514,295 bytes). Prior [initial non-factorized pilot #43866414867](https://github.com/zuizui0223/island/actions/runs/43866414867) also SUCCESS. [Original runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_source_only_joint_congener_prediction.py) · [reusable models](../../src/fd/source_only_real_joint_congener_prediction.py) · [synthetic tests](../../tests/test_source_only_real_joint_congener_prediction.py).

This is an **exploratory, AFTER-previous-results model comparison**, not preregistration, trait-provenance independence, a validated pollination syndrome, real bird/butterfly/moth visitors, within-species trait evolution or causal source recolonization. The original accepted species, private island/source names, coordinates and observer evidence never leave protected island CI.

## Design: what was actually held out

Preserve the original **4,505 accepted island botanical floras** as the population (includes all insufficiently Direct-characterized islands). In each island, **exclude every accepted species name anywhere on the entire focal island flora** from each of TWO independently frozen native GIFT potential mainland candidate rosters: five geographically closest GIFT mainland regions (`geo_k5`) and ten climate-nearest of fifty geographical candidate regions (`geo50_climate10`). Both remain *potential checklist donors*, not observed seed rain, historical ancestors or actual independent plant phenotype providers.

Only original real SAME-accepted-species paired Direct phenotypes were used for **actinomorphic×deep tube**, **zygomorphic×deep tube**, and **tubular/bell×deep tube**, for three NESTED provenance scopes (all-recorded, known-native, known-native nonendemic). Keep each original target island accepted-species set fixed across both source modes. A pair can only be predicted on target plants whose genus has at least one separately accepted source-only Direct trait-pair donor in **BOTH** modes (>=2 donors per genus both modes in the stricter nested test). Main supported regional panel requires >=5 original target directly scored accepted plants, >=2 actual donor-matched botanical genera, >=10 original candidate source-only Direct donors under each source mode, >=30 original islands and >=10 independent physical 10° blocks. This is **90 region×scope×pair×donor-gate reporting panels**, not independent evolutionary events.

At the exact same ORIGINAL target accepted plants compare:

- Source-global `pAB`: overall original candidate species-specific real joint-state fraction;
- Source shared-genus-POOL `pAB`: all source-only candidate accepted species in eligible donor genera, pooled without using the focal island's genus mix;
- **Source genus-WEIGHTED real joint `pAB`**: independently recorded actual joint-state source-only donor fraction for each accepted botanical genus, weighted by the NUMBER of eligible focal island directly documented accepted species in that genus.
- **Source genus-WEIGHTED FACTORIZED prediction `pA_g * pB_g`**: uses only separately documented trait-A and trait-B donor frequencies within each genus, NOT their genuine observed same-plant joint pairing.

Primary prediction comparison uses original known focal species' binary joint trait status scored by Brier MSE, lower is better, and produces `Brier_pooled - Brier_genus_joint` (gain of genus assembly knowledge relative to same donor species pooled). The **syndrome-specific, additional predictive gain** is `Brier_genus_factorized - Brier_genus_joint`. It is only independently informative when source donor genera have variation in BOTH original traits; singleton or fixed-trait source genera mechanically make these two predictions the same. We explicitly count per-island genus opportunity. A parallel source-only congener prediction of coupling `P(A,B)-P(A)P(B)` is a description, not a verified visitor system.

Original 999 whole-physical-block bootstrap and independent first GIFT source-group sensitivity. Include an additional STRICT documentation-balance sensitivity requiring island source-supported-genera proportion Direct-known and source-only both-mode congener Direct coverage each within **10 percentage points**; do not fill missing trait labels. The known-native scopes and the documentation subset are nested and nonrandom.

## Main tropical all-recorded source-only prediction: taxonomy YES, extra pair synergy NO

Exact COMMON-source and COMMON-target sample using **at least 1 source donor species/genus in both modes**, original ten-degree-block equal island means:

| Real source-only paired structural response | Original qualified island floras | Original median target Direct plant species / genera | Brier gain genus JOINT vs source matched-genus POOL, geo5 [95%] | Additional Brier gain genuine same-donor-species JOINT vs independent within-genus A×B marginals, geo5 [95%] |
|---|---:|---:|---|---|
| **Actinomorphic × deepest tube** | **331** | **12 / 8** | **+0.11872 [+0.10286,+0.13290]** | **+0.000054 [−0.000016,+0.000156]** |
| Zygomorphic × deepest tube | 331 | 12 / 8 | +0.02712 [+0.01804,+0.03639] | +0.000005 [−0.000193,+0.000135] |
| Tubular/bell × deepest tube | 260 | 11 / 7 | **−0.03975 [−0.05418,−0.02659]** | +0.00120 [−0.00194,+0.00300] |

Original alternative `climate10` candidate source mode: radial×deep genus-versus-pooled Brier gain **+0.12567** (geo5 +0.11872); **joint-versus-factorized marginal gain only +0.000043**, interval spans zero. The source-group conditional bootstrap also gives positive genus-versus-pooled radial×deep improvement +0.14369 (95% [+0.12429,+0.16357]), but does not demonstrate taxonomically independent island evolution. Tropical tube/bell×deep has a negative genus-versus-pooled gain under BOTH original source modes, despite the source genus trait library clearly being useful for radial×deep. These are MORPHOLOGY-specific associations, not proof that pollinators choose plant genera.

**Stricter >=2 actual accepted mainland source donors/genus in EACH mode**: n225 tropical radial/deep islands, genus-vs-pool Brier gain geo5 **+0.12987 [+0.11365,+0.14553]**, but true joint-vs-factorized gain only **+0.000025** [−0.000037,+0.000102]. Tropical tube/bell/deep n152: genus-vs-pool Brier gain geo5 **−0.01553 [−0.03068,−0.00115]**, with a small joint-vs-factorized +0.00362 [+0.00264,+0.00461]. The latter is exploratory, often selected, not a new named pollinator syndrome or proven causal floral pairing.

## Most important original data support and missingness limits

| Tropical original all-recorded, ≥1 congener/genus in both modes | Radial×deep | Bilateral×deep | Tube/bell×deep |
|---|---:|---:|---:|
| Pair-comparable original island floras | **331** | **331** | **260** |
| Original island floras with ≥1 source genus variable in BOTH traits under BOTH source modes | **13** | **16** | **202** |
| Original island floras with source-vs-island recorded trait documentation gaps <=10pp under BOTH modes | **6** | **6** | **4** |
| Original islands satisfying BOTH the variation and ≤10pp documentation criteria | **1** | **2** | **4** |

The separate source-genus variation subset for tropical radial×deep (13 islands, fewer than 30) is **NOT_IDENTIFIABLE**, and the independent trait-documentation balance subset has only SIX ISLANDS, so *none* of the seemingly strong general radial×deep predictive gains can be claimed robust to this important source recording-bias limitation. For tropical tube/bell×deep, while 202 islands have at least some source-genus variation, only four also pass source-vs-island documentation balance, leaving its small extra predictive synergy scientifically unresolved.

Known-native tropical radial×deep (n49, 18 blocks) has original genus-vs-pool Brier gain **+0.10392** geo5 and +0.10934 climate10, but balanced Direct documentation in just **2 islands**; native-nonendemic n49 has **zero** balanced cases. Neither guarantees external generalization under independent literature/observer sources or among island endemics.

In northern midlatitude all-recorded n1227 radial/deep, positive source genus-vs-pool Brier gain geo5 **+0.04583** (95% [+0.03013,+0.06269]), but documentation-balanced n9. Southern-extratropical has n19 radial/deep congener-eligible original floras across 8 blocks, below the fixed >=30 island /10 block reporting gate. High latitude n98 radial/deep has no BOTH source-side variable genera or documentation-balanced floras, so it is not evidence for lack of functional integration.

## Reconciliation with previous in-sample genus-fixed null and scientific claim ceiling

The prior original in-sample within-genus fixed-marginal label exchange on all tropical n666 real symmetry/deep source-qualified islands reconstructs the source-relative coupling `C`: observed `actinomorphic×deep` +4.018 percentage points versus genus-fixed expected +4.057, apparent residual −0.039, with **zero island/source BOTH-variable-genus-supported cases**. This reference had the focal island's trait label margins IN its null. The NEW, different accepted-name SOURCE-ONLY donor results demonstrate that genus membership can contain practical predictive information about island real pairs on a selected n331 common-donor subset **without directly fitting focal island phenotype labels**, but they DO NOT show that the *co-occurrence linkage* adds skill beyond marginal genus traits. The sample/population estimands (n666 vs n331) must not be subtracted or called replicated.

**Defensible distinction**: source candidate taxa plus the observed island's botanical genus membership help predict island radial/deep architecture. However the additional predictive value of the *source-only real two-trait pairing* above separate genus-level morphology traits is negligible/inconclusive in actual source-supported comparisons, especially after support and documentation caveats. This favors a careful description of **trait-specific botanical composition and genus-associated floral source predictability**, not identified bird/butterfly/sphingid pollinator shifts, island within-species evolution or a universal floral syndrome collapse.

**Next empirical dependency:** Acquire truly independent provider and locality-resolved island and mainland flower morphology AND effective plant-specific visitor observations. Original GloBI exact visitor evidence may be explored only with species/date/reference/geography qualifications; don't relabel non-detections as pollinator extinction. Keep [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and upstream original source branch DRAFT.
