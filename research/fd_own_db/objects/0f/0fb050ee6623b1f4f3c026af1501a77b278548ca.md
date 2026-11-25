# Global island floral reassembly — original mainland-source leave-one-donor-out negative control (2026-11-17)

**New original-source numerical test completed (SUCCESS).** Source-protected original GitHub Actions [#44074512032](https://github.com/zuizui0223/island/actions/runs/44074512032), anonymous artifact ID **14324571092**, extracted anonymous_original_triple_matched.json **SHA256 16acf3440b2e5d4ccbe94d0cff96879f7fdf52c79395842c42dcfc8488376885** (286,603 bytes). The original botanical population remained **4,505 GIFT islands**. The source provenance, native-status tags, 748 original mainland source areas, Direct species-label hashes and independent GSHHG geographic source were checked in CI. The analysis uses exactly the same accepted focal island plant species and exact same different-accepted-species candidate mainland donor species for BOTH flower-architecture responses. The original full island/species/genus records stay in protected source CI; only anonymous aggregated results are released.

## Question: does island floral architecture have EXTRA genotype-taxonomic predictability beyond the mainland baseline?

An earlier original same-plant analysis found a favorable Brier forecast gain for actinomorphic x deep tube relative to tubular/bell x deep tube when using the island's observed genus proportions and potential mainland species-Direct source donors. The old result is real but did not establish that the asymmetry is **island-specific**. It could be a generic botanical genus-associated trait pattern inherited from the mainland.

The new negative control uses a subset in which every focal genus has **at least two different accepted triply Direct source donor species in BOTH source definitions**. For each donor plant in a target genus, remove that one donor from the trained genus-pool and entire matched source pool. Evaluate (1) the SAME real focal island plants and (2) the held-out mainland donor plant. For each genus, average over the SAME omitted donor plant for the island prediction and the mainland test. Use the original island target genus weights for both, and the same botanical species for the actinomorphic x deep and tubular/bell x deep outcomes.

Within each island define:

- D_island = [Brier(pool) minus Brier(genus)] for actinomorphic/deep minus the same gain for tubular/deep, scoring **island focal plants** with donor-LOO source predictions.
- D_mainland = the same two-trait gain contrast, scoring **held-out mainland donor plants** under the matched donor-LOO source predictors, and reweighting genera by actual target island genus frequencies.
- The potential island-specific component **E = D_island minus D_mainland**. Positive E means the *difference in genus predictability between flower endpoints* is larger for island plants than matched source plants.

The LOO control equalizes the number of source plants available for the genus/pooled predictors for both the island and pseudo-mainland endpoints; ordinary full-donor genus predictions should not be subtracted from source-LOO to claim an island effect. The original 4,505 island population is not replaced with this narrower diagnostic cohort. Regions × nativeness × true-triple/deep-only × mode × documentation-support comprise **120 descriptive or non-identifiable reporting cells**, not 120 independent discoveries. Minimum support: 30 independent island floras / 10 geographical 10-degree blocks / 10 independently grouped first GIFT source entities. Conditional whole-group 999 bootstrap intervals are exploratory (no pre-registration or multiple testing).

## Original primary all-recorded, same triply Direct accepted plant / donor cohorts

All values below are **percentage-independent Brier MSE score differences**, not flower frequency percentage-points or evolutionary rates.

| Region and original source definition | Matched botanical islands / physical blocks / first GIFT source groups | D_island, equal physical block | D_mainland, equal physical block | **E=additional island excess [group bootstrap 95%]** |
|---|---:|---:|---:|---|
| **Tropical, geo5** | **69 / 27 / 21** | **+0.122057** | **+0.117545** | **+0.004512 [-0.035620,+0.045835]** |
| Tropical, climate10; **same island cohort** | 69 / 27 / 21 | (mode-specific) | (mode-specific) | **+0.015278 [-0.028225,+0.060567]** |
| **Northern midlatitude, geo5** | **873 / 51 / 73** | **+0.062565** | **+0.059538** | **+0.003028 [-0.010922,+0.015136]** |
| Northern midlatitude, climate10; same island cohort | 873 / 51 / 73 | (mode-specific) | (mode-specific) | **+0.000545 [-0.011608,+0.012414]** |

**First-source-group-equal sensitivity**, original geo5: tropical E **+0.009043 [-0.030328,+0.052749]**; northern-midlatitude E **+0.003959 [-0.008611,+0.017219]**. Exact shared accepted target and source names remain in the same respective region/source mode. Thus conditional *endpoint-dependent botanical genus predictability* is present in both the island floral assemblage and the source-only mainland negative control. **There is NO detected additional island-specific difference in the all-recorded cohorts**; the tropical interval still permits moderate effects in either direction, so this is NOT equivalence testing or proof of zero ecological change.

**Documentation-balanced control:** northern-midlatitude geo5 n**36** original eligible floras / 11 geographical blocks / 18 first source entities had within-10 percentage-point Direct documentation differences in BOTH mainland source definitions. Even here E **-0.008394 [-0.037164,+0.013303]** (block equal; source-group-equal +0.000518 [-0.019595,+0.016648]). Tropical all-recorded n69 had **zero** source-island original triple-Direct documentation-balanced cases. This severe original observation bias prohibits an unqualified global claim.

**Known-native exploratory exception**, northern-midlatitude original accepted documented natives n**90** /19 physical blocks/25 first source groups have geo5 E **+0.028569 [+0.006229,+0.055948]**, but climate10 E **+0.017340 [-0.004719,+0.044555]**, and only **29** islands meet source/island documentation balance (below the pre-fixed 30 island gate). Native nonendemic n85 geo5 E +0.009192 [+0.000853,+0.023770], climate10 +0.002903 [-0.009148,+0.017147], with 27 balanced. These are overlapping selectively documented botanically native subsets; **do NOT relabel as independent replicated island adaptation** or a confirmatory effect. Tropical native scope n10, native-nonendemic n9 are **NOT_IDENTIFIABLE** and should not be assigned zero.

Deep-tube-only all-recorded control: tropical n11 islands, northern midlatitude n5, both **NOT_IDENTIFIABLE**. Northern highlatitude all-recorded n72 /9 physical blocks /5 first sources also fails. Southern extratropical all-recorded zero supported; no inference about those climates from this diagnostic.

## What this does to the one-manuscript interpretation

**Positive original findings still stand separately**: global accepted island floras differ taxonomically from potential GIFT mainland candidates beyond naive finite source-name draws; observed floral-category repertoires overlap across taxonomic replacement; same-species Direct radial/deep categories show regional source-relative directionality; the original island genus weights help predict certain joint floral architectures. However:

1. The observed floral *architecture-dependent genus predictability* itself appears to be a **general source-flora genus-associated signature**, not a newly demonstrated island-specific restructuring beyond the mainland source donor baseline.
2. Original mainland candidate GIFT checklists are potential source opportunities and are NOT observed seed rain, true colonization, past ancestral floras or measured genus establishment rates. The genus predictor uses **known island genus membership**.
3. The difference between island vs mainland LOO control outcomes is not proof of equal real flower morphology or absence of selection. Current Direct phenotypes are species-wide categories with major provider/status missingness; provider-disjoint independent island morphology and true guild-specific service data remain outstanding.
4. Do NOT call this a confirmed guild/flower pollination syndrome shift, all-island adaptive narrowing, 17D KDE hypervolume reduction, or actual same-species field-level island evolution.

**Recommended main claim:** island floras undergo taxonomic reassembly, and documented functional profiles partly reflect a genus-associated structural template available in the mainland donor species, with genuine regional variation in floral frequency. This is an *inheritance/assembly-versus-shift identifiability* result: not a novel adaptive shift proven from those Direct source documents.

The manuscript retains **all 4,505 GIFT island floras** as its study population. This strongly selective strict triply Direct source-LOO test belongs as a falsification/sensitivity panel rather than the global gating criterion.

## Reproducibility

- [Source-only donor-LOO matched botanical predictor](../../src/fd/source_loo_flower_placebo.py), [source-null unit tests](../../tests/test_source_loo_flower_placebo.py), and [true three-character same-species source integration](../../src/fd/triple_matched_congener.py).
- [Original protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_triple_matched_congener.py), [original source-protected workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-triple-matched-congener.yml). Passed exact source hash checks, original population/anonymous-output invariants and tests. Original Actions #44074512032 SUCCESS, anonymous artifact #14324571092, exact file SHA above.
- [Earlier same-target-source Brier skill analysis](../20261117_original_same_species_triple_genus_skill/README.md): retain numerical content, but supersede any island-*specific* interpretation with the present negative control.
