# 2026-11-13 · Are globally rare flower-profile neighborhoods under-covered on islands?

**Original-data 928-island macroecological analysis — source-backed exploratory, NOT evidence of observed pollinator decline or evolutionary floral trait loss.** Scientific branch [fd PR #16](https://github.com/zuizui0223/fd/pull/16), protected original data [island PR #326](https://github.com/zuizui0223/island/pull/326), independent original source run [Actions #43378963646](https://github.com/zuizui0223/island/actions/runs/43378963646) **SUCCESS**; [anonymous original result artifact #14016260897](https://github.com/zuizui0223/island/actions/runs/43378963646/artifacts/14016260897). Before reporting, one source run initially failed because unsupported source-equal groups emitted NaN; fixed to publish JSON **null = unestimable**, never zero. The final job, synthetic tests and source-restricted original input validation passed.

## Target question

The earlier categorical finite-metric volume analysis showed a stronger source-relative deficit when all 100 recorded flower-profile categories receive equal global reference mass, but no supported deficit under weighting those categories by the number of directly characterized plant species. **Where in the recorded three-trait categorical reference is the difference concentrated?**

Here **rarity means original directly documented worldwide Direct species carriers per flower colour×form×symmetry profile**, NOT ecological function, phylogenetic rarity, pollinator specialization, threatened species or true island absence.

### Original source-verified cohort and frozen bins

Original 2,234 distinct source high/medium directly measured species, 100 distinct *actually recorded* three-trait categorical profiles. Original 928 islands with at least 30 Direct trait-complete plants; candidate five-nearest GIFT mainland **native source flora** are geographic comparators, not independently demonstrated ancestral sources.

Frequency thresholds were defined from the 100-profile original global source counts **before any island effect** and preserve ties:
- **Rare:** globally only **1–2 directly recorded plant species/profile**, **29 profiles**.
- **Common:** globally at least **26 Direct species/profile**, **25 profiles**.
- **Intermediate:** **46** remaining profiles, a descriptive reference stratum.

These are disjoint, exhaustively cover the 100 original profiles and use no estimated GloBI visitor guild labels. Values are fractions WITHIN each reference bin, so a one-profile difference counts 1/29 in rare and 1/25 in common.

For each original island and potential five-nearest GIFT source, draw **30 original directly characterized plant species** at random 40 times (same n and original source state). The rare or common reference-coverage fraction is the share of original source profiles in the bin that lie **within one categorical trait mismatch** of at least one recorded local profile. This is *categorical metric-ball coverage* and no Euclidean/KDE smoothed volume.

Two independent geometry nulls condition on **the exact number k of DISTINCT observed profiles per thinned island/source**: 399 original-reference sets of k uniformly sampled profile centers, and another 399 k-sized distinct profile draws with probability proportional to number of original worldwide Direct species carriers. The rare-vs-common difference after each k-specific expectation is therefore not a re-description of unique observed profile count, but neither k null controls genus/family/source provider.

Primary strict sensitivity removes **every accepted plant species observed on the focal island from all five candidate mainland source flora**. This excludes exact species overlap only; it does NOT match genus composition. Conditional equal-weight whole physical 10° block and first-nearest-GIFT-source-signature bootstrap distributions have 1,999 resamples. Their dependence is not equivalent to independent botanical evolutionary origins.

## Original-source aggregate results (no protected records released)

**Primary source exclusion** (island minus source; each rare/common fraction is within its fixed globally defined bin):

| Direct original macro region | Source-disjoint islands | Original physical blocks / source signatures | Rare-profile neighborhood coverage difference | Common-profile neighborhood coverage difference | Rare-minus-common gap difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| **Northern midlatitude** | **736** | 53 / 66 | **−0.02643 [−0.04883,−0.00552]** | +0.00385 [−0.00007,+0.00798] | **−0.03027 [−0.05195,−0.01060]** |
| **Tropical** | **106** | 42 / 36 | **−0.02845 [−0.04670,−0.00976]** | **+0.00305 [+0.00020,+0.00604]** | **−0.03150 [−0.04859,−0.01295]** |
| Northern high latitudes | 53 | 7 / **3** | −0.10026 [−0.13790,−0.05830] | +0.00153 [−0.01016,+0.01304] | −0.10179 [−0.13909,−0.05777] |
| Southern extratropical | 24 | 13 / 11 | +0.02080 [−0.01728,+0.06050] | +0.00969 [+0.00281,+0.01735] | +0.01111 [−0.02898,+0.05021] |

Main table numbers are equal-physical-spatial-block mean and conditional bootstrap 95% intervals, **not effect sizes of verified pollinator extinction**. North high has only **3** first-source GIFT signatures; not source-independent generality. Southern rare-minus-common contrast is indeterminate.

### The key question: is the rare-end gap simply k original recorded profile richness?

No, **under the prespecified geometric reference nulls**. With exact source species exclusion:

| Region | Rare-minus-common *excess* gap after uniform exact-k profile reference | 95% whole-block conditional interval | Excess gap after global species-frequency-weighted exact-k reference | 95% whole-block conditional interval |
| --- | ---: | --- | ---: | --- |
| **Tropical (106 islands)** | **−0.03061** | **[−0.04520,−0.01377]** | **−0.02870** | **[−0.04236,−0.01487]** |
| **Northern midlatitude (736)** | **−0.02633** | [−0.04557,−0.00749] | **−0.01962** | [−0.03621,−0.00464] |
| Southern extratropical (24) | +0.00896 | [−0.02823,+0.04641] | +0.00450 | [−0.03505,+0.04336] |

For tropical source-signature equal weighting, the stronger frequency-aware exact-k residual remains **−0.03053 [−0.04525,−0.01579]** across **36** candidate mainland signatures. Corresponding northmid source-equal result: **−0.02490 [−0.03870,−0.01205]** across 66 source signatures.

**Descriptive result:** geographic/source-nested original tropical and northmid floras lie farther from *globally infrequently represented recorded profiles* than source counterparts, beyond original profile count in the two operational reference null models. This result is most informative about the rarer end of a documented **categorical phenotype** repertoire, not demonstrated ecologically specialized or successful pollination functions.

### Critical ceiling/ascertainment issue

On tropical islands the 30-species-rarefied **common-profile neighborhood coverage has median 0.997**, almost complete, whereas rare reference coverage median ≈**0.691**. Thus a positive common-profile coverage difference has very little headroom. Comparisons involving the common side can partly reflect a **mathematical ceiling**, even though within-bin normalization and exact-k nulls constrain profile-number artifacts. A change in the categorical metric or radius could change these numbers. Do not call this a universal rare pollination-function loss.

### Critical genus and native-provenance confounding

The original source-excluding-island-species GIFT lists can replace every one of a complete **30-species** floral sample with a **different source species in the same genus** for only:
- **2/106 tropical original islands**;
- **416/743 northern midlatitude** islands;
- 27/53 northern high and 0/26 southern extratropical.

Relaxing to a **separate predeclared 15-species restricted matched-genus estimand** yields only **14/106 tropical** source-matchable islands, 655/743 northmid. That restriction selects a nonrandom subset of island floral records and mainland source genera. A genus-conditional comparison must be reported as a *restricted separate population*, not an adjustment to 106 tropical islands. Because identical species are excluded and genus-count matching is assessed without reading their floral profile states, the support gate is outcome-blind, but it does not repair taxonomic or source-quality ascertainment.

Strict source-native original trait-complete island subset: **17 tropical islands**. Their rare-minus-common k-controlled source-disjoint block estimate **−0.00630 [−0.03867,+0.02416]** has no supported directional contrast. A selected original **67 northern midlat** native-qualified cohort has +0.01127 [−0.01726,+0.04345], also indeterminate. The broad all-recorded result is **not** confirmed for strict directly recorded native-only flora.

## Result status, limitations and defensible paper claim

**Confirmed computationally:** source-backed original 928 islands, 2,234 Direct source species/100 triple profiles, real source exclusion, explicit equal-30 floristic sampling, k-specific uniform and species-frequency reference nulls, block and original first-source bootstrap; CI protected input and public aggregate-only output.

**Supported exploratory biological description:** tropical and northmid documented floras cover fewer **globally rare categorical flower-profile neighborhoods** than five-nearest GIFT mainland candidate flora after controlling recorded species number and exact profile-type number k.

**NOT identified:** change within a species' floral phenotype; phylogenetically independent flowering-plant lineages; source genus/family-conditioned difference on all tropical islands; true pollinator guild scarcity/shifts; missing pollinator functions, pollen deposition, seed set, or causal island floral trait adaptation. The choice of equal weights for colour, partly definitional form/symmetry and global study/provider Direct labels may matter.

A genus-matched 15-species restricted-cohort follow-up is recorded separately; it must NOT be projected to the original 106 tropical islands if inadequate.

## Source audit and files

- [fd reusable frequency-decomposition code](../../src/fd/rare_profile_decomposition.py); [synthetic tests](../../tests/test_rare_profile_decomposition.py).
- [Protected original source run](https://github.com/zuizui0223/island/actions/runs/43378963646) (SUCCESS), [rights-safe original result artifact](https://github.com/zuizui0223/island/actions/runs/43378963646/artifacts/14016260897), source [original island code](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_rare_profile_support.py).
- [Original 928-island metric support volume report](README.md), [predeclared study contract](../../docs/DISCRETE_CATEGORICAL_FUNCTIONAL_HYPERVOLUME.md). The two GitHub PRs stay **Draft** and must not be merged as a causal evolutionary finding.
