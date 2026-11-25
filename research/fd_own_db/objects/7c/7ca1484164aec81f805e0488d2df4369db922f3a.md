# 14 November 2026 — Phylogenetic beta diversity versus floral functional beta diversity across original islands

**Original-source DONE, TPL/LCVP/WP three dated true-tip backbones, with family-preserving whole-triple source-controlled null.** This is a descriptive **floristic assembly association**, not demonstrated adaptive flower evolution or pollinator-guild loss. Scientific work resides in **Draft** [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and source-rights-safe [island PR #326](https://github.com/zuizui0223/island/pull/326). Original CI [three true-tree macro data #43461760449](https://github.com/zuizui0223/island/actions/runs/43461760449) **SUCCESS 3/3**; stronger independent family-preserving 99-permutation complete-phenotype null [#43462900657](https://github.com/zuizui0223/island/actions/runs/43462900657) **SUCCESS 3/3**.

## Original cohort and interpretation constraints

- Original source high/medium single-state directly documented accepted species **2,234**, direct flower primary colour×form×symmetry profiles **100**, original GSHHG eligible islands **928**. Each originally represented plant has ONE fixed morphology label globally; within-species island floral evolution is NOT measured.
- Source-pinned independent CoL XR alpha2 exact-name family candidacy, `taxonomy_crosswalk.csv.gz` SHA256 `079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3`, remains **unreviewed automatic taxonomic candidate**. Do not call the family assignments independently sequenced phylogeny.
- Published original Jin & Qian *V.PhyloMaker2* GBOTB `TPL`/`LCVP`/`WP` dated rooted trees (source commit `7af3fb5152f691af2e4ec9d5e2e467d1b50505e9`), with strict original Direct source genus and candidate family agreement: **1,187 / 1,162 / 1,160 original TREE TIP names**. Genus/family inserted phylogenetic tips are excluded, not quietly imputed.
- Each original island and its original five-nearest GIFT geographically plausible mainland **native source** region flora is thinned to exactly **15 different actual published true-tip plant species on each side**, 20 independent random draws. All accepted species names observed in focal island are removed from its GIFT candidate source pool **BEFORE sampling**. Original physical island is ONE statistical row after averaging its source draws; 20 thinning repeats are NOT biological independent replications.
- Original full panel included 928 islands; the true-tip/source-positive selected all-recorded subset is **103 tropical** (42 original 10-degree spatial blocks, 36 independent first-nearest GIFT donor signatures), **741 northern midlatitude** for TPL/LCVP (51 blocks, 64 donors), **739 northern midlatitude** for WP, and northern high/south regions that fail the source-group reportability gate. The strictly native tropical cohort **16 islands/9 blocks/8 donor groups** fails the 20 island, 10 block, 10 source gate; northern-midlatitude strictly native **67** island original subset retains 16 blocks/23 donors.

## Distinct source ecological-beta measures

**Phylogenetic beta (`x`)** = `1 - 2 × (original published shared root-to-tip subtree branch length) / (original published island subtree + original published GIFT donor subtree lengths)` = 1−**PhyloSor**, with actual original dated tree edges. Original branch source is not newly inferred genetic divergence by this study. **Floral functional beta (`y`)** = Jaccard dissimilarity of the ORIGINAL globally observed 100 categorical colour×form×symmetry flower profiles reachable within one original Hamming categorical mismatch from each 15-species flora. THIS is finite categorical support dissimilarity, not original 17D Gaussian KDE density volume. Secondary responses: exact (unsmoothed) profile Jaccard and source-relative rare-minus-common functional reference coverage adjusted for each sample's actual number k of distinct observed phenotypes.

Main slope = change in functional Jaccard per 1 SD of **original region's PhyloSor beta** among independent physical islands, using each whole physical 10-degree block equal weight, with **999 conditional bootstrap block resamples**; separately equal weight of first-nearest original GIFT mainland source signatures and 999 source bootstraps. Model-adjusted slope adds original sampled GENUS and CoL candidate FAMILY Jaccard turnover, original GIFT donor exact-tip pool availability, and original island log area and original log geographic isolation. No confounder is silently dropped when the source design is rank deficient. Source/geography weights are overlapping sensitivities, NOT separate adaptive histories.

## Main observed raw and adjusted phylogenetic–functional beta correspondence

**Tropical: 103 islands / 42 geographic blocks / 36 GIFT candidate first sources.** Whole block-equal slope and 95% conditional interval:

| Original true-tip tree | Raw functional β slope per SD phylogenetic β | Adjusted for genus/family turnover, GIFT source opportunity, island area and isolation |
| --- | --- | --- |
| **TPL** | **+0.00897 [+0.00237,+0.01602]** | +0.00833 [−0.01077,+0.03494] |
| **LCVP** | **+0.01050 [+0.00326,+0.01721]** | +0.00534 [−0.00930,+0.02227] |
| **World Plants** | **+0.00947 [+0.00195,+0.01690]** | +0.00512 [−0.01474,+0.02615] |

All three **RAW** source-signature-equal bootstrap intervals also exceed zero (TPL +0.01063 [+0.00269,+0.01813]; LCVP +0.01189 [+0.00397,+0.01896]; WP +0.01086 [+0.00259,+0.01837]). **Every adjusted interval under both geographic blocks and source signatures INCLUDES zero.** Therefore the observed tropical phylogenetic versus flower-category-space turnover association has NO supported independent remainder after the five taxonomic/environment/source opportunity controls.

**Northern midlatitude:** true-tip source-comparable 741 TPL/LCVP or 739 WP. Whole-block estimates:

| Original tree | Raw slope [95%] | Adjusted slope [95%] |
| --- | --- | --- |
| TPL | +0.00425 [−0.00052,+0.00982] | +0.01664 [−0.00111,+0.02991] |
| LCVP | +0.00172 [−0.00260,+0.00573] | **+0.01603 [+0.00385,+0.02506]** |
| WP | **+0.00799 [+0.00344,+0.01308]** | **+0.02350 [+0.01026,+0.03585]** |

LCVP and WP adjusted source-signature-equal effects also have conditional intervals excluding zero (LCVP +0.01529 [+0.00230,+0.02659], WP +0.01615 [+0.00352,+0.02774]). TPL adjusted source interval **[−0.00403,+0.02559]** crosses zero. This northmid pattern is therefore **tree-source and provenance sensitive**; the full 67-island strictly native subset has adjusted 95% intervals spanning zero in all three published trees. It is NOT a general independent 928-island flower evolution effect.

**Rare-end directional response:** Under each true-tip tree, the **TROPICAL source-relative signed rare-minus-common support gap** does NOT show a phylogenetic beta association whose conditional block/source bootstrap CI excludes zero, even before multivariable controls. Northmid raw negative associations occur for some tree backbones but adj/strict-native/source effects are not reliably interchangeable. Thus documented rare-end footprint differences cannot be causally explained by island phylogenetic turnover from this model.

## Independent extra test: 99 original family-preserving whole-trait shuffles

This is a separate **post-hoc diagnostic**, NOT a preregistered independent lineage-selection test. Using EXACT same true-tip island and donor accepted species memberships and first **5 fixed random thinnings per island** (rather than 20 main-model thinnings), permute the **WHOLE** original 3-trait categorical profile as one indivisible unit among the original accepted source tree tips **within each verified candidate FAMILY**. Keep original family-specific counts of all profiles exactly unchanged, original rooted phylogenetic edges, geography/source, genus/family composition and source species ascertainment fixed. Each original GBOTB tree has **~1,120–1,146** genuinely permutable source tips across **62–63** multiply represented family groups (104 TPL, 104 LCVP, 103 WP original candidate families total). Compute block/source equal raw and multivariate-adjusted original slope for each of **99 independent whole-profile family shuffles**; compare the observed first-five-thinning slope with the expected family-preserving null, never with the 20-draw main slope as if they were identical estimands.

**Tropical original observed vs family-shuffled raw BLOCK slope (first 5 fixed plant draws):**

| Source tree | Observed 5-draw slope | Mean under 99 family-preserving trait shuffles [central 95%] | Empirical upper-tail P for observed excess |
| --- | --- | --- | ---: |
| TPL | +0.01241 | +0.01094 [+0.00050,+0.01985] | 0.40 |
| LCVP | +0.00505 | +0.00962 [+0.00068,+0.02034] | 0.81 |
| World Plants | +0.00108 | +0.00626 [−0.00275,+0.01356] | 0.90 |

All tropical block and source-equal adjustments also fall INSIDE their family-shuffle null distributions. The alternative family-preserving null therefore **can reproduce the tropical positive phylogenetic–functional beta coupling**, without any new island-specific floral evolution. This is consistent with family-level taxonomic composition and original source phenotype distribution sufficiently organizing the observational relationship. It is not a *proof* of family-mediated causality.

**Northern midlatitude** family-shuffle result is mixed: TPL raw observed family-excess p_upper 0.11 (block) /0.12 (GIFT source), LCVP 0.06/0.06, WP 0.01/0.01 (source-specific exploratory). Importantly **all three MULTIVARIABLE adjusted family-excess comparisons have p_upper >0.1** under both block and source weighting: TPL 0.37/0.71, LCVP 0.24/0.36, WP 0.13/0.24. With 99 draws these are low-resolution post-hoc empirical probabilities, not confirmatory multiple-comparison-controlled phylogenetic selection tests. Native midlat n67 has inconsistent observed signs among source backbones and cannot rescue a universal biological syndrome.

## Inferential ceiling and next decisions

**Supported:** in the original documented, source-qualified flora, the degree of true published-tree-branch phylogenetic turnover has a region-sensitive **OBSERVATIONAL** positive link with categorical floral repertoire turnover, especially tropical raw data. This relationship is not independent of family/genus/geography/source opportunity in tropical all-recorded floras. The original source-bound family-preserving phenotype null fails to show extra tropical coupling beyond documented FAMILY-level morphology distributions. Northern-midlatitude selected all-recorded floras may retain residual associations under two backbones, but not across TPL or native-only cohorts and not consistently exceeding a family-trait-preserving null after adjustment.

**Not supported:** actual historical pollinator decreases, pollinator channel switching, flower-level pollination efficacy, selection-driven flower colour/form evolution, species-population adaptation or globally universal shrinking functional hypervolume. Even true original published GBOTB tree-tip species may not be independently newly sequenced; CoL taxonomy is unreviewed candidate; GIFT native mainland source assignment is geographic hypothesis, not observed ancestor; globally rare categorical profiles are not necessarily rare ecological pollination functions. Overlap of geographic/source clusters and multiple exploratory models constrains p-value interpretation.

## Reproducible original-source code

- [PhyloSor branch-beta and category Jaccard research module](../../src/fd/phylo_beta_functional_beta_macro.py), [unit/synthetic tests](../../tests/test_phylo_beta_functional_beta_macro.py), [frozen original study design](../../docs/PHYLO_FUNCTIONAL_BETA_MACRO_SOURCE_CONTRACT.md).
- [Original source-protected 3-tree executable](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_phylo_functional_beta.py), [original source CI workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-phylo-functional-beta.yml).
- Pure original source completions [#43461760449](https://github.com/zuizui0223/island/actions/runs/43461760449) and [#43462900657](https://github.com/zuizui0223/island/actions/runs/43462900657), all three TPL/LCVP/WP jobs successful. Latest anonymous original artifacts: [TPL #14056308924](https://github.com/zuizui0223/island/actions/runs/43462900657/artifacts/14056308924), [LCVP #14057272300](https://github.com/zuizui0223/island/actions/runs/43462900657/artifacts/14057272300), [WP #14057511835](https://github.com/zuizui0223/island/actions/runs/43462900657/artifacts/14057511835). Original private taxon, island and donor source rows plus dated published tree files remain within the authorized CI workspace.

**Status:** 928-island original-source macro phylogenetic-functional beta OBSERVATIONAL ASSEMBLY association — computed, three-tree replicated numerically, but phylogenetically independent causal flower-functional evolution/pollination shortage explicitly NOT established; parent pull requests remain Draft.

## 14 Nov 2026 · CRITICAL SECOND falsification: control the number of distinct source flower profiles

The preceding raw functional beta Jaccard still allowed the exact number of occupied flower profiles to vary between two *equally species-rich* n=15 floras. This reopens the very record-type-richness confound that invalidated the previous 17D Gaussian KDE as a general island functional hypervolume measure. A separately frozen **original EXACT-(k_island,k_mainland) finite-category geometry null** was therefore computed on all actual original source comparisons. Original source CI [#43464021295](https://github.com/zuizui0223/island/actions/runs/43464021295) **SUCCESS, 3/3**; private original source artifactual identities never published.

For every original island and independently geographically assigned GIFT donor species draw, count the number of DISTINCT original colour×form×symmetry phenotypes on each side: `(k_i,k_d)` (each ≤15, species n=15 on each side). Draw **199 random sets of exactly k_i distinct original global phenotype profiles** and independently **199 random sets of exactly k_d global original profiles** from the FROZEN 100 actually observed Direct profile states. Measure each random pair's one-category-mismatch reference-support Jaccard β. The main *excess functional β* is observed β minus mean reference β for that EXACT pair of recorded profile richness numbers. Independently, repeat with profile-center draw probabilities proportional to the number of original global 2,234 Direct plant species carrying each original profile (still sampling unique centers without replacement). Neither reference null preserves genus, family, sampling provider or actual tree topology. Bootstrap the per-island mean of 20 fixed n15 draw residuals; use original 999 whole-block and first-nearest GIFT source-equal replicates. These are POST-HOC profile-count artifact checks, not a new confirmatory ecological causal hypothesis.

### Tropical 103 original islands / 42 original blocks / 36 independent source signature groups

Each column is change in functional-reference β residual per 1 ORIGINAL region SD of phylogenetic PhyloSor β, **island minus conditional expectation given exact k_i and k_d**:

| True published dated backbone | Uniform reference exact-k residual, raw [95% block] | Uniform residual adjusted for genus/family turnover, GIFT source, area, isolation [95% block] | Species-frequency-calibrated exact-k residual raw [95% block] |
| --- | --- | --- | --- |
| **TPL** | +0.00181 **[−0.00394,+0.00872]** | +0.01616 [−0.00362,+0.04172] | +0.00395 [−0.00152,+0.00984] |
| **LCVP** | +0.00150 **[−0.00407,+0.00740]** | +0.00342 [−0.01376,+0.02401] | +0.00386 [−0.00146,+0.00972] |
| **WP** | +0.00068 **[−0.00563,+0.00700]** | +0.00231 [−0.01975,+0.02310] | +0.00314 [−0.00243,+0.00902] |

**ALL three uniform-reference and species-frequency-conditioned functional Jaccard β slopes have conditional block and candidate-source 95% intervals crossing zero**, RAW or adjusted. Hence the prior raw tropical phylogenetic-functional β positive association is **not distinguishable from a geometry reference in which both floras have the same number of distinct recorded phenotype states**. This is *compatible with* unequal original profile richness and taxonomic source composition contributing to the raw correlation. It does NOT prove these are the sole mechanisms or that real pollination functionality is invariant.

### Northern-midlatitude source-qualified true-tip floras

Original TPL/LCVP n=741 physical islands, WP n=739 (51 original physical blocks, 64 candidate source signatures). The same exact island/source distinct-profile-k-conditioned **uniform** β residual produces:

| Tree | Raw conditional block-equal slope [95%] | Genus+family+area+isolation+source adjusted block-equal slope [95%] | Adjusted source-equal slope [95%] |
| --- | --- | --- | --- |
| TPL | **+0.02147 [+0.01665,+0.02517]** | −0.01358 [−0.02643,−0.00408] | −0.01051 [−0.02423,+0.00219] |
| LCVP | **+0.02209 [+0.01617,+0.02666]** | −0.00822 [−0.01959,+0.00224] | −0.00467 [−0.01445,+0.00391] |
| WP | **+0.02586 [+0.02138,+0.02985]** | +0.00241 [−0.00872,+0.01201] | −0.00508 [−0.01615,+0.00678] |

Here the **RAW** source-relative phylogenetic-functional **geometry** β coupling remains strongly positive under both block- and GIFT source-equal estimates (all three trees), even at exact recorded-profile-count k. However the five taxonomy/geography/source controls completely alter its adjusted magnitude, and **there is NO stable positive adjusted effect across tree versions, original block and GIFT donor group bootstraps**. TPL adjusted block alone is negative but its source-equal interval crosses zero; LCVP and WP intervals all cross zero. Repeat the same original profile-k null with global source species-prevalence-weighted profile center choices: raw northern midlat remains positive (+0.01735 TPL, +0.01746 LCVP, +0.02150 WP), but **all adjusted block and source 95% intervals cross zero**.

Original strictly native flora subset in northern midlat n=67/16 physical blocks/23 source groups has nonstable tree-sensitive adjusted patterns; original tropical strictly-native n=16/9 blocks/8 source groups fails original supported-coverage gate. The global original 928 physical islands should NOT be treated as a 928 independent, fully source-tip, strictly native phylogenetic replication.

**Final revised inference ceiling:** The original 100-state **categorical floral support** data do show regional differences in the extent that documented plant lineage replacement is accompanied by category-space replacement. But direct tests of actual source tree phylogenetic beta, full-family trait shuffling, exact floral profile-count-k reference geometry, and geography/source/taxonomy adjustments **do not jointly establish a phylogeny-independent floral functional niche shift**. Tropical raw correspondence is largely consistent with recorded profile richness and FAMILY-level trait composition; northern-midlatitude k-excess raw correspondence is strongly positive but source-genus/family/geographic conditional dependent. This is not a confirmation of real historical island pollinator decline, effective pollination service or within-species island flower phenotype evolution.

Original executed revised source artifact receipts: [TPL #14057508397](https://github.com/zuizui0223/island/actions/runs/43464021295/artifacts/14057508397), [LCVP #14056604565](https://github.com/zuizui0223/island/actions/runs/43464021295/artifacts/14056604565), [WP #14056524796](https://github.com/zuizui0223/island/actions/runs/43464021295/artifacts/14056524796). Original real-data source success **3/3**, source methods and accepted-plant/tip/group identities protected; this experiment stays in **Draft** and explicitly carries multiple-endpoint posthoc limitations.
