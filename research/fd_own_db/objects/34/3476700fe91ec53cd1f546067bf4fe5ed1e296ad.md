# 14 November 2026 — Regional native-provenance sensitivity of 928-island floral categorical space

**ORIGINAL DATA ANALYSIS DONE — two independently seeded original-source GitHub Actions successes; only anonymous regional outputs released.** Source CI [island Actions #43530200679](https://github.com/zuizui0223/island/actions/runs/43530200679) **SUCCESS, both jobs** (seeds 20261114 and 20261115), original aggregate [artifact #14086254224](https://github.com/zuizui0223/island/actions/runs/43530200679/artifacts/14086254224) and [artifact #14087231552](https://github.com/zuizui0223/island/actions/runs/43530200679/artifacts/14087231552). Scientific changes in draft [fd PR #16](https://github.com/zuizui0223/fd/pull/16), original protected source executable in [island PR #326](https://github.com/zuizui0223/island/pull/326).

## Why this sensitivity changes the macroecological conclusion

Previously the original all-recorded 928-island Direct floral dataset revealed a tropical and northern-midlatitude shortfall relative to five geographically plausible GIFT native mainland sources in **globally rare recorded colour×form×symmetry category-space neighborhoods**. However, native-only evidence was restricted to **87 islands (17 tropical)** and source taxonomy, family and published-tree controls did not support a source-independent evolutionary morphology shift.

The new independent goal is to test how the original floristic **origin/provenance filter itself** changes this result, rather than implicitly calling unresolved or introduced plants local natives.

Two frozen ORIGINAL source archives are cross-linked at accepted plant taxon and original island ID:
1. [Original progressive island source run #34232450884](https://github.com/zuizui0223/island/actions/runs/34232450884), artifact **#10058653212**, original Direct species×flower trait and native/introduced/unresolved GIFT flora plus original five-nearest GIFT native mainland candidate sources. Original exact trait CSV SHA256 **9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745**, source original island flora SHA256 **6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a**.
2. Independently archived official WCVP status [original run #32559322028](https://github.com/zuizui0223/island/actions/runs/32559322028), artifact **#9472549694**, original resolved WCVP species native TDWG-L3 region lists and 8,265 GSHHG original island TDWG-L3 mappings. Source archived WCVP status SHA256 **10e7c56f1cc98c1fabc43c24cb8c456476f2ae057dacaa61aae9e635e2420ab0** and original mapping SHA256 **e28351c297ab1ca7416b50b55a99fafd2c5d4ab6c8b984761fe987e0a57953b3**.

Critically, WCVP accepted-species occurrence in a TDWG-L3 native REGION is **not proof that this individual island population is native**. TDWG-L3 often covers a mainland and several islands; native plants may be introduced in individual islands within the same region. Existing GIFT source-native labels are preserved, and GIFT source-introduced taxa are **never promoted**, even with WCVP regional native distribution. Unresolved taxa lacking an unambiguous original TDWG mapping or available original WCVP native regions are retained as UNKNOWN, never nonnative or extinct by default.

The archived WCVP native-range lists originate from a prior source-status reconciliation containing **55,137** unique accepted plant species, of which only **1,379** of the **2,135** distinct phenotype-known species in the 928 islands have available native regional range codes. The full original official WCVP master may cover more species; this is a **source-snapshot restricted, nonrandom and potentially taxonomically biased subset**. Original TDWG assignment exists for 8,265 islands but was accepted for 3,451; missing/ambiguous original island geography prevents WCVP promotion.

## Original native eligibility improvement

Frozen source high/medium **2,234 original global directly documented flower species**, 100 truly recorded three-categorical-trait profiles. Original eligibility: **928 islands**, 81,635 independent original accepted island×plant Direct flower triple records across 2,135 species. Their original GIFT provenance is **7,536 GIFT-native**, **2,354 introduced**, **71,745 unresolved**. Original archived WCVP native TDWG region compatibility applies to **23,837** of the unresolved rows, without changing any of the 2,354 GIFT-introduced rows.

Among the original 928 islands, require at least **30 species** known by high/medium Direct (color/form/symmetry) record within the specified provenance tier:

| Regional group | All observed source floral islands | Original GIFT strict-native island floras >=30 | GIFT-native + WCVP regionally native-compatible island floras >=30 | Newly eligible due to independent WCVP regional compatibility |
| --- | ---: | ---: | ---: | ---: |
| **Northern midlatitude** | 743 | **67** | **287** | **+220** |
| **Tropical** | 106 | **17** | **36** | **+19** |
| Northern high latitude | 53 | 0 | 30 | +30 |
| Southern extratropical | 26 | 3 | 4 | +1 |
| **All** | **928** | **87** | **357** | **+270** |

Requiring >=30 source-known Direct plant species in the independent GIFT native mainland source donor pool after excluding ALL originally recorded focal island plant names retains **all 36 tropical regionally compatible islands / 20 original physical 10-degree blocks / 16 first-nearest GIFT source signatures**, **285/287 northern midlat / 33 blocks / 40 first sources**, 30 northern high islands but only four blocks and TWO source signatures (NOT independently reportable), four southern extratropical islands (NOT reportable).

## Fix the sample-number artifact analytically, not by turning another KDE bandwidth

The earlier original 17D Gaussian KDE hypervolume was dominated by the number of observed profile states. The latest finite categorical metric is explicitly different: for 100 ACTUALLY observed global Direct phenotypes r and an observed island/source pool of N distinct plant species, the expected one-trait-radius reference coverage under a uniform **30-plant without-replacement** draw is

\[
E[V_{\le1}(r)] = 1 - \frac{\binom{N-m_r}{30}}{\binom{N}{30}},
\]
where m_r is the number of actual plants in that source pool whose original fixed flower profile is within ONE mismatch of reference category r. Summing over the **29 original rare profiles** (source-carried by 1–2 global Direct species each) and **25 original common profiles** (source-carried by >=26 species each) gives **EXACT finite-population expected category support**, with no artificial random species-sampling uncertainty. As in all previous source analyses, the focal island species' accepted names are removed from candidate native GIFT source pools, and both sides are standardized to 30 species.

The second reference-geometry control subtracts the expected rare-minus-common category support given the *recorded number k of distinct original floral profiles* under 1,999 fixed global-uniform original profile-count reference draws. The distribution of k under a true 30-species finite-population draw is evaluated with **2,000 multivariate-hypergeometric draws per source pool**; this is a documented MC approximation, not analytically exact. Two completely independent numerical RNG seeds 20261114/20261115 reproduce original published source aggregate effect estimates to ~0.0001. **1,999 source-cluster/physical-block resamples** provide conditional original region estimates, not independent archipelago historical replicates.

## Original effect results after controlling 30 original plant species and profile number k

Below, positive = focal island flora has greater rare-minus-common one-trait categorical neighborhood support than its candidate GIFT native mainland flora, *relative to expected geometry given exactly k original profile types*. Values are original physical spatial-block-equal conditional mean and 95% bootstrap CI, based on original seed 20261114; independent 20261115 confirmation is near identical.

| Original region and plant provenance stratum | Source-comparable original islands | Original 10° spatial blocks / first GIFT source signatures | Island−source k-conditioned rare-minus-common support [conditional 95% block interval] |
| --- | ---: | ---: | --- |
| **Tropical, all recorded** | **106** | 42 / 36 | **−0.03109 [−0.04705,−0.01369]** |
| **Tropical, GIFT native + WCVP regional-compatible** | **36** | 20 / 16 | **+0.00230 [−0.02209,+0.02650]** |
| Tropical, strict GIFT-native | 17 | 10 / **9** | −0.00439 [−0.03927,+0.02974] |
| **Northern midlat, all recorded** | **736** | 51 / 64 | **−0.02627 [−0.04468,−0.00892]** |
| **Northern midlat, native+WCVP regional compatible** | **285** | 33 / 40 | **+0.01851 [−0.00320,+0.03918]** |
| Northern midlat, strict GIFT-native | 67 | 16 / 23 | +0.01086 [−0.01886,+0.03971] |

**Neither tropical nor northern-midlat WCVP regionally native-compatible subset has a block-equal interval below zero.** Independent globally fixed random seed 20261115: tropical all −0.03115 [−0.04687,−0.01500], regional native compatible +0.00243 [−0.02185,+0.02739]; northern midlat all −0.02628 [−0.04441,−0.00794], regional compatible +0.01855 [−0.00285,+0.03889]. Original all-recorded negative source contrast replicates. The regional-compatible estimate is a different, selection-restricted botanical assemblage, not proof the omitted taxa are all introduced or the native island flora has no deficit.

### Crucial exactly SAME ISLAND robustness

To distinguish **cohort composition** from merely comparing different geographic islands, directly compare the all-recorded versus WCVP regionally native-compatible flora on the SAME original physical island IDs. The original five GIFT source regions and the source flora remain identical, so the donor source-relative term cancels; only original focal botanical provenance strata change.

| Exact same original islands, native-compatible minus all-recorded | Number of original islands / blocks / GIFT source signatures | k-adjusted rare-minus-common change, **BLOCK equal** [95%] | Same comparison, **GIFT SOURCE equal** [95%] |
| --- | --- | --- | --- |
| **Northern midlat** | **287 / 33 / 40** (285 have valid GIFT donor) | **+0.02256 [+0.01186,+0.03462]** | +0.01971 [+0.01173,+0.02814] |
| **Tropical** | **36 / 20 / 16** | +0.01603 [−0.00013,+0.03547] | +0.01438 [−0.00224,+0.03881] |

For a second original independent numerical seed, midlat paired +0.02251 [0.01148,0.03397]; tropical paired +0.01612 [0.00081,0.03569]. Thus the northern-midlatitude **positive within-exact-island provenance sensitivity** is numerically and spatial/source-conditionally supported. The tropical same-island paired difference is positive in point estimate but **CI lower bound is near zero and crosses it under some valid resamples/seeds**: do NOT call the tropical paired provenance *difference* independently established.

A source-known strictly native-only cohort still has 17 tropical island floras, just nine distinct first GIFT source signatures, and lacks enough independent taxon provenance for a decisive native evolution inference. Higher northern-midlatitude sample support improves range of ecological sensitivity, not the accuracy of assigning every WCVP candidate to an individual island.

## Interpretation, biological claim ceiling and next action

**New positive, source-verified macro result:** botanical **provenance filtering materially alters observed rare-end floral morphology category support**, especially in the northern midlat. The original all-recorded floral category deficit relative to candidate mainland flora is NOT a robust negative universal effect after selecting original directly native-verified or source WCVP regional-native-compatible candidates. This is stronger evidence for the importance of **flora composition / taxonomic-source availability and regional native status representation** than for actual repeat adaptive floral simplification caused by lost island pollinators.

**NOT supported:** actual native establishment of WCVP region-compatible taxa on a specific island; historical pollinator decline, actual pollinator guild shift, measured pollination effectiveness, native-only floral trait evolution, independent floral functional niche collapse, independent island radiations, or causal evolutionary filtering. The same globally fixed accepted species has the same Direct phenotype across island and source. The three GBOTB source backbones and earlier GloBI/GloPL evidence do not repair these provenance limitations. GBIF Bombus environmental suitability in the existing island covariate file is an **environmental niche-compatibility estimate**, not real observed pollinator abundance or temporally measured decline, so it is not treated as an ecological outcome here.

**Next falsification:** if independent original taxon-local native statuses can be assigned at the actual individual island (preferably vetted indigenous checklists or repeated, field-verified native populations), re-run the source-normalized finite categorical reference on that heldout cohort. Maintain source-group and physical block uncertainty. A separate **real local guild-composition exposure** is necessary before attributing the functional floral repertoire difference to pollinator changes.

## Source and reproducibility

- [Exact finite analytical reference / regional native source method](../../src/fd/regional_native_provenance_volume.py)
- [Source-independent WCVP/TDWG + fixed-30 finite-support tests](../../tests/test_regional_native_provenance_volume.py)
- [Original SHA-verified island runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_wcvp_regional_native_macro.py)
- [Original two-seed source CI](https://github.com/zuizui0223/island/actions/runs/43530200679), original source/artifact ID and input SHA locked; official WCVP artifact from August 2026 is unchanged; anonymous original [seed 20261114 #14086254224](https://github.com/zuizui0223/island/actions/runs/43530200679/artifacts/14086254224), [seed 20261115 #14087231552](https://github.com/zuizui0223/island/actions/runs/43530200679/artifacts/14087231552). Both CI matrix jobs SUCCESS. Original botany, source species, source island IDs and exact island georeferences are not distributed.
