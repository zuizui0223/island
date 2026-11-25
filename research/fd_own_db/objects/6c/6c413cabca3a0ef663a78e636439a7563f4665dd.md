# 15 November 2026 — Island floral composition similarity versus taxonomic turnover

## Question and source execution

The user explicitly identified the measurement limit: traits were not measured separately on island and mainland populations of the SAME species. This study therefore concerns **floral composition of different island floras**, not intraspecific trait shifts, adaptation or historical pollinator replacement.

Two protected original-data executions completed:

- Main beta analysis: [island Actions 43662732106](https://github.com/zuizui0223/island/actions/runs/43662732106), code pinned to fd `a678df2c70daddcdfa7695351a698b1edd287a86`, artifact **14140842992**, ZIP SHA256 `f7d082215dab75635d353ffd123e01bacf507dbe52f7edd6d4228b672ebd720c`.
- Same-mask marginal/joint diagnostic: [island Actions 43663415375](https://github.com/zuizui0223/island/actions/runs/43663415375), code pinned to fd `9f33dcdc057ff91b1112776ce6aa5f516890b125`, artifact **14141008485**, ZIP SHA256 `3c0b2adfa7ece19cff1dd1e464197aa1ce20f8c93580121477a1022e604a93c7`.

All **4505 original GIFT floras** remain in coverage denominators. Ten independent Direct trait masks plus three actual two-trait masks are used separately in all-recorded and original-source-native flora. No photograph, common-genus, shared-phylogenetic-tip or third-trait eligibility criterion was added. All **26** original endpoint/scope eligible counts and raw finite-five Rao means reconcile with the independently archived previous analysis; maximum mean difference **1.11e-16**.

The main run passed 12 targeted tests and the full fd suite **228 tests**. The same-mask follow-up passed 11 targeted tests and the expanded full suite **232 tests**; four existing numerical/deprecation warnings were recorded, not suppressed. Local independent checks of 200 random compositional matrices matched brute-force all-pair affinities to **4.44e-16**. Actual category-beta additive identities agree to **2.22e-16**.

## Metrics and population

For an island, p is the proportion of observed species in each recorded trait category/profile, with species equally weighted. The composition affinity is **B(i,j)=sum_c sqrt(p_ic p_jc)**, range 0–1, equivalent to 1 minus normalized squared Hellinger distance. It is NOT percent species shared or percent matching flowers. The exact mean B over ALL unordered eligible pairs is computed algebraically.

For taxonomic overlap and conditional pair summaries, every eligible island is a focal and receives **16 uniformly selected other eligible islands**, selected deterministically before inspecting category values. Directed comparisons are averaged within each focal island, then equally across focal islands. Three geographical-nearest comparators form a separate descriptive graph. Complete pair distributions beyond mean B are thus sampled, not claimed exhaustive. A pair can share an endpoint with other pairs; **pair counts are not independent biological replicates**.

Taxonomic Jaccard is measured both on the full recorded flora and on the exact trait-known species mask. A secondary contrast retains selected pairs with **zero shared recorded accepted-species labels in the full flora of that scope**. This is not a requirement for the main population and is not proof of zero biological overlap under unreviewed synonym concepts.

Category/profile-set Sorensen beta is also split into replacement and nestedness-resultant components. These are comparisons of observed category sets, not reconstructed loss/extinction. This is not KDE hypervolume overlap or a continuous functional niche-volume estimate.

## Global observations

All-recorded flora; each row has its own observed species mask. Cross-row differences are therefore not controlled tests of which trait is more conserved.

| Endpoint | Eligible islands | Exact all-pair mean B | Mean B on selected pairs with no shared recorded species | Geographical-nearest-three mean B |
| --- | ---: | ---: | ---: | ---: |
| Flower colour | 3132 | 0.886732 | 0.849031 | 0.936792 |
| Floral form | 2806 | 0.805056 | 0.723757 | 0.899362 |
| Floral symmetry | 3145 | 0.973045 | 0.963927 | 0.982142 |
| Colour × form | 2324 | 0.539566 | 0.393297 | 0.750020 |

The uniform graph's mean full-flora taxonomic Jaccard is **0.02174** for the colour endpoint and **0.03596** for colour × form. These are a different metric from B and must not be compared as an effect-size ratio.

The no-shared-recorded-species colour comparison has **24127 distinct unordered sampled pairs**, 3108 focal islands and 3129 participating islands; colour × form has **11759 distinct unordered sampled pairs**, 2211 focal islands and 2277 participating islands. Every 3132/2324 eligible island remains a focal in the unrestricted uniform graph. The narrower disjoint contrast does not replace the main study.

The recorded native-only colour × form sample has **327** eligible islands. Its selected no-shared-native-species pairs have B **0.392965**, broadly similar as a descriptive value to all-recorded 0.393297. These are different flora/status samples, not a test of introduced-species effects.

## Similarity does NOT imply greater-than-reference convergence

For three fixed primary endpoints (colour, form, colour × form), the SAME disjoint-flora graph was evaluated under **199** global species-label permutations and **199** within-genus permutations. Whole two-trait profiles are moved intact, not synthesized from marginal categories. Each species gets ONE label per replicate across all islands; species incidence, sampling masks, global profile frequencies and pair graph remain fixed. Within-genus permutations also preserve each genus's profile inventory.

| Endpoint | Observed disjoint-flora B | Global-permutation mean | Within-genus-permutation mean | Observed minus within-genus mean |
| --- | ---: | ---: | ---: | ---: |
| Colour | 0.849031 | 0.908271 | 0.871196 | −0.022165 |
| Form | 0.723757 | 0.825280 | 0.746218 | −0.022462 |
| Colour × form | 0.393297 | 0.492857 | 0.418376 | −0.025079 |

The predefined **excess-similarity** direction is not supported: all six exploratory BH q-values are **1**. The lower-tail diagnostics are 0.005 for the first five contrasts and 0.020 for within-genus colour × form; these lower-tail outcomes are exploratory diagnostics, not the original upper-tail family or evidence of adaptive divergence. The 199-draw Monte Carlo resolution is 0.005. The reference is the global documented species/profile pool, not ancestral populations, a climate-matched dispersal source, a complete phylogeny or a model of missing traits.

Hence substantial absolute resemblance of coarse floral categories is real in these records, but does not establish universal island convergence above the compositional/taxonomic reference. The observed differences are consistent with geographically structured floral assemblages, with observation coverage and biogeographic history still unresolved explanations.

## Same species masks: marginal versus joint composition

A follow-up recomputed colour-only and form-only proportions using EXACTLY the same two-trait-recorded species, islands and sampled pairs as the joint endpoint. This prevents an apparent joint-vs-marginal contrast caused merely by different eligibility masks.

The product reference is computed per pair as B_colour(i,j) × B_form(i,j), then focal-weighted. It is the analytical affinity of independent marginal-product compositions. It is NOT the product of already averaged affinities, an observed flower-profile reconstruction, a finite-sample permutation null, or a significance test. Its gap is signed and can be positive or negative.

| Colour × form cohort | Colour-only B, same mask | Form-only B, same mask | Independent marginal-product reference | Actual joint B | Joint minus product reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| All recorded, uniform graph | 0.896566 | 0.756597 | 0.684749 | 0.538444 | −0.146306 |
| All recorded, no shared flora | 0.859788 | 0.643459 | 0.556222 | 0.393297 | −0.162925 |
| Original native, uniform graph | 0.879351 | 0.755917 | 0.669339 | 0.490107 | −0.179233 |
| Original native, no shared flora | 0.849004 | 0.688980 | 0.587407 | 0.392965 | −0.194442 |

This documents that knowing the separate marginal category distributions does not recover the observed inter-island joint-profile similarity. Sparse category sampling, the number of possible profiles and actual trait association all contribute; the gap alone does not isolate a mechanism or prove a significant extra biological effect. All six endpoint/scope joint affinities and their disjoint-flora masks exactly reproduce the first original run.

## Which floral combinations contribute to the shared component?

For all-recorded disjoint-flora colour × form pairs, the largest contributions to B=0.393297 are:

| Actual observed profile | Contribution to affinity |
| --- | ---: |
| White × open/radial form | 0.132586 |
| White × bell/campanulate form | 0.056098 |
| Yellow/orange × open/radial form | 0.054802 |
| White × tubular form | 0.030167 |
| Yellow/orange × tubular form | 0.024773 |

The first three sum to **61.91% of the measured shared-affinity statistic**, NOT 61.91% of species, flowers, islands, visits or pollination services. The category contribution plus the unlisted remainder exactly reconstructs B. Each released profile occurs on >=10 participating observed islands; no species/island labels are exported.

These provide concrete *pollination-syndrome-related trait combinations* to discuss. They do not justify assigning the islands to bee-, bird- or moth-pollination systems. Effective visitor identities, rewards, phenology, pollen deposition and seed set were not measured here.

## Oceanic versus continental: broad regional descriptions, not treatment effects

Within the tropical independently typed subset, all unordered eligible island pairs were enumerated and summarized symmetrically with focal-island weights.

| Pair type, tropical colour × form | All-recorded B (120 oceanic /36 continental islands) | Original-native B (63 oceanic /26 continental islands) |
| --- | ---: | ---: |
| Oceanic–oceanic | 0.662120 | 0.604754 |
| Continental–continental | 0.755612 | 0.678325 |
| Oceanic–continental | 0.693583 | 0.613974 |

All-recorded unordered counts are 7140 OO, 630 CC and 4320 OC. Native counts are 1953, 325 and 1638. These pair numbers are not independent replicates. The observed pattern is not that oceanic floras are uniformly more similar to one another; continental-island pairs are more similar descriptively in this tropical sample. Spatial distance, climate, source region, documented species richness and trait coverage are not balanced here, so a geological-origin causal effect is NOT claimed.

## Paper framing

A supported descriptive framing is **recurrent floral trait repertoires combined with geographical differentiation in their frequencies and combinations across taxonomically different island floras**. The appropriate title-level question is taxonomic turnover versus similarity/differentiation of pollination-related floral composition. The present data cannot establish island-induced within-species size/colour evolution, local adaptation, historical loss, or a realized pollinator shift.

Code: [island_floral_beta.py](../../src/fd/island_floral_beta.py), [same-mask pair structure](../../src/fd/island_floral_pair_structure.py). All individual taxa, original island IDs and pair keys remain inside protected source computation. The original all-flora population and endpoint-specific observed masks are unchanged.

Methodological references: Legendre & Gallagher (2001), doi:10.1007/s004420100716; Baselga (2010), doi:10.1111/j.1466-8238.2009.00490.x; Ollerton et al. (2009), doi:10.1093/aob/mcp031. The finite-category implementations here were independently checked against direct calculations and boundary cases.
