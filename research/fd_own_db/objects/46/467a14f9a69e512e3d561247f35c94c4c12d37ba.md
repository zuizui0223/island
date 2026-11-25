# 2026-11-15 — Island-to-island functional similarity of floral assemblages without within-species trait measurements

**Protected original botanical-data execution SUCCESS:** [island original all-unordered-pairs Actions #43667950263](https://github.com/zuizui0223/island/actions/runs/43667950263), anonymous result artifact #14142607215 (archive SHA256: 0bf6042388ddaa186f063ae80e1ae29bd32edc49fbf8ad23b5767f4ae5075b37). Every analysis uses independently accepted *species-wide* categorical floral/reproductive traits. It is NOT a measurement of island-mainland change within the same species, an ancestral-state reconstruction, or evidence of actual effective pollinator identity.

A complementary earlier independent-source implementation [original Actions #43662732106](https://github.com/zuizui0223/island/actions/runs/43662732106) SUCCESS, artifact #14140842992, uses exact-all-pairs Hellinger affinity plus fixed outcome-blind 16 uniform comparators per focal island and original source/genera-preserving label permutations. Same-mask marginal-versus-joint structure is independently calculated in [Actions #43663415375](https://github.com/zuizui0223/island/actions/runs/43663415375), artifact #14141008485. Do not conflate their Hellinger affinity with the current exact pairwise **total-variation overlap**.

## 1. Full original global sampling frame, WITHOUT three-trait completion

- **4,505** original botanical island floras, one original flora missing independently joined GSHHG geography (kept as geography-unresolved).
- **10** individual high/medium directly documented accepted-species category traits. Colour×form, colour×symmetry and form×symmetry are three separate two-trait endpoints requiring exactly two Direct labels on the SAME original accepted plant species; no arbitrary third-trait completion.
- Original all-recorded flora and independently source-labelled native flora stay separate. Introduced and unresolved-status plants are not silently called indigenous. Empty native strata remain missing, NOT zero botanical diversity.
- Original 8,264 independent geographic islands are a broader physical-polygon frame. Only **434** have independent published geological oceanic/continental classification; **425** of these intersect original botanical floras. Neither mainland isolation nor island centroid classification supplies evolutionary adaptation.
- Global source baseline unchanged: all **26** original trait×scope eligibility counts exactly match the previously saved source-Direct original frame. Original 928 three-trait/30-species subset stays a sensitivity, not the main sample.
- Original 10° physical blocks, geographical coordinates, species identities and original GIFT source references remain within original-source CI; public results are anonymous region/endpoints/aggregate only.

For island i, let p_i(c) be the fraction of its DIRECT-recorded species in categorical state c, after requiring >=5 directly recorded species FOR THAT ENDPOINT. **Functional compositional overlap** is the exact total-variation affinity O(i,j)=sum_c min[p_i(c),p_j(c)] =1-0.5 sum_c |p_i(c)-p_j(c)|. It ranges from 0 to1. Taxonomic similarity is independently calculated from the full flora and within the trait-recorded species mask; these estimates are never merged. Category-set Sørensen beta is reported with Simpson turnover and nestedness-resultant components, which DO NOT imply historical loss or ecological replacement.

A separately computed exact E[sum_c(p_i(c)-p_j(c))²] under a globally fixed source species-trait codebook preserves EACH island's original species mask, shared accepted species and observed species count, and alternatively permutes entire observed categorical profiles only within original botanical genus. Neither genus null is an ancestral genetic/null evolution model. Distances/residuals are descriptive and have no independent-island-pair p-values.

## 2. All actual original island pairs (not only 16 random neighbors)

| Categorical trait endpoint, all-recorded flora | Eligible ORIGINAL islands, >=5 Direct species | ALL unordered island pairs | ALL-pair mean TV overlap | Pairs with ZERO shared original full-flora species | TV overlap on zero-shared-species pairs |
| --- | ---: | ---: | ---: | ---: | ---: |
| Flower colour | **3,132** | **4,903,146** | **0.7003** | **2,363,622** | **0.6347** |
| Floral form | **2,806** | **3,935,415** | **0.6092** | **1,601,543** | **0.4962** |
| **Colour × form** (genuine observed combination) | **2,324** | **2,699,326** | **0.3869** | **842,058** | **0.2472** |
| Colour × symmetry | 2,699 | 3,640,951 | 0.6532 | 1,451,027 | 0.5795 |
| Form × symmetry | 2,611 | 3,407,355 | 0.5573 | 1,241,081 | 0.4277 |

**Main interpretation:** some flower functions are repeatedly represented in plant floras with NO species in common, but whole floral **combinations** are substantially more differentiated than broad individual trait categories. This alone is NOT a formal claim that one multivariate endpoint is ecologically more selected than a univariate endpoint: their state-space cardinality, native-status documentation, species masks and expected-overlap ceilings differ.

For the complete all-recorded colour×form pairwise comparison, the original recorded-flora taxonomic Sørensen category set decomposes into mean **turnover=0.3422** and **nestedness-resultant=0.2372**. For the zero-shared-species comparison, **turnover=0.4758**, nestedness-resultant **0.2115**. They describe recorded category sets, NOT pollinator guild loss, plant extirpation or the adaptation of a floral shape.

The original **GENUS-conditional squared-frequency distance** among the 842,058 zero-shared-species colour×form island pairs is expected **0.2078** under source species mask-preserving global genus shuffling, versus observed **0.2297**, descriptive excess **+0.0218**. For all 2,699,326 original pairs, corresponding values are expected **0.1365** versus observed **0.1479**, excess **+0.0114**. This is nonzero genuine category-frequency assembly difference beyond one specific genus-exchange reference, NOT independently proven adaptive divergence or a multiple-testing-supported geographical result. Region, distance, climate and documented effort can still explain the effect.

**Independent original comparison:** the earlier uniform fixed-16-focal-neighbors source audit [#43662732106](https://github.com/zuizui0223/island/actions/runs/43662732106) gives all-recorded colour×form *total-variation* distribution overlap **0.3860** compared with this full 2,699,326-pair value **0.3869**; all 26 original endpoint×status pooled-overlap means differ by at most **0.00291** between exhaustive and uniformly sampled calculations. For zero-species-shared comparisons the uniform-candidate frame is a different sampling distribution, with 0.2526 vs full-pair 0.2472; do not call those independent biological replicates.

## 3. Taxonomic turnover and functional similarity are two distinct outcomes

For colour×form all-recorded global floras, mean full-flora species taxonomic Sørensen–Dice similarity is low, **0.0613**, while the corresponding mean *functional frequency overlap* remains **0.387**; these are different categorical measurements and not proportions of the same total. Restricting to floras with zero accepted species in common retains a positive mean functional overlap 0.247. This makes **shared recorded floral-category repertoires despite different species composition** a genuine quantifiable research subject, even while the observed commonality is not unusually high relative to genus/global distributional null models.

A separate original source audit [#43662732106](https://github.com/zuizui0223/island/actions/runs/43662732106) held species incidence and observed profile frequencies fixed and permuted ORIGINAL species labels globally or within botanical genera. In its pre-sampled all-recorded zero-shared-flora graph for colour×form, observed Hellinger affinity **0.39330** vs global-profile-null **0.49286** and within-genus-null **0.41838**; 199 numerical permutations gave exploratory lower-tail p=0.005/0.020. This is a genus- and source-conditional **sampling-model** diagnostic; it is NOT a causal, independent origin, or pollination-syndrome convergence test. The direction is towards **less, not more, affinity** than those model references.

### 15 November 2026 — independent arithmetic and metric-label audit

A direct inspection of `scripts/fd_island_functional_beta.py` confirms that `incidence_partition` calculates species-set **Sørensen dissimilarity** `1−2a/(n_i+n_j)`, so `taxonomic_similarity_full` is its **Sørensen–Dice similarity** `2a/(n_i+n_j)`, **not Jaccard** `a/(n_i+n_j−a)`. The reported mean **0.0613** must not be relabeled Jaccard. A true mean Jaccard cannot be reconstructed from a mean Sørensen–Dice alone because the per-pair transformation is nonlinear. The primary floral frequency TV overlaps are unaffected by this report correction.

The published rounded all-pair means, zero-full-flora-shared-species means and their exact pair counts independently imply the following overlap for pairs sharing **at least one** accepted plant species: `([all pairs]×[all mean]−[zero-shared pairs]×[zero-shared mean])/([all pairs]−[zero-shared pairs])`. On the original all-recorded floral colour×form endpoint this gives **0.4502** on **1,857,268** shared-species pairs (versus 0.2472 on 842,058 zero-shared pairs). Other derived shared-species category-frequency overlaps are colour ≈0.7614, form ≈0.6867, colour×symmetry ≈0.7020 and form×symmetry ≈0.6315. Values are approximate because source means were rounded to four decimal places. This is a pure consistency audit of published aggregates, **not** a new original-data run, independent biological sample, causal estimate of species sharing, or ecological convergence result.

The genus-preserving expected squared frequency distance is 0.2078 versus observed 0.2297 on the zero-shared colour×form pairs (descriptive excess ≈0.0219 using rounded values; source reports +0.0218 at full precision). Pair reuse, region, spatial blocks, trait-documentation completeness, native status and source opportunity remain uncontrolled in this difference. The same tropical geology subset has **no** continental–continental species-disjoint pairs, so it cannot identify a symmetric species-disjoint geological-origin comparison.

## 4. Oceanic/continental geological context: descriptive only

The main full botanical sample is **not reduced** to the 425 geology-typed islands. Geological origin is a **secondary conditional lens** using the independently typed original botanical cohort. Tropical all-recorded colour×form requires only the two Direct traits:

| Tropical pair class | Eligible botanical islands | ALL original pairwise comparisons | Mean colour×form compositional overlap | Source Direct two-trait photo/documentation coverage averaged across the pair |
| --- | ---: | ---: | ---: | ---: |
| Oceanic–oceanic | 120 | **7,140** | **0.5040** | 0.0842 |
| Continental–continental | 36 | **630** | **0.5852** | 0.0673 |
| Oceanic–continental | 120 oceanic + 36 continental | **4,320** | **0.5243** | 0.0758 |

Oceanic–oceanic pairs are not universally more functionally similar than continental–continental pairs; this broad descriptive result is compatible with heterogeneity in island flora, phylogenetic/genus turnover, study coverage, climate and geographic distance. **Do not test 7,140 against 630 as independent sample sizes**; all pairs reuse islands. The oceanic 120 versus continental 36 islands are different geographical collections and native/introduced source proportions. No independent geographical/archipelago matched selection, standard-error, permutation or causal geology contrast is claimed.

A stricter **no-shared-accepted-species** subset here has **only 71 oceanic–oceanic** and **27 oceanic–continental** pairs, with **zero** continental–continental species-disjoint pairs. Thus the all-species-disjoint origin comparison is structurally **NOT supportable** in this tropical original geology subset. This is not an excuse to remove the whole global sample; it is a support limitation on one geological question.

## 5. Syndrome-relevant flower configurations, WITHOUT assigning effective pollinators

Global composition and regionally stratified original flower colour × floral form are reported from observed individual plants' species-wide directly measured two-trait states; categories such as white+tubular/open_radial are morphological configurations, not actual moth/bee/bird visitors. A combination can be functionally shared across species and genera without implying within-species evolution.

| Observed two-trait configuration, all-recorded source | Tropical (586 islands with >=5 genuinely paired plant traits) | Northern-midlat (1,464 islands) |
| --- | ---: | ---: |
| **white + tubular** | **14.13%** | **3.39%** |
| **white + funnel_trumpet** | **10.78%** | **1.23%** |
| **white + open_radial** | **9.27%** | **25.06%** |
| **yellow_orange + open_radial** | **5.53%** | **11.94%** |

Each cell is the island-equal average fraction of DIRECT-recorded flower-colour+form plant species, *not* all source flora, not abundance, and not the fraction of real pollen delivery. The contrasting descriptive tendencies are qualitatively reflected in independent original-native-labelled subsets (91 tropical vs215 northern-midlat eligible islands): white+tubular 18.24% vs4.29%, white+open_radial 10.40% vs16.65%. **Native-labelled source coverage is sparse**; these numbers do not establish true adaptation or active syndrome evolution.

The independent same-trait-mask source audit [#43663415375](https://github.com/zuizui0223/island/actions/runs/43663415375) reports observed colour×form Hellinger affinity minus the product of marginal same-mask affinities **−0.1629** for its sampled zero-shared-flora pairs. This is a purely *structural independence-reference comparison*, not a randomized null or a formal result about pollinator syndromes. Leading same-species-observed profiles contributing to shared affinities include white/open_radial, white/bell_campanulate and yellow_orange/open_radial. The ecological discussion can compare these configurations to published pollination-syndrome hypotheses but must explicitly abstain from assigning a guild without pollen-effective visitor data.

## 6. Methodological guardrails and next executable inference

- No biological effect from comparing millions of correlated island pairs as if they were independent independent replicates.
- The main botanical inventory includes introduced and unreviewed-status species. Native-labelled outcomes are kept as a separate original provenance-verified sensitivity; they must not be substituted silently for all flora.
- The original Direct state is FIXED for a species across all islands; this measures **floristic functional assembly and composition**, not island–mainland intraspecific change, adaptive plasticity or genetic divergence.
- The regional/coverage distribution itself is heterogeneous. Sampled 30+ species sensitivity uses a different set of islands, and the published original 928 strict triple-case only has an accessory role.
- Current exact global and genus null expectations **do not** control climate, habitat, nearest candidate GIFT-native mainland source, spatial/archipelago dependency, native-status misclassification, collector effort or chronological pollinator replacement. A conditional distance/climate/archipelago comparison is a later, separate test; no inference is fabricated now.
- Recorded profile category Sorensen turnover/nestedness are not actual local species extinctions or ecological transitions.
- Historical GIFT donor regions are not ancestors, recorded GBIF guild presence is not an effective pollinator visit, and non-detection is not historical loss.

Reproducible source: [island protected full original script](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/fd_island_functional_beta.py), [island full test](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/tests/test_fd_island_functional_beta.py), [fd reusable model](../../src/fd/island_functional_beta.py), [fd tests](../../tests/test_island_functional_beta.py), [island original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-island-functional-beta.yml). The experimental work stays in Draft PRs.