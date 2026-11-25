# Geographically disjoint origin-matched island floral repertoire test

**9 November 2026 — `fd` independent exploratory follow-up.** Public frozen source from [island PR #321](https://github.com/zuizui0223/island/pull/321). [Original PR4 matched-repertoire report](../20261108_matched_flower_repertoire/README.md), [PR9 spatial-independence correction](../20261109_coverage_rarefaction/README.md), [matching protocol](../../docs/GEOGRAPHIC_INDEPENDENCE.md).

## 1. Question and primary result

The prior **37 oceanic–continental pairs** are not 37 independent geographical contexts; they share only five transitive spatial-block and source-archipelago components. What happens if one **rematches the whole independently geological-typed island population** while **forbidding two different pairs from sharing a spatial block and/or an archipelago**?

The result is a strong negative *inference limit*: a significant **geology-specific floral repertoire shift is not established in these geographically disjoint samples**, even though the original 37-pair observational category composition difference remains visible.

## 2. Source and outcome-blind matching

Start with both the All and Direct complete-floral support populations, independently audited geology labels from Roeble/Zizka, and restrict to the original northern-midlatitude and tropical regimes. Use the EXACT original PR4 matching calipers and standardized cost from area, mainland distance, complete name richness, joint colour/form trait coverage and climate PC1–4. **No flower colour/form outcomes enter candidate selection, assignment or island-group restrictions.**

The pinned geography yields **219 candidate continental–oceanic edges** satisfying the original criteria. A global 0/1 mixed-integer optimization maximizes the number of selected pairs, while penalizing smaller geographic covariate differences. Four **post-hoc** conditions:

| Matching constraint | Eligible edges | Nonoverlapping pairs | North-midlatitude / tropical | Independent-source evidence quality |
|---|---:|---:|---|---|
| No common `spatial_block` among different pairs | 219 | 13 | 7 / 6 | Block-disjoint, but may share Weigelt archipelago |
| No common labelled Weigelt archipelago | 219 | 11 | 6 / 5 | Archipelago-disjoint, but may share geographic block |
| Neither geographic block nor known archipelago reused | 219 | 11 | 6 / 5 | Two endpoint source-archipelago labels are *unknown* |
| **Neither reused; both endpoints have verified archipelago labels** | **210** | **10** | **5 / 5** | **Both constraints satisfied and no missing source-archipelago labels** |

No physical island is reused in a group. The **10-pair verified** matching is the strictest geographic sensitivity but does not establish phylogenetically independent ancestry, true mainland source identity, or geological treatment randomization. The cohorts differ from the old fixed 37 matches and are all **post-hoc re-matched**.

## 3. Exact pair signs and complete multiplicity family

For each matched sample, calculate the squared difference of mean normalized *distinct* flower colour × form profile incidence over the **same 496 recorded state-set profiles**, and separately the full original colour/form set-Jaccard similarity-kernel contrast. All/Direct are **nested** evidence scopes, not independent phylogenetic replicates.

Enumerate **all 2^n** independent matched-pair sign assignments (8,192 for 13 pairs, 2,048 for 11, 1,024 for strict 10). Each squared statistic is invariant to overall sign reversal, so the minimum attainable exact p is respectively 0.000244, 0.000977 and 0.001953. One fixed **16-comparison BH/BY family** (four matching choices × two nested evidence scopes × two multivariate endpoints).

| Independence restriction | Pairs | Exact-profile L2 p All | Exact-profile L2 p Direct | Full categorical-kernel p All / Direct |
|---|---:|---:|---:|---:|
| Disjoint spatial block | 13 | 0.1685 | 0.0408 | 0.1011 / 0.1191 |
| Disjoint source archipelago | 11 | 0.0566 | 0.0684 | 0.0635 / 0.1592 |
| Disjoint on both, unknown archipelagos possible | 11 | 0.2549 | 0.1777 | 0.2402 / 0.3242 |
| **Disjoint on both, all source archipelagos verified** | **10** | **0.1680** | **0.2188** | **0.2070 / 0.3320** |

**Zero of the 16 tests has BH or BY q below 0.05**. The smallest BH q is **0.2734**, attained by multiple overlapping choices, including the raw p=0.0408 in the 13-pair Direct spatial-block-only case. No single raw p in a cherry-picked geography condition should be treated as a discovery.

The results do **not** demonstrate biological equivalence or absence of a geological effect; even the strictest test has only **10 independent matched pairs** (five per broad region), its own detection power and sample-selection boundaries, and imperfect source/ancestry independence.

## 4. Interpretation of the complete fd evidence chain

- **Observed:** The 37 geography-matched oceanic/continental floras show different *exact recorded* colour × form combinations. That association is common under random equal-profile and equal-complete-name thinning, and in a 21-pair gross-GBIF-effort-balanced subset.
- **Source/genus conditional:** When actual taxon→island incidence and genus × four-regime availability are fixed while shuffling whole floral profiles, **no global origin excess** passes multiplicity controls; lineage composition remains a viable explanation.
- **Spatial dependence:** The original 37 pairwise sign flips violate a prudent geographic-dependence check; only five connected components are available for conservative whole-group permutations. Those exact tests did **not** reject the category origin null.
- **Geographically disjoint alternative:** Searching the *full* geological population for independent matched comparisons yields only 10 fully source-labelled/region/block-disjoint pairs, and **no multiplicity-supported repertoire difference** across the entire 16-test family.
- **Functional geometry:** Broad original categorical-trait similarity kernel separation is not confirmed, and the projected 17-dimensional Gaussian KDE volume difference loses clear separation after controlling recorded number of category profiles.

Therefore the paper cannot presently claim universal **oceanic-island floral simplification**, geological selection, source-independent trait filtering or pollinator niche substitution. A more defensible proposition is:

> *Recorded island floral phenotypes show geographically contingent composition differences; genus and probable regional source assemblage provide useful explanatory/predictive structure, but independently established geological origin effects on floral functional dispersion remain unresolved.*

That is an observational, lineage-assembly/predictive framing. More independent source archipelagos or a tractable hierarchical spatial covariance model—and independent trait/detection histories—are needed for conclusive geological effect estimation.

## 5. Reproducibility and source rights

- Implementation: [global MILP matcher and exact sign-swap test](../../src/fd/geography_independent.py).
- [Frozen study protocol](../../docs/GEOGRAPHIC_INDEPENDENCE.md) and [source-code replay #42576004635](https://github.com/zuizui0223/fd/actions/runs/42576004635) including whole upstream G0, origin, G1, 2D/17D KDE, rarefaction, spatial-cluster and source-geographic test steps.
- **60 synthetic/integration tests passed; all four public workflows successful** in [the final source-code CI #42576004635](https://github.com/zuizui0223/fd/actions/runs/42576004635).
- Source output: aggregate matching support, full 16 exact-comparison null statistics with unadjusted p/BH/BY, and a source SHA manifest in GitHub Actions `fd-origin-strata-analysis` artifact. Results are also summarized in [support table](independent_matching_support.csv) and [16-test table](independent_origin_repertoire_tests.csv).
- Neither committed source nor uploaded workflow artifacts contain exact island IDs, taxa names, raw source-archipelago labels, source entity IDs, species-to-profile maps, or restricted original flora checklists. Uncertain archipelago labels are never silently treated as verified independent evidence.
