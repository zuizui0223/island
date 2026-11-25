# Geological origin and floral-profile repertoire: matched global follow-up

> **2026-11-09 spatial-dependence reappraisal (supersedes inferential significance below):** The original 37 matched pairs share broad geographic blocks, forming only **five transitive spatial components (23/11/1/1/1)**. Exact whole-component label sign flips produce **p=0.2500 (All exact L2)** and **p=0.3750 (Direct)**, NOT the individual-pair p=0.0001 quoted below. A separate Weigelt archipelago grouping also leaves only five connected components. The older pairwise BY q≈0.00042 therefore does **NOT** establish a spatially robust origin-specific floral trait repertoire difference. The observed descriptive composition contrast remains real, and sample-support rarefaction retains differences, but stronger geographic and source/genus controls are needed. [Updated audit](../20261109_coverage_rarefaction/README.md).

**8 November 2026 — independent exploratory `fd` analysis.** Main paper line remains the origin-aware floral phenotype assemblage study. This is **not** a claim of geological origin causing evolution or of proven pollination-functional equivalence.

**Verified execution:** [fd Actions #42445694772](https://github.com/zuizui0223/fd/actions/runs/42445694772), source-pinned public aggregates from [island PR #321](https://github.com/zuizui0223/island/pull/321). [Matching protocol](../../docs/MATCHED_REPERTOIRE_PROTOCOL.md). Geological origin: Roeble et al. (2024) and Zizka et al. (2022) crosswalk. No source-protected taxon×trait observations or geological crosswalk rows were committed.

## 1. The biological question and the novelty boundary

Prior published work has studied island functional trait spaces, source-pool filtering and geological differences in species richness. The separable question here is whether **continental and oceanic island floras have different recorded floral colour×form repertoires even when individual island richness and broad geographic/observation conditions are comparable**, and whether those repertoire changes are large in the *actual categorical distance space* or mainly exchange equivalent-looking categorical labels.

This is a new, inspectable analysis target; it is **not** evidence that we discovered the first island floral trait-space difference. The current measures represent recorded colour/form category architecture, not visitor interaction or performance.

## 2. Cohort

The independent geology crosswalk classifies 434/8,264 geographic islands; the prior All complete-floral threshold admits 260 (130 oceanic, 130 continental), while Direct admits 243. For the main test, require **both All and Direct eligible and complete covariates** (241 islands). Then match **37 unique continental–oceanic island pairs**, including 23 northern-midlatitude pairs and 14 tropical pairs, each selected with no floral-profile outcome access.

The paired selection controls regional class, island area, continent distance, complete name richness, joint floral trait coverage, and climate PCs 1–4 using calipers and one-to-one minimum-cost maximum-cardinality assignment. But balance is imperfect: maximum absolute standardized covariate mean differences are approximately **0.571** (north-mid climate PC3) and **0.569** (tropical coverage). These are an observational matched sample, **not** an experiment or perfectly matched historical source pool.

All and Direct are evaluated on the **same 37 pairs** and are nested evidence scopes, **not independent biological replications**. Global larger descriptive 260/243-island tables and the matched subset must not be conflated.

## 3. Primary multivariate repertoire tests

Profiles are the **exact recorded sets of colours × exact recorded sets of flower forms**; multistate records are never expanded into a hypothetical Cartesian product of flower morphs. Each island is normalized by its count of distinct observed profiles. A fixed reference of 496 All profiles defines the family and the similarity kernel.

Each comparison uses 9,999 within-pair origin-label sign swaps and BH/BY correction across the fixed four global tests (two endpoints × All/Direct).

| Endpoint | All observed statistic | All p / BH q / BY q | Direct p / BH q / BY q | Interpretation |
|---|---:|---|---|---|
| Exact profile repertoire difference (L2) | 0.0007564 | **0.0001 / 0.0002 / 0.00042** | **0.0001 / 0.0002 / 0.00042** | Strong difference in which exact recorded categories represent island floras |
| Original categorical-similarity-kernel difference | 0.0006964 | 0.0828 / 0.1104 / 0.2300 | 0.1609 / 0.1609 / 0.3352 | **Not supported after adjustment**; no demonstrated broad functional-distance separation |

The multivariate test is about **recorded floral profile composition**, not total Gaussian hypervolume. Jaccard profile similarity is based on original colour/form state-set distances, NOT an impoverished two-PCoA-axis projection.

In the same matched cohort, the observed *mean* oceanic–continental equal-profile Rao's Q difference is **−0.00822** (All; paired-sign-flip raw p≈0.0018) and **−0.00703** (Direct; p≈0.0090). These raw secondary tests are exploratory and do not replace prior group-specific fixed-margin nulls, which found no supported *excess* within-island packing beyond repertoire frequencies.

## 4. Which exact profiles account for the descriptive contrast?

None of the **496 individual exact-profile** tests (All or Direct) remains supported after BH/BY adjustment; minimum All BH q≈0.257 and Direct BH q≈0.369. The signal is **distributed across profiles** rather than attributable to a single independently confirmed colour×form combination.

For the matched 37 pairs, example prevalence differences (NOT FDR-supported as individual profiles):

| Exact recorded state-set profile | Continental | Oceanic |
|---|---:|---:|
| blue_purple × bell_campanulate | 73.0% | 32.4% |
| yellow_orange × bilabiate | 45.9% | 8.1% |
| blue_purple\|white × bilabiate | 40.5% | 67.6% |

The full eligible All dataset (130/130 geological classes, unmatched) likewise reports yellow_orange × bilabiate at **53.8% continental versus 12.3% oceanic**, and blue_purple × composite_head at **67.7% versus 29.2%**; these are descriptive and not source-conditioned effects.

Importantly, simple presence of each broad colour state (e.g. white, red_pink, blue_purple) is saturated across all 37 matched islands. None of the simple component colour/form states has BH-adjusted paired support. **Combinations and how categorical profiles are distributed**, not binary absence of entire broad flower-colour groups, define the current distinction.

## 5. Residual-geography and climate-balance sensitivity

In addition to the frozen initial 37-pair comparison, we re-matched under post-hoc stricter overlap rules. These sensitivity analyses are **not independent confirmatory discoveries**, do not receive new primary FDR claims and cannot reconstruct a historical source pool.

| Sensitivity | Pairs | Exact-profile L2 p All / Direct | Categorical-kernel p All / Direct | Rao Q mean contrast All / Direct |
|---|---:|---|---|---|
| Stricter coverage and climate-PC3 balance | 24 | 0.0282 / 0.0402 | 0.1260 / 0.3171 | −0.00744 / −0.00655 |
| Island centroids ≤5,000 km apart | 25 | **0.0001 / 0.0045** | 0.0297 / 0.1839 | −0.01327 / −0.01196 |
| Island centroids ≤3,000 km apart (north-midlat only) | 14 | 0.0062 / 0.0078 | 0.0338 / 0.0869 | −0.01187 / −0.01120 |

The original tropical subset included some geographically distant islands from distinct likely source regions; these proximity sensitivities limit that confounding. The exact-category composition difference **persists** under both limits, but the trait-distance-aware signal remains sensitive to evidence scope, selection and location. Even nearby islands do not necessarily draw on the *same colonist source pool*.

## 6. What the evidence says and does NOT say

**Supported as an exploratory association:** Oceanic and continental islands differ in **the exact repertoire of recorded flower-colour/form state-set profiles** on a geographically and observationally matched subset. Within-island dispersion also tends to be smaller for oceanic islands in that subset.

**Not yet supported:** (a) all colours or broad forms differ; (b) a particular exact profile reliably differs after 496-way FDR; (c) the full functional trait-space *distance* or pollination niches diverge; (d) geological origin **causes** the repertoire change; (e) adaptive convergent evolution or pollinator-channel replacement. Importantly, the original source- and genus-conditioned `island` PR #321 sharing null had no BH/BY-supported result, and the origin-specific fixed binary profile-margin null in `fd` did not support additional local packing. These are **different estimands**, and none overrides the other.

The leading alternative explanations remain **different plant source pools, genus/lineage entry and species loading, mixed native/introduced status, trait observation bias, category compression**, geographic structure, and original geological crosswalk uncertainty. A stronger origin inference requires species-to-profile/genus assignments accessible under the existing source-rights restrictions; the present public island×profile aggregates cannot reconstruct them.

## 7. Recommended paper positioning

Avoid `Oceanic island floral hypervolumes are smaller due to ecological filtering`.

Instead:

> **Geological island histories are associated with different floral trait repertoires, even where broad categorical floral-space positions remain overlapping; identifying source-pool assembly versus interaction-mediated filtering is the next causal problem.**

A possible standalone title is **"Geological contingency in the assembly of island floral phenotype repertoires"**, but the evidence supports **observational associations**, not geological causation.

## 8. Reproducibility

- Original method: `src/fd/repertoire.py`, `src/fd/repertoire_sensitivity.py`.
- Aggregate results in Actions artifact: `results/origin_repertoire/global_repertoire_tests.csv`, `paired_rao_q.csv`, `profile_specific_results.csv`, `all_eligible_profile_prevalence.csv`, `flower_state_results.csv`, and `results/origin_repertoire_sensitivity/sensitivity_repertoire.csv`.
- Source output checksums from verified run:
  - global tests SHA256 `97a0ff5c4556a56195d38ac17a47af37d04ca695ce1a1754b0a9f53bb76ebf20`
  - individual profile tests SHA256 `134f977400678586c306075d763c4b2360f6aa0f41c7f8a1a1a2e0e80c2ab3ce`
  - full eligible prevalence SHA256 `da6c058a362ee768462f7b2993e3884eb0918c23d0f33581a59e7dd87c4b3574`
  - sensitivity SHA256 `978020a49e2cfa38a082d820cbfb2a4132d052d79bab454b71c477a8078df910`
- All source taxonomy/geology matches remain local to reproducible processing; published repository reports aggregate analyses only. `island` Chapter 1 submission and PR #321 remain unchanged.
