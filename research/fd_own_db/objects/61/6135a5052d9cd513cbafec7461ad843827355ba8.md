# Equal-support rarefaction of continental and oceanic island floral repertoires

> **Critical same-day spatial correction:** The original 37 paired islands are **not demonstrably 37 independent geographic comparisons**. Pair connected-components sharing a broad `spatial_block` collapse to **five groups (23, 11, 1, 1, 1)**, with *all 23 northern-midlatitude pairs* in one group. Under **exact group-level sign flips** the original All exact-profile L2 has **p=0.2500**, Direct **p=0.3750**; in the external GBIF-effort-balanced 21 pairs, All **p=0.2500**, Direct **p=0.3125**. An alternative source-archipelago graph also has only five groups. Thus the small individual-pair p values reproduced below and the frequent small p values under category/name thinning **must not be presented as spatially independent discoveries**. The rarefaction identifies robustness of the *descriptive category contrast* to observed-support trimming; it does not fix geographic pseudoreplication. This is a deliberately severe block-dependence stress test, not proof that all islands inside a connected group have perfectly correlated errors. See the [component method](../../docs/COVERAGE_RAREFACTION_SENSITIVITY.md) and [source CI artifact](https://github.com/zuizui0223/fd/actions/runs/42555534972).

**2026-11-09 | Independent `fd` exploratory result.** [Protocol](../../docs/COVERAGE_RAREFACTION_SENSITIVITY.md) · [Source workflow](https://github.com/zuizui0223/fd/actions/runs/42554306259) (50 tests passed; complete G0/17D KDE, origin and effort runs also successful).

## Spatial dependence: exact nulls supersede pair-independent significance

The complete CI rebuilt original `spatial_block` and separately source-defined Weigelt archipelago connected components using *either* endpoint of each origin pair, preserving exact original category profile statistics. Source missing archipelago codes (e.g. literal `NA`) are explicitly **unresolved**, not a shared region. Under either geographical grouping, 37 original pairs collapse to only **five** components:

| Dependency graph | Original 37-pair component sizes | All exact-profile L2 exact p | Direct exact-profile L2 exact p |
|---|---|---:|---:|
| Broad `spatial_block` | 23 / 11 / 1 / 1 / 1 | **0.2500** | **0.3750** |
| Weigelt archipelago groups | 18 / 13 / 4 / 1 / 1 | **0.1875** | **0.5000** |

The effort-balanced 21-pair subset likewise yields five components (spatial blocks 12/5/2/1/1; Weigelt groups 8/6/4/2/1), with exact L2 p **0.2500 All** and **0.3125 Direct** under *either* graph. The original categorical similarity-kernel statistics have exact p≥0.8125 across the graph tests. There are only **2^5=32** equally weighted whole-component sign configurations; because the quadratic metric is invariant to globally reversing every pair, the smallest attainable p for these five-component data is **0.0625**. That is a resolution floor, not the observed p.

**Leaving out the largest geographic component** does not remove the descriptive colour/form exact-category shift: removing the 23-pair `spatial_block` component leaves 14 pairs with exact L2 norm **1.346×** the full 37-pair norm, but the direction in 496-profile space changes substantially (cosine to original displacement **0.511**). Removing the next 11-pair component leaves 26 pairs, norm 1.257× and directional cosine 0.905. A large observed category difference is thus spread unevenly across regions, and the two major components can partially cancel one another. The remaining singletons cannot supply broad geographic replication.

**Scientific revision:** The previously published PR #4 and PR #8 small p/q values from independent-pair sign swaps are **conditional on a stronger independence/exchangeability assumption than this five-component stress test supports**. The headline inference must be demoted to an **observed and repeatedly sampled category repertoire association, not spatially established geological-origin filtering**. This grouping is deliberately conservative and transitive; coarse-block membership does not prove perfect residual dependence, and source archipelago labels are not fully independent phylogenetic source histories.

[Inspect all sixteen exact spatial component tests](spatial_component_tests.csv) or the complete [GitHub Actions artifact and leave-component-out data](https://github.com/zuizui0223/fd/actions/runs/42556434876). **No island identifiers or private matching/trait rows are published.**

## Main conclusion: exact categorical repertoire differences are not solely unequal observed support

The same fixed 37 origin-blind geography/climate/area/coverage matched continental–oceanic pairs (23 northern-midlatitude; 14 tropical) were reconstructed from the original pinned island floral aggregates. Instead of comparing all reported profiles directly, we repeated the analysis after **standardizing observed categories or complete taxon-name carriers per island**.

The 74 original matched islands have at least **44 (All) / 36 (Direct) distinct recorded colour×form profiles**, and at least **64 (All) / 51 (Direct) complete-floral taxon-name carriers**. This is why uniform **30-profile** and **40-name** rarefaction can retain **all 37 original physical pairs** and avoids outcome-dependent island replacement.

The statistic is the square norm of the *mean* oceanic-minus-continental island vector over the 496 **exact recorded floral state-set profiles** (L2). Separate categorical-similarity-kernel comparisons use the *original* full colour/form state-set Jaccard distance, not any 2D or 17D KDE projection. Pair labels are swapped within the same 37 island pairs to construct an approximate conditional null. **No true missing species or flower traits are imputed.**

### Exact profile repertoire results

| Post-hoc support restriction | Eligible pairs | All p median | Direct p median | Fraction of 99 draws with raw p below .05 (All / Direct) |
|---|---:|---:|---:|---:|
| Original full observed support | 37 | 0.0001 | 0.0001 | one deterministic calculation |
| **30 distinct profiles per island** | **37** | **0.0070** | **0.0035** | **83.8% / 93.9%** |
| **40 complete taxon-name carriers per island** | **37** | **0.0070** | **0.0040** | **82.8% / 93.9%** |
| Within-pair equalization to the poorer island's number of profiles | **37** | **0.0005** | **0.0005** | **100% / 100%** |
| Strict single-state floral phenotype sets, ≥20 profile types per partner | **11 All / 10 Direct** | 0.0517 | 0.0930 | one deterministic calculation |

For the 30-profile restricted All comparisons, 95% central **variation over the 99 nested profile-thinning replicates** in approximate p is [0.0005, 0.1763]; Direct [0.0005, 0.0945]. For the 40-name restricted comparisons, All [0.0005, 0.2328], Direct [0.0005, 0.0640]. Some severely thinned realizations are not individually significant; **this is not a 95% confidence interval on a new population effect**.

The 99 repeats are **not 99 independent islands, datasets or experiments**; they are correlated thinnings of the same observed 74 island supports. The `p<.05` fraction above measures conditional robustness to how the finite observed support is sampled, not a hypothesis confirmation probability or external replication frequency.

The single-state-only restriction retains 11 All pairs (8 northern midlat/3 tropical) and 10 Direct (8/2). These are **small and regionally selected** subsets, not the same original population. Their approximate raw p values are not below .05. Thus **the data cannot establish that the original complete-repertoire difference independently persists in single-state-only phenotypes**; the lack of support could reflect genuinely changed encoding, severe attrition and loss of power.

### Original profile-similarity kernel: still no robust functional distance separation

Across 30-profile equalized records, kernel p median is **0.255 All / 0.334 Direct**; under 40-name equalization it is **0.341 All / 0.441 Direct**. Pair-minimum category equalization gives median p **0.078 All / 0.174 Direct**; original full support was 0.081 All / 0.162 Direct.

This is a powerful *negative limit*: the exact **identity/combinations** of recorded colour×form profiles can differ, without a supported separation when similar flower profiles are assigned similar locations in the full original categorical trait space. It is **not** proof that two flowers are pollination-functionally interchangeable.

## Interpretation relative to prior fd results

1. [Original outcome-blind geological origin comparison](../20261108_matched_flower_repertoire/README.md): exact-profile L2 All/Direct BY≈.00042, categorical similarity-kernel difference not BH-supported.
2. [Effort-balanced external GBIF screen](../20261108_gbif_effort_robustness/README.md): 21 origin pairs retain an All exact-profile difference (BY q=.0193) after a *post-hoc* gross record-count and dataset threshold screen. This does not control taxon-specific trait availability.
3. **This follow-up**: exact-category association remains common in equal-profile and equal-name **thinnings** across all 37 pairs; it cannot be attributed solely to one category having more observed carriers or distinct profiles under these screens.
4. **Most important identification caveat:** the [globally coherent genus × four-regime availability conditional randomization](../20261108_genus_source_control/README.md) found *zero* BH-supported global differences after fixing taxon incidence and genus profile repertoires. The present simple support-thinning study does **not** override that null; origin-specific source/lineage assembly remains a viable explanation.
5. [17-dimensional projected KDE](../20261108_source_opportunity_hv17/README.md): raw smoothed volume differences largely shrink toward zero once profile richness is included. This is not a contradiction because categorical profile **composition** and 17D continuous projected smoothed **volume** are different quantities.

## What remains unknown

These calculations operate entirely **inside the recorded observations**. They cannot estimate how many unrecorded species or genera are missing, their unmeasured flowers, actual establishment and colonization, trait-specific taxon detectability, full phylogenetic diversity or ecology of pollinator interactions. In particular, uniformly thinning the observed data is not equivalent to an actual detection-bias correction. A genuine species × effort × survey-replicate detection model cannot be identified from a single list and a broad historical GBIF record-count summary.

**Working paper ceiling:**

> *On geographically matched continental and oceanic islands, the recorded combinations of floral colour and form differ even after standardizing observed profile and complete-name sample support. However, these differences do not establish large categorical distance-space separation or an origin-specific effect beyond genus/source-pool assemblage structure.*

## Reproducibility

- Source-pinned [coverage rarefaction code](../../src/fd/coverage_rarefaction.py).
- [All twenty public aggregate records](support_rarefaction_results.csv), including test sizes, source support, profile L2 and categorical-kernel statistic medians, p percentiles, and Monte Carlo reproducibility fields.
- [Exact CI run #42554306259](https://github.com/zuizui0223/fd/actions/runs/42554306259) (50 unit/integration tests passed), result artifact ID 13646732666.
- All/Direct single-state data are from the same frozen [island PR #321](https://github.com/zuizui0223/island/pull/321) as the full reference. No source-protected taxon-level traits, GIFT source rows, species IDs or island-geological crosswalk records are committed or included in public results.
- Input hashes and original pair support validation are retained in the output `results/coverage_rarefaction/manifest.json` inside the GitHub Actions artifact. Full output CSV SHA256: `daba2f3c1b4c3e3a878e47b994af7c92f179616c7b080541c43b4a3d01c462f9`.
