# Independent geographic support: globally optimized origin pair design

2026-11-09. **Post-hoc alternative cohort**, independent of the original 37-pair matching algorithm's results. Original scientific analysis: [37-pair profiles](../results/20261108_matched_flower_repertoire/README.md); subsequent [54-test spatial reappraisal](../results/20261109_coverage_rarefaction/README.md).

## New inferential question

When geographically and historically related matched island pairs are counted as **one rather than many repeated comparison units**, can a genuinely broader set of independently paired continental/oceanic islands still exhibit different *recorded colour × flower form combinations*?

This is an ecological-assemblage **observational** contrast, not geological treatment randomization, plant establishment or pollinator function.

## Why re-matching rather than slicing the existing 37 pairs?

The original 37 geography/caliper-matched pairs form only **five transitive spatial components (23/11/1/1/1)**, or five **Weigelt source archipelago components (18/13/4/1/1)** after correctly parsing missing archipelago `NA` values. Exact component-wide sign flips yielded All p=.25 / Direct p=.375 under spatial blocks and All p=.1875 / Direct p=.50 under the archipelago grouping, not the earlier pair-independent p=.0001. Simply retaining one pair per original component gives little geographic replication.

This additional test searches the **entire geographically origin-classified All+Direct eligible pool**, not the original 37 matches, for a larger set of **geographically disjoint** pairs under unchanged **outcome-blind** area, isolation, trait-coverage, complete-name count and climate constraints. It is explicitly a *different*, post-hoc subset, not external biological replication.

## Frozen candidate and matching design

- Independently audited Roeble 2024 / Zizka 2022 geological origin labels (`continental`, `oceanic`; mixed/unresolved excluded).
- Both All and Direct have ≥50 complete-floral name records; same two originally supported biogeographic regimes **northern_midlatitude and tropical**. Northern-high/southern-extratropical regimes are not added after outcome inspection (they have too few paired opportunities for the original PR4 design).
- Build all possible within-regime continental × oceanic island candidate pairs. Freeze the **same original PR4 calipers**: log area ratio ≤ln3; log(1+continent distance) difference ≤ln3; complete-name count log ratio ≤ln3; joint trait-coverage absolute difference ≤0.25; standardized Euclidean four-climate-PC distance ≤3; raw climate PC3 difference ≤2; original standardized covariate cost gives 1.6× PC3 weight. No flower colour/form metric is consulted before matching.
- **Global binary integer set-packing**: maximize number of matched pairs, with covariate cost a much smaller objective contribution; each physical island may be in at most one pair. Restrictions are global across both regimes, not separate per region.
- Matching sensitivity definitions: `spatial_block` (no 2 selected pairs use a common spatial block via either endpoint); `weigelt_archipelago` (same condition on source-published archipelago, but literal missing labels create **no shared** archipelago); `both` (no reused block **or** labelled source archipelago, with missing labels explicitly unknown); and **`both_verified`**, which requires **both islands in every pair have a known archipelago** and enforces both block and archipelago exclusivity. The **10-pair `both_verified`** group is the strictest geographical inference sensitivity.
- Grouping is based entirely on **public geographic source context**, not trait dissimilarity or endpoint outcome.
- Under preliminary source-pinned geography support, 219 eligible within-regime candidate edges yield 13/11/11 non-overlapping pairs across the first three rules. **`both_verified`** excludes 9 edges involving an unknown archipelago, leaving 210 eligible candidate edges and **10 maximally disjoint pairs (5 northern-midlatitude/5 tropical)**. No taxon→profile information is used to determine these matching totals.

## Exact tests, multiplicity and precision

On the fixed matching of each rule, construct the **496-dimensional normalized distinct-profile incidence vector** from each original All or Direct island; observed difference = mean(oceanic profile − continental profile) across selected pairs. Test two multivariate endpoints:
- **exact-profile L2 squared difference**: `mean_delta^T mean_delta`;
- **original colour/form set-Jaccard similarity kernel**: `mean_delta^T (1-D) mean_delta`.

Enumerate ALL `2^n` sign assignments for each selected pair cohort (n expected 10–13, 1,024–8,192 configurations). The quadratic statistic is invariant to globally reversing every sign, so minimum possible p is `2/2^n`. No Monte Carlo approximation or artificially independent within-archipelago resampling. **One post-hoc BH/BY family of **sixteen** test slots** (4 different groupings × 2 nested All/Direct evidence scopes × 2 multivariate endpoints). Report exact upper-tail p and adjusted q; 2D/17D KDE is NOT substituted for categorical trait distance.

Geographic nonoverlap does not prove full biological independence; e.g. a genus can appear in several remote island floras, true source lineages can span archipelagos, and the two genotype/profile scopes overlap. A strong geographical result would still require source-/genus-conditioned and phylogenetically informed analysis. A nonsignificant result in 11 pairs is not absence of a true biological effect: matched support is small.

## Required release and checks

- Fail if the source-pinned geography or matching support changes without explicit audit, or if MILP cannot verify optimality.
- Unit tests with *synthetic* candidate islands establish literal `NA` handling, no reused island/region resources, matching-cardinality dominance over cost and exact p-value resolution.
- CI recomputes public frozen aggregate floral profiles, geology and origin, and writes **aggregate** `independent_matching_support.csv`, `independent_origin_repertoire_tests.csv` and `manifest.json`, without raw matched island identifiers or source archipelago names.
- Restricted original taxon×trait rows, GIFT source data and privately joined geography must never be committed.

```bash
python -m fd.geography_independent \
  --input work/island_public \
  --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/geography_independent
```

**Manuscript claim ceiling pending the 11-pair exact test:** recorded floral category associations, not confirmed geological filtering, naturalization, evolutionary simplification or pollination niche changes.
