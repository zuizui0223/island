# Geography-disjoint continental/oceanic island re-matching: methods and inference

**2026-11-09.** Independent `fd` study of recorded floral colour×form state-set profiles. The core scientific issue arose in [PR #9](https://github.com/zuizui0223/fd/pull/9): 37 geographically matched island pairs were joined into **five** broad connected spatial-source groups (23/11/1/1/1), so their independent-pair sign-swap p=.0001 cannot establish geological-origin effects after severe geographic dependence is acknowledged.

## Question

If we **select new continental–oceanic island pairs using only origin, environment and geographic source labels**, never floral traits, and ensure different pairs do not reuse a spatial block or a published archipelago, can the recorded exact flower-colour × flower-form profile repertoire difference still be distinguished from a paired sign-swap null?

This is a **post-hoc, more geographically dispersed observation-based comparison**, not independent external replication, random geological treatment, verified phylogenetic independence or known source histories.

## Candidate population and independently fixed calipers

Start from the same Roeble/Zizka independently classified geological oceanic/continental island labels and frozen `island` PR #321 public island×profile aggregates. Require eligibility in both All and Direct with at least 50 complete colour/form name units. Match within two originally specified regions, northern midlatitude and tropical. Preserve PR #4's outcome-blind constraints and all of its covariates:

- Difference in `log(area_km2)`, `log1p(distance_to_continent_km)`, `log(n_names)`: each at most `log(3)`.
- Difference in complete trait coverage `joint_coverage`: at most 0.25.
- Euclidean standardized climate PC1–4 separation: at most 3.
- Absolute raw climate PC3 difference: at most 2.
- Geographical match cost: original normalized covariate differences with PC3 1.6× weight.

The real-data public support is **222 candidate origin-labelled eligible islands and 219 eligible continental/oceanic edges**, including **125** north-midlatitude and **94** tropical candidate edges. The identities of pairs and geographic source associations are private intermediate values, never published.

## Geographic set-packing

Use a binary integer-programming set-packing model on original candidate edges. Choose as many origin-contrasting pairs as possible before optimizing environmental covariate similarity, subject to each physical island used once. In addition to no shared island ID, apply one of four outcome-blind source-resource exclusion policies to **different matched pairs**:

1. `spatial_block`: none of either pair's continental or oceanic endpoints may reuse a frozen broad spatial block used by another selected pair.
2. `weigelt_archipelago`: none may reuse a known published Weigelt archipelago label. A missing/unknown source label (`NA`, blank, etc.) does **not** create a shared block, but neither does it certify independent source ancestry.
3. `both`: enforce both recorded block and archipelago exclusions. Unknown archipelago labels do not create a shared group; their island IDs remain protected.
4. `both_verified`: require all four source archipelago labels **for both endpoints** of every pair; drop the nine edges involving unknown labels **before matching**; enforce both source-exclusion policies.

The binary MILP must return optimal status with no reused island or geographic resource. This protects against the original obvious transitive-block pseudoreplication, **not** against spatially adjacent but differently labelled blocks, cross-archipelago lineage dispersal or broad geographically mixed source groups.

## Outcome and exact inference

Each island is a vector of distinct recorded colour/form profile incidence among the same original **496** state-set profiles, normalized to sum one. For selected matched pairs, `delta_j = v_oceanic - v_continental` and `mu = mean(delta_j)`. Test two endpoints:

- Exact recorded profile composition: `mu' I mu`, a squared categorical identity shift.
- Original full categorical geometry: `mu' (1-D) mu`, with `D` the average original state-set Jaccard distance across colour and form, **not** 2D/17D PCoA KDE.

All All/Direct scopes share selected physical island pairs. Enumerate **all 2^n origin-label reversals**, with n up to 13 in this experiment, no Monte Carlo correction, and exact upper-tail fraction of null `T` at least the observed `T`. The squared statistic is invariant under a global sign reversal, yielding a minimum attainable two-sided-sign-symmetric upper-tail p of `2/2^n`. Report effect size and median/95% interval of the finite reference distribution.

Declare ONE post-hoc **16-test** BH/BY family (four alternative matching policies × two nested evidence scopes × two endpoints). No individual test p should be called a new discovery without accounting for the overlapping search and matching masks. The **observed sample population changes** between schemes; test p values assess their own geographically constrained matched samples, not a common cohort.

## Critical limitations

**Even pairs with different published archipelago IDs and grid cells are not guaranteed independent**: the source labels are operational categories, not reconstructed phylogenetic colonist sources, dispersal corridors or spatial covariance estimates. The published Weigelt archipelago categories may themselves be broad or geographically heterogeneous. Multiple origin groups across climate zones need not share plant-source genus distributions; genus×regional-availability conditional randomizations from PR #5 did **not** support a pooled geological effect.

No recorded taxon/generic zero is declared biological absence; the analysis is not an establishment probability model. No trait-specific recording detection, experimental pollination service, within-lineage evolution or true colonist ancestry is recovered. The 496 exact combinations refer to recorded **state sets**, not measured individual-plant flower morph frequencies.

## Reproducibility and privacy

Implementation: [`src/fd/geography_independent.py`](../src/fd/geography_independent.py), tested through synthetic set-packing and exact finite-sign enumerations. Run from the pinned public aggregate bundle, previously validated G0 and independently crosswalked published geology:

```bash
python -m fd.geography_independent \
  --input work/island_public \
  --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/geography_independent
```

The public workflow uploads only aggregate counts, exact tests and source file hashes. Never commit or release named matching islands, physical island IDs, source geographies, original taxon×flower traits or upstream restricted evidence source records.
