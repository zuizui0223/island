# Matched geological-origin floral repertoire: analysis protocol

2026-11-08. **Scope:** independent `fd` exploratory extension of the [geological-origin report](../results/20261108_geological_origin/README.md) and [island PR #321](https://github.com/zuizui0223/island/pull/321).

## Main question

**Which recorded colour-state-set × flower-form-state-set combinations are differently represented on oceanic and continental islands, and is the *whole* floral repertoire different after conditioning on comparable island geography and recording coverage?**

This asks about the **composition of recorded phenotypes**, not pollinator activity, selection, species evolution or species-independent colonization.

## Cohort and matching — locked before new floral-profile tests

Start from the independently classified Roeble/Zizka 2022–2024 geological origin crosswalk, All complete records with ≥50 complete taxon-name units, and the **identical All + Direct eligible island population**. Exclude `mixed` and `unresolved`. Match continental/oceanic 1:1 **within northern-midlatitude and tropical regimes**. Northern high latitude has insufficient oceanic support, southern extratropical insufficient pairs.

Each pair must satisfy: log(area ratio) ≤ ln(3), log(1+mainland distance ratio) difference ≤ ln(3), log(recorded complete-name ratio) ≤ ln(3), absolute complete-trait coverage difference ≤0.25, standardized four-climate-PC Euclidean distance ≤3.0 and absolute climate-PC3 difference ≤2.0. Among eligible edges, maximize matching cardinality and minimize standardized covariate difference without replacement. Climate PC3 receives 1.6× cost weight to mitigate balance. These thresholds were selected based on **covariate support diagnostics, not outcomes** and are exploratory, not prospective external preregistration.

Retain all regional matching-balance summaries. The number of pairs, representation by regime, and max absolute pooled standardized mean difference must be reported. The analysis cannot claim perfect matching, source-pool equivalence or species independent phylogenetic control.

## Repertoire outcomes and tests

Construct a binary incidence vector for each matched island over the same 496 **recorded** colour×form state-set profiles. Never expand multistate descriptions into a Cartesian product of hypothetical flower morphs.

1. **Primary multivariate position**: normalize each island's vector to sum 1 across observed profiles, take the oceanic minus continental vector for each matched pair, then compute `mean(delta)^T K mean(delta)`, where `K=1−D` and `D` is the equal-colour/form Jaccard dissimilarity on recorded state sets. Use an outcome-blind 9,999 paired sign-swap reference to test whether a systematic origin difference exists. Verify the profile similarity kernel is positive semidefinite; no forced projection to PC1/PC2.
2. **Independent exact-profile-space comparator**: identical test but `K=I`, detecting exact-category repertoire shifts regardless of similarity. BH and BY adjustment is one four-test family (two tests × All/Direct).
3. **Individual recorded profiles**: paired binary discordance tests for each of the **496 original All reference profiles** (even if absent in both groups, p=1), McNemar-style two-sided exact binomial; BH/BY across the full fixed 496 family in each evidence scope. Report differences and detection support even if no test survives adjustment. All/Direct share the same matched pairs and are sensitivity analyses, not independent replications.
4. **Secondary component states**: record whether each of the original colour or form states appears among recorded profiles on each island; analogous paired exact tests. States are recorded trait descriptions that might not always represent realized morphs on the specific island, so these are **hypothesis-generating**.
5. **Local functional dispersion**: compute the observed paired mean difference in equal-profile Rao Q and its two-sided sign-swap p as a parallel check, not as a new causal inference.

## Post-hoc source-geography and tighter-balance sensitivities

The initial botanical-source audit found that some tropical pairs satisfying the main environmental matching criteria were thousands of kilometres apart, potentially linking different source floras. Therefore report three additional **post-hoc diagnostic** cohorts without changing the original primary result: (i) maximum absolute complete-trait coverage 0.12 and PC3 difference 1.0; (ii) a strict spherical geographic distance between paired island centroids ≤5,000 km; (iii) ≤3,000 km. These are **proxies for shared biogeographic source opportunity**, not actual genus-level or historical plant source pools. If the 5,000-km restriction removes the apparent functional repertoire difference, the geographic-source alternative has stronger support than the geological origin mechanism. Pair sample sizes and regional balance must be disclosed for each diagnostic.

## Null interpretation and uncertainty

The sign-swap test assumes **conditional exchangeability of the two origin labels within geographically matched pairs**. Geological origin was not assigned at random, and residual climate imbalance, regional source pools and intra-archipelago dependence remain. Thus p/q are exploratory and **not** a causal test of geological origin. The paired unit is islands, not species. Any global significance establishes only a recorded assemblage composition association at the matched support; it does not demonstrate historical trait evolution or changed pollination services.

Only public, rights-cleared island×profile aggregate carriers and profile reference definitions are read. The inferred island-to-geology crosswalk must stay **private/local to CI**; no row-level source protected geology join or species×trait ledger should enter git or public action artifacts. The uploaded outputs contain cohort counts, aggregate effect estimates, category summaries and source hashes, not matched island identities.

**Next decisive gate if a composition difference emerges:** source-pool and genus/lineage reassignment using authorized taxon identity-level data; split native vs introduced taxa and replicate with independent geologic-origin records. A strong profile-level differential pattern without this gate is still an **assembly association**, not selection.

## Replay

```bash
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output work/fd_g0
python -m fd.origin_analysis --g0 work/fd_g0 --hypervolume work/fd_hv \
  --crosswalk work/fd_crosswalk --output results/origin_strata
python -m fd.repertoire --input work/island_public --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/origin_repertoire --permutations 9999
```

The GitHub Actions workflow also generates the projected 2-D KDE sensitivity files that the origin analysis reads. Do **not** label a 2-D KDE area as the entire functional hypervolume (the first two axes retain only 15.35% positive inertia).
