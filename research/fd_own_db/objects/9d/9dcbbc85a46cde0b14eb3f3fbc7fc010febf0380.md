# fd G5: observation-support standardized flower repertoire contrasts

**Prepared 2026-11-09, post-hoc robustness work.** Independent [fd floral phenotype repertoire](../results/20261108_matched_flower_repertoire/README.md), [geological origin groups](../results/20261108_geological_origin/README.md), [external GBIF effort screens](../results/20261108_gbif_effort_robustness/README.md).

## Motivation and estimand

The same 37 geographically outcome-blind, geology-labelled oceanic vs continental matched island pairs show an **exact recorded profile repertoire** difference, but no consistent categorical-trait-distance separation. The stronger genus×geographic-availability profile randomization does not show a global excess (BH q ≥0.1847). The old GBIF historical-effort balancing retains an exact-profile association in 21 pairs, but still does not identify trait-specific recording or missed species.

The present question is narrower and observable: **Does the exact colour × form repertoire difference persist when each island contributes the same number of recorded distinct profiles or recorded complete taxon-name units?** No missing species are inferred and no selection/colonization model is fitted. These are support/recording-resolution sensitivities, not unobserved-species corrections.

## Frozen matched islands and phenotype geometry

Reconstruct **exactly 37** origin-blind pairs from the All/Direct G0 ≥50 eligible complete-floral taxon-name islands. The common profile reference comprises 496 recorded colour state-set × form state-set combinations; multistate descriptions remain whole observed sets.

Retain original category identities and original mean Jaccard distances `D` in the complete 496-profile reference. Per island `v_i` is **uniform over distinct observed profiles**, not their name-carrier counts. Each island pair contributes `delta_j=v_o-v_c`. Test the squared category-repertoire mean contrast `mean(delta)^T I mean(delta)`, and the original similarity geometry statistic `mean(delta)^T(1−D)mean(delta)`. Null is **paired geological-origin sign flips** on the SAME physical pairs; exchangeability is only approximate and does not imply geological labels randomized. Pair-level taxa and profiles are not independent extra island observations.

## Sampling rules, kept distinct

- **Original full**: one deterministic full-support contrast, 9,999 paired flips (baseline original fd PR4).
- **Fixed 30 profiles per island**: every one of the original 74 islands has ≥30 profiles (observed All minimum 44; Direct minimum 36). Uniformly draw 30 complete *profile identities* from each island's original complete profile set without replacement; 99 repeated support draws, each with 1,999 sign flips. This standardizes the number of distinct recorded categories but **not number of individual taxon carriers**.
- **Fixed 40 observed names per island**: each original matched island has ≥50 complete-name carriers. Draw 40 original name-unit carriers **without replacement** from observed per-profile counts using exact multivariate hypergeometric sampling; convert drawn names to the resulting set of distinct observed profiles. Ninety-nine repeated name draws, 1,999 flips each. This standardizes complete observed **taxon-name support**, not unobserved species detection.
- **Within-pair minimum profile richness**: downsample each partner's distinct profile identities to `min(n_profiles_continental,n_profiles_oceanic)` without replacement. All 37 original pairs remain. Ninety-nine repeated profile draws, 1,999 flips. Retains more information than the global fixed-30 rule but sample size varies by pair.
- **Single-state-only**: use the original `all_single_state` and `direct_single_state` 50-name G0 eligibility masks, then require both members of an existing geology-matched pair to have **≥20 distinct single-state profiles**. No rematching or opportunistic pair replacement; preliminary support is small (about 11 All pairs and 10 Direct pairs). One deterministic 9,999 sign-flip contrast for each scope.

All five checks use the same original matched pair identities, original profile geometry and fixed seed `20261109`. For each post-hoc sampled check, report median observed statistic, its 2.5–97.5 percentile variation across **subsamples**, the median paired sign-swap p and 2.5–97.5 percentiles across subsamples, and *fraction of draws* with raw p<0.05. **Neither the distribution of p-values across related subsamples nor its fraction below .05 is an independent replication rate**. All/Direct share taxa and are overlapping confidence sensitivities.

If the single-state cohort has fewer than eight matched pairs, report an explicit insufficiency status and calculate no p-values.

## Spatial dependence audit: conservative connected components

A critical new validity risk emerged from inspecting the 37 matched island endpoints' frozen `spatial_block` labels. Link **two matched pairs if either partner's island shares a spatial block**; then apply transitive closure. The 37 original pairs collapse into **five connected components of size 23, 11, 1, 1 and 1**. All 23 northern-midlatitude pairs form a single connected component; 11 of the 14 tropical pairs form another, plus three tropical singletons.

Under the most conservative null that all paired origin-label flips *inside each connected component must be coupled*, there are **2^5 = 32 exact sign configurations**. The original multivariate L2 and categorical-kernel statistics are **quadratic** and are invariant if every sign is reversed. Therefore **even the most extreme observed result cannot yield an exact upper-tail p below 2/32=0.0625** for all 37 pairs under this null. This is a **finite-sample resolution limit**, not an empirical p-value until the real category data are evaluated. More broadly, prior p=.0001 results from **37 independent pair sign flips** must never be described as verified spatially independent.

The code in `src/fd/spatial_clusters.py` computes full-profile and effort-balanced group tests with **exhaustive** spatial-component flips, plus a **leave-one-connected-component-out** sensitivity to assess whether one geographic source cluster dominates the observed category displacement. Report component sizes and numerical results without exporting island IDs or source blocks. No mixed taxon/genus source variation can be inferred as biological adaptation from this test.

A second, independently coded source-region definition also yields the same risk: connecting pairs with a common **Weigelt published archipelago** label (missing labels treated as separate) yields **five components of size 18, 13, 4, 1 and 1** for the same 37 pairs. The archipelago names are coarse (e.g. a broad Japan region may include biologically distinct island groups), so this also provides a stringent *sensitivity*, not a new empirically estimated covariance. Separate exact component swaps and leave-one-component-out metrics must be reported for both graph definitions.

**This is deliberately a worst-case geographic correlation stress test**. Sharing a coarse 10°-like grid block does not imply biological dependence equals 1, and transitive graph linkage can join distant islands. It is also an observational sign-swap null, not geological treatment randomization. Thus a failure under five components does not prove no real category turnover, but it **does prohibit high-confidence independent replication claims** without a more defensible spatial covariance model or genuinely separate archipelagos/source lineages.

## What a difference can mean

Survival after fixed-profile and fixed-name thinning reduces the hypothesis that the observed exact-category difference is *only* a consequence of one group's larger observed category/name counts. It does NOT demonstrate absence of genus-specific differential detectability, source-pool and phylogenetic history, native/alien difference, or actual floral functional convergence.

Failure to reject on a reduced 30-profile or 40-name support set also does **not** prove the original relationship was artificial: deliberate downsampling reduces information and can reduce statistical power.

Broad single-state exclusion changes both the trait measure and eligible species and islands; it cannot be used to deny multistate phenotypes in nature. Multiple descriptive runs are **exploratory, correlated analyses**, not prospectively independent hypothesis tests. No absence is inferred from missing GBIF records or traits.

## Data rights and replay

Read only the same pinned **public** island×profile carrier aggregates, common original-state-set reference and independently published geological-origin crosswalk used by earlier fd steps. Reconstructed matching and island-level carriers stay in the CI transient scratch area. Public artifact includes only aggregate CSV, audit and SHA256 of pinned public source inputs; **never** original pair IDs, GIFT entity labels, taxon names, private trait provenance, or raw island-geology joins.

```bash
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output work/fd_g0
python -m fd.origin_analysis --g0 work/fd_g0 --hypervolume work/fd_hv \
  --crosswalk work/fd_crosswalk --output results/origin_strata
python -m fd.coverage_rarefaction \
  --input work/island_public \
  --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/coverage_rarefaction \
  --replications 99 --permutations 1999
```

The GitHub workflow executes the original 2D and 17D KDE stages before this sensitivity; neither KDE is the underlying category similarity test. The 17D projected-volume geography association was itself sensitive to original profile counts. This additional test focuses on **exact category repertoire composition**, not volume.
