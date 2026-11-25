# fd G1 + hypervolume extension — methods and reporting contract

**Status: exploratory, code-backed.** Source data are frozen at [island PR #321](https://github.com/zuizui0223/island/pull/321), commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`. This contract makes a genuine hypervolume calculation while preventing a two-axis projection from being mislabeled total functional trait space.

## Two **different** hypervolume notions

1. **Discrete recorded-category occupancy**: the number of reference colour-state-set × form-state-set profiles observed on an island, divided by the fixed global 496-profile reference. This is a finite categorical occupancy measure, **equivalent to observed profile richness** in this representation, and not an independent Gaussian or geometric-volume discovery.
2. **Smoothed projected KDE hypervolume**: a 90%-mass highest-density region (HDR) of the mixture of isotropic 2D Gaussians centered on the common reference PCoA1/PCoA2 coordinates for all observed island profiles. Both axes are scaled by their global 496-reference standard deviations. Fixed bandwidths `h=0.35` and `h=0.55` are tested; equal profile weights are primary and taxon-name-carrier weights are sensitivity. Each HDR area uses a common 88×88 midpoint integration grid spanning all 496 profile coordinates ±4 bandwidths. Grid probability-mass capture must be 0.994–1.006; no missing mass is silently normalized. HDR area has units of squared standardized PCoA axes.

**Nonclaims:** Two PCoA axes retain **15.35% of positive inertia** (80% requires 17 axes). Any KDE area reported here is a **projected 2D smoothed phenotype-density area**, not a 17-dimensional Gaussian hypervolume or volume of realized pollination functions. Continuous Gaussian smoothing invents mathematical locations between categorical states, so bandwidth sensitivity is mandatory. Total pollen-transfer niche volume, colonization filtering, and evolution cannot be inferred.

Source [geometry feasibility and PCoA audit](https://github.com/zuizui0223/island/blob/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/floral_space_20261106/geometry_notes.md).

## G1 fixed-margin reference

The estimands intentionally differ from the carrier-weighted G0 values:

- **Local:** equal-profile Rao's Q (pairwise average of original colour/form state-set Jaccard distances), median across a fixed **128-island ID-selected subset** per scope/regime (all available if fewer).
- **Between islands:** equal-profile symmetrized minimum-distance profile matching, median of **128 ID-selected upstream matched pairs** per scope/regime (all available if fewer).
- **Fixed population:** all islands appearing in the original PR321 matched-pair cohort within a regime/scope. The Curveball swaps operate on the complete observed binary island×profile matrix of that population, *not just the scored subset*.
- **Fixed row margin:** number of observed profiles per island.
- **Fixed column margin:** number of those islands where each profile appears.
- **MCMC:** two separately seeded chains, each with 12×number_of_islands burn-in trades and 3×number_of_islands trades between 200 retained draws; both binary margins are checked against their originals.
- **Mixing gate:** R-hat ≤ 1.05, absolute lag-1 autocorrelation ≤ 0.30, and ≥8 distinct rounded statistics; no p reported when this fails.
- **Tail and multiplicity:** two-sided inclusive tails with 1/(400+1) continuity correction. BH and BY adjustments across one frozen **32-test family** (2 endpoints × 4 regions × 4 evidence scopes). These are **approximate MCMC tails**, not exact i.i.d. randomization p-values. A larger independent repeat and changed seeds are necessary before confirmatory publication.

The same pair and local-island selections are used for All/Direct within a given encoding and region. All, Direct and their single-state restrictions are overlapping sensitivity populations, **not independent replications**. Pairs reuse islands and 128 sampled pairs are not 128 independent ecological observations. The selected subset changes the estimand compared with full PR321 sharing metrics. Selection depends on prior pair IDs only, not on functional distances.

**Why equal-profile weights for G1?** Fixed binary row/column margins do **not** preserve how many taxon-name carriers each profile has. Testing carrier-weighted Rao's Q or carrier-weighted NN against this null would confound the count and profile-incidence processes. G0 retains species-carrier-weighted metrics, but the current G1 significance claim (if any) pertains strictly to equally weighted recorded profiles. KDE hypervolumes are not subjected to G1 significance tests in this extension.

**What is not held fixed:** source-pool membership, genus/lineage, plant species identity, taxon sharing, real colonization opportunity, provenance quality, native status, spatial dependence, and ecological interactions. Both G1 and projected KDE are exploratory. An observed G1 departure is not evidence of plant evolution or pollinator functional replacement, and a PR321 genus-conditioned non-rejection is not contradicted by a different-margins result.

## Data contract

All inputs are public aggregates, copied at a pinned upstream commit by `scripts/fetch_island_aggregates.py`. The original taxon×trait/source ledger remains private under the [island rights audit](https://github.com/zuizui0223/island/blob/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md).

```bash
python -m pip install -e '.[test]'
python -m pytest -q
python scripts/fetch_island_aggregates.py --output work/island_public
python -m fd.pipeline --input work/island_public --output results/run_fd_v1
python -m fd.g1_fixed_margins --input work/island_public --g0 results/run_fd_v1 --output results/run_fd_g1 --draws 200 --local-n 128 --pair-n 128
python -m fd.hypervolume --input work/island_public --g0 results/run_fd_v1 --output results/run_fd_hypervolume
```

The GitHub Actions CI retains aggregate summaries, raw MCMC test statistics (no taxon IDs), output hashes and projected KDE HDR areas as an artifact. `g1_summary.csv` is the only significance output; `regional_hypervolume_summary.csv` is descriptive, not an inference result.
