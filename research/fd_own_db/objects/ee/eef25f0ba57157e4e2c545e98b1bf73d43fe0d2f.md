# Original-distance packing and fixed-support SES (2026-11-09)

**Status:** post-hoc, descriptive and exploratory. This is an honest replacement for
claiming novel floral niche contraction from a 17-dimensional Gaussian KDE volume.
It does NOT test pollinator performance, trait evolution, biological absence,
source-independent geological selection or phylogenetic independence.

## What is measured, using public aggregates only

The common 496-profile reference contains recorded pairs of colour-state sets
and flower-form-state sets. Never multiply multistate values into hypothetical
individual colour × form combinations. The ORIGINAL categorical dissimilarity is

    D_ab = 0.5 * [JaccardDistance(colour_set_a, colour_set_b)
                + JaccardDistance(form_set_a, form_set_b)]

For an island's n >= 2 **distinct** observed profiles S, define:

1. **Mean distinct-pair distance**, P(S) = sum_{a!=b in S} D_ab / [n(n-1)].
   It equals the existing equal-profile Rao Q times n/(n-1).
2. **Mean within-island nearest-profile distance**,
   NN(S) = sum_{a in S} min_{b in S,b!=a} D_ab / n.

Smaller values mean recorded category types are closer under this specific
equal-weight colour/form dissimilarity. They are not a physical volume or
demonstrated functional similarity to pollinators. Both are computed on D,
without truncation or metric distortion by PCoA.

The existing **taxon-carrier-weighted Rao Q** remains a separate abundance
sensitivity. It is not a third independent confirmation because it uses
the same trait distance. Distance-Hill effective numbers with equal similarity
and equal profile weighting can degenerate to n; no such index is presented
as a new information source here.

## Sampling-support reference nulls

- **N0_global496**: uniformly sample exactly n different profiles without
  replacement from all 496 recorded states. For full support, n is the
  observed island's number of profiles. This diagnoses whether the island's
  D-dispersion is different from an unrestricted observed reference.
  It assumes even impossible geographic availability; NOT an ecological null.
- **N1_regional_observed_LOO**: from eligible recorded island floras in
  the same source-defined broad analysis regime (e.g. northern-midlatitude,
  tropical), count *island occurrences* of each profile and remove the focal
  island entirely. Sample n distinct types uniformly from the remaining
  observed region support. This is a leave-one-island-out *empirical* reference,
  not mainland lineage opportunity, colonization, native flora or an unbiased
  complete trait pool. A focal island with more types than the remaining pool
  has no N1 score; do NOT sample with replacement to force one.
- **N1b_regional_islandfreq_LOO**, a deliberately stronger
  **descriptive falsification** of N1: use the SAME leave-one-island-out
  regional support as N1, but draw n distinct profiles without replacement
  with weights proportional to the *number of other eligible islands* in
  the same regime that record each profile. This preserves unequal profile
  prevalence at the broad regional scale. It is not a source or genus null;
  trait-coverage and correlated geography still matter. N1b is excluded
  from the fixed 32 geological paired test slots: it is a post-hoc
  robustness check for whether uniform random reference profiles
  exaggerate recorded packing.
- **N2 genus × native source opportunity** is NOT computed by the public
  pipeline. Restricted GIFT/identity-level inputs and source-donor exclusions
  are necessary. A suitable future N2 must hold trait/sampling evidence,
  name leakage and genus-source opportunity fixed. N1 does NOT stand in for N2.

For each observed endpoint T and its 999 Monte Carlo null values (N0, N1 or N1b), calculate

    delta = T_observed - mean(T_null)
    SES = delta / sd(T_null)            # sample sd, ddof=1

If sd is numerically zero, mark the score nonestimable. No artificial
pseudocount or infinite SES. Record the reference and sample size on every
aggregate comparison.

**Equal-30 rarefaction sensitivity:** only for islands with >=30 different
recorded profiles, take 99 uniform without-replacement 30-profile subsets
and use their mean T as the observed equal-support endpoint. Compare this
to a new fixed-30 null for N0 or N1. The thinning 2.5/97.5 percentiles
measure the variation across highly dependent reductions of ONE recorded
island, not a population confidence interval. No missing plant species or
unrecorded trait state is imputed.

## Geographic origin is not a randomized treatment

Make two explicitly post-hoc paired sensitivity comparisons of oceanic minus
continental **SES**, separately for two nested All/Direct evidence scopes,
two metrics, two nulls and two support rules (32 fixed slots):

1. **Original cohort**: preserve original outcome-blind geographical matching,
   but flip origin signs for entire **transitive shared-spatial-block
   components**. The original 37 pairs contain only five such components.
   This coarse graph is deliberately conservative, not proof that all
   island histories in one block are identical.
2. **Verified geography-disjoint rematching**: reuse the existing globally
   optimized outcome-blind both_verified matching that requires both endpoint
   archipelago labels and no spatial-block or archipelago reuse across
   pairs (currently 10 pairs); enumerate all pair sign flips. Geographic
   disjointness still does not guarantee lineage or ecological independence.

Use two-sided exact sign p-values for the mean SES contrast of N0/N1 only. All 32
exploratory slots (including nonestimable slots assigned p=1) form one
BH/BY family; report no single favorable raw p as independently confirmed.
N0 and N1 SES are correlated, as are All/Direct and full/equal30.
This family is post-hoc and does not retroactively create a preregistered
experiment.

The island observed data/SES are held **only in process memory**. Public CSVs
contain summaries by origin/regime (small cells <5 suppress ALL effects)
and matched aggregate scores only; no island ID, taxon, genus, mainland
source identity or geography crosswalk is persisted. The Action artifact
retains public aggregates, source hashes and methods, not confidential
GIFT/checklist rows.

## Reproduction

The pinned public files and geography source are built in the existing
origin-analysis workflow, then run:

    python -m fd.distance_ses \
      --input work/island_public \
      --g0 work/fd_g0 \
      --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
      --output results/distance_ses \
      --draws 999 --thinnings 99

Outputs are (the summary also includes N1b; the paired table excludes it):
- distance_ses_origin_summary.csv: grouped descriptive observed/null
  distance, delta, SES and thinning spread (only n >=5).
- distance_ses_paired_spatial_sensitivity.csv: 32 post-hoc test slots
  including exact grouped/pair sign probabilities and BH/BY adjustment.
- manifest.json: source hashes, sample designs and interpretation limits.

**Decision gate:** A negative N0 SES alone shows no island-specific
functional contraction. Consistent residual packing after N1 is a
descriptive regional anomaly, not demonstrated filtering. Causal or
lineage-independent inference remains blocked until N2, detection audits,
source lineage structure and adequate independent geographic replication
are supported. The existing 17D KDE belongs in a diagnostic appendix:
its volume is strongly dominated by the number of distinct recorded
profiles when kernels barely overlap and does not define discrete
flower-colour/form space.
