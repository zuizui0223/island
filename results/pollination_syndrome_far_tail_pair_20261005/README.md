# Pollination-syndrome far-tail pair contrast — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

This diagnostic resolves the remaining comparison left open by the pooled region × hinge model for the
`yellow/orange × butterfly-associated deep-tube` response. The pooled deep-tube hinge model did not
converge, although the region-specific hinge models did. We therefore compare the converged northern
mid-latitude and tropical region-specific slopes directly.

## Design

The model is unchanged from the preceding nonlinearity diagnostic:

- frozen raw `architecture | colour` island counts;
- corrected mainland distance;
- adjustment for `selfing_core`, island area and climate PC1–4;
- spatial-block cluster-robust uncertainty;
- fixed hinge at **783.017 km**, the upper edge of four-region common 5–95% isolation support.

For two disjoint regional fits, the variance of the slope difference is taken as the sum of the two
regional variances. This avoids the nonconvergent pooled parameterization without changing the regional
estimands.

## Result

### Below the common-support upper edge

Northern mid-latitude and tropical slopes are indistinguishable.

- All: north-mid minus tropical = **+0.0183**, P = **0.593**.
- Direct: north-mid minus tropical = **−0.00694**, P = **0.819**.

Thus there is no evidence for a tropical-specific slope over the distance range shared well by both
regions.

### Above 783 km

The key new result is that the same response also turns positive in northern mid-latitudes.

Northern mid-latitude post-hinge slopes:

- All: **+0.5023**, P = **0.0276**;
- Direct: **+0.4589**, P = **0.0212**.

Tropical post-hinge slopes:

- All: **+0.1419**, P = **0.113**;
- Direct: **+0.2015**, P = **0.000938**.

The north-mid versus tropical post-hinge contrasts are unsupported:

- All: difference = **+0.3605**, P = **0.141**;
- Direct: difference = **+0.2574**, P = **0.217**.

The hinge-change contrasts are also unsupported:

- All: P = **0.185**;
- Direct: P = **0.229**.

The far-tail samples are smaller in northern mid-latitudes but remain represented by 62–63 islands
across 17 spatial blocks, compared with 349–376 tropical islands across 53–54 blocks.

## Interpretation

The positive yellow/orange–deep-tube response should no longer be framed as a uniquely tropical
pollination-syndrome signal. The stronger interpretation is:

> **A butterfly-like deep-tube component emerges primarily at very large mainland distances and is
> detectable in both northern mid-latitude and tropical island floras.**

This points to **scale of isolation** rather than tropical biome identity as the cleaner organizing
variable for this specific response.

This does not imply that the realized pollinators are butterflies, nor that one pollinator guild
replaced another historically. The phenotype label remains a syndrome-concordance descriptor only.

## Relation to the broader regional picture

This result does not collapse all regional biology into one isolation-scale effect.

- Northern mid-latitude large-bee-associated form still declines within common support.
- Northern high-latitude specialized blue/purple architecture declines more strongly than northern
  mid-latitudes within like-for-like support.
- Southern bird-associated patterns remain mixed and non-monotonic.
- The yellow/orange deep-tube positive response is the component that now looks most clearly
  **remote-island rather than tropical-specific**.

The combined picture is therefore:

> **pollination-syndrome concordance depends on both regional context and isolation scale, but some
> apparent tropical signals are better explained as very-remote-island responses shared across regions.**

## Submission boundary

Do **not** add these 2026-10-05 coefficients to the current Chapter 1 manuscript or SI.

The active manuscript should only avoid treating full-range regional syndrome associations as
like-for-like biome contrasts and should retain the boundary that phenotype does not identify realized
pollinators.

## Reproducibility

- workflow run: **37294802802**
- artifact: **11338470536**
- script: `scripts/audit_pollination_syndrome_nonlinearity.py`
- machine-readable contrast:
  `results/pollination_syndrome_far_tail_pair_20261005/far_tail_pair_contrasts.csv`
