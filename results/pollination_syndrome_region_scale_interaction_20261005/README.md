# Pollination-syndrome region × isolation diagnostic — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

This diagnostic follows the regional common-support and nonlinearity analyses. It directly tests the guardrail that a significant coefficient in one region and a nonsignificant coefficient in another is not itself evidence of regional heterogeneity.

## Design

The same frozen raw `architecture | colour` responses are compared among all four regions.

Two models are used:

1. **common-support linear model:** all four regions are restricted to **0.2515–783.0171 km**, the intersection of their regional 5th–95th percentile isolation ranges;
2. **full-support hinge model:** all observations above the common lower bound are retained and a fixed hinge is placed at **783.0171 km**.

Each region receives its own intercept and its own nuisance slopes for `selfing_core`, island area and climate PC1–4. Isolation slopes are on the raw `log(1 + mainland distance km)` scale. Spatial-block cluster-robust uncertainty is retained.

The common-support model is the clean like-for-like regional comparison. The hinge model is diagnostic of scale dependence outside shared support.

## Main result

### Northern mid-latitude versus tropical

For the previously highlighted **yellow/orange × butterfly-associated deep-tube** response, the isolation slope does **not** differ between northern mid-latitude and tropical islands over shared distance support:

- All: north-mid minus tropical slope = **+0.03384**, SE = 0.03694, P = **0.360**;
- Direct: difference = **+0.00975**, SE = 0.03158, P = **0.757**.

Thus the positive tropical full-range signal is not evidence of a tropical-versus-northern difference over like-for-like isolation. Together with the preceding nonlinearity diagnostic, the safer interpretation is that the positive tropical signal is **remote-island concentrated**.

The four-region common-support omnibus for this response is nominally heterogeneous but does not survive the across-combination BH diagnostic:

- All: P = **0.00917**, q = **0.05495**;
- Direct: P = **0.00879**, q = **0.08176**.

The pooled four-region hinge model for this particular response remains numerically nonconvergent even after a 5,000-iteration retry. The separate regional hinge fits converge, but a formal pooled north-mid versus tropical far-tail contrast is therefore **unresolved**, not null.

For **yellow/orange × large-bee-associated form**, the north-mid versus tropical slope difference is also unsupported within common support:

- All: difference = **−0.03011**, P = **0.183**;
- Direct: difference = **−0.01670**, P = **0.497**.

The All-only hinge fit suggests stronger far-tail divergence, but Direct does not corroborate it. This is not promoted as a robust regional mechanism.

### Northern mid-latitude versus northern high latitude

A different pattern appears within the Northern Hemisphere. The stronger high-latitude loss of specialized architecture is detectable **within the same shared isolation range**.

For **blue/purple × large-bee-associated intermediate/deep tube**:

- All common-support north-mid minus north-high slope = **+0.08732**, P = **0.0262**;
- Direct = **+0.09494**, P = **0.0119**.

For **blue/purple × butterfly-associated form**:

- All difference = **+0.06039**, P = **0.0153**;
- Direct = **+0.05768**, P = **0.0370**.

Positive contrasts here mean that the northern-high slope is more negative than the northern-mid slope. The difference is therefore not explained simply by northern mid-latitudes containing more near-mainland islands.

The four-region common-support joint tests are only borderline after the across-combination BH diagnostic (minimum q ≈ 0.055 in All and ≈ 0.082 in Direct), so these pairwise results remain **post-hoc diagnostic contrasts**, not a new confirmatory regional theorem.

### Full-support hinge heterogeneity

For several specialized blue/purple and yellow/orange architecture responses, the full-support hinge models show strong four-region heterogeneity in pre-hinge slopes, hinge changes and/or post-hinge slopes. This confirms that regional syndrome concordance is scale dependent.

However, the hinge coefficients should not be interpreted as biological thresholds. The 783 km knot is derived from support geometry, not optimized against trait responses.

## Interpretation

The original concern was that northern mid-latitude islands might be too concentrated near continents to distinguish their pollination-related response. The combined diagnostics now reject that simple explanation.

The supported picture is:

1. **Northern mid-latitude:** the recurrent reproductive/generalized H1 signal and the large-bee-associated architecture decline are already visible inside shared regional isolation support.
2. **Northern high latitude:** specialized blue/purple architecture declines more strongly than in northern mid-latitudes even over like-for-like isolation distances.
3. **Tropics:** the positive butterfly-like/deep-tube signal is not a tropical-versus-northern difference within shared support; it is concentrated in the broader remote-island regime.
4. **Southern extratropics:** responses remain mixed/non-monotonic rather than forming one coherent bird syndrome.

The best current conceptual statement is therefore:

> **Pollination-syndrome concordance depends on both biogeographic context and the scale of isolation represented. Some northern regional differences persist after matching isolation support, whereas the tropical positive long-tongued-insect signal is primarily a remote-island phenomenon.**

This is stronger than attributing every regional contrast to sampling geometry, but weaker—and more defensible—than inferring realized pollinator replacement.

## Submission boundary

Do **not** add the 2026-10-05 common-support, nonlinear or interaction coefficients to the current Chapter 1 manuscript or SI.

The active manuscript should only:
- acknowledge that regional isolation ranges differ;
- avoid treating regional syndrome patterns as like-for-like causal pollinator comparisons;
- retain the inferential boundary that phenotype does not identify realized pollinators.

## Reproducibility

- workflow run: **37290669484**
- artifact: **11335824172**
- script: `scripts/audit_pollination_syndrome_region_scale_interaction.py`
- workflow: `.github/workflows/pollination-syndrome-region-scale-interaction.yml`
- parent diagnostics:
  - `results/h1_regional_common_support_20261005/README.md`
  - `results/pollination_syndrome_nonlinearity_20261005/README.md`
