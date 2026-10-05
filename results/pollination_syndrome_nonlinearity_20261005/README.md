# Pollination-syndrome nonlinearity diagnostic — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

> **Superseded biological interpretation (2026-10-05; PR #313).**  
> The coefficients below are retained as a historical post-hoc diagnostic trail. Later
> species/status attribution showed that the northern mid-latitude positive far tail
> disappears in confirmed native flora. In tropical native flora, the positive tail
> disappears after removing literature-qualified ocean-dispersed coastal species
> (native-nonendemic **+0.9844 → -0.2330**; **-0.2294** after recomputing
> `selfing_core`), whereas tube-depth coding flags alone do not remove it. These
> outputs therefore should **not** be interpreted as a pollination-syndrome mechanism
> or as a shared very-remote-island biological response across regions. Current
> follow-up: `results/tropical_far_tail_dispersal_20261005/README.md`.


This diagnostic follows the regional common-support analysis and asks whether the focal
colour × architecture responses change slope across isolation regimes rather than following
one linear mainland-distance response.

## Design

The analysis uses the frozen raw architecture-given-colour island counts, corrected geography,
selfing-core adjustment, island area and climate PC1–4, with spatial-block cluster-robust
uncertainty.

No knot was tuned to a trait result. Both knots come from the four-region H1 support geometry:

- **127.5559 km**: upper edge of the common regional IQR;
- **783.0171 km**: upper edge of the common regional 5–95% support.

Two diagnostics were fitted:

1. one hinge at 783.0171 km;
2. three segments: <=127.5559 km, 127.5559–783.0171 km, and >783.0171 km.

The corrected full linear models were replayed before the nonlinear fits and matched the
committed corrected raw-coupling tables.

## Main results

### Northern mid-latitudes: specialization loss is not hidden by near-mainland sampling

For yellow/orange flowers with large-bee-associated form:

- All, one-hinge model: the slope below 783 km is negative
  (beta = -0.03692, P = 0.000258); the extra slope change above 783 km is unsupported
  (P = 0.179), and the hinge does not improve AIC (Delta AIC = +0.33).
- Direct: the below-783 slope is also negative (beta = -0.02791, P = 0.00342);
  the far-tail change is unsupported and the hinge model is worse (Delta AIC = +1.51).

The three-segment All model suggests scale variation: a negative <=127.6 km slope,
near-flat 127.6–783 km slope, and a negative >783 km slope. However, that extra complexity
is not supported in Direct evidence (three-segment Delta AIC = +2.33).

The robust conclusion is therefore **not** that northern mid-latitudes lack remote islands
needed to reveal specialization loss. The negative large-bee-associated architecture response
is already detectable within the shared regional isolation support.

### Northern high latitudes: strongest specialization loss at intermediate isolation

For blue/purple butterfly-associated form:

- the <=127.6 km segment is weak;
- the slope changes strongly negative above 127.6 km;
- the 127.6–783 km segment is strongly negative in All and Direct;
- there is no supported additional change above 783 km.

The three-segment model improves AIC by 10.54 (All) and 9.41 (Direct) relative to a linear
model. This places the strongest specialization loss in the **intermediate isolation regime**,
not specifically in the most remote tail.

Blue/purple large-bee-associated intermediate/deep tubes show the same qualitative pattern:
the strongest Direct decline appears after 127.6 km and before 783 km, while the far tail is
too sparse to support a further directional change.

### Tropics: the positive long-tongued-insect signal is remote-island concentrated

For yellow/orange butterfly-associated deep tubes:

- Direct, <=783 km slope: beta = +0.02156, P = 0.231;
- Direct, change above 783 km: beta = +0.17991, P = 0.00956;
- Direct, >783 km slope: beta = +0.20147, P = 0.000938;
- one-hinge Delta AIC = -16.26 relative to the linear model.

The two-knot model gives the same biological picture but does **not** identify a sharp
single threshold: <=127.6 km is essentially flat, 127.6–783 km is weakly positive,
and >783 km is clearly positive. The individual change parameters are correlated and are not
both significant in the two-knot fit, while the far-segment slope remains supported
(beta = +0.18560, P = 0.00112; three-segment Delta AIC = -16.83).

Thus the evidence supports **curvature / remote-island concentration**, not a claim that
783 km is a biological threshold.

The All layer also prefers nonlinear fits by AIC, but its cluster-robust segment coefficients
are individually less precise. The strongest inferential evidence is therefore the Direct
pattern plus the common-support collapse already documented in the preceding diagnostic.

### Southern extratropics: no single monotonic bird syndrome

The southern All bird-associated form model is non-monotonic: negative at low isolation,
positive at intermediate isolation, and negative again in the far segment. Deep-tube
architecture is most positive in the intermediate segment. Direct evidence is much less
decisive.

This reinforces the previous conclusion that southern floral restructuring does not form one
coherent monotonic bird-syndrome response.

## Biological interpretation

The regional contrast is not explained by one sampling artefact.

- Northern mid-latitude specialization loss is detectable despite the concentration of
  near-mainland islands.
- Northern high-latitude specialization loss is strongest over intermediate isolation.
- The tropical positive long-tongued-insect-associated signal is concentrated in the
  broader remote-island regime and should not be treated as a like-for-like biome contrast
  with northern regions.
- Southern responses are non-monotonic and internally mixed.

A better conceptual model is therefore:

> **pollination-syndrome concordance depends on both biogeographic region and the scale of
> isolation represented.**

This is compatible with the possibility that very remote tropical islands retain or assemble
different pollination architectures than moderately isolated northern islands, but these
phenotypes still do not identify realized pollinators or demonstrate historical pollinator
replacement.

## Submission boundary

Do not promote these 2026-10-05 nonlinear numerical results into the current Chapter 1
submission. They are post-hoc diagnostics.

The active manuscript may acknowledge that regional isolation ranges differ and avoid
like-for-like causal pollinator claims. Explicit nonlinear inference belongs in a future
analysis/paper if it is developed prospectively.

## Reproducibility

- workflow run: **37288340871**
- artifact: **11335466621**
- script: `scripts/audit_pollination_syndrome_nonlinearity.py`
- summary: `results/pollination_syndrome_nonlinearity_20261005/syndrome_hinge_summary.csv`
