# Bombus hypervolume: final position and validation

## Decision

The Bombus hypervolume is **not a direct Chapter 1 explanatory variable**. It is
retained as:

1. an exploratory climatic-opportunity diagnostic;
2. a preserved negative result from the pre-status and status/lineage analyses;
3. possible Discussion-level context for a future multivariate Northern trait
   vector; and
4. provenance for a later mechanism chapter if direct pollinator data become
   available.

It does not justify additional Bombus occurrence acquisition for Chapter 1.

## What the score measures

The locked run discovered 169 Bombus species, scored 142, and fitted nine
environmental dimensions after thinning 278,566 occurrence rows to 78,955. For
each species, the implementation winsorizes, standardizes, and fits a
ridge-regularized Mahalanobis ellipsoid. Compatibility equals 0.5 at the
empirical 95% training-distance boundary.

The correct name is **ellipsoidal climatic-compatibility score**. It is not a
complete Hutchinsonian niche hypervolume and does not measure occurrence,
absence, abundance, functional diversity, visitation, service, historical
loss, deficit, or replacement.

## Measurement limits

| Finding | Evidence | Consequence |
|---|---:|---|
| Northern row coverage | 2,164 of 2,167 islands | low simple missing-row risk |
| Candidate set | 169 total, normally the same 142 scored | maximum is not inflated by varying candidate counts here |
| Environmental extrapolation | median any-axis share `0.775` | predictions often extend beyond well-supported environments |
| Maximum provenance | max-supplying species extrapolation not retained | no defensible max-quality subset |
| Climate overlap | Spearman with climate PC1: max `-0.572`, mean `-0.743`, compatible share `-0.693` | little clearly separable channel information |

## Pre-status opportunity-gate result

The pre-status model asked whether the Northern isolation slope strengthened at
higher recipient climatic opportunity. None of the four focal maximum-score
interactions passed FDR at 0%, 10%, or 25% trait coverage. Flower-size
interactions were consistently opposite the declared direction. A
compatible-species-share sensitivity supported an opposite-direction
`very_small` interaction at 0% and 10% coverage but not 25%.

This rejected that proxy-specific opportunity-gate prediction. It did not test
the unobserved functional-deficit and replacement gates.

## Status/lineage residual result

PR #133 then asked whether compatibility added information after floristic
status and genus composition were represented:

```text
M0: distance + area + climate PC1-4
M1: M0 + climatic compatibility
M2: M1 + distance x compatibility
```

Generalized form and SC did not improve out of sample. Plain colour improved
spatial-CV RMSE by only about 1.1% for all native species and 1.7% for native
non-endemics, while the compatibility coefficient pointed opposite to the
simple `Bombus deficit -> plain colour` prediction. No distance-by-compatibility
interaction was nominally supported.

Thus the hypervolume did not recover the predicted Northern pathway after the
major structural alternatives were represented.

## Manuscript use

- **Methods/supplement:** describe the score and retained negative diagnostic.
- **Results:** report that the proxy did not add the predicted status/lineage
  residual pathway, if mechanistic diagnostics are shown.
- **Discussion:** at most, compare a supported multivariate Northern trait vector
  with published large/long-tongued pollinator-loss expectations.
- **Prohibited:** call low compatibility Bombus absence/deficit, call high
  compatibility realized availability, or use the proxy to explain a raw
  whole-flora flower-size gradient.

The direct Chapter 1 result is status- and lineage-aware biogeographic filtering.
Bombus identity, functional diversity, effective service, and replacement are
questions for Chapter 2.

## Reproducible evidence

- Preserved pre-status configuration:
  `config/chapter1_hypothesis_evaluation.yml`
- Pre-status gate results:
  `analysis/chapter1_channel_gate_evaluation_20260825/`
- Status/lineage residual checkpoint: PR #133 run `32613425082`, artifact
  `northern-bombus-lineage-residual-32613425082`
- Canonical contract: `config/canonical_analysis_mainline.yml`
