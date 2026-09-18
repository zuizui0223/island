# Chapter 1 H4 prospective temporal replication — frozen 2026-09-18

## Status

**Prospective protocol frozen before post-2015 pollen-limitation outcome extraction.**

The existing v13 H4 remains explicitly post-hoc. This protocol does not relabel it.
Instead, it creates a temporally non-overlapping validation cohort from publications
dated 2016-01-01 through 2026-09-18.

The source GloPL descriptor covers publications from 1981 through 2015. The validation
window therefore does not overlap the GloPL publication window.

## Co-primary predictions

### H4a — reproductive assurance

Species with frozen evidence for autonomous selfing are predicted to show **lower current
experimental pollen limitation**.

- predictor: binary autonomous selfing
- direction: negative
- one-sided alpha: 0.025

### H4b — floral accessibility/generalization

A frozen equal-weight score combines:

1. generalized floral form;
2. actinomorphic symmetry;
3. shallow/open tube.

The score is the arithmetic mean of available binary components when at least two of
three components are observed. Higher accessibility/generalization is predicted to show
**lower current experimental pollen limitation**.

- direction: negative
- one-sided alpha: 0.025
- shallow/open tube is retained despite its weaker v13 geographic recurrence; it cannot
  be dropped after seeing validation support or outcomes.

## Two-stage firewall

### Stage 1 — outcome-blind preflight

Only bibliographic and design metadata may be used:

- publication identity/date;
- DOI/source ID;
- species identity;
- study/site identity;
- coordinates/geographic context;
- standardized distance;
- frozen trait states.

The preflight code explicitly rejects GloPL effect fields, raw natural/supplemented
reproductive outputs, pollen-limitation values, and generic effect-size columns.

Support thresholds are evaluated without outcomes. The generated support lock contains
input SHA-256 digests, the admitted experiment-key digest, support counts, and a
`frozen_commit=REQUIRED_BEFORE_UNBLINDING` placeholder.

### Stage 2 — outcome analysis

Outcome analysis is impossible until the support lock has been committed and its
`frozen_commit` placeholder replaced by the commit SHA. The analysis then verifies
the frozen config, metadata and trait-state file digests and refuses outcome rows whose
experiment keys were not present in the locked cohort.

## Frozen support gates

H4a requires at least:

- 30 matched species total;
- 10 species in each autonomous-selfing state;
- 10 publications total;
- 5 publications represented in each state.

H4b requires at least:

- 30 matched species;
- 10 publications;
- accessibility-score SD >= 0.10;
- 8 species with score <= 1/3;
- 8 species with score >= 2/3.

Thresholds cannot be relaxed after support counts are known.

## Primary model

For each co-primary predictor:

`PL ~ predictor + standardized distance + geographic-context intercepts + measurement fixed effects`

The analysis unit is publication × site × species × measurement cell. Duplicate cells
are averaged. Each publication receives total analysis weight 1.0, and covariance is
clustered by publication.

The response keeps the GloPL direction: positive values mean stronger pollen limitation,
so the predicted trait coefficient is negative.

## Secondary analyses

The following cannot rescue a failed co-primary test:

- self-compatibility;
- atomic generalized form;
- atomic actinomorphy;
- atomic shallow/open tube;
- supplemental-only;
- no-zero-constant;
- within-publication;
- within-publication × site.

## Claim boundary

A successful temporal replication would upgrade the statement

> frozen trait states are functionally associated with current experimental pollen limitation

from post-hoc discovery alone to **post-hoc discovery followed by prospectively specified
temporal replication**.

It still would not establish that historical pollen limitation caused the global island
syndrome, statistically mediated isolation effects, reduced pollinator abundance, or
caused within-lineage floral evolution.
