# Chapter 1 H5 — independent pollination-channel evidence chain

## Priority

This is the current highest-priority extension of the Chapter 1 double geographic filter.

Chapter 1 already supports the **plant-side** filter: increasing source separation is associated
with biogeographically contingent floral/reproductive assemblage responses, and the broad
Palearctic response is compatible with source-matched genus-level assembly. H5 asks whether the
**pollinator-side** filter can be measured independently rather than inferred from floral phenotype.

## Frozen chain

```text
regional source pollination channel
        ↓
source channel available / structurally absent
        ↓
island retained / disrupted
        ↓
realized visitation under explicit effort
        ↓
single-visit effectiveness (SVD, with no-visit control)
        ↓
rate-weighted effective service
        ↓
post-H2–H4 residual plant response
```

No link is allowed to substitute for the next one.

- source distribution is not island retention;
- occurrence is not visitation;
- visitation is not per-visit effectiveness;
- per-visit effectiveness is not effective service;
- effective service is not reproductive dependency;
- floral syndrome is not pollinator identity;
- climate compatibility is not realized loss.

## Channel states

A channel has three conceptually different outcomes.

### Retained

The channel is independently supported in the source region and independently supported as retained
on the island. Retention alone does not imply use of the focal plant, so visitation remains a separate
measurement.

### Disrupted

The channel is independently supported in the source region but is independently classified as lost
or deficient on the island under a prespecified evidence rule. Adequate non-detection may contribute
to this state; climatic compatibility alone may not.

### Structurally absent

The channel is not part of the source system. This is a context control, not a channel-loss state.
It must never be coded as disruption.

## Visitation and zeroes

Visitation is measured as observed visit bouts per monitored flower-hour. A zero is interpretable only
when usable monitoring effort is positive.

An adequately observed zero-visit unit may identify **zero realized service during that observation
window**, but it does not identify zero per-visit effectiveness. The SVD link is therefore labelled
`not_applicable_zero_visitation` rather than assigned a zero.

## Single-visit effectiveness

The default common metric is conspecific pollen deposited on a previously unvisited receptive stigma
after one observed visit. A no-visit/background control is required before the metric can enter the
formal chain.

```text
background_adjusted_SVD
  = mean_single_visit_conspecific_pollen
  - mean_no_visit_conspecific_pollen
```

This mirrors the direct effective-service measurement logic already frozen in `izu-core` without
copying its field outcomes into Chapter 1.

## Effective service

For positive visitation,

```text
effective pollen delivery per flower-hour
  = visit bouts per flower-hour × background-adjusted SVD
```

If visitation is adequately observed at zero, realized service is zero for that observation window
while per-visit effectiveness remains unmeasured/not applicable.

If positive visits occur but the no-visit control or SVD channel is missing, effective service is
withheld rather than approximated from visit rate.

## Contrast gate

A channel-side H5 contrast is marked ready only when the **same channel** has:

- at least one source-available retained unit with estimable effective service;
- at least one source-available disrupted unit with estimable effective service;
- at least one positive-visit unit in which the full visitation → SVD → service chain is observed.

Structurally absent units never satisfy the retained/disrupted loss contrast.

This gate closes only the pollinator side. Full H5 promotion still requires joining the independent
channel panel to a plant response that is residual **after H2–H4 source, lineage, area, climate,
observation and spatial safeguards**.

## Current evidence audit

### `island`

Existing Bombus infrastructure already contains pieces of the first two links: source applicability,
environmental compatibility, occurrence evidence and adequate non-detection logic. These are useful
inputs but do not by themselves pass the new H5 chain because the current global Chapter 1 analysis
does not contain linked focal-plant visitation, controlled single-visit effectiveness and effective
service.

In particular, the historical continuous Bombus availability/deficit summary remains a secondary
component score. It is **not** relabelled as H5 effective service.

### `izu-core`

The direct field contract already distinguishes monitored flower-hours, visitor identity/contact,
background-controlled SVD and rate-weighted effective pollen service, and extends further to open,
bagged-autonomous and supplemental-outcross reproduction. Its current readiness state is
`implementation_ready_field_data_missing` for the linked prospective panel.

Thus the methodological bridge is now explicit, but empirical H5 remains open until admitted linked
rows exist or an independent external system supplies the same chain.

## What this closes

The new gate closes a methodological ambiguity that previously left H5 as a prose-only requirement:
there is now one machine-auditable definition of what counts as independent pollinator-side evidence.

It does **not** yet close the ecological mechanism. The next empirical closure is:

1. populate source-channel states outcome-blind;
2. populate retained/disrupted evidence independently;
3. add effort-aware visitation;
4. add controlled SVD where visits occur;
5. compute effective service;
6. join to post-H2–H4 plant residuals and test incremental explanatory value.
