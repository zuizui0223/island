# v2: when / where island floral filtering framework

## Status and purpose

This document defines the current Chapter 1 scientific scope.

The primary question is:

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

The canonical frozen result is [`chapter1_when_where_frozen_result_20260825.md`](chapter1_when_where_frozen_result_20260825.md).

Bombus and other pollinator identities do not enter the primary fitted design. Pollination-syndrome literature is a downstream interpretation layer only.

## Hypothesis structure

### H1 — universal filtering

Isolation produces a broadly shared floral/reproductive response vector across island floras.

### H2 — conditional / biogeographically structured filtering

Isolation-associated filtering is detectable under some biogeographic/floristic conditions and can be expressed as different multivariate response vectors among contexts.

### H3 — floristic-status persistence

A regional filtering response that persists in native non-endemic flora is not confined to endemic-lineage turnover.

## Primary inference

### WHERE

Within each biogeographic context and floristic-status stratum, all atomic responses meeting the support threshold are modeled jointly. Each response has its own baseline and isolation slope; spatial-block cluster-robust covariance retains dependence among categories and islands.

```text
H0: all supported category-specific isolation slopes = 0
```

This is the headline test of whether isolation-associated floral/reproductive filtering is detectable in a context.

### BETWEEN-WHERE

For two contexts, only atomic responses meeting the same support threshold in both are retained.

```text
H0: isolation-response vector in context A = isolation-response vector in context B
```

A significant result in one region and a nonsignificant result in another never counts as evidence of regional heterogeneity by itself.

### WHEN

WHERE is repeated in `all_native`, `native_nonendemic`, and `endemic` strata. Persistence in `native_nonendemic` is used as a boundary condition indicating that the signal is not confined to endemic taxa. Because strata overlap, this is not treated as a causal status contrast.

## Support tiers

- fewer than 30 islands per retained atomic response: excluded from the declared vector test;
- 30–49: pilot;
- 50 or more: confirmatory count support.

Pairwise regional comparisons require the same threshold in both contexts.

## Current frozen result

Canonical workflow run: `32837335384`.

### Confirmatory WHERE

Filtering is detected in both:

- **northern mid-latitude** island floras;
- **tropical** island floras.

It persists in `native_nonendemic` assemblages in both contexts.

The confirmatory tests are:

| stratum | context | atomic responses | unique islands | joint Wald χ² | q |
| --- | --- | ---: | ---: | ---: | ---: |
| all native | northern mid-latitude | 21 | 240 | 243.36 | 1.76e-39 |
| all native | tropical | 17 | 136 | 170.24 | 2.52e-27 |
| native non-endemic | northern mid-latitude | 21 | 240 | 238.86 | 7.43e-39 |
| native non-endemic | tropical | 17 | 136 | 227.65 | 7.43e-39 |

Northern high-latitude and southern-extratropical contexts do not currently support equivalent confirmatory vector tests. They are **unresolved**, not demonstrated null regions.

At the pilot threshold, southern-extratropical all-native and native-nonendemic vectors are nonzero, but only 5 and 4 atomic responses contribute respectively.

### Confirmatory BETWEEN-WHERE

Northern mid-latitude and tropical vectors differ using 17 responses that meet confirmatory support in both contexts:

- all native: χ² = **69.78**, df = 17, q = **2.35e-08**;
- native non-endemic: χ² = **61.02**, df = 17, q = **7.13e-07**.

Thus the current evidence supports **biogeographically structured isolation filtering between northern mid-latitude and tropical island floras**.

## Atomic decomposition is secondary

The earlier atomic M0–M4 layer remains useful to answer:

> **What trait components make the supported northern and tropical response vectors different?**

It is not the primary when/where test. Counts of individually significant traits cannot define WHERE.

The genus-composition-preserving M3 layer also shows that broad `generalized_form`, `plain_colour`, and `self_compatibility` outcomes do not combine into one coherent confirmatory classic syndrome in the main northern/tropical strata. Therefore the two supported response vectors should not automatically be named as classical pollination syndromes.

## Pollination-syndrome interpretation

Only after the when/where result is frozen may the northern-versus-tropical vector difference be compared with literature on Bombus / long-tongued bees, birds, Lepidoptera, hawkmoths, or other pollinator functional groups.

Allowed logic:

```text
frozen response vector
<-> literature-defined ecological expectation
= concordance / mismatch / unresolved mechanism
```

Not allowed:

```text
trait vector -> inferred pollinator guild -> causal historical mechanism
```

The same evidential standard applies to every region.

## Status of Bombus products

Existing Bombus applicability, climatic-compatibility and occurrence diagnostics remain provenance-preserving descriptive/sensitivity products.

Retained semantic boundaries:

- environmental compatibility != occurrence;
- occurrence != visitation;
- visitation != effective service;
- effective service != reproductive dependency;
- nondetection != historical loss;
- floral phenotype != direct pollinator observation.

Additional Bombus acquisition must not be launched solely to rescue a preferred Chapter 1 mechanism.

## Analysis sequence

1. Freeze island universe, context definition, trait evidence, status strata, support thresholds and baseline covariates.
2. Build direct-evidence status-aware atomic trait composition.
3. Run **WHERE** response-vector omnibus tests.
4. Run **BETWEEN-WHERE** pairwise response-vector tests.
5. Classify **WHEN** through floristic-status persistence.
6. Apply the M3 genus-composition lineage guardrail.
7. Use atomic M0–M4 decomposition to characterize the supported vectors.
8. Freeze the when/where result.
9. Only then perform pollination-syndrome concordance/mismatch discussion.
10. Hand mechanism to Chapter 2.

## Chapter 1 -> Chapter 2 handoff

Chapter 1 ends with:

> **Why is isolation-associated floral/reproductive filtering detectable in both northern mid-latitude and tropical island floras, yet expressed as significantly different multivariate response vectors?**

Chapter 2 (`izu-core`) distinguishes candidate mechanisms such as pollinator identity, functional diversity, trait matching, effective service, reproductive assurance, replacement, network context and geography/history.

## Final positioning

Chapter 1 is a **when/where boundary-condition test**. Its contribution is to identify where isolation-associated filtering can be detected, under which floristic-status condition it persists, and which regional response vectors demonstrably differ—without using weak pollinator proxies to answer why.
