# Thesis positioning — Chapter 1

## Role in the dissertation

This repository is the **Chapter 1 / macroecological when-and-where** component of the dissertation.

The shared dissertation-level question is:

> **How does geographic isolation alter plant reproduction through changes in ecological interactions, and why do those changes produce different floral outcomes across islands and lineages?**

Chapter 1 asks only the first-order boundary question:

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

The primary response is a multivariate vector of directly supported floral/reproductive categories. Atomic categories are used to describe the vector after the when/where result is established; pollinator identity is not used to define the result.

See [`docs/chapter1_when_where_frozen_result_20260825.md`](docs/chapter1_when_where_frozen_result_20260825.md) for the current frozen evidence.

## Chapter 1 hypotheses

### H1 — universal filtering

> **Isolation produces a broadly shared floral/reproductive filtering response across island floras.**

A universal model predicts a detectable response vector but little supported heterogeneity among biogeographic contexts.

### H2 — conditional / biogeographically structured filtering

> **Isolation-associated filtering is detectable under some biogeographic/floristic conditions, and the multivariate isolation-response vector can differ among contexts.**

This is tested in two distinct steps:

1. **WHERE:** within each context, is the response vector jointly different from zero?
2. **BETWEEN-WHERE:** are response vectors jointly different between contexts?

A significant result in one context and a nonsignificant result in another is not evidence of a regional difference.

### H3 — floristic-status persistence

> **If filtering is not merely an endemic-lineage turnover effect, the regional response should persist within native non-endemic assemblages.**

Because `all_native`, `native_nonendemic`, and `endemic` strata overlap, this is interpreted as a persistence condition rather than a causal contrast among status classes.

## Canonical when/where tests

The primary implementation is `src/island_v2/chapter1_when_where_omnibus.py` with the frozen contract in `config/chapter1_when_where_omnibus.yml`.

### WHERE

Within each context and status stratum, supported atomic categories are stacked in one grouped-binomial model with category-specific baseline and isolation effects. Spatial-block cluster-robust covariance retains dependence among categories and islands.

The formal null is:

```text
H0: all category-specific isolation slopes in this context = 0
```

### BETWEEN-WHERE

For a pair of contexts, only categories meeting the same support threshold in both contexts are retained.

The formal null is:

```text
H0: isolation-response vector in context A = isolation-response vector in context B
```

### Support tiers

- `<30` islands per atomic response: not used for the declared vector test;
- `30–49`: pilot;
- `>=50`: confirmatory count threshold.

## Frozen result

Canonical run: `32837335384`.

### WHERE — confirmatory

Isolation-associated floral/reproductive filtering is detected in:

- **northern mid-latitude** island floras;
- **tropical** island floras.

It is detected in both `all_native` and `native_nonendemic` strata.

Key confirmatory vector tests:

| stratum | context | atomic responses | unique islands | joint Wald | q |
| --- | --- | ---: | ---: | ---: | ---: |
| all native | northern mid-latitude | 21 | 240 | 243.36 | 1.76e-39 |
| all native | tropical | 17 | 136 | 170.24 | 2.52e-27 |
| native non-endemic | northern mid-latitude | 21 | 240 | 238.86 | 7.43e-39 |
| native non-endemic | tropical | 17 | 136 | 227.65 | 7.43e-39 |

Northern high-latitude and southern-extratropical floras do not currently have enough atomic response support for an equivalent **confirmatory** vector test. They are unresolved, not demonstrated null regions.

At the 30-island pilot threshold, southern-extratropical all-native and native-nonendemic vectors are nonzero, but only 5 and 4 atomic responses contribute respectively. These results remain hypothesis-generating.

### BETWEEN-WHERE — confirmatory

Northern mid-latitude and tropical response vectors differ significantly using the same 17 confirmatorily supported atomic responses:

| stratum | comparison | unique islands | joint Wald | q |
| --- | --- | ---: | ---: | ---: |
| all native | northern mid-latitude vs tropical | 376 | 69.78 | 2.35e-08 |
| native non-endemic | northern mid-latitude vs tropical | 376 | 61.02 | 7.13e-07 |

Therefore the current Chapter 1 evidence supports **biogeographically structured filtering between northern mid-latitude and tropical island floras**.

## WHEN — current boundary condition

The key persistence result is:

```text
northern mid-latitude:
  all native             supported
  native non-endemic     supported

 tropical:
  all native             supported
  native non-endemic     supported
```

Thus the detectable filtering in these two contexts is **not confined to endemic taxa**.

The endemic-only contribution remains under-supported for a confirmatory vector test, so Chapter 1 should not claim that endemicity is either necessary or irrelevant.

## What atomic categories are for

Atomic M0–M4 results answer a secondary question:

> **What trait components make the northern and tropical response vectors different?**

They do not define WHERE by counting significant traits.

The genus-composition-preserving M3 guardrail also shows that broad `generalized + plain + SC` summaries do not form one coherent classic island syndrome. This affects interpretation of the response vector, but it is not the Chapter 1 headline question.

## Pollination-syndrome boundary

Bombus-, bird-, butterfly-, moth-, hawkmoth- and other pollination-syndrome expectations enter only after the when/where result is frozen.

Allowed:

```text
frozen northern vs tropical response-vector difference
-> compare with literature-defined ecological expectations
-> concordance / mismatch
```

Not allowed:

```text
trait vector -> inferred pollinator guild -> causal mechanism
```

Existing Bombus products remain provenance/diagnostic/sensitivity products and are not primary Chapter 1 predictors.

## Handoff to Chapter 2

Chapter 1 now ends with a specific unresolved question:

> **Why is isolation-associated floral/reproductive filtering detectable in both northern mid-latitude and tropical island floras, yet expressed as significantly different multivariate response vectors?**

Chapter 2 (`izu-core`) is where pollinator identity, functional diversity, trait matching, effective service, reproductive assurance, functional replacement, network context and geography/history can be distinguished mechanistically.

## Relationship to Chapters 2 and 3

| Chapter 1 — `island` | Chapter 2 — `izu-core` | Chapter 3 — `shimahotarubukuro` |
| --- | --- | --- |
| global island-flora scale | mechanistic response architecture | one focal lineage across five Izu islands |
| asks **when and where filtering appears** | asks **how and why responses differ** | asks **what phenotype actually diverges** |
| tests within-context vectors and between-context vector differences | distinguishes candidate mechanisms | directly measures multidimensional floral divergence |
| output: supported and unresolved boundary conditions | output: branching / propagation / buffering / mechanism evidence | output: site-corrected phenotype differentiation and Pst |
| pollinators are not primary predictors | pollinator function is mechanistic | pollinator causation is not inferred from Pst |

The thesis-level handoff is:

```text
Chapter 1: WHEN / WHERE is filtering detectable, and where do response vectors differ?
                         |
                         v
Chapter 2: WHY do those contexts generate different response architectures?
                         |
                         v
Chapter 3: WHAT phenotype axes actually diverge within one focal lineage?
```

## Claim boundary

Chapter 1 should not claim that:

- northern high-latitude or southern-extratropical floras show no isolation filtering merely because confirmatory support is insufficient;
- a one-region significant effect proves regional heterogeneity;
- northern or tropical response vectors identify a causal pollinator guild;
- Bombus absence is equivalent to service loss;
- bird- or butterfly-like traits prove functional replacement;
- a broad pollination syndrome is directly observed; or
- cross-sectional assemblage filtering demonstrates within-lineage evolution.

The Chapter 1 contribution is: **to identify where isolation-associated floral/reproductive filtering is detectable, the floristic conditions under which it persists, and which biogeographic contexts show demonstrably different multivariate responses to isolation.**
