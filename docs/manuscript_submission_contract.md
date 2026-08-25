# Manuscript submission contract

## Status

This document defines the repository surface that may support the Chapter 1 manuscript. It is intentionally narrower than the full development history.

A result is manuscript-canonical only when it follows this contract, uses locked inputs, and produces an auditable output manifest.

The current frozen scientific checkpoint is [`chapter1_when_where_frozen_result_20260825.md`](chapter1_when_where_frozen_result_20260825.md).

## Canonical scientific question

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

The manuscript is a global comparative analysis of **boundary conditions of island filtering**. It does not identify a historical pollinator-loss event and does not establish causal floral evolution from cross-sectional island data.

## Canonical hypothesis structure

### H1 — universal filtering

Isolation produces a broadly shared floral/reproductive response vector across island floras.

### H2 — conditional / biogeographically structured filtering

Isolation-associated filtering is detectable under some biogeographic/floristic conditions, and the multivariate response vector differs among at least some contexts.

### H3 — floristic-status persistence

A supported regional response that remains within native non-endemic assemblages is not confined to endemic-lineage turnover.

## Primary analysis hierarchy

### 1. WHERE — within-context response-vector omnibus

For each biogeographic context and floristic-status stratum, supported atomic floral/reproductive categories are stacked in one grouped-binomial model with category-specific baseline and isolation effects. Spatial-block cluster-robust covariance retains dependence among categories and islands.

Formal null:

```text
H0: all category-specific isolation slopes in this context = 0
```

BH correction is applied across contexts within each status stratum and support tier.

### 2. BETWEEN-WHERE — pairwise response-vector difference

For each pair of contexts, only atomic responses meeting the same support threshold in both regions are retained.

Formal null:

```text
H0: isolation-response vector in context A = isolation-response vector in context B
```

A significant result in one context and a nonsignificant result in another is not evidence of a regional difference.

### 3. WHEN — floristic-status persistence

The WHERE test is repeated for:

- `all_native`
- `native_nonendemic`
- `endemic`

Persistence in `native_nonendemic` is interpreted as evidence that filtering is not confined to endemic taxa. Because status strata overlap, this is a persistence classification, not a formal causal difference among strata.

### 4. M3 lineage guardrail

Broad outcomes are compared against a genus-composition-preserving expectation. Manuscript-level interpretation must not treat inherited lineage composition as optional.

### 5. Atomic decomposition

The existing atomic M0–M4 layer is retained to describe **what creates a supported regional response vector**. It does not define WHERE by counting individually significant traits.

## Support rule

- `<30` islands per atomic response: excluded from declared vector tests;
- `30–49`: pilot;
- `>=50`: meets the count component of confirmatory support.

Every vector test reports the number of retained atomic responses, minimum per-response island support, total unique islands, and spatial-cluster count.

## Frozen canonical result

Canonical workflow run: `32837335384`.

### Confirmatory WHERE

Filtering is detected in both northern mid-latitude and tropical contexts.

| stratum | context | responses | islands | joint Wald | q |
| --- | --- | ---: | ---: | ---: | ---: |
| all native | northern mid-latitude | 21 | 240 | 243.36 | 1.76e-39 |
| all native | tropical | 17 | 136 | 170.24 | 2.52e-27 |
| native non-endemic | northern mid-latitude | 21 | 240 | 238.86 | 7.43e-39 |
| native non-endemic | tropical | 17 | 136 | 227.65 | 7.43e-39 |

Northern high-latitude and southern-extratropical contexts do not currently support equivalent confirmatory vector tests and are **unresolved**, not demonstrated null regions.

At the pilot threshold, southern-extratropical all-native and native-nonendemic vectors are nonzero but rely on only 5 and 4 atomic responses.

### Confirmatory BETWEEN-WHERE

Northern mid-latitude and tropical response vectors differ using 17 responses with confirmatory support in both contexts:

- all native: χ² = **69.78**, df = 17, q = **2.35e-08**;
- native non-endemic: χ² = **61.02**, df = 17, q = **7.13e-07**.

Thus confirmatory biogeographic heterogeneity is established **at the multivariate vector level** between northern mid-latitude and tropical floras.

### WHEN

The filtering signal persists in `native_nonendemic` flora in both confirmatorily supported contexts. Therefore the northern and tropical results are not confined to endemic taxa.

Endemic-only vector inference remains under-supported confirmatorily.

## Atomic-category interpretation rule

Atomic results are secondary decomposition after the vector-level result is frozen.

They may be used to identify which colours, forms, and reproductive states make the northern and tropical vectors differ. They may not replace the omnibus test by counting p-values.

The genus-preserving M3 layer does not recover a coherent confirmatory `generalized + plain + SC` syndrome in the main northern/tropical strata. Therefore a supported regional response vector is not automatically a named classical syndrome.

## Pollination-syndrome interpretation contract

Pollination-syndrome expectations enter only **after the when/where result is frozen**.

Allowed:

- compare the confirmed northern-versus-tropical response-vector difference with literature on Bombus / long-tongued bees, birds, Lepidoptera, hawkmoths, or other functional pollinator groups;
- describe partial concordance, mismatch, or unresolved mechanism;
- use pollinator evidence to motivate Chapter 2 tests.

Not allowed:

- infer actual pollinator identity from floral phenotype;
- claim that Bombus loss caused the northern vector;
- claim bird/Lepidoptera replacement caused the tropical vector;
- convert climatic compatibility or opportunistic nondetection into historical loss;
- use syndrome labels as fitted predictors in the Chapter 1 when/where analysis.

## Status of Bombus analyses

Existing Bombus applicability, environmental-niche, occurrence, and channel-score products remain provenance-preserving diagnostic/sensitivity products. They are not canonical primary Chapter 1 predictors.

Additional Bombus acquisition must not be launched solely to rescue or manufacture a Chapter 1 mechanism.

## Evidence-tier rule

The manuscript reports trait evidence resolution rather than hiding it behind complete fill.

- **Confirmatory trait evidence:** direct source-backed species evidence.
- **Pilot vector support:** 30–49 islands per retained atomic response.
- **Confirmatory vector support:** >=50 islands per retained atomic response.
- **Secondary robustness:** genus/family inference where declared.
- **Sensitivity only:** global fallback or other inference layers.

## Required safeguards before submission

### Spatial dependence

The vector tests use the declared spatial-block cluster-robust covariance. Alternative spatial structures remain sensitivity analyses where feasible.

### Lineage composition

The genus-composition-preserving M3 layer is required before broad syndrome claims. A globally complete source-region assignment is not required merely to duplicate this safeguard, although source-region assignments remain useful independent sensitivity analyses.

### Biogeographic-context definition

Context grouping must be frozen independently of floral/reproductive outcomes.

### Isolation functional form

The isolation response is predeclared as log distance to continent in the current canonical run. The manuscript must report the observed distribution and a prespecified leverage/functional-form sensitivity.

### Zero-distance islands

Zero or analogous distance classes must be audited and their treatment reported.

### Outcome support

Within a given omnibus test, support rules are applied before fitting. Pairwise context tests use only responses meeting the same support threshold in both contexts.

## Reproducibility and archival rule

Before submission, all manuscript-critical inputs and outputs must be archived durably with checksums.

A fresh reader should be able to identify:

- the manuscript release/tag;
- locked input artifacts and checksums;
- the canonical when/where workflow;
- software environment;
- support and evidence-tier definitions;
- attrition from the frozen island universe to each vector test; and
- files used for every manuscript figure/table.

## Chapter 1 -> Chapter 2 boundary

Chapter 1 establishes:

> **Isolation-associated floral/reproductive filtering is detectable in both northern mid-latitude and tropical island floras, persists within native non-endemic assemblages, and is expressed as significantly different multivariate response vectors in those two contexts.**

Chapter 1 does not explain why.

The unresolved question handed to Chapter 2 (`izu-core`) is:

> **Which ecological interaction states make northern mid-latitude and tropical island floras respond differently to isolation?**

Candidate explanations include pollinator identity, functional diversity, trait matching, effective service, functional replacement, reproductive assurance, network context, and non-pollination geography/history.

## Noncanonical material

The following remain development history or sensitivity material unless explicitly promoted by a later contract change:

- Bombus-primary M0–M4 pathway variants;
- Bombus mediation/coefficient-product analyses;
- old v1/v2 bridge analyses;
- superseded Bayesian/INLA variants built around a Bombus-primary confirmatory question;
- trait-source scouting and free-bulk-source pilots;
- old interpretations based only on counts of individually significant atomic slopes.
