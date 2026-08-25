# Chapter 1 when / where frozen result — 2026-08-25

## Canonical run

This checkpoint freezes the Chapter 1 result after replacing individual-category significance counting with cluster-robust response-vector omnibus tests.

- workflow run: `32837335384`
- artifact: `chapter1-when-where-main-32837335384`
- artifact id: `9559224512`
- artifact digest: `sha256:dd8caf37fe830595083c643f871145295a6487883fe065a20a4209f159dd745c`
- direct-trait input run: `32702160934`
- floristic-status input run: `32559322028`
- isolation/context input run: `29228212586`
- M3 genus-fixed null draws: `1000`

No pollinator variable enters the fitted when/where design.

## Primary question

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

The response is the vector of supported atomic floral/reproductive categories, not a preconstructed island-syndrome score.

## Tests

### WHERE — within-context vector test

For each biogeographic context and floristic-status stratum, supported atomic categories are stacked in one grouped-binomial model. Each category receives its own intercept, baseline covariate effects and isolation slope. A spatial-block cluster-robust covariance retains dependence among categories and islands within blocks.

The formal null is:

```text
H0: all category-specific isolation slopes in this context = 0
```

BH correction is applied across contexts within each floristic-status stratum and support tier.

### BETWEEN-WHERE — response-vector difference

For each pair of contexts, the analysis uses only atomic categories meeting the same support threshold in both contexts and jointly tests all category-specific `isolation × context` contrasts.

The formal null is:

```text
H0: isolation-response vector in context A = isolation-response vector in context B
```

A significant effect in one region and a nonsignificant effect in another is not treated as evidence of a regional difference.

### WHEN — floristic-status persistence

The same WHERE test is repeated for:

- `all_native`
- `native_nonendemic`
- `endemic`

Persistence in `native_nonendemic` flora is used as evidence that the detected response is not confined to endemic-lineage turnover. Because these strata overlap, this is a persistence classification, not a formal causal contrast among status classes.

## Confirmatory WHERE result

The confirmatory support threshold requires each retained atomic response to be represented on at least 50 islands in the focal context.

### All native flora

- **northern mid-latitude:** supported
  - 21 atomic responses
  - minimum response support = 55 islands
  - 240 unique islands
  - 36 spatial clusters
  - joint Wald χ² = **243.36**, df = 21
  - p = **8.80e-40**
  - BH q = **1.76e-39**
- **tropical:** supported
  - 17 atomic responses
  - minimum response support = 50 islands
  - 136 unique islands
  - 46 spatial clusters
  - joint Wald χ² = **170.24**, df = 17
  - p = q = **2.52e-27**

### Native non-endemic flora

- **northern mid-latitude:** supported
  - 21 atomic responses
  - minimum response support = 55 islands
  - 240 unique islands
  - 36 spatial clusters
  - joint Wald χ² = **238.86**, df = 21
  - p = **7.03e-39**
  - BH q = **7.43e-39**
- **tropical:** supported
  - 17 atomic responses
  - minimum response support = 50 islands
  - 136 unique islands
  - 46 spatial clusters
  - joint Wald χ² = **227.65**, df = 17
  - p = q = **7.43e-39**

Therefore isolation-associated floral/reproductive filtering is **confirmatorily detectable in both northern mid-latitude and tropical island floras** under the current direct-evidence support.

## Confirmatory BETWEEN-WHERE result

Only the northern-mid-latitude versus tropical comparison has enough common atomic response support for a confirmatory vector comparison.

### All native flora

- 17 common atomic responses
- northern minimum support = 152 islands
- tropical minimum support = 50 islands
- 376 unique islands
- 82 spatial clusters
- joint Wald χ² = **69.78**, df = 17
- p = q = **2.35e-08**

### Native non-endemic flora

- 17 common atomic responses
- northern minimum support = 152 islands
- tropical minimum support = 50 islands
- 376 unique islands
- 82 spatial clusters
- joint Wald χ² = **61.02**, df = 17
- p = q = **7.13e-07**

Thus the current evidence supports not only filtering in both regions, but also a **confirmatory difference in the multivariate isolation-response vector between northern mid-latitude and tropical island floras**.

## WHEN result

The confirmatory status-persistence classification is:

| context | all native | native non-endemic | endemic | current classification |
| --- | --- | --- | --- | --- |
| northern mid-latitude | supported | supported | not confirmatorily testable/supported | **persists in native non-endemic** |
| tropical | supported | supported | not confirmatorily testable/supported | **persists in native non-endemic** |
| northern high-latitude | no confirmatory vector test | no confirmatory vector test | no confirmatory vector test | **unresolved for current data** |
| southern extratropical | no confirmatory vector test | no confirmatory vector test | no confirmatory vector test | **unresolved confirmatorily** |

The northern-mid-latitude and tropical results therefore cannot be explained simply as patterns confined to endemic taxa.

## Pilot boundary

At the 30-island pilot threshold, southern-extratropical all-native and native-nonendemic vectors are also nonzero, based on only 5 and 4 atomic responses respectively. Pilot pairwise comparisons suggest differences from northern mid-latitude and, less consistently, tropical vectors.

These southern results are **hypothesis-generating only**. Northern high-latitude data do not currently provide enough atomic response support even for the declared pilot vector test.

## Current Chapter 1 conclusion

> **Isolation-associated floral/reproductive filtering is detectable in both northern mid-latitude and tropical island floras, persists within native non-endemic assemblages, and differs in multivariate response direction/magnitude between those two biogeographic contexts. Current data are insufficient for confirmatory conclusions in northern high-latitude or southern-extratropical floras.**

This is a when/where conclusion. It does not identify which pollinator or ecological mechanism produced either response vector.

## Role of the atomic category results

Atomic M0–M4 results are retained only to characterize **what makes the supported regional vectors different** after the when/where result is established. They do not replace the vector-level test.

The broad `generalized + plain + SC` summary still does not survive the genus-composition-preserving M3 guardrail as a coherent classic syndrome. Therefore the regional difference should be described as a difference in **floral/reproductive response vectors**, not automatically as two named pollination syndromes.

## Pollinator interpretation boundary

Bombus / long-tongued-bee and bird / Lepidoptera syndrome literature may now be used only to interpret the already-frozen northern-versus-tropical vector difference. It remains invalid to infer pollinator identity or historical pollinator replacement from the trait vectors themselves.

## Chapter 1 -> Chapter 2 handoff

Chapter 1 now ends with:

> **Why does isolation produce detectable floral/reproductive filtering in both northern mid-latitude and tropical island floras, yet generate significantly different multivariate response vectors in those two contexts?**

That mechanism question belongs to Chapter 2 (`izu-core`).
