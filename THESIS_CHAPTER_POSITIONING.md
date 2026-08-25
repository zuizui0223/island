# Thesis positioning — Chapter 1

## Role in the dissertation

This repository is the **Chapter 1 / macroecological WHEN–WHERE** component of the dissertation.

Chapter 1 asks:

> **When and where do floral/reproductive trait probabilities change along a mainland-distance/source-pool-accessibility gradient, and where do multivariate response vectors differ?**

The chapter establishes boundary conditions. It does not identify the ecological mechanism producing them.

## Data interpretation

GBIF records are not treated as a complete island-flora census or assumed to be a simple random sample. They are treated as an opportunistic, incompletely observed sample of realised floras.

The analysis therefore separates:

1. **observation process** — where flora and direct trait evidence are visible;
2. **trait-centric ecological surface** — `P(trait state | directly trait-resolved observed flora, geography)`;
3. **WHEN / WHERE response-vector inference** — where those trait probabilities jointly vary and where regional vectors differ.

An unrecorded island/species is never coded as biological trait absence.

Mainland distance is a composite geographic gradient that may represent both dispersal limitation and changing mainland/source-pool accessibility. Chapter 1 does not call it a pure causal isolation treatment.

## Hypotheses

### H1 — common geographic filtering

> Floral/reproductive trait probabilities change systematically along the geographic/source-pool gradient in multiple island contexts.

### H2 — context-dependent response vector

> The multivariate trait-probability response to the same geographic gradient differs among biogeographic contexts.

This is tested directly with a vector contrast; significance in one region and nonsignificance in another is not sufficient evidence of heterogeneity.

### H3 — floristic-status persistence

> If the geographic response is not confined to endemic-lineage turnover, it should persist within native non-endemic assemblages.

This is a persistence condition, not a causal contrast among overlapping status strata.

### H4 — observation-robustness condition

> A manuscript-level WHEN/WHERE result must persist when measured observation and trait-resolution structure is prevented from dominating inference.

Required robustness layers include equal/capped island information weighting, response-specific trait-resolution coverage adjustment, all-island observation modelling, alternative distance transformations on the same island universe, and leave-one-spatial-block influence analysis.

## Frozen observation-robust result

Canonical run: `32845980788`.

### WHERE

Confirmatory multivariate geographic filtering is supported in:

- **northern mid-latitude** island floras;
- **tropical** island floras.

The result occurs in both `all_native` and `native_nonendemic` assemblages.

### BETWEEN-WHERE

Northern-midlatitude and tropical response vectors differ directly in the confirmatory multivariate comparison.

Canonical species-count weighting:

- all native: q = **2.35e-08**;
- native non-endemic: q = **7.13e-07**.

### Observation robustness

The complete headline is reproduced:

- **10/10** times across five information-weight schemes × two status strata, including equal-island weighting;
- **2/2** status strata after response-specific direct-trait coverage adjustment;
- **6/6** across log1p, square-root and raw distance forms while retaining the full island universe;
- **84/84** leave-one-block runs for all-native and **84/84** for native-nonendemic assemblages.

The canonical checkpoint records:

```text
observation_robust_headline = true
```

This does not prove missingness is random. It shows that the headline is not readily explained by measured trait-resolution coverage, unequal information weighting, one distance transform, or one spatial block.

## Trait-centric interpretation

Atomic trait models are read as conditional probabilities, for example:

```text
P(SC | directly trait-resolved observed flora, geography)
```

Current direct evidence does **not** show a clear positive distance-gradient slope for SC in either the northern-midlatitude or tropical all-native sample. Therefore the Chapter 1 result is not a simple Baker-rule result.

Post-freeze decomposition indicates that floral architecture contributes strongly to the north-versus-tropical vector difference, while SI/SC contributes more weakly.

## Current geographic boundary

| context | confirmatory status |
| --- | --- |
| northern mid-latitude | supported |
| tropical | supported |
| southern extratropical | pilot signal; confirmatory unresolved |
| northern high-latitude | unresolved / not confirmatorily testable |

Under-support is not interpreted as a biological null.

## Next data acquisition

Additional acquisition must not be used to rescue the northern/tropical result.

The next goal is to expand regional testability. The current bottleneck is:

```text
flora recorded
→ floristic status resolved
→ direct trait evidence resolved
→ enough islands per atomic response
```

Southern extratropical has 317 flora-recorded islands but only 34 entering the current native-status trait surface. Northern high-latitude has 424 flora-recorded islands but only 12.

Therefore the outcome-blind acquisition order is:

1. floristic-status resolution on already recorded southern/high-latitude floras;
2. direct floral-form evidence;
3. direct SI/SC evidence;
4. colour only where it materially increases island-level testability.

See [`docs/chapter1_next_acquisition_priority_20260825.md`](docs/chapter1_next_acquisition_priority_20260825.md).

## Stronger source-pool extension

Where a candidate source pool can be declared independently of island trait outcomes, a secondary stronger test can ask:

```text
P(species occurs on island | trait, candidate source pool, geography)
```

This is distinct from the primary opportunistic-flora composition estimand and should be restricted to defensible source-pool systems.

## Pollination-syndrome boundary

Bombus, bird, butterfly, moth, hawkmoth and other pollinator labels do not enter Chapter 1 fitted models.

Allowed only after the WHEN/WHERE freeze:

```text
frozen regional response vectors
→ literature-defined concordance / mismatch
```

Not allowed:

```text
trait vector → inferred pollinator guild → causal mechanism
```

## Handoff to Chapter 2

Chapter 1 ends with:

> **Why does the same mainland-distance/source-pool-accessibility gradient produce detectable floral/reproductive filtering in both northern mid-latitude and tropical island floras, yet generate different multivariate response vectors?**

Chapter 2 (`izu-core`) distinguishes candidate mechanisms such as pollinator function, functional diversity, trait matching, effective service, reproductive assurance, functional replacement, network context and non-pollination geography/history.

## Relationship to Chapters 2 and 3

| Chapter 1 — `island` | Chapter 2 — `izu-core` | Chapter 3 — `shimahotarubukuro` |
| --- | --- | --- |
| global opportunistically observed island floras | mechanistic response architecture | focal lineage across Izu islands |
| asks **WHEN / WHERE** | asks **WHY / HOW** | asks **WHAT phenotype** |
| models observation and trait-probability surfaces | distinguishes candidate interaction mechanisms | measures multidimensional phenotypic divergence |
| output: supported/unresolved boundary conditions | output: branching / propagation / buffering | output: phenotype differentiation |

## Claim boundary

Chapter 1 must not claim that:

- no record means absence;
- mainland distance is a pure causal isolation treatment;
- southern/high-latitude floras are null merely because current confirmatory support is insufficient;
- SC necessarily increases with distance;
- response vectors identify pollinator guilds;
- cross-sectional assemblage filtering demonstrates within-lineage evolution.

The Chapter 1 contribution is to identify **where and under what floristic/observation conditions trait probabilities respond to a geographic/source-pool gradient, and which contexts show demonstrably different multivariate responses.**
