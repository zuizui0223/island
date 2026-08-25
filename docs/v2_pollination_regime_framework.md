# v2: component-specific island floral reorganization framework

## Status and purpose

This document defines the current Chapter 1 scientific scope.

The primary question is:

> **Does island isolation generate one coherent floral/reproductive syndrome, or does it reorganize particular floral components independently? Where are those component responses detectable, and is regional heterogeneity formally supported?**

The current frozen result is documented in [`chapter1_frozen_result_20260825.md`](chapter1_frozen_result_20260825.md).

The primary analysis does not require a global Bombus-deficit variable. Pollination-syndrome literature enters only after plant-trait results are frozen.

## Headline hypothesis contrast

### H1 — coherent island-syndrome hypothesis

Isolation produces a coordinated floral/reproductive syndrome whose broad components move together and remain after status and lineage safeguards.

### H2 — component-specific floral reorganization hypothesis

Isolation changes particular floral components without necessarily producing one coherent syndrome. Broad summaries can cancel, weaken or disappear when multistate categories and inherited lineage composition are preserved.

The current frozen evidence favors **H2 over H1**.

## Secondary boundary test — biogeographic contingency

Context remains important, but a regional difference is not inferred from significance patterns alone.

The formal rule is:

```text
M0 = baseline + context main effects
M1 = M0 + common isolation slope
M2 = M1 + isolation × context interactions
```

For each atomic category, M2 uses a cluster-robust joint Wald test of all interaction coefficients. Joint p-values are BH-corrected across atomic states within each `stratum × trait` family.

A statement such as “north is significant but tropics are not” is **descriptive only** unless the formal interaction test supports heterogeneity.

The current frozen run has one FDR-supported joint category (`endemic × yellow_orange`), but it uses only 32 northern-midlatitude versus 30 tropical islands, so it is pilot-level rather than confirmatory.

## Claims allowed by this design

Allowed:

```text
atomic isolation-associated trait effects
+ category preservation
+ status / lineage guardrail
= evidence for component-specific assemblage reorganization
```

Also allowed:

```text
joint interaction FDR + adequate context counts
= evidence for regional slope heterogeneity
```

Also allowed in Discussion:

```text
frozen trait-component direction
<-> literature-defined pollination-syndrome expectation
= ecological concordance or mismatch
```

Not allowed:

```text
significant in region A + nonsignificant in region B = regional difference
Bombus loss caused the northern trait pattern
bird/Lepidoptera replacement caused a tropical counter-pattern
floral syndrome = direct pollinator observation
```

## Analysis domains

### Reproductive composition

Keep reproductive states distinct wherever evidence permits:

- self-incompatible;
- self-compatible;
- autonomous selfing capacity;
- delayed selfing;
- cleistogamy;
- mating system;
- unresolved.

SC, autonomous selfing and realized selfing are not interchangeable.

### Floral phenotype

Primary outcomes remain category-preserving or continuous where possible:

- flower colour;
- floral form;
- symmetry;
- tube presence / depth;
- openness / accessibility;
- flower size;
- inflorescence display;
- nectar-guide or pattern evidence where directly supported.

Binary `plain/conspicuous`, `generalized/specialized`, and `SC/SI` contrasts are broad guardrail summaries, not substitutes for atomic outcomes.

## Canonical model ladder

### M0 — context-aware baseline

```text
area + climate
+ biogeographic-context main effects
+ declared status / observation structure
-> floral and reproductive composition
```

M0 absorbs regional mean-composition differences before isolation is evaluated.

### M1 — common isolation component effects

```text
M0 + isolation -> atomic trait composition
```

M1 asks which individual floral/reproductive components change with isolation after the context-aware baseline. These atomic effects are part of the primary Chapter 1 evidence.

### M2 — regional slope heterogeneity

```text
M1 + isolation × biogeographic context -> atomic trait composition
```

M2 is a **secondary boundary-condition test**. Its primary statistic is the joint cluster-robust Wald test of all interaction terms for the atomic category.

### M3 — status / lineage guardrail

Broad outcomes are tested against a genus-composition-preserving expectation.

The current confirmatory-support northern-midlatitude and tropical results do not retain FDR-supported isolation slopes for:

- `generalized_form`;
- `plain_colour`;
- `self_compatibility`

in the main all-native/native-nonendemic strata.

Therefore those broad summaries cannot be promoted into a demonstrated classical syndrome.

### M4 — category-preserving decomposition

Colour and floral-form multistates are expanded into atomic category presences. `SC|SI` is retained once as `mixed_or_variable` rather than double-counted.

The purpose is to identify which observed categories increase, decrease, oppose or cancel broad summaries.

Current repeated northern-midlatitude signals include positive `bell_campanulate`, `composite_head`, `funnel_trumpet`, and `spurred` states, and negative `other_described`, `papilionaceous`, `salverform`, and `red_pink` states.

## Current result boundary

The canonical joint-Wald run (`32833362756`) yields:

- 43 fitted atomic categories;
- 95 context-specific slope rows;
- 17 FDR-supported within-context slopes;
- 16 slopes meeting the count component of confirmatory support, all in northern-midlatitude islands;
- 1 FDR-supported joint contingency category, based only on pilot-level 32-versus-30 island support;
- no confirmatory joint contingency category.

The correct headline is therefore:

> **Isolation-associated floral change is component-specific rather than one coherent classic syndrome; strong northern-midlatitude signals are real within that context, but confirmatory regional slope heterogeneity is not yet established.**

## Pollination-syndrome concordance layer

Pollination syndromes enter **after** M0–M4 results are fixed.

### Northern-midlatitude interpretation

The frozen component vector may be compared with Bombus / large- or long-tongued-bee literature for partial concordance or mismatch.

Because the broad `generalized + plain + SC` signature is not recovered after M3, Chapter 1 must not call the current result a demonstrated Bombus-loss syndrome.

### Tropical / southern interpretation

Bird-, butterfly-, moth-, hawkmoth- or alternative-bee expectations are considered only if supported tropical/southern component vectors exist. The absence of northern-like significance is not itself evidence for an alternative syndrome.

### Symmetry rule

```text
frozen trait pattern -> syndrome concordance or mismatch -> candidate mechanism
```

not:

```text
trait pattern -> inferred pollinator -> causal claim
```

## Status of Bombus data products

Existing Bombus products are retained for source-region provenance, environmental-compatibility diagnostics, effort-aware occurrence diagnostics, sensitivity analyses and candidate-system selection.

They are not prerequisites for Chapter 1 and cannot turn climate compatibility into occurrence, opportunistic non-detection into historical loss, or occurrence into effective pollination.

No additional Bombus acquisition should be launched solely to rescue a preferred Chapter 1 mechanism.

## Analysis sequence

1. Freeze island universe, trait evidence tiers, context definition, support thresholds and baseline covariates.
2. Build polygon-exact flora and direct-evidence trait products.
3. Produce attrition and support audits.
4. Fit M0 context-aware baselines.
5. Fit M1 common isolation effects for atomic categories.
6. Fit M2 interaction models and joint Wald tests.
7. Apply M3 genus-preserving broad-outcome guardrails.
8. Apply M4 category-preserving decomposition.
9. Freeze the atomic trait vectors and formal contingency classifications.
10. Only then compare supported patterns with pollination-syndrome literature.
11. Hand unresolved component-specific mechanisms to Chapter 2.

## Chapter 1 -> Chapter 2 handoff

Chapter 1 ends with:

> **Why does isolation reorganize some floral components but not produce one coherent syndrome, and which ecological interaction states determine which components respond?**

Candidate mechanisms include pollinator identity, large/long-tongued pollinator function, pollinator functional diversity, trait matching, effective pollination service, functional replacement, reproductive assurance, network context and non-pollination geography/history.

Those alternatives belong in Chapter 2 (`izu-core`), where they can be separated mechanistically.

## Final positioning

The Chapter 1 contribution is not a global Bombus-loss test and not a confirmed biogeographic-contingency result. It is a rigorous demonstration that **isolation-associated floral change is better represented as component-specific reorganization than as one coherent floral/reproductive island syndrome**, plus a formal framework for determining when regional heterogeneity is actually supported.
