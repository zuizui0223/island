# v2: component-specific island floral reorganization — analysis architecture

## Status and boundary

`v1-freeze` preserves the exploratory draft and its first complete results. v2 is rebuilt from new data products and new scripts; it must not import v1 scripts, derived trait labels, or v1 model outputs as analysis inputs.

The current scientific scope is defined in [`v2_pollination_regime_framework.md`](v2_pollination_regime_framework.md), and the current frozen result is in [`chapter1_frozen_result_20260825.md`](chapter1_frozen_result_20260825.md).

The primary Chapter 1 question is:

> **Does island isolation generate one coherent floral/reproductive syndrome, or does it reorganize particular floral components independently?**

The architecture therefore prioritizes direct floral/reproductive composition, category preservation, floristic status, lineage structure and geography. Biogeographic context is retained as a formal boundary-condition test, not assumed to modify the isolation effect.

```text
island geography / climate / context mean composition
        ↓
floristic status / lineage composition / observation support
        ↓
common isolation-associated atomic component effects
        ↓
formal isolation × context joint test
        ↓
category-preserving decomposition + genus-fixed broad guardrail
        ↓
frozen component vector
        ↓
pollination-syndrome concordance or mismatch in Discussion
        ↓
unresolved mechanism handed to Chapter 2
```

## What stays from earlier v2

1. Island geography is the starting context.
2. Source-pool, floristic status, lineage composition, establishment/reachability and observation process must be considered before ecological interpretation.
3. Floral signal, architecture and reproductive traits remain multistate or continuous wherever possible.
4. Binary summaries such as plain/conspicuous, generalized/specialized and SC/SI are secondary.
5. Wind pollination, self-compatibility, autonomous selfing and pollinator independence remain separate domains.
6. Existing Bombus applicability, environmental and occurrence diagnostics remain provenance-preserving data products.
7. The global analysis does not claim within-lineage evolution or historical pollinator causation from cross-sectional composition.

## What changes

1. Bombus deficit remains outside the primary Chapter 1 explanatory model.
2. **Component-specific floral reorganization** replaces biogeographic contingency as the headline Chapter 1 hypothesis.
3. Context main effects are represented before isolation; regional heterogeneity is tested only by formal interaction statistics.
4. M1 common isolation effects and M4 atomic decomposition are primary evidence; M2 is a secondary boundary-condition test.
5. M3 genus-preserving broad outcomes test whether an apparent classical syndrome survives inherited lineage composition.
6. Pollination-syndrome interpretations are allowed only after the trait vector and heterogeneity classification are frozen.
7. The mechanistic question — why particular phenotype components respond — is handed to Chapter 2 (`izu-core`).

## Interpretation guardrail borrowed from Campanula/Izu

For a trait state `z`, the Campanula programme writes:

```text
W(z) = F(z) × E(z)
```

For Chapter 1 this remains an interpretation guardrail, not a fitted decomposition.

- `F(z)` — local reproductive contribution;
- `E(z)` — arrival, establishment, persistence and reachability conditional on viable reproduction;
- `W(z)` — observed island-flora membership or trait composition.

The global dataset mainly observes `W`. Therefore an isolation-associated component effect may arise from establishment filtering, lineage turnover, ecological interaction filtering, reproductive processes, or observation structure. Chapter 1 identifies which components change and which broad syndrome claims survive safeguards without pretending to identify `F` from `W` alone.

## Canonical model architecture

### M0 — context-aware baseline

```text
area + climate + biogeographic-context main effects
+ declared status / observation structure
-> floral and reproductive composition
```

M0 absorbs baseline and regional mean differences before isolation is added.

### M1 — common isolation component effects

```text
M0 + isolation -> atomic trait composition
```

M1 is part of the primary Chapter 1 evidence. It tests which individual floral/reproductive states vary with isolation after the context-aware baseline.

### M2 — biogeographic slope heterogeneity

```text
M1 + isolation × biogeographic context -> atomic trait composition
```

M2 is a secondary boundary-condition test. A regional difference requires a supported cluster-robust joint Wald test of all interaction terms; `significant in A / nonsignificant in B` is not sufficient.

### M3 — status / lineage guardrail

Broad outcomes are compared with genus-composition-preserving expectations.

```text
observed broad trait composition
-
genus-preserving expectation
=
residual broad response
```

The current confirmatory-support northern-midlatitude and tropical M3 results do not retain FDR-supported isolation slopes for `generalized_form`, `plain_colour`, or `self_compatibility` in the main all-native/native-nonendemic strata.

### M4 — category-preserving decomposition

Analyse retained atomic categories rather than forcing one syndrome score.

Current repeated northern-midlatitude component signals include:

- `bell_campanulate` ↑
- `composite_head` ↑
- `funnel_trumpet` ↑
- `spurred` ↑
- `other_described` ↓
- `papilionaceous` ↓
- `salverform` ↓
- `red_pink` ↓

These changes do not collapse into a coherent demonstrated `generalized + plain + SC` syndrome.

## Current formal regional result

The frozen joint-Wald run (`32833362756`) fitted 43 atomic categories and produced 95 context-specific slopes.

- 17 within-context slopes survive FDR;
- 16 meet the count component of confirmatory support, all in northern-midlatitude islands;
- one joint context-heterogeneity category survives FDR (`endemic × yellow_orange`);
- that comparison is 32 versus 30 islands and therefore pilot-level;
- **no confirmatory joint regional-heterogeneity result is currently established**.

The architecture must preserve the distinction between a strong within-context signal and a supported between-context difference.

## Pollination-syndrome concordance

After M0–M4 results are frozen, compare supported component vectors against pollination-syndrome literature.

```text
northern-midlatitude component vector
<-> Bombus / large- or long-tongued-bee literature
= partial concordance / mismatch / no coherent syndrome
```

```text
supported tropical/southern component vector, if present
<-> bird / butterfly / hawkmoth / alternative-bee literature
= partial concordance / mismatch
```

A single trait such as colour must not be used to infer a guild. Absence of a northern-like signal is not evidence of an alternative pollinator syndrome.

## Status of Bombus products

The existing Bombus pipeline remains useful for source-region provenance, environmental-compatibility diagnostics, effort-aware occurrence diagnostics, sensitivity analyses and candidate-system selection.

It is not used to define the primary Chapter 1 model or to turn opportunistic non-detection into historical loss.

Retained semantic boundaries:

- environmental compatibility != occurrence;
- occurrence != visitation;
- visitation != effective service;
- effective service != reproductive dependency;
- pollinator absence != historical loss;
- a floral syndrome != direct pollinator observation.

## Trait domains retained without premature binarisation

### Reproductive domain

- self-incompatibility status;
- autonomous-selfing capacity;
- mating system;
- herkogamy;
- dichogamy;
- cleistogamy;
- sex system;
- reproductive-output evidence where directly available.

### Floral signal

- white / green-brown-inconspicuous / yellow-orange / red-pink / blue-purple / variable;
- nectar guides or patterning only when directly supported;
- size / display where measured or directly described.

### Floral architecture

- open radial;
- composite head / brush-puff;
- tubular;
- funnel / trumpet;
- salverform;
- bilabiate;
- papilionaceous;
- spurred;
- symmetry;
- tube depth;
- reward accessibility where directly supported.

## Analysis sequence

1. Freeze island universe, context definition, trait evidence tiers, status/lineage rules and baseline covariates.
2. Build polygon-exact flora and provenance-tracked direct trait evidence.
3. Produce coverage and attrition audits.
4. Fit M0 context-aware baselines.
5. Fit M1 common isolation effects for atomic categories.
6. Fit M2 interaction models and joint Wald tests.
7. Apply M3 genus-preserving broad-outcome guardrails.
8. Apply M4 category-preserving decomposition.
9. Freeze component vectors and formal regional-heterogeneity classifications.
10. Compare supported vectors with literature-defined pollination syndromes only as concordance/mismatch.
11. Leave unresolved mechanism to Chapter 2 instead of filling it with weak proxy inference.

## Required outputs

- attrition from the frozen island universe to every outcome model;
- common isolation atomic effects;
- context-specific simple slopes;
- formal joint regional-heterogeneity tests;
- status/lineage-controlled broad-outcome guardrails;
- category-preserving colour/form/reproductive outputs;
- evidence-tier sensitivity analyses;
- a frozen component-vector summary suitable for discussion-level syndrome comparison;
- explicit separation between observed trait results and ecological mechanism interpretation.

## Thesis handoff

```text
Chapter 1 — island
which floral components reorganize under isolation?
and is genuine regional heterogeneity supported?
        ↓
Chapter 2 — izu-core
which ecological interaction states make responses branch / propagate / buffer?
        ↓
Chapter 3 — shimahotarubukuro
which phenotype components actually diverge within one lineage?
```

The Chapter 1 contribution is therefore **component-specific pattern discovery with explicit boundary tests**, not pollinator causal identification.
