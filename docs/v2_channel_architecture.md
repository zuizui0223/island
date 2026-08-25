# v2: context-dependent island floral syndrome — analysis architecture

## Status and boundary

`v1-freeze` preserves the exploratory draft and its first complete results. v2 is rebuilt from new data products and new scripts; it must not import v1 scripts, derived trait labels, or v1 model outputs as analysis inputs.

The current scientific scope is defined in [`v2_pollination_regime_framework.md`](v2_pollination_regime_framework.md).

The primary Chapter 1 question is now:

> **When and where does island isolation coincide with floral and reproductive filtering, and is that filtering globally uniform or biogeographically contingent?**

The architecture therefore prioritizes direct floral/reproductive composition, floristic status, lineage/source-pool structure and geography. Pollinator-channel data are retained as interpretation and sensitivity layers rather than being required to define the primary analysis universe.

```text
isolation / geography / source-pool context
        ↓
floristic status / lineage composition / observation process
        ↓
regional floral and reproductive trait response
        ↓
category-preserving decomposition
        ↓
pollination-syndrome concordance in Discussion
        ↓
unresolved mechanism handed to Chapter 2
```

## What stays from earlier v2

1. Island geography is the starting context.
2. Source-pool, floristic status, lineage composition, establishment/reachability and observation process must be considered before ecological interpretation.
3. Floral signal, architecture and reproductive traits remain multistate or continuous wherever possible.
4. Binary summaries such as plain/conspicuous, generalized/specialized and SC/SI are secondary summaries.
5. Wind pollination, self-compatibility, autonomous selfing and pollinator independence remain separate domains.
6. Existing Bombus applicability, environmental and occurrence diagnostics remain provenance-preserving data products.
7. The global analysis does not claim within-lineage evolution or historical pollinator causation from cross-sectional composition.

## What changes

1. Bombus deficit is no longer the primary Chapter 1 explanatory variable.
2. The primary confirmatory question is `isolation × biogeographic context`.
3. Regional trait vectors are fixed before any pollination-syndrome interpretation.
4. Bombus-, bird-, Lepidoptera-, hawkmoth- and other syndrome expectations are compared symmetrically at the interpretation stage.
5. The mechanistic question — which functional pollination pathway actually generated a regional pattern — is explicitly handed to Chapter 2 (`izu-core`).

## Interpretation guardrail borrowed from Campanula/Izu

For a trait state `z`, the Campanula programme writes:

```text
W(z) = F(z) × E(z)
```

For Chapter 1 this remains an interpretation guardrail, not a fitted decomposition.

- `F(z)` — local reproductive contribution;
- `E(z)` — arrival, establishment, persistence and reachability conditional on viable reproduction;
- `W(z)` — observed island-flora membership or trait composition.

The global dataset mainly observes `W`. Therefore an isolation-associated trait pattern may arise from establishment filtering, lineage turnover, ecological interaction filtering, reproductive processes, or observation structure. Chapter 1 should identify the pattern and strongest boundary conditions without pretending to identify `F` from `W` alone.

## Primary model architecture

### M0 — biogeographic baseline

```text
area + climate
+ floristic status / endemicity
+ lineage / source-pool composition
+ establishment / reachability where represented
+ observation / evidence support
-> reproductive composition
-> floral phenotype composition
```

### M1 — universal isolation baseline

```text
M0 + isolation -> trait composition
```

This tests the classical idea that isolation produces one broadly shared direction.

### M2 — biogeographic contingency

```text
M0 + isolation × biogeographic context -> trait composition
```

This is the primary Chapter 1 model. A strong result means the same geographic isolation is associated with different trait responses in different biogeographic contexts.

### M3 — status / lineage residual analysis

Repeat M1/M2 where possible in conservative native/endemicity strata and against lineage-preserving expectations.

```text
observed trait composition
-
status / lineage / source-pool expectation
=
residual trait response
```

### M4 — category-preserving decomposition

Analyse retained categories rather than forcing one syndrome score.

Key domains include:

- reproductive compatibility / assurance states;
- colour categories;
- floral form categories;
- symmetry;
- tube depth / access;
- flower size / display;
- other direct floral architecture evidence.

## Pollination-syndrome concordance

After M0–M4 results are frozen, compare each regional multivariate trait direction against expectations from pollination-syndrome literature.

### Example interpretation structure

```text
northern-temperate trait vector
<-> Bombus / large- or long-tongued-bee literature

contrasting tropical/southern trait vector
<-> bird / butterfly / hawkmoth / alternative-bee literature
```

The comparison is **concordance**, not pollinator assignment.

A single trait such as colour must not be used to infer a guild. Interpretation should rely on combinations of direct traits and acknowledge mismatches explicitly.

## Status of Bombus products

The existing Bombus pipeline remains useful for:

- source-region provenance;
- environmental-compatibility diagnostics;
- effort-aware presence / non-detection diagnostics;
- sensitivity analyses;
- identifying candidate systems for future direct testing.

It is not used to define the primary Chapter 1 model or to turn opportunistic non-detection into historical loss.

### Retained semantic boundaries

- environmental compatibility != occurrence;
- occurrence != visitation;
- visitation != effective service;
- effective service != reproductive dependency;
- pollinator absence != historical loss;
- a floral syndrome != direct pollinator observation.

These distinctions are especially important because the Chapter 1 handoff asks Chapter 2 to discriminate among those layers mechanistically.

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

1. Freeze the island universe, biogeographic context definition, trait evidence tiers, status/lineage rules and M0 covariates.
2. Build polygon-exact flora and provenance-tracked direct trait evidence.
3. Produce coverage and attrition audits.
4. Fit M0.
5. Fit M1 global isolation baselines.
6. Fit M2 `isolation × biogeographic context`.
7. Apply M3 status/lineage/source-pool safeguards.
8. Apply M4 category-preserving decomposition.
9. Freeze regional trait vectors.
10. Compare those vectors with literature-defined pollination syndromes in Discussion.
11. Leave unresolved pollinator mechanism to Chapter 2 instead of filling it with weak proxy inference.

## Required outputs

- attrition from the frozen island universe to every outcome model;
- region/context-specific isolation slopes;
- status/lineage-controlled sensitivity results;
- category-preserving colour/form/reproductive outputs;
- evidence-tier sensitivity analyses;
- a regional trait-vector summary suitable for syndrome concordance discussion;
- explicit separation between observed trait results and ecological mechanism interpretation.

## Thesis handoff

```text
Chapter 1 — island
when / where do different regional trait syndromes appear?
        ↓
Chapter 2 — izu-core
which ecological mechanisms generate branching / propagation / buffering?
        ↓
Chapter 3 — shimahotarubukuro
which phenotype components actually diverge within one lineage?
```

The Chapter 1 contribution is therefore **boundary-condition discovery**, not pollinator causal identification.
