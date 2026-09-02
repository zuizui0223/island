# Chapter 1 explanation-gap validation plan

## Decision

The baseline Chapter 1 result is already sufficient to reject a simple universal floral-island syndrome and to establish biogeographically contingent assemblage responses. The remaining uncertainty is concentrated at the **edges between H2, H3, H4 and H5**, not in the existence of a plant-side pattern itself.

These additions are therefore frozen as **post-baseline prospective falsification tests**. They do not rewrite H1-H5 or turn a favorable post hoc result into a new primary test. Their purpose is to find where the current interpretation breaks.

## Current result tree and the remaining gaps

```text
H1 universal syndrome
   mostly rejected by direct regional vector heterogeneity
          |
          v
H2 biogeographic branching
   supported North-vs-Tropical / Palearctic structure
          |
          | unresolved: region or measured climate-dependent response?
          v
H3 source / lineage assembly
   genus entry explains much of Palearctic floral accessibility
          |
          | unresolved: family-level or genus-level sorting?
          v
H4 area / capacity moderation
   smaller-island amplification in Palearctic genus entry
          |
          | unresolved: biology or trait-support/precision artifact?
          v
H5 pollination-channel mechanism
   not identified
```

A secondary branch sits between H2 and H5:

```text
observed plant response
   -> large-bee / butterfly / bird architecture concordance
   -> source-pool sensitivity survives
   -> unresolved: shared architecture? lineage sorting? guild-specific residual?
```

## 1. V5 — trait-resolution MNAR tipping point [highest priority]

### Why this is currently underexplained

The latest Wave36 coverage is strongly asymmetric across axes: flower colour and structure are much better resolved than reproductive assurance. The existing trait-resolution model adjusts for **measured predictors of whether a trait is resolved**, but cannot prove that unresolved reproductive states are missing at random conditional on those predictors.

This matters directly for tropical/Neotropical reproductive-assurance results.

### Additional test

Add a partial-identification / selection-model sensitivity in which unresolved states are allowed to be systematically more or less likely to belong to the focal state than resolved species. The missingness shift is varied on a grid fixed before execution. Record the value at which the sign, multivariate gate or regional contrast changes.

The pre-execution parameter is the odds ratio for being trait-resolved in the focal
`selfing_core` state versus the non-focal state. The frozen grid is `0.1, 0.2, 1/3,
0.5, 1, 2, 3, 5, 10`. It is applied as (i) a shared state-dependent mechanism,
(ii) a symmetric North-versus-Tropical differential, and (iii) a symmetric
Palearctic-versus-Neotropical differential. OR = 1 must reproduce the observed island
scores exactly. Regression weights retain the observed number of resolved species, so
expected completion of unresolved states does not manufacture precision. Deterministic
all-missing-focal and all-missing-nonfocal bounds are reported separately.

### What it answers

Instead of saying only “we adjusted for trait-resolution probability,” the manuscript can say either:

- the result survives substantial unmeasured state-dependent missingness; or
- only a modest MNAR shift is needed to erase it, so the reproductive branch remains provisional.

This is especially useful while trait acquisition is still ongoing because each new snapshot should push the tipping point outward if the signal is real.

## 2. V1 — climate-overlap versus biogeographic branching

### Why this is currently underexplained

H2 includes climate PCs as main-effect covariates, but a categorical North/Tropical or Palearctic/Neotropical difference in the **distance slope** could still arise because the effect of isolation itself varies continuously with climate.

This alternative is biologically credible. Malesian island analyses have found that flower size, colour and openness depend jointly on climate, island area and isolation, with predictive importance varying among traits.

### Additional tests

1. **Outcome-blind climate/area overlap weighting** between compared contexts, using only log area and frozen climate PCs to construct common-support weights. Then rerun the direct distance x context response-vector test.
2. **Distance x climate joint sensitivity**: add predeclared distance x climate-PC interactions and ask whether a residual distance x context vector remains.

### Interpretation

- If the context contrast survives both, “biogeographic branching” becomes much stronger.
- If climate interactions absorb it, the correct conclusion becomes **environment-dependent branching**, which is still scientifically useful but narrower.
- If there is little climate overlap, the two regional slopes should not be presented as directly transportable over the same environmental support.

## 3. V4 — decompose the sampled pollination-syndrome concordance

### Why this is currently underexplained

The baseline shows:

- northern mid-latitude large-bee-like decline;
- Palearctic decline of large-bee-, butterfly- and bird-like scores together;
- tropical increase of all three together.

This is clearly not a pollinator classifier. The remaining question is whether it is:

1. one shared “specialized/attractive architecture” axis;
2. a lineage-composition effect;
3. or whether any guild-specific residual remains after those two are removed.

### Additional tests

1. Give each sampled guild axis its own **family/genus composition null and source-matched entry/loading decomposition**.
2. Estimate a shared template factor using **mainland/source species only** and no island outcomes. Project islands onto the shared factor plus residual template axes.

### Interpretation

- only common factor survives -> evidence for shared floral architecture, not a specific pollination channel;
- genus null absorbs -> architecture pattern is primarily lineage assembly;
- a residual guild axis survives both -> candidate architecture for later H5 testing, still not pollinator identity.

This is the most useful extra validation for connecting the secondary concordance layer to the main H2-H3 tree.

## 4. V2 — locate the taxonomic depth of H3

### Why this is currently underexplained

The current genus-preserving null tells us that genus composition can reproduce much of the Palearctic accessibility signal. But “genus composition” could mean very different biological things:

- deep clade/family sorting;
- genus-level source filtering;
- or fine species sorting inside represented genera.

### Additional test

Run the same decomposition at nested taxonomic scales:

```text
observed
  - family-preserving expectation
  - genus-preserving expectation
  - within-genus species loading
```

Family labels are used only as grouping variables, never to infer missing traits.

### Interpretation

- family null already absorbs the signal -> deep clade sorting;
- family residual remains but genus null absorbs -> genus-level ecological/biogeographic sorting;
- genus residual remains -> within-genus or non-taxonomic filtering remains.

This would make H3 substantially more mechanistically interpretable without pretending to identify historical ancestry.

## 5. V3 — test whether the H4 area effect is a measurement-support artifact

### Why this is currently underexplained

The Palearctic genus-entry distance effect is stronger on smaller islands in all frozen primary scenarios. That is compatible with a capacity/target-size mechanism, but larger islands also tend to contain more species and may provide more stable trait means.

Species richness is itself a possible mediator of island biogeography, so simply adding richness as a causal covariate would be unsafe.

### Additional tests

Treat this explicitly as a **measurement-artifact sensitivity**, not as a new primary causal model:

- equal-island area interaction;
- capped-information-weight area interaction;
- direct-only area interaction;
- common-support diagnostics for `n_species` and mean informative traits across the area gradient;
- heteroskedastic simulations asking whether the observed support-area relationship alone can generate an interaction under a true zero interaction.

### Interpretation

If the interaction remains stable, “area/capacity moderation” becomes much more defensible. If it is precision-sensitive, area should remain a statistical modifier without a capacity interpretation.

## What should not be added yet

### Do not add another floral proxy for H5

H5 remains the only step that fundamentally needs **new independent data**:

- source-channel exposure;
- retained / disrupted / structurally absent state;
- effort-aware occurrence or abundance;
- interaction/visitation;
- ideally per-visit effectiveness.

Adding more Bombus suitability, more syndrome labels or more floral proxies does not solve that gap.

## Recommended implementation order

1. **V5 MNAR tipping point** — highest leverage because reproductive coverage is still the weakest and this directly determines how much to trust the current tropical reproductive branch.
2. **V1 climate-overlap test** — directly challenges the word “biogeographic” in H2.
3. **V4 guild common-factor + lineage decomposition** — clarifies the strongest unresolved bridge toward H5.
4. **V2 family-to-genus depth** — refines what H3 actually means biologically.
5. **V3 area-support artifact test** — strengthens or bounds the capacity interpretation of H4.

The full machine-readable specification is `config/chapter1_explanation_gap_validation.yml`.
