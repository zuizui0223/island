# Manuscript submission contract

## Status

This document defines the repository surface that may support the Chapter 1 manuscript. It is intentionally narrower than the full development history.

A result is manuscript-canonical only when it follows this contract, uses a locked input artifact, and produces an auditable output manifest.

## Canonical scientific question

> **When and where does island isolation coincide with floral and reproductive filtering, and is that filtering globally uniform or biogeographically contingent?**

The manuscript is a global comparative analysis of **pattern and boundary conditions**. It does not identify a historical pollinator-loss event and does not establish causal floral evolution from cross-sectional island data.

## Canonical hypothesis

### Biogeographically contingent floral island syndrome

Island isolation should not generate one universal floral and reproductive syndrome. Instead, the magnitude and direction of isolation-associated trait filtering may differ among biogeographic contexts after floristic status, lineage/source-pool composition, climate, geography and observation support are represented.

## Canonical analysis hierarchy

1. **M0 — biogeographic baseline**
   - area, climate, floristic status/endemicity, lineage/source-pool composition, establishment/reachability where represented, and observation/evidence support;
   - estimates how much trait composition is explained before isolation is added.

2. **M1 — universal isolation baseline**
   - adds isolation to M0;
   - tests whether one broadly shared isolation-associated direction remains.

3. **M2 — biogeographic contingency**
   - adds `isolation × biogeographic context`;
   - this is the primary Chapter 1 model;
   - the context definition must be frozen independently of floral/reproductive outcomes.

4. **M3 — status / lineage residual analysis**
   - tests M1/M2 within conservative floristic-status subsets and/or against genus/family/source-pool-preserving expectations;
   - manuscript-level inference must not treat inherited lineage composition as optional.

5. **M4 — category-preserving decomposition**
   - retains observed colour, floral-form, architecture and reproductive categories;
   - identifies which categories create, oppose or cancel a broad regional response;
   - composite syndrome scores are secondary only.

## Pollination-syndrome interpretation contract

Pollination-syndrome expectations enter only **after M0–M4 results are frozen**.

Allowed:

- compare a northern-temperate multivariate trait direction with literature on Bombus / large- or long-tongued-bee associated systems;
- compare tropical/southern counter-patterns with literature on bird, butterfly, moth, hawkmoth or alternative-bee syndromes;
- discuss concordance and mismatch explicitly.

Not allowed:

- infer actual pollinator identity from floral phenotype;
- claim that Bombus loss caused a northern pattern;
- claim that bird/Lepidoptera replacement caused a contrasting regional pattern;
- use syndrome concordance as a fitted predictor or as evidence of historical functional replacement.

## Status of Bombus analyses

Existing Bombus applicability, environmental-niche, occurrence and channel-score products remain valid provenance-preserving development products, but they are **not canonical primary predictors for Chapter 1**.

They may be used only as:

- descriptive regional context;
- exploratory diagnostics;
- declared sensitivity analyses;
- candidate-system selection for mechanistic follow-up.

They must not be used to:

- define the primary analysis universe;
- convert climatic compatibility into realized occurrence;
- convert opportunistic non-detection into historical absence/loss;
- equate occurrence with visitation, effectiveness or reproductive dependency;
- rescue a weak/null Chapter 1 pattern by additional mechanism-targeted acquisition.

## Evidence-tier rule

The manuscript must report trait evidence resolution rather than hide it behind complete fill.

- **Confirmatory:** direct source-backed species evidence, with any broader direct-evidence track declared separately.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback or explicitly labelled inference layers.

Every primary figure/table must state the evidence tier and the number of islands, species, and trait-resolved trials contributing to the estimate.

## Required model safeguards before submission

### Spatial and lineage structure

A manuscript-canonical global or regional inference must include, or explicitly test sensitivity to:

- spatial dependence among islands; and
- inherited floristic/lineage or source-pool composition.

A model containing only geography and climate does not satisfy the Chapter 1 safeguard.

### Biogeographic-context definition

The primary context grouping must be:

- frozen before primary outcome interpretation;
- based on defensible geography/source-pool/biogeographic information;
- independent of floral/reproductive trait outcomes;
- accompanied by at least one reasonable sensitivity grouping where feasible.

### Isolation functional form

The isolation response must be predeclared and checked for leverage. The manuscript must report the observed isolation distribution and compare the chosen functional form with a prespecified robust alternative such as log-distance or a constrained smooth.

### Zero-distance islands

`distance_to_continent_km <= 0` or analogous zero-distance classes must be audited rather than silently dropped. Their count, geography and treatment must be reported.

### Outcome support

Each response model may use its maximum defensible support. A common complete-case set is required only for comparisons that mathematically require identical support. Attrition must be reported per equation/outcome.

### Category preservation

Multistate/category outcomes precede broad binary summaries. A binary result cannot override an opposing or heterogeneous category decomposition.

## Reproducibility and archival rule

Before submission, all manuscript-critical inputs and outputs must be archived durably with checksums.

A fresh reader should be able to identify from the repository root:

- the exact manuscript release/tag;
- locked input artifacts and checksums;
- one canonical command/workflow per main analysis;
- the software environment;
- evidence-tier definitions;
- attrition from the frozen island universe to every fitted model; and
- the files used for each manuscript figure and table.

## Chapter 1 -> Chapter 2 boundary

Chapter 1 may establish that regional floral/reproductive trait syndromes differ. It should not use weak proxy data to decide why.

The unresolved mechanism is handed to Chapter 2 (`izu-core`), where candidate explanations such as pollinator identity, functional diversity, trait matching, effective service, functional replacement, reproductive assurance and network context are evaluated under a stricter mechanistic standard.

## Noncanonical material

The following remain development history or sensitivity material unless explicitly promoted by a later contract change:

- trait-source scouting and free-bulk-source pilots;
- validation/core-pilot acquisition workflows;
- EOL/TraitBank inventory experiments;
- Bombus-primary M0–M4 pathway variants;
- Bombus mediation/coefficient-product analyses;
- old v1/v2 bridge analyses;
- superseded Bayesian/INLA variants built around a Bombus-primary confirmatory question.

Historical code may remain in git for provenance. It should not be presented beside the current canonical Chapter 1 workflow as if it were still the primary manuscript method.
