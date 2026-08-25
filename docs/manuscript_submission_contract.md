# Manuscript submission contract

## Status

This document defines the repository surface that may support the Chapter 1 manuscript. It is intentionally narrower than the full development history.

A result is manuscript-canonical only when it follows this contract, uses locked inputs, and produces an auditable output manifest.

The current frozen scientific checkpoint is [`chapter1_frozen_result_20260825.md`](chapter1_frozen_result_20260825.md).

## Canonical scientific question

> **Does island isolation generate one coherent floral/reproductive syndrome, or does it reorganize particular floral components independently? Where are those component responses detectable, and is formal regional heterogeneity actually supported?**

The manuscript is a global comparative analysis of **trait-component reorganization and its boundary conditions**. It does not identify a historical pollinator-loss event and does not establish causal floral evolution from cross-sectional island data.

## Canonical hypothesis contrast

### H1 — coherent island-syndrome hypothesis

Isolation produces a coordinated floral/reproductive syndrome whose broad components move together and remain after floristic-status and lineage safeguards.

### H2 — component-specific floral reorganization hypothesis

Isolation is associated with changes in particular floral components, while broad syndrome summaries can cancel, weaken or disappear after category preservation and inherited lineage composition are constrained.

The current frozen evidence favors **H2 over H1**.

### Secondary boundary hypothesis — biogeographic contingency

Isolation slopes may differ among biogeographic contexts. This is retained as a formal boundary-condition test, but it is **not** the current headline result. A regional difference is claimed only from the joint interaction test, never from comparing one significant regional slope with one nonsignificant slope.

## Canonical nested analysis hierarchy

The model ladder must distinguish regional mean composition differences from regional differences in the isolation slope.

1. **M0 — context-aware baseline**
   - area and climate plus frozen biogeographic-context main effects;
   - captures baseline structure and regional mean-composition differences before isolation is added.

2. **M1 — common isolation effect**
   - adds one common isolation slope to M0;
   - `M0 -> M1` tests whether isolation adds information beyond the context-aware baseline;
   - atomic M1 effects are part of the primary evidence for component-specific reorganization.

3. **M2 — regional slope heterogeneity / boundary test**
   - adds only `isolation × biogeographic context` interactions to M1;
   - the primary heterogeneity statistic is a cluster-robust **joint Wald test** of all interaction coefficients for each atomic category;
   - joint p-values are BH-corrected across atomic states within each `stratum × trait` family;
   - individual interaction coefficients are secondary decomposition, not the headline contingency test.

4. **M3 — status / lineage guardrail**
   - broad outcomes are compared with a genus-composition-preserving expectation;
   - manuscript-level inference must not treat inherited lineage composition as optional;
   - broad `generalized_form`, `plain_colour`, and `self_compatibility` outcomes cannot be promoted if M3 does not retain them at adequate support.

5. **M4 — category-preserving decomposition**
   - colour and floral-form multistates are represented as atomic category presences; `SC|SI` is retained once as `mixed_or_variable`;
   - identifies which categories create, oppose or cancel broad summaries;
   - composite syndrome scores are secondary only.

## Current frozen result boundary

The canonical joint-Wald run (`32833362756`) fitted 43 atomic categories and produced 95 context-specific slopes.

- 17 within-context isolation slopes survived FDR;
- 16 met the count component of confirmatory support and all 16 occurred in northern-midlatitude islands;
- only one atomic category survived the formal joint contingency FDR (`endemic × yellow_orange`);
- that contrast used 32 northern-midlatitude versus 30 tropical islands and is therefore pilot-level;
- there is currently **no confirmatory biogeographic-contingency result**;
- M3 does not retain confirmatory northern/tropical isolation slopes for the broad `generalized_form`, `plain_colour`, or `self_compatibility` outcomes in the main all-native/native-nonendemic strata.

Therefore manuscript text must describe the current result as **component-specific floral reorganization**, not as a demonstrated classic island syndrome and not as a confirmed region-specific syndrome.

## Context-support rule

Small regional subsets cannot define a primary regional-difference result.

- `<30` islands for a context/category: exploratory only;
- `30–49`: pilot-level support;
- `>=50`: meets the count component of confirmatory support.

All effect tables retain the actual number of contributing islands. Statistical significance does not override an inadequate context count.

## Pollination-syndrome interpretation contract

Pollination-syndrome expectations enter only **after M0–M4 results and the atomic trait vectors are frozen**.

Allowed:

- compare the northern-midlatitude component vector with literature on Bombus / large- or long-tongued-bee associated systems;
- state explicitly whether the vector shows partial concordance, mismatch, or no coherent syndrome;
- compare tropical/southern vectors with bird, butterfly, moth, hawkmoth or alternative-bee expectations only where supported component patterns exist.

Not allowed:

- infer actual pollinator identity from floral phenotype;
- claim that Bombus loss caused the northern pattern;
- call the current northern result a demonstrated `generalized + plain + SC` Bombus-loss syndrome;
- claim bird/Lepidoptera replacement from the absence of a northern-like pattern;
- treat “significant here, nonsignificant there” as evidence of regional heterogeneity;
- use syndrome concordance as a fitted predictor or as evidence of historical functional replacement.

## Status of Bombus analyses

Existing Bombus applicability, environmental-niche, occurrence and channel-score products remain valid provenance-preserving development products, but they are **not canonical primary predictors for Chapter 1**.

They may be used only as descriptive regional context, exploratory diagnostics, declared sensitivity analyses, or candidate-system selection for mechanistic follow-up.

They must not be used to convert climatic compatibility into realized occurrence, convert opportunistic non-detection into historical loss, equate occurrence with visitation/effectiveness/dependency, or rescue a weak/null Chapter 1 mechanism by further mechanism-targeted acquisition.

## Evidence-tier rule

The manuscript reports trait evidence resolution rather than hiding it behind complete fill.

- **Confirmatory:** direct source-backed species evidence.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback or explicitly labelled inference layers.

Every primary figure/table must state the evidence tier and the number of islands, species, and trait-resolved trials contributing to the estimate.

## Required model safeguards before submission

### Spatial and lineage structure

A manuscript-canonical inference must include or explicitly test sensitivity to spatial dependence among islands and inherited floristic/lineage or source-pool composition.

The current primary lineage safeguard is the genus-composition-preserving M3 null. A globally complete source-region assignment is not required merely to duplicate that safeguard; source-region assignments remain a useful independent sensitivity where defensible.

### Biogeographic-context definition

The context grouping must be frozen before primary outcome interpretation, based on geography/biogeographic information independent of floral/reproductive outcomes, and accompanied by a reasonable sensitivity grouping where feasible.

### Isolation functional form

The isolation response must be predeclared and checked for leverage. The manuscript must report the observed isolation distribution and compare the chosen functional form with a prespecified robust alternative such as log-distance or a constrained smooth.

### Zero-distance islands

`distance_to_continent_km <= 0` or analogous zero-distance classes must be audited rather than silently dropped. Their count, geography and treatment must be reported.

### Outcome support

Each response model may use its maximum defensible support. M0/M1/M2 for a given category must use identical islands because the comparisons are nested. Attrition is reported per category.

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

Chapter 1 currently establishes that isolation-associated floral change is **component-specific rather than one coherent syndrome**, while formal confirmatory regional heterogeneity remains weak.

The unresolved question handed to Chapter 2 (`izu-core`) is:

> **Which ecological interaction states make particular floral/reproductive response components branch, propagate or buffer under isolation?**

Candidate explanations include pollinator identity, functional diversity, trait matching, effective service, functional replacement, reproductive assurance, network context and non-pollination geography/history.

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
