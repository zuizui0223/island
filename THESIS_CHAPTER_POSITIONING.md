# Thesis positioning — Chapter 1

> **Current scientific surface: corrected geography + final traitwise H1 (4 October 2026).**  
> Active selectors: config/chapter1_submission_current.json and config/chapter1_h1_final_traitwise_t_20261004.yml.  
> Earlier composite-score, North-versus-Tropical, lineage-first and Bombus-centred Chapter 1 designs are historical provenance.

## Role in the dissertation

This repository is the **Chapter 1 macroecological evidence layer**.

Chapter 1 now asks:

> **Which components of plant reproductive response recur with island isolation, which floral phenotypes remain region dependent, and is the same geographic gradient independently associated with stronger pollen limitation?**

The chapter therefore identifies a broad ecological pattern and its functional alignment. It does **not** identify the historical causal mechanism that generated contemporary island assemblages.

## Current empirical spine

### H1 — recurrent reproductive function, contingent floral phenotype

Seven binary traits are fitted separately in four geographic regions with corrected source-matched isolation, island area and climate covariates.

The clearest recurrent result is **self-compatibility**, which increases with isolation in all four regions in both broad All and WCVP regional-native-compatible All analyses.

Other components are less uniform:

- selfing mating system increases in several regions;
- autonomous selfing increases in northern high latitudes and the tropics;
- generalized form and actinomorphy increase in several regions;
- shallow/open tubes increase in northern high latitudes but decrease in southern extratropical broad flora;
- plain colour is concentrated in the southern response and appears in tropical WCVP-compatible flora.

The chapter therefore does not claim a universal floral checklist. Its primary interpretation is **functional recurrence with phenotypic contingency**.

### H2 — reproductive assurance does not absorb all floral response

The reproductive-assurance core is separated from floral accessibility and colour.

After conditioning on measured reproductive assurance, generalized accessibility remains positively associated with isolation in all four primary strata, with FDR support concentrated in northern high latitudes and tropical all-analysis. Tropical Direct-only remains positive but is not FDR-supported.

H2 therefore supports a partially separable response structure rather than a compulsory serial model in which isolation acts only through selfing and all floral change follows from it. This is conditional decomposition, not causal mediation.

### H3 — independent reproductive constraint

Independent GloPL pollen-supplementation experiments show increasing pollen limitation with corrected geographic isolation:

- beta = 0.09191
- SE = 0.03806
- finite-publication two-sided p = 0.01594

The offshore-only continuous-gradient sensitivity is also positive.

This is evidence for an isolation-associated **reproductive-service constraint**, not proof of a global decline in pollinator abundance or visitation.

### H4 — functional compatibility

In exact-species post-hoc overlap:

- reproductive-assurance score: beta = -0.29830, p = 0.00417;
- generalized-accessibility score: beta = -0.29566, p = 0.02334.

Thus the two response families highlighted by H2 are associated with lower current pollen limitation. H4 strengthens functional interpretation but does not establish historical mediation.

## Chapter 1 synthesis

The current inferential structure is:

geographic isolation
  -> recurrent reproductive-assurance traits
  -> regionally contingent floral phenotype

and independently:

geographic isolation
  -> stronger experimental pollen limitation

while:

reproductive assurance / accessibility
  -> lower current pollen limitation in exact-species overlap

The missing historical edge is:

past isolation-associated reproductive constraint
  -> selection / sorting / persistence
  -> contemporary island trait composition

Chapter 1 does not claim to have identified that edge.

## Handoff to Chapter 2 — izu-core

Chapter 1 ends with a macroecological result:

> **The reproductive problem associated with isolation is more repeatable than the floral phenotype through which plants respond to it.**

Chapter 2 asks the mechanism question directly:

> **How do realized changes in pollination channels alter effective pollen transfer, reproductive success and plant-specific floral responses, and why do species respond differently to the same deterioration in pollination service?**

The dissertation handoff is therefore:

Chapter 1 — WHERE / WHAT RECURS
global isolation gradients
-> recurrent reproductive function
-> regionally contingent floral response
-> independent pollen-limitation correlate

Chapter 2 — HOW THE RESPONSE IS GENERATED
realized pollination channels
-> effective pollen transfer
-> reproductive outcome
-> plant-specific response

Bombus remains a concrete local mechanism only where independent visitation/effectiveness evidence supports it. It is not inferred from Chapter 1 floral architecture.

## Claim ceiling

Chapter 1 may claim:

- recurrent self-compatibility enrichment with island isolation;
- regionally contingent expression of selfing, colour and floral structure;
- selected selfing-adjusted accessibility responses;
- increasing experimental pollen limitation with isolation;
- functional compatibility between reproductive-assurance/accessibility traits and lower current pollen limitation.

Chapter 1 must not claim:

- a family-wise significant universal floral syndrome;
- historical causal mediation from pollen limitation to trait evolution;
- global pollinator abundance or visitation decline;
- realized pollinator identity from phenotype;
- within-lineage evolution rather than assemblage filtering;
- corrected geography as prospective confirmation.

## Sources of truth

Use, in order:

1. config/chapter1_submission_current.json
2. config/chapter1_h1_final_traitwise_t_20261004.yml
3. submission/chapter1_current/MANUSCRIPT.md
4. results/h1_final_traitwise_t_20261004/
5. results/geography_20260924/
6. docs/PAPER_PIPELINE.md
