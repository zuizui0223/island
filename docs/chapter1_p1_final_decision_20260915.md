# Chapter 1 P1 final decision — 2026-09-15

Status: post-baseline robustness decision. This note integrates the frozen P1a, P1c and P1d results. It does not replace the historical H3 result or alter any preregistered threshold.

## Question

Does the strong Palearctic floral-island response disappear after genus adjustment because of a trivial analysis artifact, or do true genus boundaries carry biological structure relevant to the attenuation?

## P1a — paired support

P1a passed. Observed, family-adjusted and genus-adjusted stages use the same frozen species/island support and the same `n_species` information weights in the focal audit cells. Family and genus are grouping variables, not trait-imputation devices. Stage-specific sample loss is therefore not a viable explanation for the attenuation.

Formal source: run `34935183075`, artifact `10382749051`, digest `sha256:4c4689d5d8e0eaca7e3ac66e1e866a48942b2f7220e962f7b1bb7e7ac1d43af2`.

## P1c — matched-complexity pseudo-genus null

P1c tested the strongest model-flexibility objection. Species were randomly repartitioned *within family* while preserving the exact real genus-count and genus group-size multiset in every family. Thus pseudo-genera had the same grouping complexity as true genera but different biological membership.

All 2,000 frozen permutations were valid. The observed median conditional genus attenuation was `0.7206615`; the matched pseudo-genus null median was `0.2243375`. Only 57/2,000 null permutations equalled or exceeded the observed statistic, giving a one-sided randomization p-value of `0.0289855`.

Decision: `true_genus_exceeds_matched_complexity_null`.

This rejects the simple explanation that any equally fine within-family grouping would absorb the same response merely because the model is flexible.

Permutation source: run `34936193944`; aggregate-only repaired run `34941827774`, artifact `10385820775`, digest `sha256:861d18fe87b9190619f925a2446be5fd4d460b818825578930883257c0a6ed16`. The aggregate repair did not regenerate permutations; it replaced bitwise equality of independently recomputed observed statistics with a fail-closed `max-min <= 1e-12` invariant. The observed-statistic spread was `3.33e-16`.

## P1d — paired spatial-block uncertainty

P1d reproduced the frozen direct-only slopes to numerical precision and then propagated spatial uncertainty with 2,000 paired spatial-block bootstrap draws.

The observed total genus-adjusted attenuation remained large (`0.788–0.802` across the eight primary direct-only profiles), with bootstrap lower bounds still positive (`0.277–0.402`). However, the *incremental* family-to-genus attenuation was not precisely localized: all eight 95% intervals for the additional family-to-genus attenuation included zero.

Decision: the total genus structuring is strong, but the amount specifically attributable to the family-to-genus increment is spatially imprecise.

Formal source: run `34939113182`, artifact `10383754545`, digest `sha256:822c69e9a1be4eb2c340391c1fb84dfba91b45b8a765549bb8ad1ec6379ef480`.

## Integrated decision

The evidence supports the following claim:

> **The Palearctic floral-island response is strongly structured by genus-level taxonomic composition. True genus boundaries absorb more of the response than arbitrary within-family partitions of identical grouping complexity, while the exact incremental family-to-genus attenuation is not precisely estimated across spatial blocks.**

This is stronger than saying only that “adding genus terms makes significance disappear,” because sample loss and generic fine-grouping complexity have been explicitly challenged. It is narrower than saying that the syndrome has a precisely measured breakpoint at the family-to-genus transition.

## Claim ceiling

P1 may support:

- genus-specific taxonomic structure beyond matched arbitrary fine grouping;
- strong total attenuation after source-matched genus adjustment;
- the interpretation that the broad Palearctic syndrome is not a robust beyond-genus response repeated uniformly across the flora;
- `assembly depth` as a useful localization concept, provided it is described as taxonomic localization rather than a causal mechanism.

P1 may not support:

- a precisely estimated amount of incremental family-to-genus attenuation;
- dispersal-only assembly;
- absence of within-lineage evolution;
- genus itself as a causal mechanism;
- a pollinator-specific cause.

## Manuscript consequence

Replace strong phrases such as “the response is concentrated at the family-to-genus transition” or “most attenuation occurs at the family-to-genus transition” with wording that separates two findings:

1. **true genus boundaries carry non-random structure relevant to attenuation** (`P1c`), and
2. **the exact family-to-genus increment is spatially imprecise** (`P1d`).

The island-first story therefore survives, but in a more defensible form: the classic floral-island response is strongly genus-structured, and that structure is not a trivial consequence of sample loss or arbitrary fine grouping.
