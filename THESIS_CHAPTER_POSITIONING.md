# Thesis positioning — Chapter 1

> **Current scientific surface: corrected geography baseline (24 September 2026).**
> Machine-readable selector: `config/chapter1_submission_current.json`.
> Older WHEN/WHERE, branching, lineage-first and Bombus-centered Chapter 1 designs are historical provenance, not the current paper claim.

## Role in the dissertation

This repository is the **Chapter 1 macroecological evidence layer** of the dissertation.

Chapter 1 now asks:

> **Does geographic isolation repeatedly reorganize floral/reproductive function across island floras, is that response separable into reproductive-assurance and pollinator-facing accessibility components, and is the same geographic gradient independently associated with stronger pollen limitation?**

The chapter establishes a global pattern, an independent ecological-pressure correlate, and post-hoc functional compatibility. It does **not** identify the historical causal mechanism that generated the contemporary assemblages.

## Current corrected baseline

- analysis universe: **8,264 island units**;
- broad H1 union: **4,379 islands**;
- plant species: **106,295**;
- resolved trait cells: **222,688 / 318,885 = 69.83%**;
- GloPL: **2,969 experiments / 1,248 sites / 919 publications**;
- geographic exposure: source-matched GSHHG 2.3.7 coastline separation.

The 24 September 2026 geography correction repaired 1,113 spurious island zero distances and excluded one continental split component. It is the primary measurement baseline, but remains a post-hoc measurement correction rather than prospective confirmation.

## Current H1–H4 scientific spine

### H1 — recurrent multivariate island response

The seven-response floral/reproductive isolation vector is jointly supported in all four predeclared geographic strata in both evidence scopes.

The result is **recurrent, not uniform**. Individual traits can be weak or move in the opposite direction; in particular, southern shallow/open tube is negative in the corrected analysis.

### H2 — two partially separable plant-response components

Reproductive assurance is separated from additional floral responses.

After conditioning on measured reproductive assurance, generalized/accessibility responses remain positive in all four primary strata and are FDR-supported in northern high latitudes and the tropics. Tropical Direct-only accessibility is nominally positive but not FDR-supported (`q=0.1196`).

Therefore floral reorganization is not reducible to a compulsory serial pathway of:

```text
isolation -> selfing -> floral simplification
```

The supported interpretation is a **partially separable reproductive-assurance route plus an additional pollinator-facing accessibility route**. H2 is conditional decomposition, not causal mediation.

### H3 — independent ecological-pressure correlate

Independent GloPL pollen-supplementation experiments show increasing pollen limitation with corrected geographic isolation:

```text
beta = 0.09191
SE   = 0.03806
p    = 0.01575
```

This is evidence for an isolation-associated **pollination-service constraint**, not proof of a global decline in pollinator abundance or visitation.

### H4 — functional compatibility

In exact-species post-hoc triangulation, the two H2 trait families are associated with lower current pollen limitation:

- reproductive assurance: `beta=-0.29830`, `p=0.00396`;
- generalized accessibility: `beta=-0.29566`, `p=0.02187`.

These associations are functionally compatible with reduced dependence on external pollen delivery. They do not establish that historical pollen limitation mediated the contemporary island trait pattern.

## What Chapter 1 no longer claims

The current paper is **not** organized around:

- rejection of a universal syndrome in favour of North-vs-Tropical branching;
- a Palearctic-only floral architecture result;
- source/lineage decomposition as the main paper spine;
- Bombus loss as a global or primary Chapter 1 mechanism;
- pollination-syndrome scores as realized visitor identities;
- within-lineage evolution rather than assemblage filtering.

Those analyses remain useful provenance and sensitivity history, but they are not the current submission-level result hierarchy.

## Current inferential hierarchy

```text
geographic isolation
      |
      +--> recurrent multivariate floral/reproductive response (H1)
      |        |
      |        +--> reproductive assurance
      |        |
      |        +--> additional accessibility/generalization response (H2)
      |
      +--> stronger experimental pollen limitation (H3)

H2 trait states
      |
      +--> lower current pollen limitation in exact-species overlap (H4)
```

The dashed historical causal bridge remains unresolved:

```text
past isolation-associated pollination constraint
      -> selection / sorting / persistence
      -> contemporary island trait composition
```

Chapter 1 does not claim this sequence has been directly identified.

## Handoff to Chapter 2 — `izu-core`

Chapter 1 ends at the macroecological claim ceiling:

> **Isolation is repeatedly associated with a functional shift toward reproductive assurance and floral accessibility, and with stronger pollen limitation, but the causal route from interaction change to plant response remains unresolved.**

Chapter 2 asks the mechanism question directly:

> **How do changes in realized pollination channels alter effective pollen transfer, reproductive success and plant-specific floral responses, and why do species respond differently to the same deterioration in pollination service?**

The dissertation handoff is therefore:

```text
Chapter 1 — GLOBAL PATTERN + ECOLOGICAL PRESSURE + FUNCTIONAL COMPATIBILITY
                |
                v
Chapter 2 — REALIZED CHANNELS + EFFECTIVE SERVICE + PLANT-SPECIFIC MECHANISM
```

Bombus can remain a concrete local interaction mechanism in Chapter 2 where independently measured visitation/effectiveness supports it. It is not inferred from Chapter 1 floral architecture.

## Claim ceiling

Chapter 1 may claim:

- recurrent multivariate floral/reproductive response to geographic isolation;
- partially separable reproductive-assurance and accessibility components;
- increasing experimental pollen limitation with geographic isolation;
- exact-species functional compatibility between the two response families and lower current pollen limitation.

Chapter 1 must not claim:

- historical causal mediation from pollen limitation to trait evolution;
- global pollinator abundance or visitation decline;
- one universal named pollinator mechanism;
- uniform positive change in every floral trait;
- realized pollinator identity from flower colour or syndrome templates;
- species sorting versus within-lineage evolution as already identified;
- the corrected geography analysis as prospective confirmation.

## Sources of truth

Use, in order:

1. `config/chapter1_submission_current.json`
2. `submission/chapter1_current/MANUSCRIPT.md`
3. `docs/chapter1_corrected_submission_20260924.md`
4. `results/geography_20260924/`
5. `docs/PAPER_PIPELINE.md`

Historical design documents are indexed in `docs/CHAPTER1_HISTORY.md`.
