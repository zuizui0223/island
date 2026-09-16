# Chapter 1 H5 — global pollination integration v3

## Current H5 decision

H5 now separates three different empirical layers that should not be conflated:

1. **plant response:** reproductive assurance and pollination-associated floral architecture;
2. **interaction structure:** GloBI source/genus interaction breadth and exact-island functional-channel occurrence;
3. **experimental pollination service:** GloPL pollen-supplementation experiments.

The new full-scale GloPL extension changes one part of the H5 conclusion. A broad isolation-associated pollen-limitation relationship is now supported across the global GloPL sampling frame. It does **not** reproduce the predeclared North–Tropical service-gradient branching expected if one context-specific pollination mechanism generated the Chapter 1 plant H2 response.

## 1. Plant-side response remains two-route

The plant data continue to support partially separable responses:

- reproductive assurance;
- pollination-associated floral architecture.

The floral response is not reducible to measured `selfing_core`, so the data reject an obligatory serial `isolation -> selfing -> floral simplification` explanation.

Named large-bee-like, butterfly-like and bird-like templates remain highly shared (>86% common architecture). The orthogonalized named-syndrome analysis is support-limited because the strict complete-template subset does not robustly retain the original six-atomic H2 across all-analysis and Direct evidence. Syndrome labels therefore remain phenotype descriptions rather than observed visitor identities.

## 2. GloBI covers the broad island scale but does not close mechanism

Global GloBI v3 uses **3,252 islands** in the all-observed primary layer, spanning northern mid-latitude, northern high-latitude, tropical and southern extratropical contexts.

The four-context distance-slope heterogeneity test is FDR-supported in **3/4** frozen source definitions. Because one source definition fails, the result remains `source_definition_sensitive_context_heterogeneity`, not a universal pollinator filter. Distance × area heterogeneity is supported in **0/4** definitions.

Thus GloBI shows that interaction structure changes geographically on almost the same island scale as plant H2, but it does not identify effective pollen delivery or a unique upstream filter.

## 3. Exact-island channel occurrence remains weak for identity-specific mechanism

No single functional channel has balanced retained/disrupted support in both primary contexts. Pooling five channels makes total attrition estimable but does not support the two plant pathways. Retaining channel identity produces a suggestive tropical reproductive-assurance signal but fails the full promotion gate.

The strict same-island identity-preserving comparison remains structurally non-identifiable in the northern context: only five northern mixed islands exist and disruption is confounded with channel identity after island and channel×context fixed effects.

These occurrence analyses therefore remain evidence about channel presence/disruption, not effective pollination service.

## 4. Exact-island GloPL: strict but small

The first GloPL route required experiments to lie inside the frozen Chapter 1 island polygons.

Outcome-blind overlap:

- northern mid-latitude: **15 islands / 55 studies**;
- tropical: **14 islands / 23 studies**.

The frozen adjusted model retained **28 islands / 77 publications / 25 spatial blocks** and controlled island area plus climate PC1–4.

Results:

- North isolation slope = **+0.3762**, one-sided p=`0.1209`;
- Tropical−North interaction = **−0.0311**, one-sided negative p=`0.4577`.

The North point estimate is directionally diffuse under leave-one-island and leave-one-publication audits, but the frozen primary result remains unsupported. This route is the stricter island-specific adjusted subset, not the full information available in GloPL.

## 5. Full-scale GloPL: all 2,969 experiments

The global extension was frozen before the distance exposure was merged to outcomes at commit `6f570b3aa97b66a8012f55953fdf8c1dd4aba953`.

Every georeferenced GloPL site was placed on the same continuous distance-to-six-seeded-major-continents definition already used by Chapter 1. No near/remote threshold was introduced.

### Outcome-blind support

All **2,969 / 2,969** GloPL rows have valid coordinates.

- **1,248** unique sites;
- **919** publications;
- **988** zero-distance mainland sites;
- **260** positive-distance offshore sites.

Context coverage:

| context | sites | studies | mainland sites | offshore sites |
|---|---:|---:|---:|---:|
| northern mid-latitude | 737 | 526 | 609 | 128 |
| tropical | 263 | 252 | 218 | 45 |
| southern extratropical | 248 | 146 | 161 | 87 |
| northern high-latitude | 0 | 0 | 0 | 0 |

The global and North–Tropical gates passed before outcome inspection. The four-context gate failed structurally because no northern-high-latitude GloPL sites exist under the frozen latitude definition.

### Frozen global model

After publication × coordinate × measurement-cell aggregation:

- **1,408** cells;
- **1,248** sites;
- **919** publications;
- each publication has total weight 1;
- measurement type/treatment/constant/supplementation category are fixed effects;
- uncertainty is publication-cluster robust.

Primary global distance result:

- slope = **+0.07937**;
- SE = `0.03773`;
- two-sided p=`0.03543`;
- predeclared one-sided positive p=`0.01772`.

Under the frozen rule, this is a **supported global association**.

Sensitivities:

- supplemental-only: slope **+0.03349**, one-sided p=`0.2149`;
- no-zero-constant: slope **+0.07789**, one-sided p=`0.02095`.

Both retain the positive direction required by the frozen global promotion rule, although the supplemental-only sensitivity is not independently statistically supported.

Canonical provenance: run **35090599662**, artifact **10444156159**, digest `sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623`.

## 6. The global association is not only a mainland/offshore step

Because most sites lie on major continents, a post-hoc shape audit asked whether the parent signal was just a discontinuity between distance zero and positive distance.

The diagnostic reused the parent **1,408 cells**, their exact weights, outcome values, measurement controls and publication clustering. It cannot alter the parent p-value or promotion status.

Support:

- mainland cells: **1,126** / 988 sites;
- offshore cells: **282** / 260 sites / 158 publications.

Two-part decomposition:

- mainland → average offshore step: **+0.13745**, two-sided p=`0.1084`;
- within-offshore distance slope: **+0.19983**, p=`0.01110`.

Offshore-only model:

- distance slope: **+0.23545**, p=`0.005865`.

Classification: `within_offshore_gradient_present`.

This is post-hoc descriptive evidence, but it shows that the supported parent association is not plausibly summarized as only a mainland-versus-offshore level shift. Pollen limitation continues to increase across the sampled offshore distance gradient.

Provenance: run **35091385147**, artifact **10444198356**, digest `sha256:ee7700a03dcbe9b0d6184195c482e3ec7915a5e7a9e8f1c6a85d0e4fe2de8e44`.

## 7. North–Tropical branching does not match plant H2

The full-scale global result does **not** establish the specific H5 route that would explain the plant H2 branching.

Frozen North–Tropical pair:

- North slope = **+0.06170**, one-sided positive p=`0.1073`;
- Tropical−North interaction = **+0.14300**;
- predeclared negative-interaction p=`0.9119`;
- implied tropical slope = **+0.20470**.

The interaction therefore points in the opposite direction to the frozen H5 branching prediction: the tropical pollen-limitation gradient is point-estimated as stronger rather than weaker.

The supplemental-only sensitivity makes this mismatch clearer: the North slope is approximately zero/slightly negative while the tropical-minus-North interaction remains positive.

Therefore the correct distinction is:

> **General isolation-associated pollen limitation: supported at the full global GloPL scale.**
>
> **The predeclared North–Tropical pollen-service branching proposed as the common upstream generator of plant H2: not supported.**

## Updated H5 hierarchy

```text
isolation-associated reproductive-assurance response
        supported
        +
isolation-associated floral-architecture response
        supported and partially separable from selfing
        ↓
named visitor identity from syndrome phenotype
        not identified
        ↓
GloBI global interaction-structure heterogeneity
        present but source-definition-sensitive
        ↓
functional-channel occurrence/disruption
        no robust global identity-specific mechanism
        ↓
strict same-island channel identity test
        northern target non-identifiable
        ↓
experimental pollen limitation — exact-island adjusted subset
        North direction positive but imprecise / unsupported
        ↓
experimental pollen limitation — full global GloPL frame
        GLOBAL DISTANCE ASSOCIATION SUPPORTED
        and post-hoc offshore-only gradient present
        ↓
North–Tropical service-gradient branching matching plant H2
        NOT SUPPORTED; interaction direction is opposite prediction
```

## Final H5 statement

> **Isolation is associated globally with stronger experimental pollen limitation across the full GloPL sampling frame, and the association continues across positive-distance offshore sites rather than being only a mainland/offshore step. This supplies independent experimental support for a broad isolation–pollination-service relationship. It does not, however, reproduce the North–Tropical branching of the floral/reproductive response: the tropical service gradient is not weaker than the northern one. Consequently, pollination-service limitation becomes a supported general correlate of geographic isolation, but it is not identified as the context-specific common upstream generator of the Chapter 1 plant H2 pattern.**

## Claim ceiling

Do not convert this result into any of the following claims:

- distance causes pollen limitation;
- pollinator abundance or visitation decline was measured;
- historical pollinator loss was identified;
- the full global model controls island area or climate;
- the post-hoc offshore shape audit is confirmatory;
- GloPL mediates the six-atomic plant H2 response;
- the North–Tropical plant branching is explained by pollen limitation.
