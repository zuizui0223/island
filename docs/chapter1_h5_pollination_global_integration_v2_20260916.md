# Chapter 1 H5 — global pollination integration v2

## Current H5 decision

H5 now separates plant-side response structure from four increasingly independent mechanism layers. The plant data continue to support two partially separable response components — reproductive assurance and pollination-associated floral architecture — but none of the independent global mechanism tests identifies a common upstream pollination-service explanation.

The most direct new layer is GloPL experimental pollen supplementation. Unlike occurrence-based channel proxies, this directly measures reproductive pollen limitation. Its exact-island North/Tropical analysis is testable, but the frozen isolation-gradient prediction is not supported.

## Evidence hierarchy

1. **Plant-side response structure** — reproductive assurance and floral architecture are partially separable; an obligatory `isolation -> selfing -> floral simplification` pathway is rejected.
2. **Named floral architecture** — large-bee-like, butterfly-like and bird-like templates share >86% of their covariance. Orthogonalized named-architecture tests are support-limited because the same strict complete-template support does not robustly retain the original six-atomic H2 in Direct evidence.
3. **Global interaction structure** — GloBI breadth heterogeneity appears in 3/4 source definitions but is source-definition-sensitive and does not explain the area result.
4. **Exact-island functional channels** — no individual channel has adequate retained/disrupted overlap in both primary contexts; pooled attrition is testable but unsupported; identity-aware tropical signals are suggestive only; strict same-island identity is structurally non-identifiable in the North.
5. **Experimental pollen limitation (GloPL)** — exact-island overlap passes an outcome-blind admission gate, but experimental pollen limitation does not show the frozen northern isolation gradient or the predicted North–Tropical slope difference.

## Named-architecture boundary

The source-trained PCA1 explains 86.94% of all-analysis and 86.44% of Direct covariance among the three named floral templates. Under the strict >=50 scored-species/island gate:

- rank-2 template-specific residual H2: all p=.6143; Direct p=.3038;
- shared-factor H2 interaction: all p=.4324; Direct p=.2782.

However, the original six-atomic H2 on exactly that complete-template support is all-analysis p=.0431 but Direct p=.2719. Therefore the named-architecture null is **support-limited**, not a clean demonstration that the broad H2 lies outside pollination-syndrome geometry.

Provenance:

- structural correction run **35065912097**, artifact **10433982600**;
- common-support gate run **35083454858**, artifact **10440734711**.

## Functional-channel evidence

### Pooled channel attrition

The pooled five-channel exposure passes overlap but does not explain either plant route:

- North `selfing_core`: beta=-0.0113, p=.738;
- Tropical `selfing_core`: beta=+0.0588, p=.536;
- North `generalized_accessible | selfing_core`: beta=-0.00496, p=.834;
- Tropical: beta=-0.00863, p=.897.

### Identity-aware channel rows

Direct evidence gives a nominal tropical reproductive-assurance coefficient (beta=+0.21983, p=.0153), but q=.0611 and the all-analysis sensitivity does not reproduce it. The matched floral response has the expected tropical sign but remains unsupported.

### Strict same-island identity

North has only five mixed islands and, after island and channel-by-context fixed effects, the two context-specific disruption targets add only one independent design dimension. The North target is therefore structurally non-estimable. Lowering the support threshold cannot recover the missing identity-crossing variation.

## Experimental pollen limitation — GloPL

### Outcome-blind admission gate

The GloPL preflight read only latitude, longitude, DOI, author, year and species name. It did not read `PL_Effect_Size` or treatment outcomes.

Exact matching to the frozen 8,265 GSHHG islands found:

- 37 islands overall;
- North: **15 islands / 55 studies**;
- Tropical: **14 islands / 23 studies**.

Both contexts passed the frozen minimum of 10 islands and 5 studies, so the effect-size stage was opened without threshold relaxation.

Preflight run **35085378171**, artifact **10442126584**.

### Frozen effect-size test

The model was frozen before effect-size inspection at commit `0f7479c3607f86333c8d904a4786c0b938912b0b`. The response is the published GloPL master log response ratio of hand-pollen treatment versus natural reproductive output. Positive values mean a supplementation advantage and therefore stronger experimental pollen limitation at the measured reproductive endpoint.

The primary analysis contains:

- 262 exact-island finite-effect rows;
- 80 publication × island cells;
- 28 islands;
- 77 publications;
- 25 spatial blocks.

Primary estimates:

- North distance slope = **+0.3762**, one-sided positive p=.1209;
- Tropical-minus-North interaction = **-0.0311**, one-sided negative p=.4577;
- Tropical distance slope = +0.3451, two-sided p=.4774.

The northern point estimate is positive but unsupported. The context interaction is near zero and unsupported.

Frozen sensitivities do not rescue the claim:

- supplemental-only: North +0.4318; interaction **+0.3006** (opposite predicted context sign);
- no-zero-constant: North +0.3231; interaction -0.0893;
- equal-island: North +0.7473; interaction -0.2284.

None passes the frozen northern directional test, and the North–Tropical contrast is not directionally robust.

Canonical run **35086524129**, artifact **10442910503**, digest `sha256:5d39d3658100edabcbb9a7f6b5e6d7130ffc7d32e643e7900fcebfbd137efaa0`.

Classification: **`pollen_limitation_gradient_not_supported`**.

## Integrated interpretation

The H5 evidence now rules out several overly simple explanations without ruling out pollination biology itself.

- Floral architecture is not reducible to measured reproductive selfing.
- Raw named syndrome labels cannot be treated as realized visitor identity.
- Total functional-channel attrition does not explain both plant response routes.
- Identity-aware occurrence evidence is either suggestive but non-robust or structurally non-identifiable.
- Most importantly, independent experimental pollen limitation does **not** recover the frozen context-specific isolation gradient on the available exact-island subset.

Therefore the current defensible H5 statement is:

> **Island isolation is associated with partially separable reproductive-assurance and floral-architecture responses, but current independent global evidence does not support increasing pollen limitation with isolation as their common upstream generator.**

This is stronger than saying the mechanism is merely unmeasured, because GloPL directly measures experimental pollen limitation. It is still not evidence that pollinators are irrelevant: the GloPL subset is sparse, cross-sectional and does not measure abundance, visitation, historical loss or local selection gradients.

## Claim ceilings

Do not infer that:

- pollinators are irrelevant to island floral evolution;
- pollen limitation is absent on islands;
- pollinator abundance or visitation never declines with isolation;
- historical pollinator loss is absent;
- local pollination-mediated selection is absent;
- the Chapter 1 plant H2 response has no pollination contribution in any lineage or island system.

The global upstream claim that is specifically unsupported is the frozen proposition that experimental pollen limitation increases with mainland isolation in northern mid-latitude islands and does so more strongly than in tropical islands.

## Remaining identification target

A stronger H5 test would require island-level abundance, visitation or pollen-delivery measurements with repeated functional identities across isolation gradients and enough within-context identity crossing to separate service intensity from geography, lineage assembly and observation effort. The present global public data do not provide that design at scale.
