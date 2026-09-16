# Chapter 1 v12 two-panel figure architecture — 2026-09-16

## Decision

The main result is not one H1→H5 causal ladder. The manuscript/figure surface is organized as two nested evidence layers, with the taxonomic-depth contrast shown explicitly side-by-side.

The numbering is fixed as:

- H1 — universal response rival
- H2 — biogeographic branching
- H3 — taxonomic representation depth
- H4 — area moderation
- H5 — independent mechanism gate

## Figure 1 — Global contemporary-flora response

### Panel A — one universal response is inadequate

Show the four-context six-atomic isolation-response vectors from the all-observed beta-binomial model.

Headline evidence:

- global four-context heterogeneity: rejected universal response;
- predeclared North–Tropical contrast, all evidence: 3,626 islands / 187 blocks, p = 0.000763;
- direct High/Medium: 3,569 islands / 187 blocks, p = 0.004502;
- matched grouped-binomial sensitivity: p = 0.001385 / 0.004501.

Visual rule: show six atomic coefficients, not the historical two species-level composite scores.

Claim: contemporary observed island floras show biogeographically contingent floral/reproductive responses to isolation.

Do not claim native assembly from this panel.

### Panel B — floristic-status defensibility gate

Show immediately why the broad response cannot be read as a native-assembly result.

On the same 375 islands:

- all observed species: p = 0.0005607 all evidence / 0.0003959 direct;
- native species only: p = 0.2491 / 0.3018.

Expanded WCVP regional-native-compatible sensitivity:

- 2,372 islands total;
- North–Tropical p = 0.07954 all / 0.17472 direct.

Claim: the global contemporary-flora response is broader than the portion currently defensible as native assembly.

## Figure 2 — H3 taxonomic-depth contrast: the two panels behave differently

This is the critical comparison figure. Put H3A and H3B on the same page, on the same visual grammar, so that the asymmetry is visible without reading the caption.

### Panel A — H3A broad all-observed response is not erased by source-free genus residualization

Show three aligned stages for the North–Tropical vector on the exact common taxonomic support:

1. observed;
2. after source-free LOO family residualization;
3. after source-free LOO genus residualization.

First show the exact-model gate:

- same beta-binomial H2 model on H3 common support:
  - all evidence: 3,560 islands / 187 blocks, p = 8.39e-05;
  - direct: 3,487 islands / 187 blocks, p = 1.21e-04.

Then show the equal-island decomposition:

- all evidence: after-family p = 0.01043; after-genus p = 0.01150;
- direct: after-family p = 0.02497; after-genus p = 0.004593.

Do not headline the point attenuation percentages because paired-block intervals are broad. If attenuation is shown, display interval uncertainty explicitly.

Key retained post-genus components:

- generalized_form: p = 0.00363 all / 0.00327 direct;
- self_compatibility: p = 0.0212 all / 0.00768 direct.

Preferred wording: **not erased by source-free genus residualization**.

Avoid using `below genus` as shorthand in the panel title because that can be misread as a within-lineage evolutionary claim.

Claim: broad context dependence is not simply a raw family/genus-composition effect and is not the same signal as the defended native genus-assembly result.

### Panel B — H3B defended Palearctic native response collapses after source-matched genus adjustment

Use the same three-stage visual grammar as Panel A:

1. observed;
2. family-adjusted;
3. source-matched genus-adjusted.

Retain frozen v11/P1 evidence:

- observed → family-adjusted → genus-adjusted support: 4/4 → 4/4 → 0/4;
- broader frozen genus attenuation: approximately 78.8–85.9%;
- matched-complexity pseudo-genus randomization: p = 0.0289855;
- paired spatial-block uncertainty keeps the exact family→genus increment imprecise.

Claim: the defended Palearctic native response is strongly genus-structured beyond arbitrary grouping complexity.

Do not use this panel as the explanation of Figure 1.

### Panel C — explicit comparison statement

Put a compact comparison strip or bracket under Panels A and B:

- **all-observed / broad H2:** residual vector remains supported after source-free genus residualization;
- **defended native / Palearctic:** response loses support after source-matched genus adjustment.

Caption statement:

> Taxonomy matters in both layers, but it matters at different depths and under different inferential conditions. The broad contemporary-flora response is not explained by the defended native genus-assembly result.

This panel is the visual reason the manuscript needs two evidence layers rather than one causal ladder.

## Figure 3 — H4 area moderation is secondary and measurement-sensitive

Place area after H3 rather than between H2 and taxonomic decomposition.

Show the supported beta-binomial distance×area patterns together with equal-island/common-support guardrails. The interpretation remains:

> area is a measurement-sensitive modifier, not established evidence for island capacity, founder filtering, or pollinator-population persistence.

This panel must not be drawn as the causal bridge from global H2 to native H3B.

## Figure 4 — H5: two plant pathways, then progressively stricter upstream tests

Figure 4 should distinguish what is supported on the plant side from what remains unidentified upstream. Do not draw `pollinator decline -> plant response` as a solid causal arrow.

### Panel A — two partially separable plant response components

Show two branches from isolation/source separation:

1. reproductive assurance / `selfing_core`;
2. pollination-associated floral architecture / attraction-accessibility.

Use the Palearctic conditional decomposition as the key separation evidence:

- attraction/access response remains positive after conditioning on `selfing_core` across four frozen source definitions;
- all-native conditional distance estimates approximately 0.091–0.101, q=0.0079–0.034.

Add the tropical counterexample to a compulsory serial selfing-syndrome model:

- reproductive assurance increases while accessibility/generalization declines in direct all-native data.

Panel claim:

> reproductive assurance and floral architecture are partially separable plant response components; `isolation -> selfing -> floral simplification` is not an obligatory serial pathway.

Do not label the floral branch `direct pollinator selection`; the upstream driver is not identified.

### Panel B — syndrome consistency without visitor identity

- Northern direct large-bee-like concordance decreases with isolation;
- Tropical butterfly-like concordance increases with isolation;
- tropical warm-colour × tubular architecture is strongly enriched relative to northern mid-latitudes;
- named-template common architecture factor explains >86% of covariance.

Panel claim: the regional floral responses are compatible with different pollination-associated architectures, but named templates do not identify realized visitor guilds.

### Panel C — total channel attrition vs identity-aware disruption

This panel should explicitly compare two upstream models rather than stopping at the pooled null.

#### C1. Pooled five-channel attrition

Pool strict exact-island states for Bombus, non-Bombus bees, Lepidoptera, flower-visiting birds and Diptera. Retained and disrupted states use the same frozen effort gate.

Primary exposure: any documented disruption among at least two evaluable channels.

Support:

| context | islands | disrupted | no documented disruption |
|---|---:|---:|---:|
| northern mid-latitude | 303 | **15** | 288 |
| tropical | 99 | **21** | 78 |

Route A (`selfing_core`):

- North: beta = -0.0113, p=0.738;
- Tropical: beta = +0.0588, p=0.536.

Route B (`generalized_accessible | selfing_core`):

- North: beta = -0.00496, p=0.834;
- Tropical: beta = -0.00863, p=0.897.

Four-test FDR q ≈ **0.897**.

Visual label: **total functional-channel attrition — not supported**.

#### C2. Identity-aware channel stacking

Now keep island × channel rows rather than collapsing all visitor guilds. Give each island equal total weight and absorb channel × context identity with fixed effects.

Route A uses all five channels and predicts higher `selfing_core` with disruption:

- North: beta = +0.06735, p=0.356, q=0.356;
- Tropical: beta = **+0.21983, p=0.0153, q=0.0611**.

Route B uses predeclared channel-matched architectures and conditions on `selfing_core`:

- Bombus -> large-bee-like;
- Lepidoptera -> butterfly-like;
- flower-visiting birds -> bird-like.

Results:

- North: beta = +0.06469, p=0.125, q=0.1668;
- Tropical: beta = **-0.13820, p=0.1215, q=0.1668**.

The tropical Route A effect is nominal/direct-only and disappears in all-analysis sensitivity. The tropical Route B sign matches the identity-specific prediction but is unsupported in both evidence scopes.

Visual label: **identity-aware turnover — suggestive tropical pattern, not promoted**.

Panel-C claim:

> a simple count of lost pollinator channels is too coarse. Preserving functional identity reveals a possible tropical reproductive-assurance signal and the expected direction of matched floral change, but neither survives the full inferential gate.

### Panel D — remaining falsification wall and the unmeasured quantity

Use a compact cross-examination strip:

- individual-channel overlap: **0/5** channels have adequate retained/disrupted overlap in both primary contexts;
- Bombus-only bridge: overlap failure + unsupported adjusted associations;
- GloBI global context heterogeneity: 3/4 source definitions, not robustly promoted;
- GloBI distance×area heterogeneity: 0/4;
- N1 channel heterogeneity: p = 0.65516;
- H5c biotic vs wind: p = 0.41221;
- H5d distributed threshold: 0/8 qualified.

Beside the strip, show the remaining unmeasured target:

`pollination-service limitation = abundance × visitation × pollen delivery × effectiveness`

with modifiers:

`functional identity / turnover | compensation among channels | lineage assembly | local selection`

Figure-4 synthesis:

> The data support two partially separable plant-side response components and regional pollination-syndrome consistency. Total channel attrition does not explain them. Identity-aware disruption gives a suggestive tropical signal but still does not identify causation. The remaining mechanistic target is effective pollination-service limitation, not coarse channel presence alone.

## One-sentence synthesis

> **Global contemporary island floras show context-dependent floral and reproductive responses that are not erased by source-free genus residualization, whereas the narrower native Palearctic response loses support after source-matched genus adjustment; reproductive assurance and floral architecture are partially separable, total pollinator-channel attrition is unsupported, and identity-aware disruption is suggestive in the tropics but does not yet identify pollination-service limitation as the common upstream cause.**

## Forbidden shortcuts

Do not write:

- the global H2 response is explained by native genus assembly;
- all-observed H3A residual proves within-lineage evolution;
- post-genus residual means a process literally operating below the genus level;
- unresolved/non-native-status records are introduced species;
- area proves island capacity or founder filtering;
- butterfly-like/bird-like/large-bee-like scores identify realized visitors;
- channel detection is effective pollination service or abundance;
- adequate non-detection proves historical extinction or temporal decline;
- the pooled or identity-aware bridge proves pollinator decline causes selfing or floral change;
- the nominal tropical identity-aware result is confirmatory evidence;
- failure of pooled-channel/GloBI gates means pollinators are irrelevant.
