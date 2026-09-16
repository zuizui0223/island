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

## Figure 4 — H5: two plant pathways, then an independent upstream falsification wall

Figure 4 should distinguish what is supported on the plant side from what is not identified upstream. Do not draw `pollinator decline -> plant response` as a solid causal arrow.

### Panel A — two partially separable plant response components

Show two branches from isolation/source separation:

1. reproductive assurance / `selfing_core`;
2. pollination-associated floral architecture / `attraction_shift`.

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

### Panel C — exact-island Bombus upstream bridge fails identification

Use the canonical independent exact-island Bombus Search, not floral phenotype, to test the upstream arrow.

Canonical Bombus state counts among 7,154 source-available islands:

- detected/retained: 824;
- adequate non-detection/disrupted: 23;
- insufficient effort: 5,983;
- unresolved: 324.

Show the primary-context overlap matrix prominently:

| context | retained | disrupted |
|---|---:|---:|
| northern mid-latitude | 751 | **5** |
| tropical | **5** | 17 |

Frozen overlap gate: >=10 of each state in both contexts. **0/2 contexts pass.**

Beside it, show adjusted direct-evidence estimates only as diagnostics:

- `selfing_core`: beta(disrupted)=+0.047, p=0.722;
- `attraction_shift`: +0.047, p=0.292;
- `attraction_shift | selfing_core`: +0.038, p=0.231.

Post-genus residual checks are also unsupported:

- `generalized_form`: p=0.197;
- `self_compatibility`: p=0.428.

Panel claim:

> canonical occurrence data do not identify Bombus disruption as the upstream cause; reliable retained/disrupted states have inadequate within-context overlap, and adjusted associations are unsupported.

This is an identifiability failure plus negative adjusted evidence, not evidence that Bombus is irrelevant.

### Panel D — broader independent cross-examination

Global GloBI extension:

- 3,252 islands across all four contexts;
- 3/4 source definitions support four-context breadth-slope heterogeneity, but the fourth misses the frozen FDR gate;
- global distance×area heterogeneity: 0/4 supported;
- classification: source-definition-sensitive context heterogeneity, not a promoted global mechanism.

Other frozen checks:

- N1 channel heterogeneity: p = 0.65516;
- H5c biotic vs wind: p = 0.41221;
- H5d distributed threshold: 0/8 qualified.

Figure-4 synthesis:

> The data support two partially separable plant-side response components and regional pollination-syndrome consistency, but neither coarse interaction structure nor strict exact-island Bombus occurrence evidence identifies pollinator decline as the upstream cause.

## One-sentence synthesis

> **Global contemporary island floras show context-dependent floral and reproductive responses that are not erased by source-free genus residualization, whereas the narrower native Palearctic response loses support after source-matched genus adjustment; reproductive assurance and floral architecture are partially separable, but current independent interaction and exact-island Bombus data do not identify pollinator decline as their upstream cause.**

## Forbidden shortcuts

Do not write:

- the global H2 response is explained by native genus assembly;
- all-observed H3A residual proves within-lineage evolution;
- post-genus residual means a process literally operating below the genus level;
- unresolved/non-native-status records are introduced species;
- area proves island capacity or founder filtering;
- butterfly-like/bird-like/large-bee-like scores identify realized visitors;
- Bombus detection is effective pollination service or abundance;
- Bombus adequate non-detection proves historical extinction;
- the Bombus bridge proves pollinator decline causes selfing or floral change;
- failure of the Bombus/GloBI gates means pollinators are irrelevant.
