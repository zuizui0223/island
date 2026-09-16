# Chapter 1 H5 — global pollination synthesis v5

This document supersedes `chapter1_h5_pollination_global_integration_v4_20260916.md` by adding the frozen GloPL × atomic floral-architecture functional-moderation test.

## Current evidence hierarchy

### 1. Plant response

Chapter 1 retains two partially separable response components:

1. reproductive assurance;
2. pollination-associated floral architecture.

The floral response is not reducible to measured `selfing_core`, so a compulsory `isolation → selfing → floral simplification` chain is already rejected.

### 2. Interaction and occurrence evidence

GloBI spans 3,252 islands and shows source-definition-sensitive geographic heterogeneity in interaction breadth, but it does not directly measure effective pollination service. Exact-island individual-channel occurrence lacks sufficient identity-crossing overlap for a global North–Tropical interpretation; the strict same-island northern target is structurally non-identifiable.

### 3. Experimental pollen limitation is globally associated with isolation

The full georeferenced GloPL frame provides a direct experimental measure of reproductive pollen limitation.

Frozen global-distance model:

- 2,969 experiments;
- 1,248 sites;
- 919 publications;
- global standardized distance slope = **+0.07937**;
- predeclared one-sided positive p = **0.01772**.

Classification: `global_pollen_limitation_gradient_supported_context_specificity_not_established`.

A post-hoc shape audit shows that the association is not well summarized as only a mainland/offshore step:

- mainland → offshore level step: +0.13745, p=.1084;
- two-part within-offshore slope: **+0.19983, p=.01110**;
- offshore-only slope: **+0.23545, p=.005865**.

The shape audit is descriptive post-hoc evidence and does not alter the parent promotion.

### 4. The plant H2 regional branching is not reproduced by GloPL

Frozen North–Tropical full-scale model:

- North slope = +0.06170, one-sided p=.1073;
- Tropical−North distance interaction = **+0.14300**;
- predeclared negative-interaction p=.9119;
- implied tropical slope = +0.20470.

Thus general isolation-associated pollen limitation is supported, but the service gradient does not reproduce the proposed stronger northern pattern.

### 5. Route A — reproductive-assurance buffering is not supported

Frozen species-level test: `chapter1_h5_glopl_reproductive_assurance_moderation_v1`.

- design freeze: `baf5b7994171e31c01d91e9a4c936408bc4db6cf`;
- canonical run: **35093274622**;
- artifact: **10444749163**;
- exact matched GloPL species across the three traits: 642.

Results:

- self compatibility: distance × SC = **+0.0351**, one-sided buffering p=.6790 — opposite the prediction;
- selfing mating system: not evaluable; selfing class has only 11 matched species and one offshore site;
- autonomous selfing: distance × autonomous/delayed = **−0.0736**, p=.1543; both frozen sensitivity interactions remain negative but the primary effect is unsupported.

Frozen family rule: at least two supported traits. Observed support: **0/3**.

Classification: `reproductive_assurance_buffering_not_supported`.

### 6. Route B — atomic floral-architecture buffering is not supported

Frozen species-level test: `chapter1_h5_glopl_floral_architecture_moderation_v1`.

- design freeze: `4f478c76b4ac0ed963a0b7d0ae92a8b143e72283`;
- RED run: `35094224454`;
- canonical run: **35094588521**;
- artifact: **10445189257**;
- digest: `sha256:ef81551d41cc747b07e927b9a5c4f13eaf0fc013a5cf8ac3fb1a5491aadec8c4`;
- exact matched architecture overlap: 624 GloPL species.

#### Generalized form

Outcome-blind support:

- restricted architecture: 169 species / 245 sites / 174 publications / 56 offshore sites;
- generalized architecture: 77 species / 114 sites / 79 publications / 24 offshore sites.

Primary model:

- restricted distance slope = +0.07409;
- generalized distance slope = +0.00735;
- distance × generalized interaction = **−0.06674**;
- one-sided buffering p=.2441.

The primary point estimate is directionally compatible with buffering, but it is imprecise and the supplemental-only sensitivity reverses sign to **+0.06455**. No-zero-constant remains negative at −0.07565.

Decision: unsupported; sensitivity direction not retained.

#### Actinomorphic symmetry

Outcome-blind support:

- zygomorphic: 213 species / 263 sites / 180 publications / 89 offshore sites;
- actinomorphic: 369 species / 434 sites / 315 publications / 117 offshore sites.

Primary model:

- zygomorphic distance slope = +0.08518;
- actinomorphic distance slope = +0.12013;
- distance × actinomorphic interaction = **+0.03496**;
- one-sided buffering p=.6690.

This is opposite the buffering prediction. The supplemental-only interaction is slightly negative, but the no-zero-constant interaction is again positive.

The admitted North–Tropical secondary diagnostic is also unsupported:

- North interaction = +0.04181;
- Tropical interaction = +0.12071;
- Tropical−North = +0.07891, p=.7130.

Decision: unsupported; sensitivity direction not retained.

#### Shallow/open tube

Only 52 exact-matched species were available. The deep class had 4 offshore sites and the open/shallow class had 7, so the frozen global support gate failed. This is non-evaluable, not a biological null.

Frozen Route-B family rule: at least two supported atomic contrasts. Observed support: **0/3**.

Classification: `floral_architecture_buffering_not_supported`.

## Integrated H5 result

The current evidence now separates four claims:

1. **General experimental service limitation:** supported. Pollen limitation increases with geographic separation from major continental landmasses at full GloPL scale.
2. **North–Tropical branching of service limitation:** not supported. The tropical point-estimated distance gradient is stronger, not weaker.
3. **Functional buffering by reproductive assurance:** not supported under the frozen Route-A family rule.
4. **Functional buffering by atomic generally accessible floral architecture:** not supported under the frozen Route-B family rule.

The evidence hierarchy is therefore:

```text
isolation
  |
  +--> experimental pollen limitation increases globally
          |
          +--> offshore-only gradient also positive (post-hoc)
          |
          +--> North-Tropical service branching: not supported
          |
          +--> Route A reproductive-assurance buffering: 0/3 supported
          |
          +--> Route B atomic architecture buffering: 0/3 supported
          |
          `--> no identified functional bridge to six-atomic H2
```

## Current H5 statement

> **Geographic isolation is associated with increasing experimental pollen limitation at global scale, including across sampled offshore sites, but the current independent data do not connect that general service gradient to the context-dependent Chapter 1 plant response through either measured reproductive-assurance buffering or the three predeclared atomic generally accessible floral-architecture contrasts. The predeclared North–Tropical service branching is also unsupported. Pollination-service limitation is therefore a supported general correlate of geographic isolation, not an identified common upstream generator of the Chapter 1 H2 branching.**

## What this does not mean

The result does not show that pollinators are irrelevant, that floral architecture has no local fitness effect, that reproductive assurance has no evolutionary role, or that island pollination-mediated selection is absent. The functional tests are specific to the current exact-matched GloPL species and predeclared trait contrasts; shallow/open tube depth and explicit selfing mating system are notably support-limited.

## Stop rule

Further post-hoc GloPL trait-moderation searches should not be used to hunt for a bridge after both predeclared functional families failed promotion. With current public data, the remaining H5 uncertainty concerns region-specific pollinator-community composition, assembly history, or unmeasured service dimensions rather than another unfrozen trait interaction.

## Claim ceiling

The current evidence does not establish:

- causal island effects on pollen limitation;
- pollinator abundance or visitation decline;
- historical pollinator loss;
- causal evolution of reproductive or floral traits;
- causal selection by pollen limitation;
- mediation of Chapter 1 H2 by GloPL pollen limitation;
- absence of local pollinator-mediated selection;
- biological nulls for support-limited trait contrasts.
