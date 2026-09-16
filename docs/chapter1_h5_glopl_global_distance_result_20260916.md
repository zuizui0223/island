# Chapter 1 H5 — all-site GloPL global distance result

## Why this analysis was added

The exact-island GloPL test was deliberately strict: it asked only about GloPL experiments whose coordinates fall inside the frozen Chapter 1 island polygons. That produced a defensible but small adjusted subset (28 model islands). The present extension instead uses every georeferenced GloPL experiment and places mainland and offshore sites on one continuous geographic-isolation axis.

This extension was frozen before the global distance exposure was merged to pollen-limitation outcomes at commit `6f570b3aa97b66a8012f55953fdf8c1dd4aba953`.

## Geography and outcome-blind support

The distance definition reuses the six seeded major continental landmasses already used by the Chapter 1 `purpose_shortest_analysis` route. Mainland sites on one of those landmasses have distance zero; offshore sites have positive distance. No post-hoc near/remote threshold is used.

Outcome-blind preflight:

- GloPL rows with coordinates: **2,969 / 2,969**
- unique coordinate sites: **1,248**
- unique publications/study keys: **919**
- zero-distance sites: **988**
- positive-distance sites: **260**

Context support:

| Context | Sites | Studies | Zero-distance | Positive-distance |
|---|---:|---:|---:|---:|
| Northern mid-latitude | 737 | 526 | 609 | 128 |
| Tropical | 263 | 252 | 218 | 45 |
| Southern extratropical | 248 | 146 | 161 | 87 |
| Northern high-latitude | 0 | 0 | 0 | 0 |

The frozen global and North–Tropical gates passed. The frozen four-context gate failed before effect-size inspection because GloPL has no northern-high-latitude sites under the Chapter 1 latitude definition.

## Primary all-site result

All 2,969 master pollen-limitation effect rows were finite. After the frozen publication × coordinate × measurement-cell aggregation:

- **1,408** measurement cells
- **1,248** sites
- **919** publications
- each publication has total analysis weight 1
- uncertainty is publication-cluster robust
- measurement type/treatment/constant/supplementation categories enter as fixed effects.

The frozen continuous distance model gives:

- global standardized distance slope = **+0.07937**
- SE = **0.03773**
- two-sided p = **0.03543**
- predeclared one-sided positive p = **0.01772**

Under the frozen promotion rule this is a supported global association between distance from the major continents and experimental pollen limitation.

### Frozen sensitivities

**Supplemental-hand-pollination only**

- slope = **+0.03349**
- one-sided p = **0.2149**

The direction remains positive, satisfying the frozen directional sensitivity rule, but this sensitivity is not independently statistically supported.

**No zero-constant cases**

- slope = **+0.07789**
- one-sided p = **0.02095**

This closely reproduces the primary global result.

## North–Tropical mechanism prediction

The broad global association does **not** reproduce the Chapter 1 H2 mechanism prediction.

Primary North–Tropical model:

- Northern mid-latitude slope = **+0.06170**, one-sided p = **0.1073**
- Tropical − Northern interaction = **+0.14300**
- predeclared one-sided negative interaction p = **0.9119**
- implied tropical slope = **+0.20470**

Thus the tropical gradient is point-estimated as stronger, not weaker, than the northern-midlatitude gradient. Supplemental-only strengthens this mismatch: North is approximately zero/slightly negative while the tropical-minus-North interaction is strongly positive in direction.

Therefore the all-site result is classified as:

`global_pollen_limitation_gradient_supported_context_specificity_not_established`

## Relationship to the exact-island result

The two analyses answer different questions and should both remain visible.

- **Exact-island subset:** asks whether pollen limitation varies with isolation among identified Chapter 1 islands while controlling island area and climate. It had only 28 model islands and the North slope was positive but statistically unsupported.
- **All-site global analysis:** uses the full GloPL geographic frame and gains much more information, but does not have harmonized island-area or climate controls across mainland and offshore sites.

The all-site analysis therefore supplies evidence for a broad geographic association, while the exact-island analysis remains the stricter island-specific adjusted test.

## What this changes in H5

The previous H5 statement that independent global experimental data do not support an isolation-associated pollen-limitation gradient is no longer correct at the full GloPL sampling scale.

The defensible updated statement is:

> **Across the full GloPL sampling frame, experimental pollen limitation increases with geographic distance from the major continental landmasses. This broad association does not reproduce the predeclared North–Tropical branching expected if one context-specific pollination-service mechanism generated the Chapter 1 plant H2 response.**

This result therefore strengthens the plausibility of a general isolation–pollination-service link while leaving the specific upstream explanation of the floral/reproductive H2 branching unresolved.

## Remaining shape question

Because **988/1,248 sites have distance zero** and only **260** have positive distance, the supported continuous coefficient could arise primarily from a mainland-versus-offshore step rather than progressively stronger pollen limitation among increasingly remote offshore sites.

A separate post-hoc shape diagnostic will decompose:

1. mainland versus offshore step;
2. continuous distance slope within positive-distance offshore sites.

That diagnostic is not allowed to promote the parent result or alter the frozen parent p-value.

## Provenance

- RED run: **35090344592**
- canonical GREEN run: **35090599662**
- job: **104775735096**
- artifact: **10444156159**
- artifact digest: `sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623`
- GloPL SHA-256: `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`
- Natural Earth SHA-256: `9e0729ee253ca7d7a5c4ae9395fb1902264c5377c52e224d13dd85010e2835d9`

## Claim ceiling

This analysis does not establish pollinator abundance decline, visitation decline, a historical island effect, floral selection, or mediation of the Chapter 1 plant H2 response. It also does not establish the predeclared North–Tropical pollination-service mechanism.
