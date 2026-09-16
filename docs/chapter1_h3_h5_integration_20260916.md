# Chapter 1 — H3 / H5 integration map (2026-09-16)

## Purpose

This note fixes the relationship between H3 and H5 after the global all-observed taxonomic-depth decomposition and the multichannel/identity-aware pollinator diagnostics.

The central distinction is:

- **H3 asks where the observed response is represented taxonomically.**
- **H5 asks which biological mechanism, if any, can explain the response.**

H3 is therefore a localization/decomposition question. H5 is a mechanism-identification question. H5 must not be used to retroactively relabel an H3 result as causal, and H3 must not be treated as evidence for a specific pollinator mechanism.

## H3 has two inferential layers

### H3A — broad contemporary flora

For the global all-observed North–Tropical H2 response, the exact common-support beta-binomial signal remains strong before decomposition:

- all evidence: 3,560 islands / 187 blocks, p = 8.39e-05;
- direct High/Medium: 3,487 islands / 187 blocks, p = 1.21e-04.

After source-free LOO taxonomic residualization, the North–Tropical difference remains supported:

- after family: p = 0.01043 all / 0.02497 direct;
- after genus: p = 0.01150 all / 0.004593 direct.

The clearest retained post-genus components are `generalized_form` and `self_compatibility`.

**Interpretation:** the broad contemporary-flora response is not erased by raw family/genus composition. This does not identify within-lineage evolution, because introduced/unresolved floristic history, deeper or shallower phylogenetic structure, environmental filtering and other unmeasured composition remain possible.

### H3B — defended native Palearctic assembly

The defended native Palearctic response behaves differently:

- observed -> family-adjusted -> source-matched genus-adjusted support: 4/4 -> 4/4 -> 0/4;
- broader frozen genus attenuation: approximately 78.8–85.9%;
- true genera outperform matched pseudo-genera: p = 0.0289855.

**Interpretation:** the defended native response is strongly represented in real genus composition. The exact family-to-genus increment remains spatially imprecise, but the qualitative contrast with H3A is clear.

## H5 contains a plant-side decomposition and an upstream mechanism test

### H5A — two partially separable plant response routes

The plant data support two response components:

1. **Route A — reproductive assurance**
   - `selfing_core`, self-compatibility and related selfing/autonomous-assurance traits;
2. **Route B — pollination-associated floral architecture**
   - accessibility/generalization and named syndrome-compatible floral architecture.

The key separation evidence is that, in the defended Palearctic analyses, the attraction/access response remains associated with isolation after conditioning on `selfing_core` across four frozen source definitions (conditional distance approximately 0.091–0.101, q approximately 0.0079–0.034).

Tropical assemblages also show increasing reproductive assurance without a compulsory shift toward floral generalization. Therefore `isolation -> selfing -> floral simplification` is not an obligatory serial pathway.

### H5B — syndrome consistency

Plant-side syndrome templates provide biological interpretation but not realized visitor identity:

- northern large-bee-like architecture declines with isolation;
- tropical butterfly-like architecture increases with isolation;
- tropical warm-colour x tubular architecture is enriched relative to northern mid-latitudes;
- >86% of covariance among large-bee-like, butterfly-like and bird-like templates is shared architecture.

Thus syndrome scores identify pollination-associated architecture, not actual visitors.

### H5C — independent upstream pollinator tests

The independent tests ask whether reduced pollination opportunity explains either Route A or Route B.

1. **GloBI source interaction breadth** — context-sensitive but not a robust global mechanism.
2. **single exact-island channels** — 0/5 channels have adequate retained/disrupted overlap in both primary contexts.
3. **pooled five-channel disruption** — overlap becomes adequate, but the simple two-route mechanism is unsupported:
   - Route A `selfing_core`: North beta = -0.0113, p = 0.738; Tropical beta = +0.0588, p = 0.536;
   - Route B `generalized_accessible | selfing_core`: North beta = -0.00496, p = 0.834; Tropical beta = -0.00863, p = 0.897;
   - four-test FDR q approximately 0.897.
4. **identity-aware channel diagnostic** — retains functional identity rather than treating all channels as interchangeable:
   - Tropical Route A: beta = +0.21983, p = 0.0153, q = 0.0611, direct-only and not reproduced in all-analysis;
   - Tropical Route B matched architecture after selfing: beta = -0.13820, p = 0.1215, q = 0.1668, expected direction but unsupported.

**Interpretation:** simple total pollinator-channel attrition is not supported as the common upstream cause. Identity/turnover and service intensity remain plausible targets, but no causal pollinator mechanism is identified.

## How H3 and H5 connect

The relationship is not `H3 -> H5` or `H5 -> H3` as a proved causal chain. It is a two-dimensional decomposition:

```text
                     H3: WHERE is the response represented?
                     ----------------------------------------
                     broad contemporary        native Palearctic
                     genus residual remains    genus adjustment removes support
                     (H3A)                     (H3B)

H5: WHY?      Route A reproductive assurance
mechanism     Route B floral architecture
candidate     independent pollinator tests
```

### Broad contemporary flora: H3A x H5

H3A shows that the global context contrast, including `self_compatibility` and `generalized_form`, is not erased by source-free genus residualization.

Therefore H5 cannot explain the broad H2 signal merely by saying that different pollinator environments select different genera. A pollination mechanism for Panel A, if real, would have to operate through one or more of:

- finer-than-raw-genus taxonomic/phylogenetic structure;
- within-lineage change;
- introduced/non-native assembly;
- unmeasured lineage composition;
- direct ecological filtering not captured by raw genus means.

Current H5 occurrence data do not distinguish among these possibilities.

### Defended native Palearctic: H3B x H5

H3B shows that the defended native Palearctic response is strongly genus-structured. This changes the interpretation of H5.

If pollination contributes to the Palearctic response, the present data are more consistent with an **assembly-mediated route** than with a demonstrated repeated within-species floral adaptation:

```text
pollination environment (candidate, not identified)
        -> differential persistence / colonisation / establishment of genera
        -> genus-structured floral + reproductive composition
```

However, H5 does not identify pollination as the upstream generator of that genus structure. The genus pattern could also reflect other lineage-associated ecological or historical filters.

## Route-specific H3 x H5 reading

| Plant component | H3A broad contemporary flora | H3B defended native Palearctic | H5 interpretation |
|---|---|---|---|
| Reproductive assurance / selfing | `self_compatibility` remains in the post-genus residual contrast | native response overall is strongly genus-structured | Route A exists, but current channel-disruption data do not identify pollinator reduction as its common cause |
| Floral accessibility / architecture | `generalized_form` remains in the post-genus residual contrast | native attraction/access response is strongly genus-structured | Route B is partly independent of measured selfing; syndrome-compatible, but no independent pollinator cause is identified |
| Pollinator identity / turnover | not resolved by H3 | not resolved by H3 | identity-aware tropical diagnostic is suggestive only; effective service remains unmeasured |

## Manuscript logic

The manuscript should present H3 and H5 in this order:

1. **H3 localizes the response first.**
   - Panel A: broad response is not erased by source-free genus residualization.
   - Panel B: defended native Palearctic response loses support after source-matched genus adjustment.
2. **H5 then asks whether pollination can explain either localization pattern.**
   - two plant routes are supported as partially separable;
   - syndrome consistency is present;
   - simple total channel loss fails;
   - identity-aware tropical evidence is suggestive but not confirmatory.
3. **No solid causal arrow is drawn from pollinator disruption to H3 taxonomic structure.**
   - at most use a dashed candidate arrow labelled `pollination-service / identity filter?`.

## Compact interpretation

> **H3 tells us that the broad contemporary response and the defended native response live at different taxonomic depths; H5 shows that both contain reproductive-assurance and pollination-associated floral components, but current independent pollinator data do not identify a common upstream pollinator-loss mechanism. In the native Palearctic, any future pollination mechanism would have to explain why the response is carried largely by genus composition; in the broad contemporary flora, it would additionally have to explain why a context contrast remains after raw genus composition is removed.**

## Claim ceilings

Do not write:

- H5 explains H3B;
- genus structuring proves pollinator filtering;
- H3A post-genus residual proves within-lineage evolution;
- the two H5 routes are two independently proven causal pathways;
- total channel attrition is the measured pollinator decline process;
- identity-aware tropical signals identify effective pollination service.
