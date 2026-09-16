# Chapter 1 — H3 / H5 integration map (2026-09-16)

## Purpose

This note fixes the relationship between H3 and H5 after the global all-observed taxonomic-depth decomposition, the exact-island pollinator diagnostics, and the full-scale GloPL pollen-limitation analysis.

The central distinction is:

- **H3 asks where the observed plant response is represented taxonomically.**
- **H5 asks which ecological pressure or biological mechanism, if any, can explain the response.**

The new GloPL result changes H5 in an important way. H5 is no longer only a sequence of failed pollinator proxies: **experimental pollen limitation itself increases with geographic isolation at global scale.** What remains unresolved is how that general service limitation is translated into the context-dependent plant response, and whether that translation occurs through reproductive assurance, floral architecture, lineage assembly, or some combination.

H3 and H5 therefore form two orthogonal dimensions rather than one proved causal ladder.

## H3 — where is the plant response represented?

### H3A — broad contemporary flora

For the global all-observed North–Tropical H2 response, the exact common-support beta-binomial signal remains strong before decomposition:

- all evidence: 3,560 islands / 187 blocks, p = 8.39e-05;
- direct High/Medium: 3,487 islands / 187 blocks, p = 1.21e-04.

After source-free LOO taxonomic residualization, the North–Tropical difference remains supported:

- after family: p = 0.01043 all / 0.02497 direct;
- after genus: p = 0.01150 all / 0.004593 direct.

The clearest retained post-genus components are `generalized_form` and `self_compatibility`.

**Interpretation:** the broad contemporary-flora response is not erased by raw family/genus composition. This does not identify within-lineage evolution because introduced/unresolved floristic history, finer phylogenetic structure, genus-internal species sorting, environmental filtering and other unmeasured composition remain possible.

### H3B — defended native Palearctic assembly

The defended native Palearctic response behaves differently:

- observed -> family-adjusted -> source-matched genus-adjusted support: 4/4 -> 4/4 -> 0/4;
- broader frozen genus attenuation: approximately 78.8–85.9%;
- true genera outperform matched pseudo-genera: p = 0.0289855.

**Interpretation:** the defended native response is strongly represented in real genus composition. The exact family-to-genus increment remains spatially imprecise, but the qualitative contrast with H3A is clear.

## H5 — what ecological pressure is supported, and what bridge remains missing?

### H5A — two partially separable plant response routes

The plant data support two response components:

1. **Route A — reproductive assurance**
   - `selfing_core`, self-compatibility and related selfing/autonomous-assurance traits;
2. **Route B — pollination-associated floral architecture**
   - accessibility/generalization and pollination-syndrome-compatible floral architecture.

In the defended Palearctic analyses, attraction/access response remains associated with isolation after conditioning on `selfing_core` across four frozen source definitions (conditional distance approximately 0.091–0.101, q approximately 0.0079–0.034).

Tropical assemblages also show increasing reproductive assurance without a compulsory shift toward floral generalization. Therefore `isolation -> selfing -> floral simplification` is not an obligatory serial pathway.

### H5B — experimental pollen limitation provides a positive global ecological signal

The full georeferenced GloPL frame directly measures reproductive pollen limitation rather than visitor occurrence or floral phenotype.

Frozen full-scale support:

- 2,969 experiments;
- 1,248 sites;
- 919 publications;
- global standardized distance slope = **+0.07937**;
- frozen one-sided positive p = **0.01772**.

Classification: `global_pollen_limitation_gradient_supported_context_specificity_not_established`.

A post-hoc shape audit further shows a positive sampled offshore gradient:

- within-offshore slope = +0.19983, p = 0.01110;
- offshore-only slope = +0.23545, p = 0.005865.

These shape results are post hoc and do not upgrade the parent claim, but they make a pure mainland/offshore step an incomplete description of the sampled pattern.

**Interpretation:** geographic isolation is associated with increasing experimental pollen limitation at global scale. This is the strongest current H5 evidence that an ecologically relevant pollination-service pressure covaries with isolation.

### H5C — the service gradient does not reproduce H2 biogeographic branching

The predeclared North–Tropical GloPL comparison does not match the plant H2 pattern:

- northern slope = +0.06170, one-sided positive p = 0.1073;
- Tropical−North distance interaction = **+0.14300**;
- predeclared negative-interaction p = 0.9119;
- implied tropical slope = +0.20470.

Thus the point-estimated pollen-limitation gradient is stronger in the tropics, not weaker.

**Interpretation:** H2 cannot be reduced to `the North experiences stronger pollen limitation with isolation`. The global service pressure is real, but its biological translation into floral/reproductive composition is context dependent.

### H5D — direct functional bridges from GloPL to the two plant routes are not established

#### Route A — reproductive assurance buffering

Exact species-matched GloPL × trait tests:

- self compatibility: distance × SC = +0.0351, buffering p = 0.6790 — opposite the predicted buffering direction;
- selfing mating system: non-evaluable under the frozen support gate;
- autonomous selfing: distance × autonomous/delayed = −0.0736, p = 0.1543; direction compatible but unsupported.

Frozen family rule: **0/3 supported**.

Classification: `reproductive_assurance_buffering_not_supported`.

#### Route B — atomic floral-architecture buffering

Exact species-matched GloPL × atomic architecture tests:

- `generalized_form`: interaction = −0.06674, buffering p = 0.2441; direction compatible but imprecise and supplemental-only sensitivity reverses sign;
- `actinomorphic_symmetry`: interaction = +0.03496, p = 0.6690 — opposite the buffering prediction;
- `shallow_open_tube`: non-evaluable under the frozen support gate.

Frozen family rule: **0/3 supported**.

Classification: `floral_architecture_buffering_not_supported`.

**Interpretation:** pollen limitation is a supported general isolation-associated pressure, but current exact species-matched data do not establish either proposed functional bridge from that pressure to the six-atomic H2 response.

### H5E — interaction/occurrence data remain supporting cross-examination, not the main positive result

- GloBI source interaction breadth: context-sensitive, but not a robust global mechanism;
- individual exact-island channels: 0/5 have adequate retained/disrupted overlap in both primary contexts;
- pooled five-channel disruption: estimable but does not support the simple two-route model;
- identity-aware exact-island diagnostics: tropical hints are suggestive only and do not survive the full gate;
- strict same-island northern identity target is structurally non-identifiable.

These results narrow what GloPL means: the supported pollen-limitation gradient should not be translated into a claim of simple pollinator-channel loss, abundance decline, or a specific visitor identity.

## How H3 and H5 connect

The strongest current conceptual map is **common pressure, divergent biological translation, different taxonomic depth**.

```text
                         geographic isolation
                                |
                                v
                 experimental pollen limitation
                  increases globally (H5; supported)
                                |
                     translation not identified
                 ...............+...............
                 :                              :
                 v                              v
      Route A reproductive assurance    Route B floral architecture
                 \                              /
                  \                            /
                   v                          v
                  H2 context-dependent plant response
                                |
                 +--------------+--------------+
                 |                             |
                 v                             v
        H3A broad contemporary         H3B native Palearctic
        genus residual remains         genus adjustment removes support
```

The dotted middle connection is the missing bridge. GloPL supports the upstream ecological pressure; H3 localizes the downstream plant response; the current data do not identify the transformation connecting them.

## Broad contemporary flora: H3A x H5

H3A shows that the global context contrast, including `self_compatibility` and `generalized_form`, is not erased by source-free genus residualization.

GloPL independently shows that pollen limitation increases with isolation globally. But GloPL does **not** reproduce the North–Tropical branching, and species-matched buffering through either reproductive assurance or atomic accessible architecture is not supported.

Therefore the broad contemporary-flora result is best read as:

> **a broadly shared isolation-associated service pressure is translated into different plant responses across biogeographic contexts, and those differences cannot be explained solely by raw genus composition or by the two currently measured trait-buffering bridges.**

A complete mechanism for Panel A would still have to involve one or more of:

- finer taxonomic/phylogenetic structure;
- genus-internal species sorting;
- within-lineage change;
- introduced/non-native assembly;
- region-specific pollinator identity/turnover or compensation;
- service dimensions not captured by the present GloPL effect-size model;
- other ecological filters.

## Defended native Palearctic: H3B x H5

H3B shows that the defended native Palearctic response is strongly genus-structured.

The new GloPL result makes an assembly-mediated pollination interpretation biologically more plausible than before because a general isolation-associated pollen-limitation pressure is independently observed. But it still does not identify pollination as the generator of the genus structure.

The candidate chain is therefore:

```text
isolation-associated pollen limitation (general pressure; supported globally)
        -> lineage-specific demographic / reproductive consequences ?
        -> differential persistence / colonisation / establishment of genera ?
        -> genus-structured floral + reproductive composition (H3B; supported)
```

Only the first and last boxes are presently supported. The middle arrows remain mechanistic hypotheses.

## Route-specific H3 x H5 reading

| Plant component | H3A broad contemporary flora | H3B defended native Palearctic | H5 / GloPL interpretation |
|---|---|---|---|
| Reproductive assurance / selfing | `self_compatibility` remains in the post-genus residual contrast | native response overall is strongly genus-structured | Global pollen limitation increases with isolation, but reproductive-assurance buffering is 0/3 under the frozen family rule |
| Floral accessibility / architecture | `generalized_form` remains in the post-genus residual contrast | native attraction/access response is strongly genus-structured | Route B is partly independent of measured selfing, but atomic floral-architecture buffering is 0/3 under the frozen GloPL family rule |
| Regional translation | North–Tropical plant response is robust on broad support | Palearctic is the clearest defended native branch | GloPL North–Tropical branching is not supported; the tropical service slope is point-estimated stronger |
| Pollinator identity / turnover | not resolved by H3 | not resolved by H3 | occurrence/GloBI identity evidence remains insufficient or non-identifiable; GloPL measures service outcome, not visitor identity |

## Manuscript logic

The manuscript should now present H3 and H5 in this order:

1. **H3 localizes the response first.**
   - broad contemporary response is not erased by source-free genus residualization;
   - defended native Palearctic response loses support after source-matched genus adjustment.
2. **H5 establishes one general ecological pressure.**
   - experimental pollen limitation increases with geographic isolation globally.
3. **H5 then asks whether that pressure explains H2.**
   - North–Tropical service branching does not match H2;
   - reproductive-assurance buffering fails its frozen family rule;
   - atomic floral-architecture buffering fails its frozen family rule.
4. **Independent visitor/interactions data constrain interpretation.**
   - coarse channel attrition, source breadth and strict identity tests do not supply the missing bridge.
5. **No solid causal arrow is drawn from pollen limitation to H3 taxonomic structure.**
   - use a dashed candidate arrow labelled `service limitation -> lineage/trait translation?`.

## Compact interpretation

> **H3 shows that island floral/reproductive responses are represented at different taxonomic depths, whereas H5 now shows that pollen limitation itself increases with geographic isolation globally. The crucial unresolved step is the translation between that general ecological pressure and the context-dependent plant response: the service gradient does not reproduce North–Tropical H2 branching, and neither reproductive-assurance nor atomic floral-architecture buffering provides a supported species-level bridge. In the native Palearctic, any future mechanism must additionally explain why the response is carried largely by genus composition.**

## Claim ceilings

Do not write:

- GloPL proves pollinator abundance or visitation decline;
- the global GloPL slope is a causal island effect;
- the global GloPL gradient mediates H2;
- H5 explains H3B;
- genus structuring proves pollinator filtering;
- H3A post-genus residual proves within-lineage evolution;
- the two H5 routes are two independently proven causal pathways;
- unsupported functional buffering means reproductive assurance or floral architecture are biologically irrelevant;
- support-limited trait contrasts are biological nulls;
- post-hoc offshore shape results are confirmatory.
