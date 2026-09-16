# Chapter 1 H5 — global pollination integration

## Current decision

H5 now has four evidence levels.

1. **Plant-side consistency:** floral trait combinations are compatible with changes in pollination-associated architecture.
2. **Independent interaction structure:** GloBI breadth and coarse channel tests ask whether interaction structure covaries with the plant pattern.
3. **Independent exact-island single-channel states:** Bombus and the other functional channels test whether a particular source-available channel is retained or strictly disrupted.
4. **Independent pooled channel-disruption bridge:** Bombus, non-Bombus bees, Lepidoptera, flower-visiting birds and Diptera are pooled into a conservative source-to-island functional-channel disruption exposure and tested against the two plant response routes.

The fourth level directly addresses the broader `pollinator reduction/disruption -> two plant pathways` hypothesis without privileging Bombus. It is estimable in both primary contexts, but the predicted two-route associations are not supported.

## Plant-side consistency retained

The predeclared direct-evidence predictions remain supported on native data:

- northern mid-latitude: large-bee-like floral architecture declines with isolation;
- tropical: butterfly-like floral architecture increases/does not decline with isolation.

More than 86% of covariance among large-bee-like, butterfly-like and bird-like templates is shared plant architecture, so these named templates cannot be treated as realized pollinator identities.

The reproductive and floral components are partially separable. In the defended Palearctic analyses, attraction/access response remains associated with distance after conditioning on `selfing_core` across four source definitions. Tropical assemblages also show increasing reproductive assurance without the compulsory floral-generalization direction. Thus the plant data reject an obligatory serial `isolation -> selfing -> floral simplification` model.

## Independent GloBI evidence is global but not mechanism-identifying

The frozen GloBI genus-breadth predictor was applied to northern mid-latitude, northern high-latitude, tropical and southern extratropical island floras.

Primary all-observed coverage: **3,252 islands**.

The four-context distance-slope heterogeneity test was FDR-supported in **3/4** frozen source definitions and narrowly missed in the fourth. The result is therefore **source-definition-sensitive context heterogeneity**, not a promoted mechanism.

Descriptively:

- northern mid-latitude: weak positive breadth slopes;
- northern high-latitude: mixed near zero;
- tropical: consistently negative breadth slopes;
- southern extratropical: positive breadth slopes, with one source definition supported.

The global GloBI `distance x area` heterogeneity test was supported in **0/4** source definitions, so sampled source-genus interaction breadth does not explain the plant-side area pattern.

## Exact-island single-channel states

The canonical exact-island Search campaign provides outcome-blind occurrence evidence for five pollination functional channels on the frozen global island universe. Missing records are never automatically coded as absence.

For Bombus, among 7,154 source-available islands the canonical Search recorded:

- detected: **824**;
- adequate non-detection: **23**;
- insufficient effort: **5,983**;
- unresolved: **324**.

The Bombus-specific upstream bridge (`chapter1_h5_bombus_upstream_bridge_v1`; run **35049673581**, artifact **10428646713**) failed its within-context overlap gate and showed no supported adjusted association with either `selfing_core` or the floral attraction/access pathway.

A subsequent five-channel overlap audit showed that **0/5 individual channels** had at least 10 retained and 10 disrupted islands in both northern-midlatitude and tropical contexts. Thus no single channel provides a globally balanced retained-versus-disrupted comparison for the primary context contrast.

## New pooled five-channel two-route bridge

The broader mechanism hypothesis concerns pollinator opportunity as a whole rather than Bombus specifically. We therefore pooled:

- Bombus;
- non-Bombus bees;
- Lepidoptera;
- flower-visiting birds;
- Diptera.

A channel counts as `disrupted` only when it is source-available and has an effort-qualified adequate non-detection. To make the retained state comparable, a detected channel enters the retained-versus-disrupted comparison only when the **same background-effort gate** is satisfied. Insufficient-effort and unresolved states remain missing.

The primary exposure is whether an island has **any documented channel disruption among at least two effort-qualified evaluable source-available channels**. Models additionally control the evaluable-channel composition, source separation, island area, climate PC1-4 and spatial block.

This pooled exposure solves the main identifiability problem of the single-channel analyses:

| context | islands | disrupted | no documented disruption |
|---|---:|---:|---:|
| northern mid-latitude | 303 | **15** | 288 |
| tropical | 99 | **21** | 78 |

The frozen overlap reference is therefore passed in both primary contexts.

### Route A — reproductive assurance

The simple pollinator-reduction model predicts higher `selfing_core` where at least one source-available functional channel is disrupted.

Direct High/Medium primary estimates are:

- northern mid-latitude: beta = **-0.0113**, p=`0.738`;
- tropical: beta = **+0.0588**, p=`0.536`.

The predicted positive direction is not reproduced in both contexts and neither coefficient is supported.

A threshold-sensitive tropical signal appears when the minimum number of evaluable channels is changed: `selfing_core` is +0.1895 (p=`0.0257`) at minimum one channel and +0.1799 (p=`0.00213`) at minimum three channels. Because the frozen primary minimum-two estimate is unsupported and the second pathway is not simultaneously recovered, this remains suggestive sensitivity rather than promoted evidence.

### Route B — pollinator-facing floral architecture independent of measured selfing

If reduced pollinator opportunity directly relaxes or reorganizes selection/filtering on pollinator-facing floral architecture, `generalized_accessible` should increase with disruption after conditioning on `selfing_core`.

Primary estimates are instead:

- northern mid-latitude: beta = **-0.00496**, p=`0.834`;
- tropical: beta = **-0.00863**, p=`0.897`.

The secondary shared named-architecture score and the legacy attraction-shift diagnostic are likewise unsupported in the primary minimum-two analysis.

The four primary context × pathway tests have FDR q approximately **0.897**. The pooled overlap gate passes, but the directional and FDR gates fail.

Therefore the strongest available global test now says:

> **Pooling all five functional pollinator channels makes the broad upstream hypothesis testable, but does not support a simple global mechanism in which source-to-island pollinator-channel disruption jointly drives reproductive assurance and a selfing-independent floral-accessibility shift.**

Workflow provenance: run **35051029608**, artifact **10429046082**, digest `sha256:ec56ad97899bc90774c05fd08c7384bebc853c703a4fed4746f2741aa2b0f065`.

## Why this does not close the pollination question

The pooled exposure is stricter and broader than the Bombus-only test, but it still measures **functional-channel retention/disruption**, not the biological quantity most likely to matter directly to plant fitness. Several alternatives remain open:

- pollinator abundance may decline while a functional channel remains detectable;
- visitation frequency and pollen delivery may decline without channel loss;
- different channels may compensate for one another;
- floral architecture may track **which** channel is retained rather than total channel attrition;
- native lineage assembly may absorb interaction effects before residual plant responses are measured;
- global cross-sectional occurrence states cannot identify local selection gradients.

Thus `pollinator decline` should not be used as a literal measured exposure in the manuscript. The defensible upstream term is **reduced or disrupted pollination-channel opportunity**, with temporal decline and effective-service decline retained as biological hypotheses requiring stronger data.

## Other independent gates remain closed

- N1 independent channel heterogeneity: `p=0.65516`;
- H5c independent biotic-vs-wind specificity: `p=0.41221`;
- H5d distributed-threshold identifiability: `0/8` cells qualified;
- original source-breadth promotion: `0/4` context x stratum cells.

## H5 synthesis

The evidence hierarchy is now:

```text
isolation-associated reproductive-assurance response
        supported
        +
isolation-associated pollination-architecture response
        supported and not reducible to measured selfing core
        ↓
named visitor identity inferred from floral phenotype
        no — >86% shared architecture
        ↓
GloBI source breadth identifies one global filter
        no
        ↓
individual exact-island channel retained/disrupted state identifies mechanism
        no — each channel has inadequate cross-context overlap
        ↓
pooled five-channel disruption explains both plant pathways
        no — overlap passes, but route A/B coefficients are unsupported
        ↓
coarse channel / biotic-wind / threshold tests identify mechanism
        no
```

Final H5 claim:

> **Island isolation is associated with partially separable reproductive-assurance and pollination-associated floral responses, but current independent interaction and exact-island occurrence data do not identify reduced pollinator-channel opportunity as their common upstream cause. Importantly, this conclusion is no longer based only on Bombus: a pooled five-channel exposure is estimable in both primary contexts and still fails the predicted two-route test.**

The remaining identification target is more specific: within-context estimates of pollinator abundance, visitation, pollen delivery and effective service, together with channel identity/turnover and plant lineage assembly.
