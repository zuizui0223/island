# Chapter 1 H5 — global pollination integration

## Current decision

H5 now has five evidence levels.

1. **Plant-side consistency:** floral trait combinations are compatible with changes in pollination-associated architecture.
2. **Independent interaction structure:** GloBI breadth and coarse channel tests ask whether interaction structure covaries with the plant pattern.
3. **Independent exact-island single-channel states:** Bombus and the other functional channels test whether a particular source-available channel is retained or strictly disrupted.
4. **Independent pooled channel-disruption bridge:** Bombus, non-Bombus bees, Lepidoptera, flower-visiting birds and Diptera are pooled into a conservative source-to-island functional-channel disruption exposure and tested against the two plant response routes.
5. **Identity-aware stacked bridge:** island × channel rows retain channel identity, equalize island weight, and test reproductive assurance plus channel-matched floral architecture separately by biogeographic context.

The fourth level rejects a simple total-channel-attrition explanation. The fifth shows that this pooled null is not the end of the pollination question: once channel identity is retained, there is a nominal direct-only tropical reproductive-assurance signal and a tropical matched-architecture effect in the predicted direction, but neither passes the full FDR/sensitivity gate. No upstream pollinator-causal mechanism is promoted.

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

## Pooled five-channel two-route bridge

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

Direct High/Medium primary estimates are:

- northern mid-latitude: beta = **-0.0113**, p=`0.738`;
- tropical: beta = **+0.0588**, p=`0.536`.

The predicted positive direction is not reproduced in both contexts and neither coefficient is supported.

A threshold-sensitive tropical signal appears when the minimum number of evaluable channels is changed: `selfing_core` is +0.1895 (p=`0.0257`) at minimum one channel and +0.1799 (p=`0.00213`) at minimum three channels. Because the frozen primary minimum-two estimate is unsupported and the second pathway is not simultaneously recovered, this remains suggestive sensitivity rather than promoted evidence.

### Route B — generic accessibility after reproductive assurance

The pooled model predicted that `generalized_accessible` should increase with disruption after conditioning on `selfing_core`.

Primary estimates are instead:

- northern mid-latitude: beta = **-0.00496**, p=`0.834`;
- tropical: beta = **-0.00863**, p=`0.897`.

The four primary context × pathway tests have FDR q approximately **0.897**. Thus the pooled overlap gate passes, but the directional and FDR gates fail.

> **Pooling all five functional pollinator channels makes the broad upstream hypothesis testable, but does not support a simple global mechanism in which total source-to-island pollinator-channel disruption jointly drives reproductive assurance and a selfing-independent generic accessibility shift.**

Workflow provenance: run **35051029608**, artifact **10429046082**, digest `sha256:ec56ad97899bc90774c05fd08c7384bebc853c703a4fed4746f2741aa2b0f065`.

## Identity-aware two-route diagnostic

The pooled exposure assumes that losing different pollinator channels has the same biological meaning. The next diagnostic removes that assumption.

The analysis unit is **island × pollinator channel**. It keeps only strict effort-qualified retained/disrupted states, gives every island equal total weight across its channel rows, absorbs static channel and regional differences with **channel × context fixed effects**, and controls distance, area and climate PC1-4.

### Identity-aware Route A — reproductive assurance

All five channels contribute to a common within-channel disruption test of `selfing_core`.

Direct High/Medium results:

- northern mid-latitude: beta = **+0.06735**, p=`0.3563`, q=`0.3563`;
- tropical: beta = **+0.21983**, p=`0.01528`, q=`0.06113`.

The tropical estimate is sizeable and nominally supported, but it narrowly misses the four-test FDR family and is not reproduced in the all-analysis sensitivity (beta = +0.07396, p=`0.5416`). It is therefore **suggestive only**.

### Identity-aware Route B — channel-matched architecture after selfing

Rather than assuming all disruptions should produce more generic flowers, each channel is matched to its predeclared named architecture:

- Bombus -> `large_bee_like`;
- Lepidoptera -> `butterfly_like`;
- flower-visiting birds -> `bird_like`.

Non-Bombus bees and Diptera are excluded from this identity-matched Route B because no equally specific named template was frozen for them. The model conditions on `selfing_core`.

Direct High/Medium results:

- northern mid-latitude: beta = **+0.06469**, p=`0.1251`, q=`0.1668`;
- tropical: beta = **-0.13820**, p=`0.1215`, q=`0.1668`.

The tropical sign is the predicted one: disruption is associated with lower concordance to the corresponding pollination-associated architecture after measured selfing is held constant. The all-analysis sensitivity keeps the negative sign (beta = -0.12066, p=`0.1624`), but neither analysis is supported.

Therefore:

> **Channel identity changes the biological interpretation of the pooled null, but does not close the mechanism. The direct-only tropical data are compatible with a reproductive-assurance response to channel disruption and with reduced channel-matched floral architecture, yet the evidence is not robust enough for promotion.**

Workflow provenance: run **35053988220**, artifact **10430156411**, digest `sha256:9f845fbb1b61e862e5b9956c43d132322304f93c82ad47c538828afa7105ec30`.

## Why this still does not close the pollination question

The sequence of tests now separates three increasingly biological hypotheses:

1. **total channel attrition** — not supported;
2. **channel-identity-aware disruption** — context-specific suggestive signal, not robustly supported;
3. **effective pollination-service limitation** — not yet measured.

Several mechanisms therefore remain open:

- pollinator abundance may decline while a functional channel remains detectable;
- visitation frequency and pollen delivery may decline without channel loss;
- different channels may compensate for one another;
- the identity and relative contribution of retained channels may matter more than channel count;
- native lineage assembly may absorb interaction effects before a residual plant response is measured;
- global cross-sectional occurrence states cannot identify local selection gradients.

Thus `pollinator decline` should not be used as a literal measured exposure in the manuscript. The defensible terms are **pollinator-channel disruption** for the present independent occurrence analysis and **pollination-service limitation** for the upstream biological hypothesis.

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
individual exact-island channels identify a global retained/disrupted mechanism
        no — cross-context overlap is inadequate per channel
        ↓
pooled five-channel attrition explains both plant pathways
        no — overlap passes, coefficients are unsupported
        ↓
identity-aware channel disruption
        tropical Route A nominal/direct-only; tropical Route B expected sign but unsupported
        ↓
effective pollination-service limitation
        not measured
```

Final H5 claim:

> **Island isolation is associated with partially separable reproductive-assurance and pollination-associated floral responses. Simple total pollinator-channel attrition does not explain both routes. Retaining channel identity reveals a suggestive tropical reproductive-assurance association and a directionally compatible matched-floral response, but neither survives the full inferential gate. Current independent data therefore do not identify pollinator decline or pollination-service limitation as the common upstream cause.**

The remaining identification target is now precise: within-context estimates of pollinator abundance, visitation and effective pollen delivery, linked to functional channel identity/turnover and plant lineage assembly.
