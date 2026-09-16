# Chapter 1 H5 — global pollination integration

## Current decision

H5 now has three evidence levels.

1. **Plant-side consistency:** floral trait combinations are compatible with changes in pollination-associated architecture.
2. **Independent interaction structure:** GloBI breadth and coarse channel tests ask whether interaction structure covaries with the plant pattern.
3. **Independent exact-island Bombus state:** a new post-hoc bridge asks whether source-available Bombus `retained` versus strictly `disrupted` states predict the two partially separable plant response components.

The third level is the closest current test of the proposed upstream `pollinator decline -> plant response` arrow, but it fails its overlap/identifiability gate and does not support adjusted associations.

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

## New exact-island Bombus upstream bridge

A separate canonical exact-island Search campaign provides outcome-blind Bombus occurrence evidence on the frozen global island universe. For Bombus, among 7,154 source-available islands the canonical Search recorded:

- detected: **824**;
- adequate non-detection: **23**;
- insufficient effort: **5,983**;
- unresolved: **324**.

The frozen observation policy maps `detected -> retained` and `adequate_non_detection -> disrupted`. Insufficient-effort and unresolved islands are missing, never absence. Detection means potential channel presence, not visitation or effective service; disruption means adequate non-detection under the frozen policy, not proven historical extinction.

The new bridge is frozen as `chapter1_h5_bombus_upstream_bridge_v1` and was run successfully as **35049673581**, artifact **10428646713**.

### Overlap gate fails

Direct-evidence evaluable support after joining Chapter 1 outcomes is:

| context | retained | disrupted |
|---|---:|---:|
| northern high-latitude | 47 | 0 |
| northern mid-latitude | 751 | **5** |
| southern extratropical | 21 | 1 |
| tropical | **5** | 17 |

The frozen primary gate required at least 10 retained and 10 disrupted islands in both northern-midlatitude and tropical contexts. **0/2 primary contexts pass.** Strict Bombus state is therefore nearly confounded with biogeographic context, preventing a defensible upstream effect estimate.

### Adjusted associations are unsupported

North-midlatitude plus tropical, direct High/Medium evidence, adjusting for context, distance, area, climate PC1-4 and spatial block:

- `selfing_core`: beta(disrupted)=**+0.04725**, p=`0.7224`;
- `attraction_shift`: **+0.04749**, p=`0.2924`;
- `attraction_shift | selfing_core`: **+0.03799**, p=`0.2308`;
- `large_bee_like`: **-0.01809**, p=`0.5948`;
- `generalized_accessible`: **+0.07771**, p=`0.2242`.

The two primary pathway coefficients have the expected positive sign but are imprecise and unsupported. All-analysis sensitivity is also unsupported (`selfing_core` p=`0.3251`; conditional attraction p=`0.9171`).

The same Bombus-state bridge does not explain the two clearest H3A post-genus residual components either:

- `generalized_form`: p=`0.1965`;
- `self_compatibility`: p=`0.4277`.

This result does not show that Bombus is irrelevant. It shows that the current independent exact-island occurrence data do not identify the proposed upstream Bombus-disruption mechanism because reliable within-context state overlap is inadequate and the adjusted associations are not supported.

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
exact-island Bombus retained/disrupted state explains both pathways
        no — overlap gate failed; adjusted associations unsupported
        ↓
coarse channel / biotic-wind / threshold tests identify mechanism
        no
```

Final H5 claim:

> **Island isolation is associated with partially separable reproductive-assurance and pollination-associated floral responses, but current independent interaction and exact-island Bombus occurrence data do not identify pollinator decline as the upstream cause. In particular, strict Bombus retained/disrupted states have insufficient within-context overlap and do not show supported adjusted associations with either pathway.**

The remaining identification target is therefore more specific: repeated within-context measurements of visitor continuity/turnover and effective pollination service, rather than floral phenotype or coarse source interaction breadth alone.
