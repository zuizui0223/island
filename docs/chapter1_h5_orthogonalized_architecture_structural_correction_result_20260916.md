# Chapter 1 H5 — orthogonalized named-architecture structural correction

## Result

The post-hoc structural correction does **not** recover a robust North–Tropical H2 signal in either the dominant shared floral-architecture dimension or the template-specific residual subspace.

Canonical provenance:

- workflow run: **35065912097**
- workflow head: `dcc46e8a5e0439819d083fd7098cbbe055d60875`
- artifact: **10433982600**
- digest: `sha256:a1b259d7657757704d3437701e93b2a9f3daa88366ca4dc95f6d4fb1f938ac59`
- validation: **12 tests passed**; Ruff passed.

This analysis is explicitly a **post-hoc structural correction**, not a new confirmatory pollinator-identity test.

## Why a structural correction was required

The earlier frozen v1 correctly recognized that `large_bee_like`, `butterfly_like`, and `bird_like` share most of their floral-trait geometry, and it reused the existing source-trained V4 PCA decomposition. Two interface defects were then exposed before interpretation.

First, removing one shared PCA dimension from three template dimensions leaves a residual subspace of maximum rank **two**, not three. The v1 contract had required `df=3` merely because three residual labels remained. That requirement was mathematically impossible.

Second, the generic multivariate branching helper deliberately refuses a one-response vector, so the one-dimensional shared PCA factor was mechanically returned as `not_testable` even though a scalar distance × context interaction is estimable.

A third implementation error was also found and corrected before the final analysis: the frozen support condition was **at least 50 scored species per island per component**, but the first wrapper had accidentally connected that value to the helper's minimum number of islands. The final analyses apply the intended `n_species >= 50` filter before model fitting.

The v1 result therefore remains locked as a failed design rather than being overwritten. V2 changes only the linear-algebra/interface defects while preserving the biological estimand and support rule.

## Shared architecture dominates the named templates

The source-trained first principal component explains:

- all-analysis: **86.94%** of variance;
- Direct: **86.44%**.

Thus raw named-template slopes are mostly contrasts along a common floral-architecture dimension. They cannot be read as three independent pollinator-guild responses.

## Rank-two template-specific residual H2

After applying the frozen `>=50` scored species per island/component gate, all three residual labels are retained, but their estimable joint covariance rank is two.

| evidence scope | islands | spatial blocks | chi-square | df | North–Tropical p |
|---|---:|---:|---:|---:|---:|
| all-analysis | **154** | 57 | 0.975 | 2 | **0.6143** |
| Direct | **144** | 54 | 2.383 | 2 | **0.3038** |

Therefore the broad H2 contrast is **not supported in the template-specific residual subspace** in either evidence scope.

This matters because it removes a tempting interpretation of the earlier raw syndrome-concordance results: the North–Tropical contrast is not robustly recovered as a guild-labelled residual geometry after the dominant shared component is removed.

It does **not** show that pollinator identity is biologically irrelevant. The residual scores are still plant-trait contrasts, and the analysis contains no realized visitation, abundance or pollen-delivery measurement.

## Shared one-dimensional architecture H2

The shared PCA factor is tested separately with a scalar cluster-robust distance × tropical interaction while controlling island area and climate PC1–4.

### All-analysis

- North: 82 islands
- Tropical: 72 islands
- total: **154 islands / 57 blocks**
- North slope: +0.0322, p=`0.435`
- Tropical slope: +0.0716, p=`0.0251`
- direct North–Tropical interaction: **p=`0.4324`**

### Direct

- North: 75 islands
- Tropical: 69 islands
- total: **144 islands / 54 blocks**
- North slope: +0.0311, p=`0.445`
- Tropical slope: +0.0825, p=`0.00794`
- direct North–Tropical interaction: **p=`0.2782`**

The tropical shared factor has a positive distance slope in both evidence scopes, but the formal H2 question is the **difference in distance response between North and Tropical**, and that interaction is unsupported in both scopes. A significant slope in one context and a nonsignificant slope in another is not itself heterogeneity.

## H5 interpretation

This creates a useful decomposition of the earlier syndrome-consistency layer:

```text
raw named-template concordance
        ↓
~86–87% shared floral architecture
        ↓
shared architecture North–Tropical H2
        unsupported (p=.432 / .278)
        +
rank-two identity-labelled residual North–Tropical H2
        unsupported (p=.614 / .304)
```

The conclusion is therefore narrower but cleaner:

> **The current all-observed North–Tropical floral/reproductive H2 signal is not robustly reproduced by either the dominant shared dimension or the template-specific residual geometry of the three sampled named pollination-syndrome templates under the strict component-support gate. Raw syndrome concordance should therefore remain a compatibility layer, not evidence of pollinator identity or the mechanism generating H2.**

This result strengthens the separation between the broad plant response and the mechanism tests. It is consistent with the independent H5 results showing that total channel attrition is unsupported, identity-aware channel disruption is suggestive only, and the strict same-island channel test is structurally non-identifiable in the northern context.

## Claim ceiling

The result does not establish any of the following:

- pollinators are irrelevant to island floral evolution;
- large bees, butterflies or birds are absent or ineffective;
- shared architecture is not pollination-associated;
- identity-specific selection does not occur;
- the post-hoc structural correction is confirmatory evidence.

The remaining mechanism target is still independent pollination-service variation — abundance, visitation and effective pollen delivery — with enough within-context identity crossing to separate channel identity from service state.
