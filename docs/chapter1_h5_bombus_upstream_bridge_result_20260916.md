# Chapter 1 H5 Bombus upstream bridge — 2026-09-16

Status: **completed; canonical exact-island Bombus occurrence evidence does not identify Bombus-channel disruption as the upstream cause of either plant response component.**

## Provenance

- branch: `ch1-all-data-primary`
- contract: `chapter1_h5_bombus_upstream_bridge_v1`
- workflow: `Run Chapter 1 H5 Bombus upstream bridge`
- run: **35049673581**
- artifact: `chapter1-h5-bombus-upstream-bridge-35049673581`
- artifact ID: **10428646713**
- artifact digest: `sha256:3dfdd96837e99970a2e70e858fb2dd1a0e344909f3ac328fd6bdc4bcda40d734`

Independent Bombus evidence came from the canonical exact-island Search run **34731405488**, Bombus artifact **10313926043**. That Search was frozen independently of focal plant traits and used the predeclared effort-aware channel observation policy.

## Independent Bombus state definition

Only source-available islands with one of two evaluable observation states enter the bridge:

- `detected` -> `retained`;
- `adequate_non_detection` -> `disrupted`.

`insufficient_effort` and `unresolved` are missing, never absence. A detection means potential Bombus-channel presence, not realized visitation or effective pollination service. An adequate non-detection means source-available and adequately non-detected under the frozen policy, not proven historical extinction.

Raw Bombus record counts are not used as abundance because the canonical Search may stop after a confirmatory detection.

## Canonical global Search support

Among 7,154 source-available islands in the canonical Bombus Search:

- detected: **824**;
- adequate non-detection: **23**;
- insufficient effort: **5,983**;
- unresolved: **324**.

Thus only 847 islands provide a strict retained/disrupted state before joining plant outcomes.

## Critical overlap result

After joining the Chapter 1 covariates and direct-evidence plant scores, strict evaluable state support is:

| context | retained | disrupted | evaluable |
|---|---:|---:|---:|
| northern high-latitude | 47 | 0 | 47 |
| northern mid-latitude | 751 | **5** | 756 |
| southern extratropical | 21 | 1 | 22 |
| tropical | **5** | 17 | 22 |

The frozen overlap gate required at least 10 retained and 10 disrupted islands in both the northern-midlatitude and tropical primary contexts. **Neither context passes.** The bridge classification is therefore `overlap_gate_failed` before coefficient interpretation.

This is the main inferential result: strict Bombus retained/disrupted states are almost confounded with biogeographic context. The available global occurrence evidence cannot cleanly separate a Bombus disruption effect from regional geography.

## Direct High/Medium primary diagnostics

Models use equal island weight and adjust for biogeographic context, log distance, log island area, climate PC1-4, and spatial-block clustered covariance.

North-midlatitude + tropical:

| response | Bombus disrupted coefficient | SE | p |
|---|---:|---:|---:|
| `selfing_core` | +0.04725 | 0.13301 | 0.7224 |
| `attraction_shift` | +0.04749 | 0.04511 | 0.2924 |
| `attraction_shift` conditional on `selfing_core` | +0.03799 | 0.03170 | 0.2308 |
| `large_bee_like` | -0.01809 | 0.03401 | 0.5948 |
| `generalized_accessible` | +0.07771 | 0.06394 | 0.2242 |

The two pathway coefficients (`selfing_core` and attraction shift conditional on `selfing_core`) have the expected positive sign, but neither is statistically supported and the overlap gate has already failed.

Across all four contexts, results are likewise unsupported:

- `selfing_core`: beta=+0.06368, p=0.6080;
- `attraction_shift`: beta=+0.04838, p=0.2355;
- conditional attraction shift: beta=+0.03493, p=0.2463;
- `large_bee_like`: beta=-0.01594, p=0.6145;
- `generalized_accessible`: beta=+0.07782, p=0.1742.

## All-analysis evidence sensitivity

North-midlatitude + tropical:

- `selfing_core`: beta=+0.09835, p=0.3251;
- `attraction_shift`: beta=+0.02554, p=0.3774;
- conditional attraction shift: beta=+0.00322, p=0.9171;
- `large_bee_like`: beta=-0.02139, p=0.4667;
- `generalized_accessible`: beta=+0.02973, p=0.4129.

Thus broad evidence does not rescue the upstream Bombus association.

## H3A post-genus residual diagnostic

The bridge was also applied to the two clearest direct-evidence H3A post-genus residual components.

North-midlatitude + tropical:

- `generalized_form` post-genus residual: beta=-0.01245, p=0.1965;
- `self_compatibility` post-genus residual: beta=+0.02786, p=0.4277.

Across four contexts:

- `generalized_form`: beta=-0.01357, p=0.1125;
- `self_compatibility`: beta=+0.02873, p=0.3833.

Therefore canonical Bombus retained/disrupted status does not explain the broad post-genus residual response either.

## What this proves and what it does not

The Chapter 1 plant data already support two partially separable response components:

1. reproductive assurance / `selfing_core`;
2. pollination-associated floral architecture / `attraction_shift`.

This new independent cross-examination asks whether canonical exact-island Bombus-channel disruption closes the upstream arrow into those two components. It does **not**.

The failure has two parts:

1. **identifiability:** retained/disrupted overlap within the primary contexts is too poor for a defensible cross-context causal comparison;
2. **association:** even conditional adjusted coefficients are unsupported in both evidence scopes.

The correct conclusion is therefore:

> **Isolation is associated with two partially separable plant response components, but the current independent exact-island Bombus occurrence evidence does not identify Bombus-channel disruption as their upstream cause. Strict evaluable Bombus states have insufficient within-context overlap, and adjusted associations with both reproductive assurance and floral architecture are unsupported.**

This is not evidence that Bombus or pollinators are irrelevant. It shows that the present global occurrence data cannot turn the two-pathway plant pattern into a demonstrated pollinator-decline mechanism. Stronger identification would require within-context variation in reliable channel state and, ideally, visitation or effective-service data rather than occurrence alone.
