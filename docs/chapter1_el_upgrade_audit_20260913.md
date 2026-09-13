# Chapter 1 Ecology Letters upgrade audit — 2026-09-13

## Decision

The current Chapter 1 is strong enough for a **GEB-level macroecological paper** without further mechanism rescue. An Ecology Letters submission becomes reasonable only if one additional matched analysis supports a more general claim: the taxonomic depth at which isolation-associated floral architecture is expressed differs among biogeographic contexts.

The current evidence does **not** yet justify that sentence as a headline. It does justify running the matched test.

## 1. Endemic versus native-nonendemic cannot yet be promoted to the main axis

The latest 69.83% trait checkpoint was inspected directly from workflow run `34210169033`, artifact `10049591876`.

Tropical endemic support in the status-stratified lineage layer is:

| outcome | islands | support |
|---|---:|---|
| generalized form | 30 | pilot |
| plain colour | 46 | pilot |
| self compatibility | 29 | below pilot |

The endemic WHEN/WHERE omnibus retains only two colour responses at the pilot threshold and is unsupported (`p = 0.204481`). There is no confirmatory endemic response family at the frozen >=50-island threshold.

Therefore endemicity remains a useful **distribution-history / evolutionary-opportunity stratum**, but not a time axis and not a confirmatory assembly-versus-evolution discriminator.

This prevents a high-profile headline from depending on underpowered endemic strata.

## 2. Replace `4/4 -> 0/4` with effect sizes and intervals

The Wave52 progressive artifact (`33587356935`, artifact `9830619059`) was inspected directly.

For the Palearctic all-analysis, all-native result under `geo50_climate10`:

| response | observed slope [95% CI] | after family [95% CI] | after genus [95% CI] |
|---|---|---|---|
| generalized accessibility | 0.0369 [0.0084, 0.0655] | 0.0081 [0.0007, 0.0156] | 0.0035 [-0.0043, 0.0113] |
| selfing core | 0.0405 [0.0206, 0.0604] | 0.0338 [0.0101, 0.0575] | 0.0071 [-0.0018, 0.0160] |

Across all four source modes, the absolute all-analysis all-native slopes are reduced by about 82–91% after genus adjustment for these two primary axes. In native non-endemics the corresponding reduction is about 77–83%.

Direct-only evidence is less clean for generalized accessibility, so the correct interpretation is **not** “genus explains 100%”. The correct statement is:

> the broad Palearctic primary vector is strongly attenuated by source-matched genus composition and no longer passes the predeclared post-genus vector gate.

This wording preserves the evidence-scope asymmetry instead of selecting the sharper all-analysis result.

## 3. The tropical beyond-genus obstacle is partly outdated

Wave52 V4 contains stronger evidence than the earlier summary suggested.

The tropical `large_bee_like_residual` is positive after genus adjustment under all four source modes in both evidence scopes and both `all_native` and `native_nonendemic` strata. This remains a **floral-architecture residual**, not evidence for large-bee abundance or replacement.

More importantly, the direct **northern-midlatitude versus tropical beyond-genus architecture-vector contrast** is supported:

- direct-only: all four source modes in both primary floristic strata;
- all-analysis, native non-endemic: all four source modes;
- all-analysis, all-native: three of four source modes, with `geo50_climate10` narrowly above the threshold (`q ~= 0.054`).

This is the strongest current route toward an EL-level abstraction.

## 4. Why this still does not prove context-dependent taxonomic depth

The Palearctic taxonomic-depth result is based on the **primary pollinator-name-free H2/H3 response** (`generalized_accessible + selfing_core`). The tropical beyond-genus result comes from the **secondary V4 shared-architecture decomposition**.

Those are not the same response family.

Therefore the manuscript must not currently say:

> Palearctic responses occur at genus level whereas tropical responses occur beyond genus.

That sentence would compare unlike estimands.

The missing analysis is a direct, matched test using the same architecture response vector before and after genus adjustment in both northern-midlatitude and tropical contexts.

## 5. Required matched test before EL routing

The rules are frozen in `config/chapter1_el_hierarchical_depth_audit.yml` before implementation.

The required statistic is a multivariate `context x taxonomic_stage` interaction using the same four V4 architecture components in both contexts and all four source modes. A significant vector in one context and a non-significant vector in another is explicitly insufficient.

Promotion requires source-mode robustness in native non-endemics, directional concordance in all natives, and no contradiction between all-analysis and direct-only evidence.

If that gate passes, the general claim becomes:

> ecological filtering can be expressed at different taxonomic depths depending on biogeographic context.

If it fails, the paper retains the stronger island-specific conclusion and routes naturally to GEB.

## 6. GloBI source-side interaction breadth is feasible but secondary

The pinned GloBI 0.9 streaming code already retains, for every accepted flower interaction, the plant taxon, channel identity and independent reference key before aggregating to the pollinator-side catalog. The current artifact discards plant-level rows after aggregation, but the source runner can be extended without changing interaction filtering.

A valid source-side extension must estimate **sampled functional-channel breadth** for source-pool plant genera, not raw partner richness. No GloBI record must remain missing rather than be interpreted as specialization.

The predictor must be computed before joining island genus-entry outcomes and must retain independent-reference count as effort information. Any result remains D3 association evidence and cannot rescue H5/N1.

The useful question is:

> are source-available genera with narrower documented pollination-channel breadth increasingly under-represented as island isolation increases?

This belongs as a secondary H3 extension, not as the flagship mechanism test.

## 7. H4 and failed N1 should be retained as negative results

H4 is not “area has no effect”. The frozen result is that apparent small-island amplification fails the heteroskedastic-null promotion gate in all 16 primary classifications. Area therefore remains a measurement-sensitive modifier rather than an established capacity mechanism.

The N1 channel test is a stronger negative result: it was prospectively frozen, executable, non-significant, deletion-robust, and stopped before N2 with no rescue. It should appear in the manuscript/supplement as the empirical reason Chapter 1 does not promote a pollinator-channel mechanism.

## 8. Submission framing

The old `when and where` label should not be used as a temporal claim. The current study has no island-age, colonisation-time or divergence-time variable.

Recommended Chapter 1 question:

> **Where, whether, and at what assembly level does floral/reproductive island filtering emerge?**

Current routing:

- **GEB:** already defensible with biogeographic branching + source/genus assembly + negative H4/H5 gates.
- **Ecology Letters:** rational first submission only if the new matched taxonomic-depth interaction passes and the paper is rewritten around a general hierarchical-response principle rather than an island-specific syndrome correction.
