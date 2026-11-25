# 14 November 2026 — Native provenance sensitivity of marginal, joint, and pairing-specific island floral functional diversity

**Original data executed, 2 independent fixed random seeds completed. Status: exploratory observational assemblage analysis, NOT identified pollinator-guild loss or floral evolutionary adaptation.**

- Original protected \`island\` Actions **#43564079792 SUCCESS**, two seeds (20261114 and 20261115), anonymous source artifacts \`#14101806132\` and \`#14101965769\`.
- Independent numeric implementation: [fd finite Rao/coupling module](../../src/fd/regional_native_flora_fd.py), [exhaustive combinatorial/fixed-margin synthetic tests](../../tests/test_regional_native_flora_fd.py), [island original-source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_wcvp_native_flora_fd.py), [original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-wcvp-native-flora-fd.yml). Historical 17D KDE not used.
- Input eligibility and original source SHA gate: **2,234** species-direct high/medium complete flower colour/form/symmetry names, **928** island floras with at least 30 source-complete recorded species, and original protected five-nearest GIFT native mainland source lists; **all exact accepted species names recorded on each focal island are excluded from its source flora before donor evaluation**.
- Native origin tiers: original **all-recorded 928**, strictly GIFT-native **87**, and independent WCVP TDWG-L3 **regionally native-compatible 357**. The last is NOT proof of local-island nativeness, and no GIFT-introduced record was upgraded.
- Source-comparable flora counts for the current stricter exact-species donor exclusion are below the corresponding eligible counts in some regions (e.g. northern midlatitude all **736/743**, WCVP-compatible **285/287**).

## Estimands — finite population exactness

For each original local plant population of N accepted trait-complete species, the expected categorical Gini–Simpson Rao Q of an **unbiased n=30 species without-replacement draw** is analytic:
\[
E[Q_{30}] = \frac{29}{30}\left(1-S\right),
\]
where S is the probability two DISTINCT species share the relevant fixed trait state. \`marginal_rao\` averages three separate trait-specific Q's, and \`joint_rao\` uses the ORIGINAL complete recorded triple state (never Cartesian-expands multistate flower labels). The open-radial and actinomorphic shares are finite-population prevalence expectations.

The new pairing-specific endpoint preserves each original FORM+SYMMETRY architecture pair while permuting only the recorded colour labels. Write S_C for distinct-species same-colour probability, S_FS for same-architecture probability, and S_CFS for same complete-triple probability. Then
\[
D_{\rm pairing} = E[Q_{30}^{\rm observed\ joint}] -
 E[Q_{30}^{\rm colour-shuffled\ joint}]
 = \frac{29}{30}(S_C S_{FS}-S_{CFS}).
\]
This is **NOT mutual information** and NOT an excess of documented effective pollination; positive values indicate less same-colour+same-architecture clustering than in the fixed-marginal colour-shuffle reference. The MAIN reported response is **island-minus-GIFT-source difference in D_pairing**, not a statement that absolute island coupling is positive.

All outcome calculations are deterministic conditional on original source records. Two independently seeded 1,999-draw physical 10° block and nearest GIFT source-signature cluster bootstraps alter intervals slightly but reproduce *identical point estimates*. The two groupings overlap; neither is phylogenetic or ancestral-source independence. No FDR discovery claim is made across these post-hoc endpoints.

## Main cross-sectional source-relative results

**Original original-data island minus its five-nearest GIFT candidate native mainland-source flora, equal 30 species.** Each entry is equal-geographic-block mean [conditional 95% block-bootstrap interval]. A separate source-signature-equal interval is used as geographic-source sensitivity.

| Region/origin tier | Source-comparable islands / 10° blocks / first GIFT sources | Marginal Rao Q Δ | Joint Rao Q Δ | Colour–architecture pairing Δ |
| --- | ---: | --- | --- | --- |
| Northern midlatitude, all-recorded | **736 / 51 / 64** | **−0.03119 [−0.04346,−0.01903]** | **−0.01777 [−0.02744,−0.00824]** | −0.00261 [−0.00647,+0.00142] |
| Northern midlatitude, WCVP regional-compatible | **285 / 33 / 40** | **−0.03043 [−0.04844,−0.01190]** | **−0.01472 [−0.02689,−0.00261]** | −0.00058 [−0.00546,+0.00424] |
| Northern midlatitude, strict GIFT-native | **67 / 16 / 23** | −0.01487 [−0.03551,+0.00655] | −0.00965 [−0.02430,+0.00520] | −0.00311 [−0.00817,+0.00117] |
| Tropical, all-recorded | **106 / 42 / 36** | −0.00503 [−0.01428,+0.00349] | +0.00639 [−0.00016,+0.01345] | **+0.00651 [+0.00258,+0.01077]** |
| Tropical, WCVP regional-compatible | **36 / 20 / 16** | −0.00994 [−0.02257,+0.00310] | +0.00783 [+0.00001,+0.01733] | **+0.00881 [+0.00269,+0.01605]** |
| Tropical, strict GIFT-native | **17 / 10 / 9** | −0.01063 [−0.02768,+0.00599] | +0.00896 [−0.00199,+0.02138] | **+0.01171 [+0.00236,+0.02203]** |
| Southern extratropical, all-recorded | **24 / 11 / 9** | +0.04519 [+0.02229,+0.07005] | +0.03984 [+0.01595,+0.06849] | +0.03291 [+0.01683,+0.04968] |

**Separate nearest-source-signature-equal findings (not replications):**

- Northern-midlat all-recorded: marginal Δ **−0.01945 [−0.03058,−0.00782]**, joint Δ **−0.01044 [−0.01925,−0.00216]**; pairing Δ **−0.00082 [−0.00415,+0.00231]**.
- Northern-midlat WCVP-compatible: marginal Δ **−0.00822 [−0.02585,+0.00831]**, joint Δ **−0.00078 [−0.01228,+0.01058]**; both now include zero; pairing still uncertain.
- Tropical all-recorded: pairing Δ **+0.00749 [+0.00325,+0.01164]**.
- Tropical WCVP-compatible: pairing Δ **+0.01005 [+0.00288,+0.01818]**.
- Tropical strict GIFT-native: pairing Δ **+0.01182 [+0.00137,+0.02245]**, but only 17 islands/9 source signatures.
- Southern has only 24 all-recorded islands and fewer independent geographic groups; just 4 regionally compatible and 3 strict-native islands, both unreportable. Northern-high has only 3 or fewer nearest GIFT source signatures.

**Critical contrast:** northern-midlatitude source-relative marginal and joint diversity losses remain visible on the geographically BLOCK-equal WCVP-compatible 285-island cohort, but are NOT directional under nearest SOURCE-signature equal weighting. In tropical floras, *marginal* Q has no supported source-relative decline, whereas the analytical pairing-specific source-relative departure is positive on the broad and both much smaller native-origin tiers. These are different diversity dimensions; a universal one-dimensional "floral syndrome contraction" is not supported. The tropical pairing endpoint is new, post-hoc, and not genus/phylogeny/provider conditioned.

## Stronger exact-SAME-island origin filter: floristic composition, not independent island cohorts

Restrict to identical physical island IDs and identical original candidate mainland donor. The source term cancels between the two tiers, isolating a change in which original records are retained by plant provenance.

| Region, WCVP compatible minus all-recorded on SAME islands | Island count / blocks / source signatures | Marginal Rao Q shift [block interval] | Joint Rao Q shift [block interval] | Open-radial flower fraction shift [block interval] | Actinomorphic fraction shift [block interval] |
| --- | ---: | --- | --- | --- | --- |
| Northern midlatitude | **285 / 33 / 40** | +0.00030 [−0.00527,+0.00636] | −0.00202 [−0.00543,+0.00136] | +0.00662 [−0.00448,+0.01830] | −0.00918 [−0.02262,+0.00350] |
| **Tropical** | **36 / 20 / 16** | −0.00653 [−0.01448,+0.00264] | +0.00010 [−0.00393,+0.00571] | **−0.03075 [−0.04959,−0.01389]** | **−0.02015 [−0.03780,−0.00345]** |

The tropical within-exact-island open-radial and actinomorphy fractions are ALSO negative under independent nearest-source-signature equal resampling (open **−0.03523 [−0.05653,−0.01509]**; symmetry **−0.02612 [−0.04267,−0.00948]**). Thus a native-regional provenance filter materially changes the recorded *morphology composition* of the **same tropical islands**, even while neither marginal nor joint Rao difference has a conditional interval excluding zero.

The contrast does **not** demonstrate that the omitted plants were introduced or that native status directly altered morphology: many excluded entries have unresolved origin, WCVP region status is coarse, and island flora traits are species-fixed. For the 17 strict-native tropical islands, the nested sample is still small.

## Result status and next falsifications

**New source-verified descriptive association:** source-relative tropical **colour–architecture pair-geometry differs** despite near-unchanged marginal diversity; source-relative northern-midlatitude recorded Rao declines are also present, but their robustness depends on native cohort and geographic/source weighting. Source provenance changes apparent tropical open-flower/symmetry composition within the same islands. This is informative about *floristic trait assembly* and the distinct roles of diversity versus trait coordination, not verified bird/bee/moth pollination syndromes.

**Not identified:** local decline or temporal replacement of pollinator channels, pollen delivery, species-level fitness, common botanical genus/phylogeny held constant across island and source, independent historical oceanic-vs-continental selection, local island nativeness of WCVP region-compatible taxa, or new reproductive phenotype evolution. GloBI/GloPL/GBIF matched-guild support remains separately insufficient for causal ecological mediation.

A stronger next test would fix genus and family identity or the three published true-tip backbone distances inside WCVP-compatible island/source floras *before* evaluating the new pairing endpoint, and repeat the same-source contrast with independent taxon-local provenance. The existing 928-island source-qualified GBIF channel labels provide **0/20** region×guild supported retained-vs-adequately-nondetected effects, so they must not be repurposed as a pollinator selection coefficient. For an actual pollinator shift, require same-plant locally observed visitation/pollen success, time/effort and matched mainland/island sites.

Reproduction: original dual-seed [Actions #43564079792](https://github.com/zuizui0223/island/actions/runs/43564079792). Branches are draft; no paper-level universal causal claim is promoted.
