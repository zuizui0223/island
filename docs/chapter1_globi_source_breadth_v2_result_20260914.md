# Chapter 1 GloBI source-side breadth V2 — frozen result

## Decision

The source-side GloBI interaction-breadth extension is **not promoted** as an explanation for the Chapter 1 H3 genus-assembly pattern.

The final V2 analysis matched candidate source genera not only on source prevalence and source species richness, but also on a broad GloBI independent-reference effort class. This was frozen before the island-outcome join to reduce the obvious risk that well-studied genera look artificially generalised.

## Outcome-blind predictor receipt

The predictor was built before island occurrence, island distance or Chapter 1 trait outcomes were loaded.

- GloBI version DOI: `10.5281/zenodo.20546682`
- raw interaction rows: **24,577,183**
- retained explicit flower-interaction rows: **714,785**
- plant-matched evidence rows: **524,922**
- unmatched evidence rows: **189,863**
- exact plant-name matches: **524,729**
- unique-binomial fallback matches: **193**
- genus breadth rows: **2,228**
- predictor SHA-256: `9f43ee056f285f8afcf92393170cf615236f7d99bd3d53e65fe998c415361386`

No GloBI record was interpreted as missing breadth, not as specialization.

## Final effort-matched analysis

- workflow run: `34790842133`
- artifact: `chapter1-globi-source-breadth-v2-34790842133`
- artifact ID: `10328696801`
- digest: `sha256:475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`
- island-enrichment rows: **38,260**
- fitted slope rows: **576**
- primary promoted context × stratum cells: **0 / 4**

The primary estimand was the isolation slope of source-matched genus `entry_enrichment` using `effective_channel_number`, three or more independent GloBI references, `prevalence + source richness + reference-effort` matching, and the four frozen source definitions.

All four context × floristic-stratum cells were classified `not_promoted`:

- northern-midlatitude / all-native;
- northern-midlatitude / native-nonendemic;
- tropical / all-native;
- tropical / native-nonendemic.

Representative estimates are near zero and their intervals overlap zero. For example:

- northern-midlatitude / all-native / `geo50_climate10`: slope `0.0060`, 95% CI `[-0.0092, 0.0213]`, q=`0.971`;
- tropical / all-native / `geo50_climate10`: slope `0.0077`, 95% CI `[-0.0052, 0.0207]`, q=`0.361`;
- northern-midlatitude / native-nonendemic / `geo_k10`: slope `-0.0067`, 95% CI `[-0.0140, 0.0006]`, q=`0.4243`;
- tropical / native-nonendemic / `geo_k5`: slope `-0.0150`, 95% CI `[-0.0306, 0.0006]`, q=`0.3568`.

The frozen robustness rule also failed: the direction did not remain stable across the 2/3/5-reference effort thresholds in the required way.

## Scientific consequence

The result does **not** mean pollination interaction dependence is irrelevant to island assembly. It means the specific independently documented source-genus breadth proxy does not recover the H3 distance-associated genus-assembly pattern after source availability and broad sampling effort are controlled.

This is useful because it removes a tempting shortcut explanation:

> **The Palearctic primary response is strongly genus-structured, but that structure is not simply recovered by island area, a global nonlinear threshold, prospective channel-specific pollinator isolation, or source-genus sampled functional-channel breadth.**

The strongest positive result therefore remains the context-dependent plant assemblage hierarchy itself rather than a promoted pollination mechanism.

## Claim ceiling

This is D3 association evidence only. It does not identify true ecological specialization, historical pollinator loss, effective service, colonization failure, extinction or causal interaction dependence. It cannot rescue N1 or open N2.
