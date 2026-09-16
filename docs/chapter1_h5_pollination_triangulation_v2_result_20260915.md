# Chapter 1 H5 pollination triangulation v2 — result

Status: **completed; plant-side pollination architecture is consistent, but no single independent pollinator mechanism is promoted.**

## Provenance

- branch: `ch1-all-data-primary`
- workflow: `Run Chapter 1 H5 pollination triangulation v2`
- run: **34977453917**
- artifact: `chapter1-h5-pollination-triangulation-v2-34977453917`
- artifact ID: **10399424739**
- digest: `sha256:86c17147942d34e400b0d4eb1ac5ca0c0467d0ed6a947a9ca294a05195157b9a`
- frozen plant input: run `34232450884`
- frozen GloBI V2 input: run `34790842133`, artifact `10328696801`

The H5 v2 contract was frozen before the new GloBI distance-by-area result was inspected.

## 1. Plant-side syndrome consistency

The primary consistency layer used High/Medium direct evidence and the two predeclared regional predictions that are testable on confirmatory native support:

1. northern mid-latitude: `sampled_large_bee_concordance` decreases with isolation;
2. tropical: `sampled_butterfly_concordance` is maintained or increases with isolation.

Both predictions were FDR-supported in both native strata under direct evidence.

### Northern mid-latitude

- all native: large-bee-like slope `-0.04317`, q=`0.02871`, n=234;
- native non-endemic: `-0.04544`, q=`0.02579`, n=234.

The all-analysis sensitivity retained the negative sign but not FDR support:

- all native: `-0.02565`, q=`0.26835`;
- native non-endemic: `-0.02756`, q=`0.20077`.

Thus the northern large-bee-like consistency is **direct-evidence supported but evidence-scope sensitive**.

### Tropical

Butterfly-like concordance increased robustly with isolation:

- direct all native: `+0.05765`, q=`0.01204`, n=123;
- direct native non-endemic: `+0.06817`, q=`0.00465`, n=122;
- all-analysis all native: `+0.05722`, q=`0.01131`;
- all-analysis native non-endemic: `+0.06632`, q=`0.00257`.

This is strong plant-side consistency with the predeclared tropical pollination-architecture prediction.

## 2. Named syndromes are not pollinator identities

The source-trained common architecture factor explains most covariance among large-bee-like, butterfly-like and bird-like templates:

- all-analysis: **86.94%**;
- direct-only: **86.44%**.

Both exceed the frozen 80% identity-specificity ceiling. Therefore the named templates are retained as **pollination-associated floral architecture concordances**, not as realized pollinator assignments.

## 3. New area-conditioned GloBI bridge

The new independent test was aligned to redesigned H3. It asked whether genera with broader source-side documented flower-interaction channels become more enriched with isolation especially on small islands.

Primary predictor:

- GloBI `effective_channel_number`;
- >=3 independent references;
- source prevalence + source richness + reference-effort matching;
- source-matched genus `entry_enrichment`;
- equal island weight;
- model: `entry_enrichment ~ distance + area + distance:area + climate PC1-4`;
- spatial-block cluster-robust covariance.

The predeclared northern prediction was a **negative distance x area** coefficient: broader-channel genera should be increasingly enriched with isolation more strongly on smaller islands.

This prediction failed decisively as a source-definition-robust mechanism.

### All native

- `geo50_climate10`: `-0.00241`, p=`0.6953`;
- `geo_k10`: `+0.00787`, p=`0.0918`;
- `geo_k20`: `+0.00325`, p=`0.3862`;
- `geo_k5`: `+0.01067`, p=`0.0473`, FDR q=`0.1836`.

### Native non-endemic

- `geo50_climate10`: `-0.00160`, p=`0.7401`;
- `geo_k10`: `+0.00955`, p=`0.00130`, q=`0.00520`;
- `geo_k20`: `+0.00425`, p=`0.1416`, q=`0.1888`;
- `geo_k5`: `+0.01159`, p=`0.00443`, q=`0.00886`.

The signs are not only non-robust; two source definitions in native non-endemics show a significant interaction in the **opposite** direction from the simple small-island generalist-filter prediction.

The direct North-Tropical distance x area x context comparison was also not source-definition robust. At the primary 3-reference effort threshold, none of the four source modes survived FDR within either native stratum.

Verdict: **simple_pollination_breadth_filter_not_promoted** in both native strata.

## 4. Independent triangulation remains negative

The new result joins the frozen independent tests:

- original GloBI distance-only breadth: `0/4` context x stratum cells promoted;
- N1 independent channel heterogeneity: `W=1.6187`, df=3, p=`0.65516`;
- H5c independent biotic-vs-wind specificity: distance x biotic `+0.06495`, 95% CI `[-0.09030, 0.22020]`, p=`0.41221`;
- H5d distributed-threshold identifiability: `0/8` cells qualified.

No independent test currently identifies one global pollinator-filter mechanism.

## Integrated interpretation

The evidence now separates three statements that should not be conflated.

1. **Plant-side consistency is real.** Northern floral architecture moves away from the predeclared large-bee-like template under direct evidence, while tropical flora retains/increases butterfly-like architecture.
2. **Named pollinator identity is not identified.** More than 86% of named-template covariance is common plant architecture.
3. **The obvious independent shortcut mechanisms fail.** Neither source-genus interaction breadth, coarse pollination-channel heterogeneity, nor biotic-vs-wind specificity explains the defended plant pattern.

The current H5 conclusion is therefore:

> **Isolation-associated floral responses are consistent with altered pollination-associated architecture, but the present independent interaction data do not support one global pollinator-filter mechanism. If pollination contributes to the genus-structured native response, it must operate through a more specific interaction dimension than sampled source-genus channel breadth or coarse biotic-versus-wind class.**

This is not evidence that pollinators are irrelevant. It narrows the mechanism search.
