# Chapter 1 global GloBI V3 — four-context result

Status: **completed; global coverage achieved, but no source-definition-robust global pollinator-breadth mechanism is promoted.**

## Provenance

- branch: `ch1-all-data-primary`
- contract: `chapter1_globi_global_v3`
- workflow: `Run Chapter 1 global GloBI v3`
- successful run: **34979722831**
- artifact: `chapter1-globi-global-v3-34979722831`
- artifact ID: **10401900158**
- digest: `sha256:c821eebc056cf7aeaad35e24682137f5eb9a50edad2247923428e68db241f1cf`
- frozen GloBI genus-breadth predictor SHA-256: `9f43ee056f285f8afcf92393170cf615236f7d99bd3d53e65fe998c415361386`
- frozen Chapter 1 input run: `34232450884`

The predictor itself was unchanged from the outcome-blind GloBI V2 build. The global extension was frozen before its outcome inspection.

## Design

The analysis extended GloBI source-genus `effective_channel_number` to the same four analysis regimes used by redesigned H1-H3:

1. northern mid-latitude;
2. northern high-latitude;
3. tropical;
4. southern extratropical.

Primary flora scope: `all_observed`.

Biological sensitivities: `all_native` and `native_nonendemic` where support permitted.

Primary source expectation:

- source prevalence + source species richness + GloBI reference-effort matching;
- minimum 3 independent GloBI references;
- four frozen source assignments: `geo_k5`, `geo_k10`, `geo_k20`, `geo50_climate10`;
- equal island weight;
- climate PC1-4, island area and distance included;
- spatial-block cluster-robust covariance.

The broad `all_observed` layer is a contemporary-flora association. It is not interpreted as native colonisation or historical source reconstruction.

## Coverage

The global primary layer contained **3,252 unique islands with estimable GloBI source-breadth enrichment**.

At the primary 3-reference threshold, approximately:

- northern mid-latitude: **1,814-1,818 islands**;
- northern high-latitude: **281-283 islands**;
- tropical: **908-922 islands**;
- southern extratropical: **180-184 islands**.

Thus all four broad contexts passed the confirmatory island-count gate.

Native support remained geographically incomplete: northern high-latitude had only 11 native-supported islands and southern extratropical about 28-29, so those native cells were not testable at the confirmatory threshold.

## 1. Within-context isolation-breadth slopes

At the primary effort-matched, 3-reference layer:

### Northern mid-latitude

All four source definitions were weakly positive (`+0.0036` to `+0.0095`) but none was FDR-supported (`q=0.244-0.544`).

### Northern high-latitude

Signs were mixed around zero (`-0.0052` to `+0.0139`) and none was supported.

### Tropical

All four source definitions were negative (`-0.0159` to `-0.0206`), but the source-mode FDR values remained just above the frozen threshold (`q=0.0768-0.0965`).

### Southern extratropical

All four source definitions were positive (`+0.0097` to `+0.0260`). One source mode (`geo_k5`) was FDR-supported (`q=0.0125`), while the other three were not (`q=0.132-0.314`).

There is therefore no source-definition-robust within-context breadth filter that can be promoted in any context under the frozen rule.

## 2. Four-context heterogeneity

The direct four-context test asked whether the isolation-breadth slopes differ among contexts.

At the primary 3-reference effort-matched layer:

- `geo50_climate10`: p=`0.03382`, q=`0.04510`;
- `geo_k10`: p=`0.05998`, q=`0.05998`;
- `geo_k20`: p=`0.03207`, q=`0.04510`;
- `geo_k5`: p=`0.01342`, q=`0.04510`.

Thus **3/4 source definitions support context heterogeneity**, but the fourth misses the frozen FDR threshold narrowly. The predeclared classification is therefore:

> **source-definition-sensitive context heterogeneity**

—not a promoted global mechanism.

The pattern behind that test is biologically suggestive: northern mid-latitude and southern extratropical estimates tend positive, tropical estimates tend negative, and northern high-latitude estimates are near zero/mixed. But these signs are descriptive because the robustness gate is not fully passed.

## 3. North-Tropical direct contrast

The direct northern-midlatitude versus tropical distance-slope contrast was also source-definition sensitive:

- `geo50_climate10`: q=`0.0206`;
- `geo_k20`: q=`0.0202`;
- `geo_k10`: q=`0.0555`;
- `geo_k5`: q=`0.0975`.

Therefore only **2/4** source definitions retain the direct contrast after FDR. This is compatible with the four-context result but does not satisfy the robustness rule.

## 4. Area moderation does not explain the global GloBI pattern

The four-context heterogeneity of `distance x area` was unsupported in **0/4** source definitions.

All primary FDR q-values were approximately `0.632`.

Within individual contexts, distance-by-area terms were also generally unsupported. The global GloBI extension therefore does not reproduce the redesigned plant-side H3 area-conditioning as a robust source-genus interaction-breadth mechanism.

## 5. Native sensitivity

Where native support was sufficient, the broad result weakened further.

Northern mid-latitude native slopes were small and source-definition inconsistent; tropical native slopes tended negative but none was FDR-supported. Northern high-latitude and southern extratropical native cells were below the confirmatory support gate.

Consequently the all-observed global signal cannot be promoted as a native island-assembly or pollinator-filter result.

## Integrated H5 interpretation

The global extension changes H5 in a useful way.

1. **The GloBI test is no longer geographically narrow.** It now covers all four broad Chapter 1 contexts on 3,252 islands.
2. **There is suggestive context dependence in independent interaction breadth.** Three of four source definitions support four-context slope heterogeneity, with tropical slopes tending negative while northern-midlatitude/southern-extratropical slopes tend positive.
3. **That heterogeneity is not source-definition robust enough to promote.** The fourth source mode fails the frozen FDR gate.
4. **Area does not rescue the mechanism.** Global `distance x area` heterogeneity is unsupported in all four source modes.
5. **Native evidence remains incomplete globally.** The broad all-observed pattern cannot be interpreted as native assembly.

The correct H5 conclusion is therefore:

> **Independent GloBI interaction breadth shows suggestive but source-definition-sensitive biogeographic heterogeneity across the global island flora, while no robust universal breadth filter or area-conditioned pollinator mechanism is identified.**

This strengthens the broader Chapter 1 conclusion that the island floral response is context dependent, but it still does not identify historical pollinator loss, effective service or a single global pollination mechanism.
