# Chapter 1 floristic-status resolution audit — 2026-09-15

Status: **completed observation-process diagnostic; not a biological response result**.

## Provenance

- branch: `ch1-all-data-primary`
- workflow run: `34963956873`
- artifact: `chapter1-status-resolution-34963956873`
- artifact ID: `10394441857`
- digest: `sha256:09d0bd8c7d0d69de8e716ac389e54d544405f837e4642f405ff11277c2921586`
- frozen upstream input: Chapter 1 progressive artifact run `34232450884`

## Question

Is floristic-status resolution itself geographically structured across the observed island flora?

For every island, the response is the count of observed island-species rows whose origin status is resolved (`native` or `introduced`) out of all observed island-species rows. The audit uses the same beta-binomial logit framework, distance, island area, climate PC1--PC4 and spatial-block cluster-robust covariance as the all-data biological analysis. No floral trait outcome enters the model.

## Raw coverage

| analysis regime | islands | resolved status rows | observed species rows | resolved fraction |
|---|---:|---:|---:|---:|
| northern mid-latitude | 2,206 | 85,752 | 611,303 | 0.1403 |
| northern high latitude | 424 | 440 | 39,219 | 0.0112 |
| tropical | 1,558 | 71,733 | 325,206 | 0.2206 |
| southern extratropical | 317 | 12,048 | 64,029 | 0.1882 |

The marginal status-resolution fraction is therefore strongly heterogeneous among broad geographic regimes.

## Distance association within regimes

Status resolution itself increases with the distance covariate within every fitted regime:

- northern mid-latitude: slope `+1.1326`, p=`2.63e-7`;
- northern high latitude: `+14.8162`, p=`0.0225`;
- tropical: `+1.4118`, p=`0.00464`;
- southern extratropical: `+1.6521`, p=`0.00598`.

These coefficients describe the observation/status-resolution process, not biological trait change.

## North--Tropical comparison

The direct North--Tropical difference in the distance slope of status resolution is not supported:

- 3,748 islands;
- 190 spatial blocks;
- interaction estimate `+0.0996`;
- cluster-robust SE `0.5231`;
- p=`0.8490`.

Thus the previously observed difference between the all-observed and native-only biological analyses cannot be reduced to a simple claim that the *distance gradient of status resolution* differs between northern mid-latitude and tropical islands.

## Interpretation

This audit sharpens, rather than removes, the floristic-status problem:

1. status resolution is highly incomplete and geographically heterogeneous in level;
2. its distance association is non-zero within broad regimes;
3. but the North--Tropical difference in that observation-process slope is not supported;
4. the same-island diagnostic already shows that changing which species are admitted as native changes the biological result even when the island frame is held fixed.

Therefore a simple island-level inverse-probability correction for status *availability* is unlikely to resolve the core issue. The remaining uncertainty concerns the biological identity of unresolved island-species records (native versus introduced), not merely which islands have more status information.

## Claim boundary

Do not interpret unresolved status as introduced, native, or random missingness. The audit supports explicit status-uncertainty analysis or improved external status data; it does not justify automatic imputation.
