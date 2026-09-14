# Chapter 1 V6 species-detection tipping result — 2026-09-14

## Decision

V6 closes the remaining species-list observation-bias gap with a **bounded falsification analysis**, not an occupancy estimate. The analysis asks how strong jointly distance-dependent flora-list incompleteness and trait-state-dependent species recording would have to be to erase or reverse the frozen Chapter 1 accessibility results.

The canonical run reproduces the frozen PR142 baseline at `OR_D = 1` to numerical precision and never increases regression information weight when hypothetical unrecorded species are added.

## Sensitivity model

For the focal generalized/accessibility state,

`OR_D = odds(recorded | focal state) / odds(recorded | nonfocal state)`.

List completeness varies with distance as

`logit C_i = logit(C0) + log(OR_C) * z_distance`.

The frozen grid was:

- median completeness `C0`: `0.50, 0.67, 0.80, 0.90, 0.95`;
- distance-completeness odds ratio `OR_C`: `1.0, 0.8, 0.667, 0.5, 0.25`;
- state-recording odds ratio `OR_D`: `1.0, 0.8, 0.667, 0.5, 0.333, 0.2, 0.1`.

`OR_C < 1` means more distant islands are assumed to have lower flora-list completeness. `OR_D < 1` means generalized/accessibility-positive species are assumed to be under-recorded relative to nonfocal species.

## Canonical receipt

- frozen source run: `34232450884`;
- frozen source artifact: `10058653212`;
- source digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`;
- V6 workflow run: `34800498716`;
- V6 artifact: `10331282464`;
- V6 artifact digest: `sha256:095506a214143a2312bcb1d115644ab21d2d130b7813261019bc386206fd3849`;
- headline scenario rows: `2100`;
- tipping-surface rows: `300`;
- maximum `OR_D=1` score discrepancy: `3.33e-16`;
- hypothetical species increase regression precision: `false`.

## Main result

### Palearctic accessibility

Across `2 evidence scopes × 2 floristic strata × 25 completeness surfaces = 100` baseline-supported conditions:

- `99 / 100` had **no break anywhere on the frozen OR_D grid down to 0.1**;
- all `80 / 80` conditions with `OR_C < 1` had no break;
- therefore, under the specific concern that distant islands are more under-surveyed and generalized/accessibility-positive species are preferentially missed, the sensitivity correction generally strengthens rather than generates the positive Palearctic slope.

The sole break occurred for `direct-only × all-native × C0=0.50 × OR_C=1`. The baseline slope remained positive but FDR support crossed 0.05 at approximately a twofold state-dependent recording difference (`OR_D=0.5`). This exception does **not** correspond to distance-dependent under-survey.

### Tropical accessibility

The tropical accessibility component is materially more sensitive.

Among the `75` completeness surfaces in which the baseline tropical accessibility result was supported:

- `40 / 75` had a break on the frozen grid;
- under strong distance-dependent completeness decline, breaks can occur at state-recording differences in the `1.5–5×` range;
- under `OR_C=1` (no distance trend in flora-list completeness), the baseline-supported tropical result survives the frozen OR_D grid.

Thus the relevant vulnerability is **joint** trait-selective recording and isolation-dependent inventory incompleteness, not trait-selective recording alone.

### Northern-midlatitude vs tropical vector difference

Among the `75` completeness surfaces in which the frozen North–Tropical vector difference was baseline-supported:

- `70 / 75` had no break;
- `5 / 75` broke;
- all breaks occurred under the strongest distance-completeness decline (`OR_C=0.25`).

The direct regional vector difference is therefore more robust than the tropical accessibility component considered alone.

## Existing observation-process diagnostic

The older observation model concerns whether an island has **any flora record**, not species-list completeness, so it is descriptive context only and is not used to estimate `C0` or `OR_C`.

- northern-midlatitude odds ratio for any flora record per distance SD: `0.7068`, p=`0.0659`;
- tropical combined odds ratio: `0.9382` (no standalone combined p-value stored).

These values do not prove inventory completeness and do not calibrate the V6 sensitivity grid.

## Scientific consequence

V6 supports an asymmetric interpretation:

1. The Palearctic accessibility response is difficult to explain away by the specified class of joint distance-dependent inventory incompleteness and focal-state under-recording. Under the biologically motivated `OR_C < 1` concern direction, the bias is conservative for the positive Palearctic response across the entire frozen grid.
2. The tropical negative accessibility component is more observation-process-sensitive and should not be presented as an equally robust opposite trajectory.
3. The broader North–Tropical response-vector difference remains substantially more stable than the tropical accessibility slope alone.

## Claim ceiling

V6 can quantify the amount of **specified** joint species-list incompleteness and state-dependent recording needed to erase or reverse the frozen accessibility result. It does **not** estimate the true number of missing species, prove GBIF/flora inventory completeness, identify the actual detection mechanism, or establish robustness to arbitrary taxon-dependent omission or arbitrary joint species-detection and trait-resolution missingness.

The family→genus attenuation result remains a separate piece of evidence about assemblage hierarchy; V6 does not prove that taxonomic attenuation cannot arise from taxon-dependent observation bias.
