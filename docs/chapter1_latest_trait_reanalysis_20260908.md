# Chapter 1 latest-trait reanalysis — 2026-09-08

## Inputs and canonical boundary

This reanalysis uses the latest public reproducible trait checkpoint from source-scale integration Run `34191508045`:

- filled species×axis cells: **222,688 / 318,885**
- reproductive-assurance cells: **48,497 / 106,295**
- cumulative reviewed source-batch gain: **+206**
- total gain from the recovered public baseline: **+313**
- integration loss: **0**

The canonical Chapter 1 analysis remains Bombus-free. The primary estimand is context-dependent isolation-associated change in multistate floral/reproductive composition, with Bombus/pollination syndromes reserved for post-freeze interpretation.

Core latest-data analysis Run: `34210169033`  
Artifact: `chapter1-latest-trait-34210169033` (`10049591876`)

Robustness Run: `34210793330`  
Artifact: `chapter1-latest-robustness-34210793330` (`10049795988`)

## Primary WHERE results

All four confirmatory WHERE tests are supported.

| stratum | context | responses | minimum islands per response | joint Wald chi-square | df | q | supported |
|---|---|---:|---:|---:|---:|---:|---|
| all_native | northern_midlatitude | 21 | 53 | 134.0014 | 21 | 1.8222e-18 | yes |
| all_native | tropical | 18 | 50 | 267.2185 | 18 | 5.0471e-46 | yes |
| native_nonendemic | northern_midlatitude | 21 | 53 | 245.8791 | 21 | 5.5133e-40 | yes |
| native_nonendemic | tropical | 18 | 50 | 194.2089 | 18 | 1.4367e-31 | yes |

Thus isolation-associated multivariate filtering is confirmatorily detectable in both northern-midlatitude and tropical contexts and persists in native non-endemic flora.

Northern high latitude remains not testable at current support. Southern extratropical retains a pilot signal but is not confirmatorily testable.

## BETWEEN-WHERE results

The northern-midlatitude and tropical multivariate isolation-response vectors differ confirmatorily in both primary strata.

| stratum | common responses | joint Wald chi-square | df | q | supported |
|---|---:|---:|---:|---:|---|
| all_native | 18 | 88.1590 | 18 | 3.0858e-11 | yes |
| native_nonendemic | 18 | 102.8081 | 18 | 6.7558e-14 | yes |

This supports biogeographically structured filtering rather than a single universal island floral syndrome.

## Trait-domain decomposition

The supported multivariate vectors are not the same biological syndrome in the two regions.

### Northern midlatitude

The strongest confirmatory domain is **floral form**:

- all_native: domain q = **0.003337**
- native_nonendemic: domain q = **0.000218**

Flower colour and self-incompatibility are not independently supported domains in the northern-midlatitude confirmatory analysis.

Representative nominal atomic slopes per SD isolation include:

- bell/campanulate form: **+0.2449** log-odds, p=0.00122 (all_native)
- tubular: **+0.1898**, p=0.00355
- spurred: **+0.3737**, p=0.00967
- salverform: **−0.1179**, p=0.0122
- red/pink colour: **−0.1011**, p=0.00242

These atomic slopes are descriptive contributors to the supported response vector; individual p-values do not replace the omnibus test.

### Tropical

The confirmatory domain signal is concentrated in **flower colour** and **self-incompatibility / compatibility composition**:

- flower colour all_native: q across domains = **4.1e-05**
- flower colour native_nonendemic: **0.000218**
- self-incompatibility all_native: **0.003414**
- self-incompatibility native_nonendemic: **0.002628**

Floral-form domain support is weaker: not supported after across-domain correction in all_native or native_nonendemic.

A representative reproductive atomic result is tropical `SC`, all_native:

- isolation slope = **+0.1130** log-odds per SD
- p = **0.00897**

Thus the latest data do **not** support a simple story in which northern temperate islands alone become self-compatible/generalized while tropical islands do not respond. Instead, both regions show isolation-associated filtering, but different components of the floral/reproductive phenotype carry the signal.

## M3 lineage guardrail

The genus-composition-preserving M3 residual does not recover a coherent broad `generalized_form + plain_colour + self_compatibility` syndrome.

For all_native flora:

- northern_midlatitude generalized_form: p=0.0677
- northern_midlatitude plain_colour: p=0.2677
- northern_midlatitude self_compatibility: p=0.8127
- tropical generalized_form: p=0.9674
- tropical plain_colour: p=0.4588
- tropical self_compatibility: positive distance slope, p=0.0412

Therefore broad syndrome claims remain inappropriate. The stronger result is multivariate, category-preserving, context-dependent filtering.

## Observation / trait-resolution support

Across the full geographic universe:

- islands: **8,265**
- flora recorded: **4,505**
- any direct trait evidence: **425 islands**
- direct evidence for all three primary traits: **406 islands**

Direct all-three support by region:

- northern_midlatitude: **235 islands**
- tropical: **133**
- southern_extratropical: **31**
- northern_high_latitude: **7**

Coverage-adjusted analyses reproduce the headline result in both all_native and native_nonendemic strata.

## Robustness

### Information weighting

All headline tests replicate under every tested information-weight scheme:

- canonical
- cap_100
- cap_50
- cap_20
- equal_island

For every mode and both primary strata, northern WHERE, tropical WHERE, and northern-vs-tropical BETWEEN-WHERE remain supported.

### Distance definition

All headline tests replicate for:

- log1p distance
- sqrt distance
- raw distance

All scenarios retain the full 8,265-island geographic universe.

### Leave-one-spatial-block deletion

For each primary stratum, 84 block-deletion runs were testable and **84/84 reproduced the complete headline result**.

- northern WHERE supported fraction: **1.0**
- tropical WHERE supported fraction: **1.0**
- BETWEEN-WHERE supported fraction: **1.0**

## Updated interpretation

The latest trait additions strengthen support without changing the central Chapter 1 result:

> Geographic isolation is associated with floral/reproductive filtering in more than one biogeographic context, but the multivariate response differs strongly among contexts. The northern-midlatitude signal is primarily floral-form structured, whereas the tropical signal is more strongly expressed in colour and compatibility composition. Consequently, the evidence supports conditional biogeographically structured filtering, not a universal classic island floral syndrome.

The main remaining data limitation is not instability of the northern/tropical headline result. It is lack of confirmatory support in northern-high-latitude and southern-extratropical contexts, plus incomplete mechanistic identification of the ecological channels generating the regional response vectors.
