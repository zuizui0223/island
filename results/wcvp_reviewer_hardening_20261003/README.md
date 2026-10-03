# WCVP reviewer hardening — 2026-10-03

## Questions

This audit addresses two limitations of the Chapter 1 WCVP regional-native compatibility sensitivity.

1. Does the complementary regionally incompatible / source-introduced flora carry its own within-region H1 signal?
2. Does the geographic resolution of TDWG Level-3 units covary with island isolation strongly enough to create the regional-native H1 result?

The WCVP native-range summary is the exact frozen artifact used in the manuscript
(run 37089905542). No WCVP update was introduced.

## 1. Complementary partition carries multivariate isolation structure

The combined `regionally_incompatible_or_introduced` partition contains source-backed introduced records plus unresolved records whose mapped TDWG Level-3 unit is not listed in the species' WCVP native range.

Within-region seven-response H1 remained supported in all four regions in the primary all-analysis scope:

| Region | Islands | q |
| --- | ---: | ---: |
| Northern mid-latitude | 818 | 3.46e-6 |
| Northern high latitude | 163 | 0.0205 |
| Tropical | 745 | 0.0110 |
| Southern extratropical | 132 | 0.0118 |

Direct-only also returned q <= 0.0430 in all four regions. The northern-high Direct-only vector is numerically stable but retains an optimizer termination warning for actinomorphic symmetry: increasing the optimizer limit from 1,000 to 5,000 iterations changed neither the fitted slopes nor q (=0.00210165).

Atomic direction is not identical across regions. In the tropical combined partition all seven isolation slopes are positive in both evidence scopes. Northern mid-latitude and northern-high partitions have five positive and two negative atomic slopes; southern extratropical composition is mixed.

The strict source-backed introduced partition is support-limited outside the tropics. In the tropics it is jointly supported, but the atomic response is mixed rather than a uniformly positive H1 direction. Therefore the combined result must not be read as proof that introduced plants alone generate H1.

### Interpretation

The WCVP sensitivity does **not** isolate a native-specific syndrome. Its defensible contribution is narrower:

- regional-native-compatible records independently reproduce H1, so the broad signal is not confined to known introduced records;
- regionally incompatible / introduced records also carry substantial isolation-associated trait structure;
- floristic-status provenance therefore does not by itself identify the biological process generating H1.

## 2. TDWG Level-3 resolution covaries with isolation

Official WGSRPD Level-3 geometry was pinned to TDWG repository commit
`52da7828aba9d461dd133c27b3bd7a4407161f54`, file git-blob SHA1
`91104e5159e31f88154833a51d3b0d1c9271083f`.

Level-3 polygon area was calculated in equal-area projection EPSG:6933. The concern is real in the two southern geographic strata:

| Region | Spearman rho(distance, log L3 area) |
| --- | ---: |
| Northern mid-latitude | +0.040 |
| Northern high latitude | +0.028 |
| Tropical | **-0.621** |
| Southern extratropical | **-0.737** |

The ratio of Level-3 area to focal-island area shows the same qualitative structure and is strongly negatively correlated with isolation in all four strata, especially tropical and southern extratropical islands.

## 3. Adjusting for Level-3 scale does not remove regional-native H1

The WCVP regional-native H1 models were re-fitted with `log_tdwg_l3_area_km2` as an additional covariate, on top of the existing island-area and climate controls.

To separate covariate adjustment from the loss of islands lacking usable Level-3 area, the unadjusted H1 was first re-fitted on the **same Level-3-area-complete island support**. All four vectors were already supported on that matched support:

| Region | Matched-support all-analysis q | Matched-support Direct-only q |
| --- | ---: | ---: |
| Northern mid-latitude | 0.00138 | 0.00226 |
| Northern high latitude | 0.04396 | 0.00195 |
| Tropical | 0.00166 | 0.00142 |
| Southern extratropical | 6.20e-11 | 8.48e-9 |

After adding Level-3 area, all four regional vectors remained supported:

| Region | L3-adjusted all-analysis q | L3-adjusted Direct-only q |
| --- | ---: | ---: |
| Northern mid-latitude | 0.00259 | 0.00257 |
| Northern high latitude | 0.00917 | 1.83e-5 |
| Tropical | 3.36e-6 | 1.14e-5 |
| Southern extratropical | 0.000388 | 0.00634 |

Thus unequal TDWG Level-3 spatial resolution is a real property of the sensitivity design, but neither complete-case support restriction nor explicit Level-3-area adjustment explains away the four-region regional-native result.

## Claim boundary

After this audit the WCVP result should be written as:

> The recurrent H1 pattern is reproducible in a large regional-native-compatible flora and remains supported after adjustment for TDWG Level-3 geographic scale. However, the complementary incompatible/introduced partition also contains strong isolation-associated trait structure, so WCVP partitioning does not establish a native-specific ecological or evolutionary mechanism.

## Provenance

- validation workflow run: **37103668310**
- artifact: **11266114760**
- artifact digest: `sha256:33e63c2e018425b79ed8bf01ef432c84487dd4dfba95de31705db73fa29a33f2`
- WCVP native ranges: frozen run **37089905542**
- TDWG WGSRPD geometry commit: `52da7828aba9d461dd133c27b3bd7a4407161f54`

Compact committed outputs:

- `incompatible_or_introduced_H1.csv`
- `resolution_matched_baseline_H1.csv`
- `resolution_adjusted_H1.csv`
- `tdwg_resolution_distance_audit.csv`
- `direct_high_convergence_audit.json`
