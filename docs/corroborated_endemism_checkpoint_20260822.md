# WCVP-corroborated endemism checkpoint — 2026-08-22

## Corrected status layer

Official Kew WCVP bulk run `32559322028` completed the full pipeline:

- official `wcvp.zip` names/distribution tables;
- exact accepted-species matching against the frozen island taxon master;
- valid native distribution records only (`introduced=0`, `extinct=0`,
  `location_doubtful=0`);
- official TDWG WGSRPD level-3 geometry;
- frozen GSHHG island geometry;
- independent GIFT endemic source claim.

WCVP coverage of the frozen master:

- target accepted species: 115,328;
- exact accepted WCVP species: 101,101;
- species with a valid native TDWG-L3 range: 100,887;
- species native to exactly one TDWG-L3 region: 47,537.

The resolved 198,559 GIFT-origin `island x species` rows are now classified as:

- native: 150,531;
- introduced: 19,442;
- origin unresolved: 28,586.

Among all rows, conservative endemism status is:

- non-endemic native: 118,429;
- TDWG-L3-corroborated endemic native: 17,444;
- endemism unresolved: 62,686.

`endemic` here means **GIFT source-reported endemic + exactly one matching WCVP
native TDWG-L3 region**. It is a regional endemism class. It must not be called
single-island endemism unless the TDWG region itself corresponds to the island.

## Status ceilings

Corrected endemic flora is present on:

- northern mid-latitude: 59 islands;
- tropical: 70 islands;
- southern extratropical: 13 islands.

Thus northern and tropical endemic analyses can in principle reach the 50-island
confirmatory-count target. Southern endemic analysis is status-limited before
trait evidence is considered.

## Direct trait support in the corrected endemic stratum

### Northern mid-latitude

| outcome | >=1 direct | >=3 | >=5 | >=10 |
|---|---:|---:|---:|---:|
| plain-colour contrast | 31 | 19 | 14 | 10 |
| generalized-form contrast | 20 | 6 | 5 | 1 |
| SI/SC | 12 | 8 | 3 | 0 |

### Tropical

| outcome | >=1 direct | >=3 | >=5 | >=10 |
|---|---:|---:|---:|---:|
| plain-colour contrast | 46 | 24 | 19 | 17 |
| generalized-form contrast | 27 | 12 | 11 | 6 |
| SI/SC | 27 | 10 | 5 | 2 |

### Southern extratropical

The endemic-status ceiling itself is 13 islands, so all three outcomes are
status-limited and should not trigger trait acquisition to reach 50 islands.

## Acquisition decision

This produces the first genuinely model-driven acquisition boundary:

- northern endemic colour: `31 / ceiling 59` -> targeted trait acquisition zone;
- tropical endemic colour: `46 / ceiling 70` -> targeted trait acquisition zone;
- northern endemic form: `20 / ceiling 59` -> below pilot support; test recoverable
  endemic species before any larger campaign;
- tropical endemic form: `27 / ceiling 70` -> same;
- northern endemic SI/SC: `12 / ceiling 59` -> same;
- tropical endemic SI/SC: `27 / ceiling 70` -> same;
- southern endemic outcomes -> status-limited, not trait-limited.

Most corroborated endemic species occur in only one focal island. Consequently,
additional acquisition does not have the old genus-amplified / widespread-species
leverage: one species usually unlocks at most one currently unsupported endemic
island. Acquisition must therefore be island-specific and recoverability-aware.

## Isolation -> endemism

Endemism is modeled among endemism-resolved native species only:

```text
endemic / (endemic + confirmed native non-endemic)
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

with cluster-robust spatial-block uncertainty. The primary support threshold is
>=30 endemism-resolved native species per island.

| regime | n islands | blocks | isolation beta | robust SE | p |
|---|---:|---:|---:|---:|---:|
| northern mid-latitude | 222 | 35 | **+0.902** | 0.191 | **2.4e-6** |
| tropical | 104 | 39 | **+1.207** | 0.499 | **0.0156** |
| southern extratropical | 20 | 13 | +0.058 | 0.915 | 0.949 |

The northern coefficient corresponds to about `exp(0.902) = 2.46` times the odds
of regional endemism per 1-SD increase in log isolation, conditional on the
specified area/climate baseline. The tropical estimate is larger but less
precise.

This is an island-level floristic association. It does not by itself distinguish
differential colonization, extinction, lineage turnover, or in-situ speciation.

## Measurement warning and resolution sensitivity

Endemism-resolution fraction decreases with isolation:

- northern mid-latitude: Spearman rho about -0.318, p about 1.4e-6;
- tropical: rho about -0.509, p about 3.4e-8;
- southern extratropical: rho about -0.528, p about 0.017.

Therefore the endemism model must retain an explicit status-resolution
sensitivity rather than treating the resolved denominator as missing at random.

### Minimum endemism-resolution fraction

With >=30 resolved native species still required:

| minimum resolved fraction | north n | north beta | north p | tropics n | tropics beta | tropics p |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 222 | +0.902 | 2.4e-6 | 104 | +1.207 | 0.0156 |
| 0.70 | 204 | +0.676 | 0.00014 | 81 | +1.385 | 0.0150 |
| 0.80 | 195 | +0.641 | 0.0010 | 72 | +1.429 | 0.0215 |
| 0.90 | 172 | +0.638 | 3.0e-6 | 58 | +1.203 | 0.0369 |
| 0.95 | 139 | **+0.514** | **0.0016** | 45 | -0.250 | 0.069 |

The northern isolation-endemism association remains positive and supported even
among islands with >=95% endemism resolution. The tropical result is not robust
to the strictest resolution filter and drops below the 50-island count at that
threshold. Therefore:

- northern `Isolation -> regional endemism` is the current robust result;
- tropical `Isolation -> regional endemism` is provisional and measurement-
  sensitivity dependent;
- southern inference remains under-supported.

## Consequence for the main decomposition

The current evidence now separates two facts:

1. isolation strongly changes the **composition/status structure** of northern
   native floras through greater regional endemicity;
2. after genus composition is held fixed, the broad direct-native plain/form/SC
   isolation slopes are approximately null.

That combination makes lineage/endemic turnover a serious explanation for raw
island-flora trait patterns and raises the evidentiary bar for any residual
pollination-channel mechanism.
