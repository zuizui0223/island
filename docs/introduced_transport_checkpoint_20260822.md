# Introduced-flora transport contrast — 2026-08-22

## Role

Introduced/naturalized flora is retained as a secondary human-transport contrast.
It is **not** assumed to be a clean negative control because propagule pressure,
horticultural selection, planting history and environmental filtering can all
be geographically structured.

Inputs are the same frozen sources used for the direct-native checkpoint:

- GIFT origin status run `32555972362`;
- direct species trait run `32551742699`;
- locked distance/area/climate/regime artifact `29228212586`.

Confirmed introduced origin contains 19,442 island-species rows across 129 frozen
status-covered islands.

## Direct support

Direct introduced-trait support is below the 50-island target in northern and
southern regimes, but exceeds it for the tropical subset:

- northern mid-latitude: colour 42, floral form 43, SI/SC 42 islands;
- tropical: colour 63, floral form 59, SI/SC 63 islands;
- southern extratropical: colour 22, floral form 21, SI/SC 21 islands.

Thus only the tropical introduced contrast has confirmatory-count support under
the current status-covered universe. Northern and southern introduced results
remain exploratory.

## Direct introduced M0 isolation slopes

The same grouped-binomial M0 used in the native checkpoint was fitted:

```text
introduced trait successes / introduced direct trait trials
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

with cluster-robust uncertainty over frozen spatial blocks.

| regime | outcome | n islands | isolation beta | robust SE | p |
|---|---|---:|---:|---:|---:|
| northern mid-latitude | plain colour | 42 | +0.0709 | 0.1086 | 0.514 |
| northern mid-latitude | generalized form | 43 | **+0.2101** | 0.0430 | **<0.001** |
| northern mid-latitude | SC | 42 | +0.0132 | 0.0778 | 0.866 |
| tropical | plain colour | 62 | **+0.1943** | 0.0381 | **<0.001** |
| tropical | generalized form | 58 | +0.0506 | 0.0675 | 0.454 |
| tropical | SC | 62 | **-0.2379** | 0.1105 | **0.031** |
| southern extratropical | plain colour | 22 | -0.0590 | 0.0679 | 0.385 |
| southern extratropical | generalized form | 21 | -0.0367 | 0.1106 | 0.740 |
| southern extratropical | SC | 21 | -0.1607 | 0.0847 | 0.058 |

## Interpretation

The introduced stratum is not uniformly flat with isolation. Therefore it cannot
be used as a simplistic `no biological island process` negative control.

The result is still useful: geographic structure among introduced taxa warns that
raw island-distance slopes can also arise through human transport, horticultural
selection and establishment filters. Native mechanism claims must therefore be
supported by status stratification, lineage controls and the predeclared
pollination-regime tests rather than by distance alone.
