# Direct-native M0 and genus-fixed lineage checkpoint — 2026-08-22

## Question

Before adding Bombus or endemicity, does isolation by itself predict floral or
self-compatibility composition within source-backed confirmed-native island
floras, and do any apparent slopes exceed what the observed genus composition
already predicts?

This checkpoint uses only:

- GIFT-origin rows from successful status run `32555972362`;
- direct species trait ledger from PR #132 run `32551742699`;
- locked distance/area/climate/regime data from run `29228212586`.

GIFT endemism is disabled in this checkpoint. Only confirmed native origin is
used.

## Binary direct-evidence contrasts

The broad binary contrasts are secondary summaries only; category-preserving
INLA remains responsible for the manuscript category decomposition.

- plain colour: `white` or `green_brown_inconspicuous` versus
  `yellow_orange`, `red_pink`, or `blue_purple`;
- generalized form: `open_radial`, `brush_puff`, or `composite_head` versus
  tubular/bell/funnel/urn/bilabiate/salverform/papilionaceous/spurred forms;
- self compatibility: `SC` versus `SI`.

A multistate record is used only when every state lies on the same side of the
contrast. Cross-side or otherwise ambiguous records are excluded.

## M0

For each regime and outcome the grouped-binomial model is:

```text
native trait successes / native direct trait trials
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

Uncertainty is cluster-robust over the frozen 10-degree spatial block.
Coefficients below are per 1 SD of log distance in the complete analysis set.

| regime | outcome | n islands | spatial blocks | isolation beta | robust SE | p |
|---|---|---:|---:|---:|---:|---:|
| northern mid-latitude | plain colour | 238 | 36 | +0.0232 | 0.0256 | 0.366 |
| northern mid-latitude | generalized form | 234 | 36 | +0.0045 | 0.0365 | 0.902 |
| northern mid-latitude | SC | 238 | 36 | +0.0162 | 0.0168 | 0.335 |
| tropical | plain colour | 136 | 46 | **+0.0380** | 0.0184 | **0.039** |
| tropical | generalized form | 132 | 43 | -0.0314 | 0.0389 | 0.420 |
| tropical | SC | 133 | 45 | +0.0429 | 0.0337 | 0.203 |
| southern extratropical | plain colour | 33 | 23 | +0.0731 | 0.0688 | 0.288 |
| southern extratropical | generalized form | 31 | 21 | +0.0546 | 0.1141 | 0.632 |
| southern extratropical | SC | 32 | 22 | -0.0311 | 0.1439 | 0.829 |

### Immediate interpretation

There is **no direct-native northern raw isolation syndrome** in these broad
binary outcomes. The northern Bombus hypothesis therefore cannot be justified
by treating distance itself as the mechanism.

The only nominal raw isolation association is greater plain-colour composition
in the tropical subset. This is not treated as a final biological result until
lineage composition is tested.

The southern subset is below the 50-island confirmatory-count target and is
status-limited in the current GIFT-covered universe.

## Genus-composition-preserving null

For each direct binary outcome:

1. native island-species memberships are kept fixed;
2. each island's observed genus composition is therefore kept fixed;
3. direct species trait states are permuted only among species of the same genus;
4. one permuted state is shared by a species across all island occurrences in a
   draw;
5. 300 draws are used to estimate the genus-fixed expected island composition;
6. `observed share - genus-fixed null mean` is then modeled with the same
   distance/area/climate structure, weighted by direct species trials and using
   cluster-robust spatial-block uncertainty.

| regime | outcome | n islands | residual isolation beta | robust SE | p |
|---|---|---:|---:|---:|---:|
| northern mid-latitude | plain colour | 238 | -0.00368 | 0.00312 | 0.238 |
| northern mid-latitude | generalized form | 234 | -0.00110 | 0.00189 | 0.560 |
| northern mid-latitude | SC | 238 | +0.00017 | 0.00214 | 0.937 |
| tropical | plain colour | 136 | -0.00062 | 0.00207 | 0.764 |
| tropical | generalized form | 132 | +0.00049 | 0.00240 | 0.838 |
| tropical | SC | 133 | +0.00520 | 0.00405 | 0.199 |
| southern extratropical | plain colour | 33 | +0.00972 | 0.01299 | 0.455 |
| southern extratropical | generalized form | 31 | -0.01285 | 0.01000 | 0.199 |
| southern extratropical | SC | 32 | +0.00096 | 0.00901 | 0.915 |

The tropical raw plain-colour association disappears completely after genus
composition is held fixed. No northern outcome has a residual isolation slope.

## Current scientific consequence

The direct-native evidence supports neither:

```text
isolation -> universal plain/generalized/SC island syndrome
```

nor a northern mechanism inferred from distance alone.

Instead, the next tests must distinguish:

1. whether the native/non-endemic versus corroborated-endemic split changes the
   pattern;
2. whether a Bombus environmental/channel variable adds information beyond the
   weak northern M0 without being relabelled as observed Bombus loss;
3. which multistate colour/form categories move under INLA even when broad
   binaries cancel;
4. whether status/missingness support changes these conclusions.

This is a useful falsification checkpoint, not a failed analysis: the simple
universal-isolation explanation has already become less plausible, and the
conditional mechanism now has to earn support beyond lineage composition.
