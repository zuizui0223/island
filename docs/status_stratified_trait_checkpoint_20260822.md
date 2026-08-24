# Native non-endemic versus corroborated-endemic trait checkpoint — 2026-08-22

## Purpose

After establishing a robust northern `Isolation -> regional endemism`
association, split direct trait composition into:

- confirmed native non-endemics; and
- WCVP-corroborated regional endemics.

This distinguishes a change in the proportion of endemic lineages from an
additional within-stratum floral/reproductive isolation slope.

Models use the same broad direct-evidence binary M0:

```text
trait successes / direct trait trials
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

with cluster-robust spatial-block uncertainty.

## Native non-endemic stratum

### Northern mid-latitude

| outcome | n islands | isolation beta | p |
|---|---:|---:|---:|
| plain colour | 239 | +0.0316 | 0.173 |
| generalized form | 235 | +0.0041 | 0.916 |
| SC | 239 | +0.0257 | 0.129 |

No broad non-endemic-native outcome shows a supported isolation slope in the
northern regime.

### Tropical

| outcome | n islands | isolation beta | p |
|---|---:|---:|---:|
| plain colour | 138 | +0.0444 | 0.0108 |
| generalized form | 132 | -0.0110 | 0.755 |
| SC | 134 | +0.0339 | 0.321 |

The raw tropical plain-colour slope is the only nominal association. The
all-native genus-composition-preserving checkpoint already showed that the
corresponding tropical plain-colour association disappears after genus turnover
is held fixed. Therefore this raw slope is not interpreted as a pollination or
within-lineage floral response without the stratum-specific lineage sensitivity.

Southern non-endemic support remains below the 50-island target.

## Corroborated endemic stratum

### Colour

Colour is the only endemic outcome currently close to the model-support target:

- northern: n=31, beta -0.0349, p=0.889;
- tropical: n=46, beta -0.0609, p=0.381;
- southern: n=10, status-limited.

There is no evidence in the current direct endemic-colour subset for a simple
plain-colour isolation slope.

### Floral form and SI/SC

Current endemic support is below the pilot threshold:

- northern form n=20; SI/SC n=12;
- tropical form n=27; SI/SC n=27;
- southern much lower.

Some low-support fitted coefficients are numerically large because a few islands
and outcome states approach or reach separation. They are **not** treated as
biological findings. The predeclared support rule routes these outcomes to a
recoverability / targeted-acquisition decision instead.

## Current decomposition

The reliable pattern is presently:

```text
Isolation
  -> regional endemic share increases strongly in northern mid-latitudes
  -> native non-endemic broad floral / SC slopes remain weak
  -> endemic colour shows no simple distance slope on current direct evidence
  -> endemic form and SI/SC are under-supported rather than positive/negative results
```

Thus the strongest current isolation signal is in **floristic status / lineage
composition**, not a demonstrated universal shift of broad floral categories.
The conditional pollination mechanism, if present, must be sought as residual
information beyond this status/lineage structure rather than inferred from the
raw distance gradient.
