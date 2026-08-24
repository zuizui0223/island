# Status-stratified genus-fixed lineage checkpoint — 2026-08-23

## Purpose

This checkpoint asks whether an isolation-associated broad trait shift remains
after two structural filters are made explicit:

1. floristic status (`native_nonendemic` versus WCVP-corroborated regional
   `endemic`); and
2. the observed genus composition of each island.

The lineage null is a measurement-preserving randomization, not trait
imputation. Direct species trait states are permuted only among species within
the same genus; island species memberships and therefore island genus
composition stay fixed.

Real-data workflow:

- run: `32613129024`
- artifact: `status-stratified-lineage-32613129024`
- status input: official WCVP-corroborated run `32559322028`
- direct trait input: PR #132 run `32551742699`
- covariates: locked run `29228212586`
- null draws: 100
- minimum direct trial per island: 1

The fitted response is:

```text
observed direct trait share - genus-fixed null mean
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

with direct species trials as weights and cluster-robust uncertainty over the
frozen spatial blocks.

## Native non-endemic stratum

This is the cleanest test of whether isolation is associated with a broad trait
shift among confirmed native colonists after genus turnover is held fixed.

| regime | outcome | n islands | blocks | residual isolation beta | robust SE | p |
|---|---|---:|---:|---:|---:|---:|
| northern mid-latitude | plain colour | 238 | 36 | -0.00329 | 0.00279 | 0.238 |
| northern mid-latitude | generalized form | 234 | 36 | -0.00176 | 0.00225 | 0.432 |
| northern mid-latitude | SC | 238 | 36 | +0.00053 | 0.00222 | 0.811 |
| tropical | plain colour | 136 | 46 | +0.00051 | 0.00210 | 0.809 |
| tropical | generalized form | 131 | 43 | +0.00093 | 0.00231 | 0.688 |
| tropical | SC | 133 | 45 | +0.00394 | 0.00399 | 0.323 |
| southern extratropical | plain colour | 33 | 23 | -0.00225 | 0.01237 | 0.856 |
| southern extratropical | generalized form | 30 | 20 | +0.00196 | 0.00560 | 0.726 |
| southern extratropical | SC | 32 | 22 | -0.00127 | 0.00827 | 0.878 |

### Interpretation

No native-nonendemic outcome has a supported residual isolation slope.
Importantly, the tropical raw plain-colour association seen before lineage
control disappears within the non-endemic stratum itself. It therefore cannot
be used as evidence for a tropical isolation-driven floral shift.

The southern estimates are also null, but their 30–33-island support remains
below the 50-island confirmatory target and is status-limited.

## All-native cross-check

The all-native lineage residuals remain null as well:

- northern plain colour: beta -0.00310, p=0.320;
- northern generalized form: beta -0.00168, p=0.390;
- northern SC: beta +0.00038, p=0.866;
- tropical plain colour: beta -0.00044, p=0.833;
- tropical generalized form: beta +0.00028, p=0.916;
- tropical SC: beta +0.00649, p=0.109.

This independently reproduces the earlier all-native lineage checkpoint with
the status-stratified workflow.

## Corroborated endemic stratum

Direct endemic trait support is much thinner. The support class must therefore
be read before any coefficient.

| regime | outcome | n islands | blocks | residual isolation beta | p | support class |
|---|---|---:|---:|---:|---:|---|
| northern mid-latitude | plain colour | 31 | 13 | +0.01487 | 0.774 | targeted acquisition zone |
| northern mid-latitude | generalized form | 20 | 12 | +0.07591 | 0.0268 | below pilot support |
| northern mid-latitude | SC | 12 | 7 | -0.14731 | 0.0266 | below pilot support |
| tropical | plain colour | 45 | 28 | -0.02050 | 0.0596 | targeted acquisition zone |
| tropical | generalized form | 26 | 18 | -0.04264 | 0.0964 | below pilot support |
| tropical | SC | 26 | 17 | -0.11026 | 0.128 | below pilot support |
| southern extratropical | plain colour | 10 | 9 | +0.07182 | 0.0600 | below pilot support |

The apparently large northern endemic form and SC coefficients are numerical
diagnostics from only 20 and 12 islands. They are **not biological findings**.
The predeclared support rule routes them to recoverability testing before any
interpretation or broad acquisition campaign.

Endemic colour is the only endemic outcome currently near a usable support
range. Neither northern nor tropical endemic colour shows a clear residual
isolation slope on current direct evidence.

## Current decomposition

The evidence now supports the following narrower structure:

```text
Isolation
  -> regional endemic share increases robustly in northern mid-latitudes
  -> broad native-nonendemic floral / SC composition does not shift with distance
     once genus composition is held fixed
  -> endemic colour does not show a clear residual distance slope on current data
  -> endemic form and SI/SC remain under-supported, not positive or negative results
```

Therefore the strongest robust isolation-associated signal is presently
**floristic status / endemic-lineage composition**, not a universal broad floral
syndrome.

This does not reject a conditional pollination-channel mechanism. It changes the
burden of proof: Bombus or another channel variable must add information beyond
area, climate, distance, floristic status, and genus composition. Distance alone
can no longer stand in for pollinator loss.

## Acquisition boundary

- northern endemic colour: continue targeted acquisition only if it adds new
  supported islands / isolation support;
- tropical endemic colour: same, with a smaller gap to the 50-island target;
- northern/tropical endemic form and SI/SC: run recoverability diagnostics first;
  do not launch a global fill campaign;
- southern endemic outcomes: status-limited, so trait acquisition cannot solve
  the current ceiling.

Probabilistic genus inference remains a separate measurement sensitivity. It is
not substituted for this direct-evidence lineage checkpoint.
