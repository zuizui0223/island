# Northern Bombus lineage-residual checkpoint — 2026-08-23

## Purpose

This checkpoint asks whether the existing species-specific Bombus environmental
hypervolume adds information about northern floral/reproductive composition
after the major structural alternatives already identified in v2 are represented.

The response is not raw island composition. It is:

```text
observed direct trait share - genus-composition-preserving null mean
```

for source-backed confirmed native floras. The same analysis is repeated for
native non-endemics, removing the corroborated-endemic stratum that carries the
strongest robust isolation signal.

Real-data workflow:

- run: `32613425082`
- artifact: `northern-bombus-lineage-residual-32613425082`
- lineage-null input: run `32613129024`
- locked covariate / Bombus compatibility input: run `29228212586`
- regime: `northern_midlatitude`
- spatial blocks: 36 in the complete analysis sets

## Predictor boundary

`environmental_compatibility_max` is the existing species-specific Bombus
hypervolume summary. It is climatic-environmental compatibility only. It is not
observed occurrence, abundance, visitation, pollination service, or historical
Bombus loss.

Because the hypervolume is climate-derived, a Bombus interpretation requires it
to add out-of-sample information beyond the flexible climate baseline.

Three nested models were therefore predeclared for this checkpoint:

```text
M0: distance + area + climate PC1-4
M1: M0 + Bombus environmental compatibility
M2: M1 + distance x compatibility
```

M2 is specifically a sensitivity test of the "remote but environmentally
compatible" arrival-limitation idea. Model value is judged with deterministic
spatial-block cross-validation, not coefficient significance alone.

## Spatial-block cross-validation

### All confirmed native

| outcome | M0 RMSE | M1 change vs M0 | M2 change vs M0 | decision |
|---|---:|---:|---:|---|
| generalized form | 0.02495 | **-0.00118** | -0.00095 | worse with Bombus terms |
| plain colour | 0.02948 | **+0.00034** | +0.00032 | very small improvement |
| SC | 0.03424 | **-0.00010** | -0.00069 | no improvement |

### Native non-endemic

| outcome | M0 RMSE | M1 change vs M0 | M2 change vs M0 | decision |
|---|---:|---:|---:|---|
| generalized form | 0.02587 | **-0.00062** | -0.00045 | worse with Bombus terms |
| plain colour | 0.02793 | **+0.00048** | +0.00031 | small improvement |
| SC | 0.03407 | **-0.00007** | -0.00054 | no improvement |

Positive change means lower RMSE than M0. The plain-colour improvement is only
about 1.1% for all native and 1.7% for native non-endemics. Generalized form and
SC do not improve out of sample.

## Focal coefficients

### M1: add compatibility

| stratum | outcome | compatibility beta | p |
|---|---|---:|---:|
| all native | generalized form | +0.00327 | 0.224 |
| all native | plain colour | **+0.00675** | **0.029** |
| all native | SC | +0.00205 | 0.655 |
| native non-endemic | generalized form | +0.00453 | 0.0785 |
| native non-endemic | plain colour | **+0.00606** | **0.0418** |
| native non-endemic | SC | +0.00264 | 0.555 |

The only nominal M1 coefficients are for plain colour, and they are **positive**:
higher Bombus environmental compatibility is associated with more plain colour
in this residual model. That is opposite the simple `Bombus deficit -> plain`
prediction if compatibility is interpreted as stronger channel availability.
The small cross-validated improvement therefore cannot be used as support for
that deficit hypothesis.

### M2: distance x compatibility

No distance-by-compatibility interaction is nominally supported:

- all-native generalized form: beta -0.00344, p=0.122;
- all-native plain colour: +0.00376, p=0.187;
- all-native SC: -0.00514, p=0.121;
- non-endemic generalized form: -0.00297, p=0.234;
- non-endemic plain colour: +0.00360, p=0.261;
- non-endemic SC: -0.00473, p=0.131.

M2 also fails to improve spatial CV relative to M0 for generalized form and SC,
and improves plain colour less than the simpler M1.

## Scientific decision

The existing Bombus environmental hypervolume **does not recover the predicted
northern Bombus-deficit pathway after floristic status and genus composition are
accounted for**.

Specifically:

- there is no incremental predictive support for generalized form;
- there is no incremental predictive support for SC;
- the only small predictive gain is plain colour, but its compatibility direction
  is opposite the simple deficit prediction;
- the `distance x compatibility` arrival-limitation sensitivity is unsupported.

Therefore this proxy must not be promoted to evidence that isolation causes a
Bombus-channel deficit. The hypervolume remains useful as an environmental
compatibility diagnostic / sensitivity layer.

This is a falsification result, not a pipeline failure. At the present evidence
level, the robust manuscript signal is stronger for:

```text
Isolation -> regional endemicity / floristic turnover
```

than for:

```text
Isolation -> Bombus environmental deficit -> broad plain/generalized/SC shift
```

## Next analysis boundary

1. Do not acquire Bombus occurrence data solely to rescue the hypervolume result;
   observation evidence is too uneven to serve as a clean direct channel measure.
2. Keep endemic colour acquisition targeted because that stratum is close to the
   50-island target; test recoverability before acquiring endemic form / SI-SC.
3. Run category-preserving colour/form decomposition only as a test for broad
   binary cancellation, not as a way to override the failed incremental Bombus
   checkpoint.
4. If no category-level, status-aware residual signal survives, position the main
   result around isolation-associated floristic/endemic turnover and the failure
   of a universal or simple Bombus-mediated floral syndrome.
