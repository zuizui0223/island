# Status-aware category decomposition checkpoint — 2026-08-23

## Purpose

This checkpoint tests whether the broad binary null results were hiding opposing
fine-category responses. It is intentionally run **after** floristic status and
lineage structure are represented.

Inputs:

- corrected GIFT + WCVP floristic status: run `32559322028`;
- direct species-trait evidence: PR #132 run `32551742699`;
- locked isolation / area / climate / spatial covariates: run `29228212586`;
- category + coverage-sensitivity workflow: run `32614027696`;
- artifact: `status-category-decomposition-32614027696`.

Only resolved species/synonym-direct evidence is used. Four colour and four
floral-form categories are retained. For each category, the response is:

```text
observed island category share
- expected category share from the direct trait distribution of the same genera
```

The model is weighted by direct species trials and uses cluster-robust spatial
blocks:

```text
category residual
~ z(log distance to continent)
+ z(log island area)
+ climate PC1 + PC2 + PC3 + PC4
```

Multiple testing is controlled with Benjamini-Hochberg FDR separately within
`stratum x regime x domain` (colour or form). The primary category decision uses
`q < 0.05`. Coefficients from fewer than 50 islands are retained as diagnostics
but are not confirmatory category results.

## Confirmatory-support results (>=50 islands)

Only two category-level isolation slopes survive FDR with at least 50 islands.

| stratum | regime | domain | category | n islands | blocks | beta | robust SE | p | BH q |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| all confirmed native | northern mid-latitude | colour | red/pink | 238 | 36 | **-0.00251** | 0.000863 | 0.00369 | **0.0147** |
| native non-endemic | tropical | colour | red/pink | 136 | 46 | **-0.00410** | 0.00107 | 0.000116 | **0.000465** |

No floral-form category with >=50-island support survives FDR.

## Measurement-coverage sensitivity

Trait resolution itself can vary with isolation, so the two FDR-supported
contrasts were stress-tested without filling missing species. For each island,
`direct_fraction = direct colour trials / species in the relevant floristic
stratum` was calculated. The analysis then used:

1. a continuous direct-fraction covariate; and
2. prespecified direct-fraction subsets at 0.25, 0.30, 0.40 and 0.50 where at
   least 30 islands remained.

### Northern all-native red/pink

| sensitivity | n islands | blocks | beta distance | robust SE | p |
|---|---:|---:|---:|---:|---:|
| baseline | 238 | 36 | -0.00251 | 0.000863 | 0.00369 |
| + continuous direct fraction | 238 | 36 | -0.00250 | 0.000872 | 0.00411 |
| direct fraction >=0.30 | 232 | 35 | -0.00253 | 0.000861 | 0.00334 |
| direct fraction >=0.40 | 188 | 33 | -0.00222 | 0.000827 | 0.00735 |
| direct fraction >=0.50 | 78 | 15 | -0.00202 | 0.00128 | 0.114 |

The coefficient is essentially unchanged after continuous coverage adjustment
and remains supported through the >=40% subset. At >=50% coverage the estimate
retains the same direction and similar magnitude but precision drops sharply
with only 78 islands / 15 spatial blocks. This is therefore a **small but
comparatively coverage-robust colour turnover**.

Its magnitude remains small: roughly 0.2–0.25 percentage points of residual
red/pink category share per 1-SD increase in standardized log isolation.

### Tropical native-nonendemic red/pink

| sensitivity | n islands | blocks | beta distance | robust SE | p |
|---|---:|---:|---:|---:|---:|
| baseline | 136 | 46 | -0.00410 | 0.00107 | 0.000116 |
| + continuous direct fraction | 136 | 46 | -0.00447 | 0.00108 | 0.000034 |
| direct fraction >=0.30 | 123 | 42 | -0.00445 | 0.00115 | 0.000113 |
| direct fraction >=0.40 | 65 | 26 | -0.00090 | 0.00422 | 0.831 |

The tropical effect survives continuous coverage adjustment and the >=30%
subset, but it collapses rather than merely losing precision in the >=40%
subset. The >=40% subset still spans a broad isolation range, so this contrast
is classified as **measurement-coverage sensitive**, not a stable manuscript
result at this stage.

### Northern interpretation

The all-native northern red/pink residual decreases slightly with isolation.
However, the corresponding native-nonendemic estimate is smaller and does not
survive FDR (beta -0.00155, p=0.0194, q=0.0776; n=238).

Therefore this is not evidence for a general native-colonist red/pink filtering
rule. The effect is instead a small all-native compositional signal that remains
after genus expectation and direct-coverage adjustment. It should not be
assigned to Bombus: the separate Bombus checkpoint (`32613425082`) found no
supported distance-by-compatibility interaction and no incremental Bombus value
for generalized form or SC.

### Tropical interpretation

The tropical native-nonendemic red/pink residual initially survives genus
control and FDR, but it is not robust to the higher-coverage subset. Thus it
cannot currently support a claim that tropical isolation drives loss of red/pink
flowers, nor can it be used to infer failure or success of bird/butterfly
functional replacement.

The broad tropical plain-colour residual remains null. At most, the current data
suggest a possible narrow colour-compositional turnover that requires stronger
measurement support before interpretation.

## Under-supported diagnostics

Several nominal/FDR signals occur below the confirmatory support target and are
not promoted to biological findings:

- northern endemic blue/purple: n=31, beta +0.0942, q=0.0073;
- northern endemic tubular/trumpet: n=20, beta -0.0984, q=0.0455;
- southern non-endemic blue/purple / yellow-orange: n=33;
- southern all-native yellow-orange / blue-purple: n=33;
- southern endemic colour signals: n=10.

The northern endemic colour stratum is already in the targeted-acquisition zone;
the northern endemic form result remains below pilot support. Southern estimates
are status-limited and cannot be rescued by trait acquisition alone.

## Scientific consequence

The combined evidence now separates three levels:

1. **Floristic status / lineage composition:** northern isolation robustly predicts
   greater regional endemicity.
2. **Broad floral syndrome:** after status and genus composition are represented,
   broad plain/generalized/SC isolation slopes are approximately null.
3. **Fine categories:** one small northern all-native red/pink decline remains
   reasonably robust to direct-trait coverage; the tropical non-endemic red/pink
   signal is coverage-sensitive. No confirmatory floral-form category remains
   after FDR.

Therefore the manuscript should not be framed as a recovered universal island
floral syndrome. The stronger current result is that isolation restructures the
floristic/endemic composition, with limited colour-specific turnover remaining
after lineage and measurement controls.

## INLA decision

The old PR #115 category-preserving INLA analysis used `sensitivity_all` and
predates the status/lineage correction, so it is not reinstated as canonical.

A new heavy all-category INLA model is not justified as a mechanism-rescue step.
If a final category model is desired, it should be narrowly focused on the
coverage-robust northern red/pink contrast and must retain the same status,
lineage and evidence boundaries.

## Acquisition boundary

- no broad native colour/form campaign;
- northern all-native red/pink already has adequate direct support and needs no
  further broad acquisition;
- tropical non-endemic red/pink is measurement-sensitive, so any additional
  acquisition should specifically improve high-isolation / lower-coverage
  tropical non-endemic colour support rather than increase global colour counts;
- targeted acquisition remains relevant for northern endemic colour (31/59
  status ceiling) and tropical endemic colour (45/70 complete-case support in
  this checkpoint; 46 direct-supported islands in the status checkpoint);
- endemic form and SI/SC remain recoverability questions before any larger
  campaign;
- southern endemic outcomes remain status-limited.
