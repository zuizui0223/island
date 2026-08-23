# Status-aware category decomposition checkpoint — 2026-08-23

## Purpose

This checkpoint tests whether the broad binary null results were hiding opposing
fine-category responses. It is intentionally run **after** floristic status and
lineage structure are represented.

Inputs:

- corrected GIFT + WCVP floristic status: run `32559322028`;
- direct species-trait evidence: PR #132 run `32551742699`;
- locked isolation / area / climate / spatial covariates: run `29228212586`;
- real category workflow: run `32613833400`;
- artifact: `status-category-decomposition-32613833400`.

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

### Northern interpretation

The all-native northern red/pink residual decreases slightly with isolation.
However, the corresponding native-nonendemic estimate is smaller and does not
survive FDR (beta -0.00155, p=0.0194, q=0.0776; n=238).

Therefore this is **not evidence for a general native-colonist red/pink filtering
rule**. The difference between all-native and confirmed non-endemic strata means
that endemic or endemism-unresolved components may contribute to the all-native
signal. Its magnitude is also small: about 0.25 percentage points of residual
category share per 1-SD increase in the standardized log-distance predictor.

This result does not rescue the simple Bombus-deficit hypothesis. The separate
Bombus checkpoint (`32613425082`) found no supported distance-by-compatibility
interaction and no incremental Bombus value for generalized form or SC.

### Tropical interpretation

The tropical native-nonendemic red/pink residual decreases with isolation even
after the genus expectation is removed. Other colour categories do not survive
FDR, and the broad plain-colour residual itself remains null.

Thus the earlier binary result was hiding a **narrow colour-compositional
turnover**, not a coherent plain-versus-showy syndrome. Importantly, this result
runs against a simple prediction that mobile bird/butterfly replacement should
maintain or increase red/pink composition on more isolated tropical islands.
Without island-level alternative-pollinator availability, it should be reported
as a compositional pattern rather than assigned to a pollination mechanism.

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

1. **Floral status / lineage composition:** northern isolation robustly predicts
   greater regional endemicity.
2. **Broad floral syndrome:** after status and genus composition are represented,
   broad plain/generalized/SC isolation slopes are approximately null.
3. **Fine categories:** a small number of colour-specific residual shifts remain,
   chiefly red/pink declines in northern all-native and tropical non-endemic
   floras. No confirmatory floral-form category remains after FDR.

Therefore the manuscript should not be framed as a recovered universal island
floral syndrome. The stronger current result is that isolation restructures the
floristic/endemic composition, with limited category-specific colour turnover
remaining after lineage control.

## INLA decision

The old PR #115 category-preserving INLA analysis used `sensitivity_all` and
predates the status/lineage correction, so it is not reinstated as canonical.

A new heavy category model is justified only for the two supported colour
contrasts above, or if targeted endemic-colour acquisition materially increases
support. There is currently no reason to run a large all-category INLA model to
rescue a failed broad syndrome.

## Acquisition boundary

- no broad native colour/form campaign;
- tropical non-endemic red/pink already has adequate direct support and requires
  modelling/robustness, not acquisition;
- northern all-native red/pink already has adequate direct support;
- targeted acquisition remains relevant for northern endemic colour (31/59
  status ceiling) and tropical endemic colour (45/70 complete-case support in
  this checkpoint; 46 direct-supported islands in the status checkpoint);
- endemic form and SI/SC remain recoverability questions before any larger
  campaign;
- southern endemic outcomes remain status-limited.
