# Final Chapter 1 directional inference — 2026-10-03

## Decision

The former high-dimensional joint Wald tests are not used as confirmatory evidence for H1.

The final H1 uses the seven frozen, pre-oriented v14 indicators and a one-dimensional
classic-island directional score. The three biological domains receive equal total weight:

- reproductive assurance: 1/3 total;
- colour dulling: 1/3 total;
- accessibility/generalization: 1/3 total.

Flower size and inflorescence display remain in the raw-state descriptive audit and are
not added to the confirmatory score after seeing results.

## H1 regional results

All four regional point estimates are positive, but northern mid-latitudes are weak.

| scope | north-mid | north-high | tropical | south-extra |
|---|---:|---:|---:|---:|
| all-analysis estimate | 0.0161 | 0.1137 | 0.0978 | 0.0698 |
| finite-cluster one-sided P | 0.134 | 0.00211 | 1.83e-8 | 0.00649 |
| wild-cluster P | 0.145 | 0.0001 | 0.0001 | 0.0013 |
| Direct-only estimate | 0.0193 | 0.1153 | 0.0816 | 0.0696 |
| finite-cluster one-sided P | 0.106 | 0.00254 | 8.31e-5 | 0.00494 |
| wild-cluster P | 0.119 | 0.0001 | 0.0001 | 0.0019 |

Therefore the strict four-region intersection-union recurrence claim is **not supported**:

- all-analysis: P_t = 0.1345; P_wild = 0.1452;
- Direct-only: P_t = 0.1062; P_wild = 0.1190.

## H1a — global-average direction

Paule-Mandel random effects with modified Hartung-Knapp inference:

- all-analysis: mean = **0.06906**, SE = 0.02191, one-sided P = **0.02560**;
- Direct-only: mean = **0.06349**, SE = 0.01957, one-sided P = **0.02384**.

Thus the data support a positive global-average classic-island direction, not a
region-invariant universal syndrome.

## H1b — regional heterogeneity

- all-analysis: Q = **16.568**, P = **0.000867**, I2 = **0.819**;
- Direct-only: Q = **9.511**, P = **0.02321**, I2 = **0.685**.

Regional realization is therefore a primary result rather than residual noise.

## Sensitivities

Changing score weights does not remove the positive global-average result. Equal-indicator
and reproductive-assurance-plus-accessibility weighting both retain one-sided P < 0.032
in both evidence scopes. These are sensitivities, not post-hoc replacements for the
primary equal-domain score.

Dropping the single spatial block contributing the largest contrast variance in every
region preserves H1a in all-analysis (mean = 0.09291, P = **0.04383**) but weakens
Direct-only slightly beyond 0.05 (mean = 0.08786, P = **0.05326**). H1a should therefore
be described as positive but leverage-sensitive.

## Raw-state role

The three-axis raw-state analysis remains useful for describing floral/reproductive
reorganization, provenance differences and regional phenotype detail. It is
direction-free and high-dimensional and therefore cannot be used to rescue a failed
classic-island directional test. In particular, raw-state reorganization in strict
native tropical flora does not convert the earlier directional null into confirmatory
support.

## H2-H4 finite-cluster audit

- H2: five of sixteen primary H2b context-by-response cells remain FDR-supported after
  finite-spatial-block t inference. Conditional decomposition is not causal mediation.
- H3: corrected global pollen-limitation slope = **0.09191**, finite-publication
  two-sided P = **0.01594**. Offshore-only robustness slope = **0.22031**, P = **0.02459**.
- H4: exact H2 reproductive-assurance score beta = **-0.29830**, finite-publication
  P = **0.00417**; accessibility score beta = **-0.29566**, P = **0.02334**.
  H4 remains explicitly post-hoc functional triangulation.

## Claim ceiling

The submission may claim:

> Geographic isolation is associated with a positive average shift toward reproductive
> assurance and accessible/generalized floral architecture, but the magnitude and
> phenotypic realization of that tendency are strongly region dependent. The same
> geographic gradient is independently associated with greater pollen limitation, and
> the island-associated reproductive-assurance and accessibility states are associated
> with lower current pollen limitation.

It may not claim:

- independent classic-syndrome support in all four regions;
- a universal trait checklist;
- causal mediation from historical pollen limitation to present traits;
- named pollinator loss or replacement from floral phenotype;
- within-lineage evolution rather than assemblage sorting/persistence;
- that raw-state omnibus significance rescues the directional H1.

## Provenance

- branch: `analysis/h1-final-directional-inference`
- validation run: **37124405765**
- artifact: **11274313934**
- artifact digest:
  `sha256:05389ffe438876c11527ecb8b6be7c43ac3f45b14fa343970826e142e36dab5b`
