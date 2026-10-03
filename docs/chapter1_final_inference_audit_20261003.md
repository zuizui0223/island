# Chapter 1 final inference audit — 2026-10-03

## Decision

The submission-facing analysis should no longer claim that a classic floral island
syndrome is independently supported in all four regions.  The high-dimensional joint
Wald tests are retained only as historical/reference outputs.

The defensible result is narrower and more informative:

1. **H1a — positive global-average direction:** supported under the predeclared
   directional alternative when the four prespecified regional effects are synthesized
   with Paule–Mandel heterogeneity and modified Hartung–Knapp inference.
2. **H1b — regional heterogeneity:** clearly supported.
3. **Strict four-region recurrence:** not supported because the northern-midlatitude
   directional score is positive but individually weak.
4. **Raw-state three-axis analysis:** useful for describing how colour, structure and
   reproductive assurance reorganize, but not a confirmatory direction-free syndrome
   test.
5. **H2–H4:** primary conclusions survive finite-cluster reference distributions.

## Why this is not a post-result hypothesis replacement

The approved v13 design dated 2026-09-17 explicitly defined H1 as:

> Is there a positive global common component in the isolation response toward the
> predeclared classic island-syndrome direction?

It also stated that regional coefficients may differ in magnitude and exact equality is
not required.  The final H1 therefore returns to that pre-result estimand.  v14 later
encoded seven pre-oriented indicators across reproductive assurance, colour dulling and
accessibility/generalization.  The present analysis collapses those already-oriented
coefficients to a one-degree-of-freedom contrast and repairs finite-cluster inference; it
does not add a result-selected trait state.

## H1 primary directional score

The seven v14 indicators are retained:

- reproductive assurance: self-compatibility, selfing mating system, autonomous selfing;
- colour: plain colour;
- accessibility/generalization: generalized form, actinomorphy, shallow/open tube.

Indicators are averaged within domain and the three domain means receive equal total
weight.  Positive coefficients point toward the predeclared classic island-syndrome
direction.

Inference uses:

- beta-binomial outcome models with the corrected geography and frozen covariates;
- spatial-block cluster-robust covariance;
- Student-t reference with G-1 degrees of freedom;
- deterministic linearized Rademacher sign-flip as a small-cluster sensitivity.

### Regional estimates

| scope | northern mid-latitude | northern high-latitude | tropical | southern extratropical |
|---|---:|---:|---:|---:|
| all-analysis | 0.0161 (p1=0.1345) | 0.1137 (p1=0.00211) | 0.0978 (p1=1.83e-8) | 0.0698 (p1=0.00649) |
| Direct-only | 0.0193 (p1=0.1062) | 0.1153 (p1=0.00254) | 0.0816 (p1=8.31e-5) | 0.0696 (p1=0.00494) |

All eight regional point estimates are positive.  The strict intersection-union test
fails because the northern-midlatitude p value is the limiting component:

- all-analysis IUT p = 0.1345; sign-flip p = 0.1452;
- Direct-only IUT p = 0.1062; sign-flip p = 0.1190.

Therefore the paper must not state that all four regions independently support the
classic syndrome.

## H1a — global-average common component

The four regional effects are treated as prespecified replication strata and synthesized
with Paule–Mandel between-region variance plus modified Hartung–Knapp uncertainty
(k=4).

- **all-analysis:** random-effects mean = 0.0691, SE = 0.0219,
  one-sided p = 0.0256, two-sided p = 0.0512;
- **Direct-only:** random-effects mean = 0.0635, SE = 0.0196,
  one-sided p = 0.0238, two-sided p = 0.0477.

Because the positive direction was specified in the v13 design before this reanalysis,
the one-sided test is the hypothesis-aligned primary reference.  The two-sided values
must also be reported because the all-analysis result is borderline under a
direction-agnostic reference.

Weighting sensitivities preserve the positive global-average conclusion:

- equal weighting of all seven indicators: p1 = 0.0288 all-analysis, 0.0277 Direct-only;
- reproductive-assurance + accessibility core only: p1 = 0.0316 in both scopes.

The core-two-domain sensitivity is consistent with the earlier v13 biological framing in
which reproductive assurance and accessibility/generalization were the two main
pathways and colour was a context-dependent display component.

## H1b — regional heterogeneity

Regional heterogeneity is clear in the primary sandwich-based analysis:

- all-analysis: Q = 16.57, df = 3, p = 0.000867, I² = 0.819;
- Direct-only: Q = 9.51, df = 3, p = 0.0232, I² = 0.685.

An exact delete-one-spatial-cluster jackknife stress test preserved H1a but made H1b less invariant. Re-synthesizing the regional point estimates with exact jackknife SEs gave:

- H1a all-analysis: mean = 0.06352, SE = 0.02243, one-sided p = 0.03305;
- H1a Direct-only: mean = 0.05481, SE = 0.01862, one-sided p = 0.03018;
- H1b all-analysis: Q = 11.86, p = 0.00789, I² = 0.747;
- H1b Direct-only: Q = 5.63, p = 0.131, I² = 0.467.

Several northern-high leave-one-cluster refits failed numerically, so the exact jackknife is a conservative stress test rather than the sole primary estimator. The appropriate interpretation is therefore **a positive global-average island-syndrome direction with regionally variable realization**, while formal evidence for heterogeneity is strongest in all-analysis and not invariant to every evidence scope/uncertainty estimator.

## Influence sensitivity

The northern-high-latitude fit has a large maximum cluster variance share (0.81
all-analysis; 0.84 Direct-only).  Removing the highest-influence spatial block from each
region does not reverse any regional estimate.

The global-average result after this leaveout is:

- all-analysis p1 = 0.0438;
- Direct-only p1 = 0.0533.

Thus H1a is not perfectly influence-robust in the Direct-only scope and should be
described as moderate evidence rather than a definitive universal law.  H1b
heterogeneity becomes stronger after leaveout.

## Raw-state three-axis analysis

The raw-state analysis remains valuable for biological description:

- reproductive-assurance composition reorganizes with isolation;
- structural composition reorganizes with isolation;
- colour is more geographically contingent;
- strict tropical native and known introduced reproductive-assurance response vectors
  are nearly opposed (cosine -0.758 all-analysis; -0.853 Direct-only).

However, high-dimensional direction-free raw-state joint Wald p values are not used as
confirmatory H1 evidence because the number of response states is large relative to the
number of spatial clusters in some regions.  A significant raw-state omnibus means
composition changed; it does not mean the change followed the classic syndrome
direction.

The strict source-backed tropical native directional score remains unsupported
(all-analysis p1 = 0.499; Direct-only p1 = 0.398).  Raw-state reorganization therefore
must not be used to relabel that directional null result as support.

## H2 finite-cluster audit

The same spatial-block cluster-robust SEs were re-evaluated with G-1 Student-t
references and the frozen eight-test BH family within evidence scope.

FDR-supported H2b cells are:

- all-analysis northern-high accessibility: q = 0.00243;
- all-analysis tropical accessibility: q = 0.0191;
- all-analysis southern plain colour: q = 0.00243;
- Direct-only northern-high accessibility: q = 0.01098;
- Direct-only southern plain colour: q = 1.8e-5.

Direct-only tropical accessibility remains positive but not FDR-supported
(q = 0.1267).  H2 therefore remains a conditional decomposition, not causal mediation.

A separate exact delete-one-spatial-cluster jackknife of the continuous H2 pathways preserves the same core pattern under one-sided Holm correction: all-analysis accessibility is supported in northern high latitudes (adjusted p = 0.00326) and the tropics (0.0181); Direct-only accessibility is supported in northern high latitudes (0.0144) but not in the tropics (0.0985).

## H3 finite-publication audit

Using publication-cluster t inference:

- primary global pollen-limitation gradient:
  beta = 0.09191, two-sided p = 0.01594;
- no-zero-constant sensitivity:
  beta = 0.09089, two-sided p = 0.01891;
- supplemental-only:
  beta = 0.04410, p = 0.311;
- post-hoc offshore-only continuous gradient:
  beta = 0.22031, p = 0.02459.

H3 remains an isolation-associated reproductive constraint.  It does not identify
historical pollinator decline or mediation.

## H4 finite-publication audit

The exact-species post-hoc bridge remains supported:

- reproductive-assurance score:
  beta = -0.29830, publication-t two-sided p = 0.00417;
- accessibility/generalization score:
  beta = -0.29566, publication-t two-sided p = 0.02334.

These results demonstrate functional compatibility with lower current pollen limitation,
not that historical pollen limitation selected the present island trait distributions.

## Reviewer-facing claim ceiling

The manuscript can support:

- a positive global-average directional island-syndrome tendency;
- strong regional heterogeneity in how that tendency is realized;
- a recurrent core centered on reproductive assurance and floral
  accessibility/generalization, with colour more context dependent;
- independent association of isolation with experimental pollen limitation;
- post-hoc exact-species functional alignment of the two core plant-response families
  with lower current pollen limitation.

The manuscript must not claim:

- a universal four-region classic phenotype;
- that direction-free raw-state reorganization proves the classic syndrome;
- direct pollinator loss, abundance decline or named pollination-syndrome mechanisms;
- causal mediation from isolation through pollen limitation to trait evolution;
- within-lineage evolution rather than species sorting, colonization, persistence or
  taxonomic assembly;
- a native-only global replication;
- that the corrected distance analysis was prospective.

Existing lineage/source, floristic-origin, trait-missingness and species-list
incompleteness audits remain important boundary analyses.  They should constrain
interpretation rather than be presented as causal identification.


## Exact-jackknife sensitivity provenance

- workflow run: **37122292456**
- artifact: **11273547268**
- artifact digest: `sha256:52329a2fbe68520417693d2572510958adf1a93d02fc06cafd39768ed0b2ad34`
- caveat: several northern-high leave-one-cluster model refits failed numerically, so successful-refit jackknife summaries are stress-test evidence rather than a replacement for the primary finite-cluster sandwich/sign-flip analysis.
