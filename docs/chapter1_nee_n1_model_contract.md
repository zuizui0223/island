# Chapter 1 NEE challenge — N1 model contract

Status: `prospective_pre_N1_fit`

This document freezes the inferential N1 model before the functional-channel catalog and island-retention results are inspected.

## Biological question

N1 asks only whether independently measured, source-available pollination channels differ in how their island retention changes with isolation.

It does **not** ask whether those channel differences explain plant genus assembly. That is N2.

Primary unit:

```text
island x pollination channel
```

Only rows with:

```text
source_state = available
channel_state in {retained, disrupted}
```

are eligible.

`retained = 1`, `disrupted = 0`. Structural absence and unresolved observation are excluded rather than recoded as zero.

## Qualification before fitting

A channel enters primary N1 only when the already-frozen channel-input receipt classifies it as `confirmatory` and `N1_gate_eligible`.

Primary N1 requires at least **three confirmatory channels**. Two-channel comparison is not sufficient for the NEE gate. Failure to reach three channels means `N1_not_evaluable_for_NEE_gate`; it does not trigger channel substitution.

## Frozen model comparison

The common-slope model is:

```text
retained
  ~ isolation
  + channel
  + island area
  + climate PC1-PC4
  + observation-effort terms
  + source-region fixed effects
```

The candidate model adds:

```text
isolation x channel
```

The isolation exposure remains the frozen Chapter 1 `log1p_distance_to_continent_km`, interpreted as the composite source-separation/connectivity gradient. A source-specific water-gap exposure is a prespecified sensitivity when available, not a replacement selected after outcomes.

Continuous predictors are standardized on the fixed common analysis support. Sparse source-region levels (<5 primary rows) are pooled to `other_source_region` without looking at channel state.

Observation effort is represented by predeclared transforms of:

- background record count;
- spatial coverage;
- temporal coverage;
- dataset count;
- recency relative to 2026.

An evaluable retained/disrupted row with missing effort fields is a hard failure.

## Primary heterogeneity test

N1 is **not** established by saying one channel has a significant distance slope while another does not.

The primary test is the joint cluster-robust Wald test:

```text
H0: all channel-specific isolation slopes are equal
```

with `(number of qualified channels - 1)` degrees of freedom.

Cluster-robust uncertainty uses the frozen 10-degree Chapter 1 spatial blocks. Individual channel slopes and pairwise slope contrasts are reported only as interpretation after the global gate.

## N1 pass rule

All of the following are required:

1. at least three confirmatory channels;
2. global `isolation x channel` Wald `P < 0.05`;
3. heterogeneous-slope model log likelihood exceeds the common-slope model on exactly the same rows;
4. deleting one spatial block at a time reverses the global conclusion in no more than 20% of evaluable deletions;
5. deleting one source region at a time reverses the global conclusion in no more than 20% of evaluable deletions.

No one channel is required to have a nominally significant slope. No channel can be relabelled, split, combined or substituted after the outcome is seen.

## Claim ceiling

If N1 passes, Chapter 1 may say:

> Independently measured pollination channels have different isolation-retention curves.

It still may **not** say:

- pollinator loss caused plant genus filtering;
- channel retention explains H3;
- one mobile guild functionally replaced another.

Those statements remain blocked until N2 tests `channel retention x lineage dependency` against source-matched genus entry/persistence.

## Stop rule

If N1 fails, the NEE double-filter challenge stops before N2. The frozen Chapter 1 result remains: the broad response is represented at genus-level assembly, while the cause of differential genus representation remains unidentified.
