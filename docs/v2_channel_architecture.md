# v2: when / where island floral filtering — analysis architecture

## Status and boundary

Current v2 Chapter 1 is a **when/where boundary-condition analysis** using direct floral/reproductive evidence.

Canonical result: [`chapter1_when_where_frozen_result_20260825.md`](chapter1_when_where_frozen_result_20260825.md).

The primary question is:

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

```text
island geography / climate / frozen context
        ↓
status-aware direct trait composition
        ↓
WHERE: within-context response-vector omnibus
        ↓
BETWEEN-WHERE: response-vector difference
        ↓
WHEN: floristic-status persistence
        ↓
lineage guardrail + atomic decomposition
        ↓
pollination interpretation only after freeze
        ↓
mechanism handed to Chapter 2
```

## Architectural rules

1. **Response vectors precede individual p-values.** WHERE is not defined by counting significant atomic categories.
2. **Regional difference requires a direct vector test.** Significant in A and nonsignificant in B is not evidence that A differs from B.
3. **Insufficient support means unresolved.** Data-poor contexts are not interpreted as ecological nulls.
4. **Status persistence operationalizes the current WHEN boundary.** Persistence in native non-endemic flora shows that the signal is not confined to endemic taxa.
5. **Lineage composition remains a guardrail.** Broad syndrome claims require genus-composition-preserving checks.
6. **Pollinator labels do not enter the Chapter 1 primary model.** They are Discussion-level interpretations after the when/where result is frozen.

## Primary implementation

### WHERE

`src/island_v2/chapter1_when_where_omnibus.py` stacks supported atomic outcomes within each context and status stratum. Each response has its own intercept, baseline effects, and isolation slope. Dependence is retained using spatial-block cluster-robust covariance.

```text
H0: all supported isolation slopes in the focal context = 0
```

### BETWEEN-WHERE

For each context pair, only outcomes meeting the same support threshold in both contexts are used.

```text
H0: isolation-response vector in context A = isolation-response vector in context B
```

### WHEN

Repeat WHERE for `all_native`, `native_nonendemic`, and `endemic`. Persistence in `native_nonendemic` is treated as a boundary condition, not a causal comparison among overlapping strata.

## Support policy

- `<30` islands per retained response: excluded;
- `30–49`: pilot;
- `>=50`: confirmatory count threshold.

## Frozen result

Formal run: `32837335384`.

### Confirmatory WHERE

- **northern mid-latitude:** supported in `all_native` and `native_nonendemic`;
- **tropical:** supported in `all_native` and `native_nonendemic`;
- **northern high-latitude:** unresolved from current response support;
- **southern extratropical:** unresolved confirmatorily; pilot vector signal exists.

### Confirmatory BETWEEN-WHERE

Northern mid-latitude and tropical response vectors differ significantly using 17 common confirmatorily supported responses:

- all native: joint Wald χ² = **69.78**, df = 17, q = **2.35e-08**;
- native non-endemic: joint Wald χ² = **61.02**, df = 17, q = **7.13e-07**.

Thus Chapter 1 has a confirmatory **northern-midlatitude versus tropical multivariate boundary**, while other regional boundaries remain data-limited.

## Lineage and broad-syndrome guardrail

M3 compares broad trait summaries against a genus-composition-preserving expectation. The broad `generalized + plain + SC` package does not survive as a coherent confirmatory classical syndrome in the main northern/tropical strata.

Therefore:

```text
supported regional response vector
!= automatically a named pollination syndrome
```

## Atomic decomposition

The atomic M0–M4 models identify which categories produce the confirmed northern-versus-tropical vector difference. They are decomposition, not the primary WHERE test.

## Interpretation guardrail from the Izu programme

For trait state `z`:

```text
W(z) = F(z) × E(z)
```

Chapter 1 mainly observes `W`, island-flora composition. A response vector can therefore reflect establishment filtering, lineage composition, ecological interactions, reproductive processes, or combinations thereof. Chapter 1 identifies **where** filtering is detectable and where responses differ, but does not identify `F` mechanistically.

## Pollination boundary

Only after the when/where result is frozen may the northern-versus-tropical vector difference be compared with pollination-syndrome literature.

Retained semantic boundaries:

- environmental compatibility != occurrence;
- occurrence != visitation;
- visitation != effective service;
- effective service != reproductive dependency;
- nondetection != historical loss;
- floral phenotype != direct pollinator observation.

## Analysis sequence

1. Freeze island universe, context definition, trait evidence, status strata, support thresholds, and baseline covariates.
2. Build direct-evidence status-aware atomic trait composition.
3. Run WHERE response-vector omnibus tests.
4. Run BETWEEN-WHERE pairwise vector tests.
5. Classify WHEN using floristic-status persistence.
6. Apply M3 lineage guardrails.
7. Use atomic M0–M4 decomposition to characterize supported vectors.
8. Freeze the when/where result.
9. Perform ecological concordance/mismatch discussion only afterward.
10. Hand mechanism to Chapter 2.

## Thesis handoff

```text
Chapter 1 — island
WHEN / WHERE is filtering detectable?
WHERE do response vectors differ?
        ↓
Chapter 2 — izu-core
WHY do supported contexts generate different response architectures?
        ↓
Chapter 3 — shimahotarubukuro
WHAT phenotype axes actually diverge within one lineage?
```

The Chapter 1 contribution is therefore **boundary-condition discovery**, not pollinator causal identification.
