# Chapter 1 analysis execution surface

## Canonical question

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

Canonical definitions live in:

- `config/chapter1_when_where_omnibus.yml`
- `config/analysis_models.yml`
- `docs/chapter1_when_where_frozen_result_20260825.md`
- `docs/manuscript_submission_contract.md`
- `THESIS_CHAPTER_POSITIONING.md`

## Canonical inference hierarchy

```text
WHERE
within-context response-vector joint Wald             [primary]

BETWEEN-WHERE
pairwise response-vector difference joint Wald        [primary]

WHEN
persistence across floristic-status strata            [primary boundary]

M3
genus-composition-preserving lineage guardrail        [required]

atomic M0-M4
category-level decomposition of supported vectors     [secondary]
```

The primary Chapter 1 result is not defined by counting individually significant traits.

## Active implementation

- `src/island_v2/chapter1_context_input.py`
  - builds Bombus-free status-aware atomic trait composition;
  - expands colour/form multistates into atomic presences;
  - retains `SC|SI` once as `mixed_or_variable`;
  - fails closed on conflicting direct evidence.
- `src/island_v2/chapter1_when_where_omnibus.py`
  - runs within-context response-vector tests;
  - runs pairwise response-vector difference tests;
  - uses spatial-block cluster-robust covariance;
  - applies declared confirmatory/pilot support gates;
  - classifies floristic-status persistence.
- `src/island_v2/chapter1_context_analysis.py`
  - retains atomic M0/M1/M2 models for decomposition and diagnostic checks.
- `src/island_v2/genus_fixed_trait_null.py`
  - builds the genus-composition-preserving M3 null.
- `src/island_v2/status_stratified_lineage_analysis.py`
  - evaluates broad lineage-controlled outcomes.
- `src/island_v2/chapter1_trait_vector_freeze.py`
  - freezes atomic directions for post-when/where interpretation.

Canonical workflow:

- `.github/workflows/run-chapter1-context-main.yml`

## Frozen run

Canonical when/where run: `32837335384`.

### WHERE

Confirmatory filtering is detected in:

- northern mid-latitude `all_native`;
- northern mid-latitude `native_nonendemic`;
- tropical `all_native`;
- tropical `native_nonendemic`.

### BETWEEN-WHERE

Northern mid-latitude and tropical response vectors differ confirmatorily using 17 common supported atomic responses:

- all native: q = `2.352889e-08`;
- native non-endemic: q = `7.125861e-07`.

### WHEN

Both regional signals persist in native non-endemic flora, so neither is confined to endemic taxa.

Northern high-latitude and southern-extratropical contexts remain unresolved at the confirmatory tier. Southern extratropical shows a pilot vector signal; northern high-latitude currently lacks enough supported responses even for the declared pilot vector test.

See `docs/chapter1_when_where_frozen_result_20260825.md` for exact response counts, island counts, Wald statistics and claim boundaries.

## Atomic interpretation rule

Atomic M0–M4 results answer:

> **Which trait categories make the confirmed northern-versus-tropical vectors different?**

They do not answer WHERE by themselves.

The genus-preserving M3 layer still prevents promotion of a broad `generalized + plain + SC` package into a demonstrated classical island syndrome.

## Legacy analysis

`run_bayesian_m0_m4_main.R` is a historical Bombus-primary analysis. It remains for provenance/reproducibility and is not the canonical Chapter 1 workflow.

## Protected boundaries

The canonical workflow must not:

- require Bombus or any pollinator label as a primary predictor;
- infer pollinator guild from flower colour/form;
- treat opportunistic nondetection as historical loss;
- define WHERE from counts of individual significant traits;
- treat one-region significance as evidence of between-region difference;
- interpret insufficient support as an ecological null;
- omit lineage safeguards from broad syndrome claims.
