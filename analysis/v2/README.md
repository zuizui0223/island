> **HISTORICAL / SUPERSEDED EXECUTION SURFACE.** This directory documents the pre-corrected WHEN/WHERE Chapter 1 program. It is retained for replay/provenance and is **not** the current paper pipeline. Use `config/chapter1_submission_current.json`, `submission/chapter1_current/MANUSCRIPT.md`, and `docs/PAPER_PIPELINE.md` for the corrected H1–H4 surface.

# Historical Chapter 1 WHEN/WHERE analysis

## Historical question

> **When and where was isolation-associated floral/reproductive filtering detectable, and where did multivariate response vectors differ?**

The retained definitions are:

- `config/chapter1_when_where_omnibus.yml`;
- `config/analysis_models.yml`;
- `docs/chapter1_when_where_frozen_result_20260825.md`;
- `docs/manuscript_submission_contract.md`.

These files describe an earlier analysis stage. They do not override the current corrected submission selector.

## Historical inference hierarchy

```text
WHERE
within-context response-vector joint Wald

BETWEEN-WHERE
pairwise response-vector difference joint Wald

WHEN
persistence across floristic-status strata

M3
genus-composition-preserving lineage guardrail

atomic M0-M4
category-level decomposition
```

The frozen WHEN/WHERE run was `32837335384`. It supported northern-midlatitude and tropical response vectors and their direct difference under the then-current geography and trait surface.

That result remains provenance. It is **not** the present Chapter 1 headline because the submission was subsequently reorganized around the corrected 8,264-unit H1–H4 surface.

## Retained implementation

The following modules remain available because they are useful for historical replay and diagnostics:

- `src/island_v2/chapter1_context_input.py`;
- `src/island_v2/chapter1_when_where_omnibus.py`;
- `src/island_v2/chapter1_context_analysis.py`;
- `src/island_v2/genus_fixed_trait_null.py`;
- `src/island_v2/status_stratified_lineage_analysis.py`;
- `src/island_v2/chapter1_trait_vector_freeze.py`.

There is **no active canonical Chapter 1 workflow selected from this directory**. Pre-corrected scientific runners are preserved as manual-only historical replay where retained.

## Historical interpretation boundary

The old WHEN/WHERE program did not permit:

- inferring pollinator guild from flower colour/form;
- treating opportunistic non-detection as historical pollinator loss;
- defining a region difference from significance in one region and non-significance in another;
- interpreting insufficient support as an ecological null;
- claiming within-lineage evolution from cross-sectional assemblage composition.

Those restrictions remain useful provenance, but the current scientific hierarchy and numerical baseline are defined elsewhere.

## Current source of truth

Use, in order:

1. `config/chapter1_submission_current.json`
2. `submission/chapter1_current/MANUSCRIPT.md`
3. `docs/chapter1_corrected_submission_20260924.md`
4. `results/geography_20260924/`
5. `docs/PAPER_PIPELINE.md`
