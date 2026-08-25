# Chapter 1 analysis execution surface

## Current canonical model contract

The current Chapter 1 hypothesis is **biogeographically contingent floral and reproductive island syndromes**, not a primary Bombus-deficit pathway.

Canonical scientific definitions live in:

- `config/analysis_models.yml`
- `docs/v2_pollination_regime_framework.md`
- `docs/manuscript_submission_contract.md`
- `THESIS_CHAPTER_POSITIONING.md`

The canonical model ladder is:

```text
M0 biogeographic baseline
M1 universal isolation baseline
M2 isolation × biogeographic context        [primary]
M3 status / lineage / source-pool residual
M4 category-preserving decomposition
```

Pollination-syndrome comparisons are performed only after M0–M4 results are frozen.

## Legacy script in this directory

`run_bayesian_m0_m4_main.R` is a **historical Bombus-primary analysis**. It remains for reproducibility and sensitivity work but is not the current Chapter 1 manuscript workflow.

The associated GitHub Actions workflow has been demoted to manual-only legacy status.

Do not cite its old M0–M4 labels as the current model ladder.

## Current implementation gap

A new canonical fitting workflow should not be written on top of the old Bombus input until the status / lineage support layer is available on current `main`.

The required order is:

1. port the clean floristic-status / genus-lineage utilities from the former PR #135 without its obsolete acquisition assumptions;
2. construct one island-level model table containing:
   - isolation;
   - area / climate;
   - frozen biogeographic context;
   - floristic status / endemicity support;
   - lineage/source-pool expectation or residual;
   - observation/evidence support;
   - category-preserving floral/reproductive outcomes;
3. fit M0–M2 with maximum defensible support per outcome;
4. fit M3 using status/lineage-constrained support;
5. fit M4 category decomposition;
6. freeze regional trait vectors;
7. only then perform literature-based pollination-syndrome concordance.

## Protected boundaries

The new workflow must not:

- require `bombus_deficit` to enter the primary analysis universe;
- infer pollinator guild from flower colour/form;
- treat opportunistic absence as historical pollinator loss;
- collapse all categories to one syndrome index before decomposition;
- omit status/lineage/source-pool safeguards from manuscript-level inference.
