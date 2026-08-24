# Wang et al. 2023 family-validated symmetry checkpoint

This checkpoint turns a large published species-row table into a conservative
Medium candidate package. It does not accept all 279,877 rows.

Selection gates are:

1. exact species-rank name in the 106,295-species island master;
2. exact source/master family agreement;
3. exclusion of Apiaceae, Asteraceae, Brassicaceae, Cyperaceae, and Poaceae,
   whose symmetry definitions the paper treats specially;
4. at least 10 independently resolved species-direct comparisons per family;
5. at least 0.95 agreement within that family.

Thirty-eight families pass. They contain 28,083 exact master species rows.
The independent audit contains 1,636 species, of which 1,601 agree
(precision 0.978606). It includes 16 binary-state conflicts plus 19 conflicts
with existing non-binary symmetry states; all 35 remain in the audit and are
rejected by the common source-package gate.

The paper describes Data S2 as species-level, but the table does not map each
row to an upstream flora. For that reason every row is Medium, not High. Rows
overlapping existing direct evidence conservatively reuse an existing lineage
instead of being counted as a second independent source. Novel rows share one
dataset lineage, which prevents one compilation from fabricating
leave-one-lineage-out genus support.

Rebuild after downloading the pinned formal direct ledger with:

```bash
python -m island_v2.wang_2023_symmetry_checkpoint \
  --source-csv data/v2/external/traits/wang_2023_symmetry/wang_2023_species_symmetry_snapshot.csv.gz \
  --baseline-direct-csv direct_species_trait_ledger.csv.gz \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/wang_2023_symmetry_checkpoint_20260811
```

No family inference, global fallback, n=2 rule, axis-only join, pollen-vector
substitution, or reward substitution is used.
