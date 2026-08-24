# BSdb species-level floral-symmetry checkpoint (2026-08-11)

This checkpoint reuses the immutable Zell et al. (2025) BSdb release at
`dirtyplants/BSdb@9e87946d1e3121d39e657b702cf9f92ccc10936e` without treating
the BSdb redistribution as a second biological source.

Only `FlrSymCheck` rows with a non-empty `trait.source` are admitted. The
source documentation defines that column as individually collected
species-level data. Rows supported only by `genus_trait_source` or
`family_trait_source` are excluded, as are infraspecific names, family
conflicts, already resolved species-traits, and `wind` (a pollination mode,
not a floral-symmetry state).

After exact species-rank TNRS/master-list and family gates, 877 records for
843 species remain: 720 actinomorphic and 157 zygomorphic. A reproducible
hash sample of 200 structured records was reviewed; all passed identity,
value, provenance, and cultivar checks. This is an extraction audit of the
published database rows, not a new biological remeasurement.

Rebuild with:

```bash
python -m island_v2.bsdb_morphology_checkpoint \
  --source-csv Zell_df_12_29_23.csv \
  --baseline-direct-csv direct_species_trait_ledger.csv.gz \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --reviewed-audit-csv bsdb_morphology_manual_audit_200_20260811.csv \
  --output-dir bsdb_morphology_checkpoint_20260811
```

The manifest pins all inputs and output hashes. Formal promotion still uses
the common source-package gate and the shared `genus x trait_name` Validated
Low rebuild.
