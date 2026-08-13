# Cached evidence recovery checkpoint (2026-08-13)

This checkpoint recovers direct trait statements that were already present in
completed, repository-pinned acquisition ledgers but were not present in the
preceding formal public-Web ledger. It does not perform another 106,295-species
search.

Selection was made against strict coverage Run `31667163418` and formal
reviewed-Web Run `31666705636`. A row must be species/synonym-direct High or
Medium evidence, have an exact excerpt, stable source URL and source lineage,
map to the same biological trait, and belong to a currently unresolved strict
axis. Existing species-trait pairs, same-trait value conflicts, cultivar text,
cross-species statements and floral-part ontology mismatches are excluded.

The output contains 235 novel species-trait records after conservative
same-clause colour, whole-flower size, inflorescence-state and reproductive
trait gates. One hundred were expressed
in the line-by-line curated schema and appended to the preceding 1,766 reviewed
records. The package audit is a deterministic 200-row stratified sample; only a
trait with at least ten audited records and precision/cultivar gates passing may
scale. The formal PR #131 implementation remains authoritative for source-lineage
deduplication, direct conflict resolution, `genus × trait_name` Validated Low,
and strict three-axis coverage. Family inference, global fallback, `n=2` formal
inference, axis-only genus joins, and cross-trait reproductive substitutions are
not used.

Rebuild with:

```text
python -m island_v2.cached_evidence_recovery_checkpoint \
  --strict-coverage-csv <run-31667163418>/strict_species_axis_coverage.csv.gz \
  --prior-formal-csv <run-31666705636>/broad_web_medium_evidence.csv.gz \
  --master-csv data/v2/staging/gbif/collected/island_taxa.csv \
  --prior-curated-evidence-csv data/v2/staging/traits/open_web_pilot/independent_lineage_recovery_checkpoint_20260813/combined_curated_evidence_20260813.csv \
  --prior-curated-audit-csv data/v2/staging/traits/open_web_pilot/independent_lineage_recovery_checkpoint_20260813/combined_curated_manual_audit_20260813.csv \
  --output-dir data/v2/staging/traits/open_web_pilot/cached_evidence_recovery_checkpoint_20260813
```
