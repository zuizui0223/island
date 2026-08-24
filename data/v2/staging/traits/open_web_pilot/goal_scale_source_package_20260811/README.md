# Goal-scale reviewed source package (2026-08-11)

This checkpoint merges two already acquired but previously separate reviewed
source packages before the all-evidence rebuild:

- Go Botany + BaseFlor synonym recovery + PROTEUS: 1,062 candidate rows;
- BSdb + CSIRO structured traits + Wang et al. (2023): 30,200 candidate rows.

The shared trait-specific review gate approves 31,227 of 31,262 rows.  Its
2,836-row independent/structured audit contains 2,801 accepted rows
(`accepted_correct / all_reviewed = 0.987659`) and no cultivar contamination.
All 35 rejected rows are retained in the audit, not silently promoted.

The package keeps `pollen_vector_mode` and `reward_type` outside the strict
three axes.  It does not use family inference, global fallback, n=2 rules, or
cross-trait substitution.  Source-lineage deduplication and the
`genus x trait_name` join remain the responsibility of the shared all-evidence
formalizer.

Rebuild with:

```powershell
$env:PYTHONPATH = "src"
py -3.13 -m island_v2.merge_reviewed_source_packages `
  --evidence data/v2/staging/traits/open_web_pilot/combined_next_source_package_20260811/combined_next_evidence_20260811.csv.gz `
  --evidence data/v2/staging/traits/open_web_pilot/bsdb_csiro_wang_checkpoint_20260811/bsdb_csiro_wang_symmetry_evidence_20260811.csv.gz `
  --audit data/v2/staging/traits/open_web_pilot/combined_next_source_package_20260811/combined_next_reviewed_audit_600_20260811.csv `
  --audit data/v2/staging/traits/open_web_pilot/bsdb_csiro_wang_checkpoint_20260811/bsdb_csiro_wang_symmetry_independent_audit_20260811.csv.gz `
  --output-dir data/v2/staging/traits/open_web_pilot/goal_scale_source_package_20260811 `
  --output-stem goal_scale_source_package_20260811 `
  --generated-at-utc 2026-08-11T10:30:00Z
```

The manifest pins every input and output hash and records the exact common-gate
summary.
