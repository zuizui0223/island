# BSdb + CSIRO + Wang reviewed source package (2026-08-11)

This package appends the family-validated Wang et al. (2023) floral-symmetry
rows to the prior BSdb + CSIRO checkpoint without changing any existing
source lineage.

- Candidate evidence: 30,200 rows.
- Independent or structured audit: 2,236 rows.
- Accepted audit rows: 2,201 (precision 0.984347).
- Explicitly rejected Wang conflicts: 35.
- Rows selected by the common trait-specific package gate: 30,165.

The Wang component contributes 28,083 Medium species rows from 38 approved
families. It is intentionally not an independent lineage when the species
already has a direct symmetry record, and all novel Wang rows share one
dataset lineage because per-row upstream sources are unavailable. This keeps
source-lineage leave-out validation conservative.

Formal promotion must use
`bsdb_csiro_wang_symmetry_evidence_20260811.csv.gz` and
`bsdb_csiro_wang_symmetry_independent_audit_20260811.csv.gz` through the
shared High > Medium > Validated Low implementation. The earlier 2,117-row
package remains committed as its own reproducible checkpoint.
