# Chapter 1 Database 1.0 public subset

Database 1.0 remains the immutable scientific analysis database. Its frozen species-axis SHA-256 is `a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`.

The public subset is a rights-filtered derivative, not a replacement analysis database. It is built from the canonical post-PR172 cell-rights audit:

- workflow run: `34566271397`
- artifact: `chapter1-database-v1-cell-rights-34566271397`
- artifact id: `10186293306`
- artifact digest: `sha256:6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a`
- resolved Database 1.0 cells evaluated: `222,688`
- redistributable cells included: `46,274`
- blocked/review-required cells omitted: `176,414`
- semantic mismatch cells omitted: `1`

The selection rule is exact: include a species × axis cell only when `CELL_RELEASE_RIGHTS.csv.gz` reports `release_status == redistributable`.

No trait value, quality tier, H1-H5 contract, or source provenance is rewritten. Missing public cells mean "not authorized for redistribution under the frozen rights audit", not biological missingness and not low scientific quality.

The public compilation is distributed under CC BY-SA 4.0 as the conservative umbrella licence introduced in PR #172. Row-level upstream provenance and source-specific licence/attribution obligations remain authoritative.

Do not use the public subset as the canonical analysis input for the Chapter 1 paper. Analyses remain pinned to immutable Database 1.0; the public subset exists only to make the currently redistributable portion citable and inspectable without silently dropping or relicensing blocked evidence.
