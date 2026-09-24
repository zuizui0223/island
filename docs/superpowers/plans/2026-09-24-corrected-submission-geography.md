# Corrected submission geography implementation plan

> Implement inline with superpowers:executing-plans; use test-driven development for the geometry and admission changes.

**Goal:** Make the verified 24 September geography correction the publication baseline, with reproducible code, original/corrected results and an immutable audit of the superseded v14 exposure.

**Architecture:** Correct GSHHG sibling admission at source. Add a separately versioned spherical coastline-distance module and a replay entry point that validates locked input hashes, retains original model specifications and records corrected tables. Promote a new submission contract and manuscript while retaining frozen v14 files as parent provenance.

**Tech stack:** Existing Python package, NumPy/SciPy, GeoPandas/Shapely, pytest, Ruff, GitHub Actions. Python >=3.11 as declared by the repository.

**Spec:** User explicitly requested on 24 September 2026 that the inflated-zero correction be the submission baseline, not merely a sensitivity. Exact original source archive and analytical inputs are retained. No outcome-dependent exclusion or threshold selection.

## Constraints and review focus

- Recombine GSHHG `sibling_id` (including split 0-E/0-W) before mainland area filtering; reject missing/invalid sibling linkage rather than admitting a continental fragment.
- Spherical radius 6371.0088 km, explicit distinction from WGS84 ellipsoid; ignore artificial dateline seams, not physical coastlines; guard degenerate and invalid segments.
- Retain true mainland point distances at zero; do not substitute island boundary distance for GloPL site coordinates.
- Exclude the one locked-universe fragment only on source identity and exact geometry evidence; record the mapping and reason.
- Preserve traits, models, strata, region definitions, controls, clustered covariance and multiplicity families. Keep weakened/null results.
- Reject missing, duplicated or unverified inputs; no silent fallback to old projected distances.
- Current publication docs/config must point to corrected results. Old locks must not be silently overwritten.

## Tasks

- [x] Add failing admission tests for sibling-linked mainland fragments and invalid linkage; implement correction and run source tests.
- [x] Add failing spherical-distance tests for crossings, endpoint minima, nearshore gaps, dateline seams, interior mainland points and exhaustive-candidate equivalence; implement reusable geometry module and validate real input replay.
- [x] Add a frozen-input correction/refit runner and manifest with all H1/H2, raw colour/architecture, H3 and H4 comparisons; test missing-input/provenance gates and verify the reported scientific outputs.
- [x] Promote corrected submission config, complete manuscript/results explanation, and README/PAPER_PIPELINE links; retain prior v14 surface as superseded provenance. Add consistency checks for current-vs-parent roles.
- [ ] Run focused tests, lint, reproducible result checks and CI; commit only task paths, push branch, open and attach PR. Report remote commit and CI state accurately.
- [ ] Finish the corresponding corrected A0 PPTX/PDF and scientific package, inspect native render and verify all numbers/claims against the corrected publication baseline.

## Validation ledger

- Original sibling tests: 5 expected failures, then 8 passing after correction. Spherical and hash-gate tests failed before implementation, then passed. Final focused suite: 48 passing.
- Full repository-module replay (Python 3.13) completed; 26 numeric tables / 10,664 rows match the published corrected tables at atol=1e-6, rtol=1e-5. Independent numerical geometry validation: 24 cases, maximum discrepancy 1.13e-8 km.
- Fresh-context whole-branch review found no remaining important issues after adding the exact locked-GPKG SHA gate.
- Ruling: preserve frozen v14 manuscript verbatim as historical provenance; select a self-contained corrected methods/results document and explicitly prohibit submitting the old manuscript unchanged. This avoids silently rewriting a frozen document while making corrected results primary.
