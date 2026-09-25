# Chapter 1 historical surfaces

This index separates the current submission from preserved replay/provenance surfaces.

## Current
- selector: `config/chapter1_submission_current.json`
- manuscript: `submission/chapter1_current/MANUSCRIPT.md`
- methods/results: `docs/chapter1_corrected_submission_20260924.md`
- corrected tables: `results/geography_20260924/`

## Superseded v14
Preserved for exact provenance only:
- `config/chapter1_v14_canonical_result_lock.json`
- `docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md`
- `docs/chapter1_submission_freeze_v14_20260918.md`
- `.github/workflows/run-chapter1-v14-reordered-hypotheses.yml`

The v14 exposure is superseded by the source-matched coastline correction. Its result lock remains immutable evidence of what was reproduced on 18 September 2026.

## Frozen v13 parent
v13 remains parent scientific provenance beneath v14. It is not an active submission surface.

## Pre-v13
Archived under `legacy/chapter1-pre-v13/`.

## Alpha1 database
The alpha1 bundle retains its original 8,265-unit contract for exact historical database replay. It is not the corrected 8,264-unit Chapter 1 analysis universe.


## Historical workflow boundary

Superseded publication-surface workflows for v8–v14 are retained only for explicit historical replay. They are **manual-only** (`workflow_dispatch`) and must not run on push, pull request or schedule.

The guard is:
- `tests/test_chapter1_historical_workflows_manual_only.py`

The active automatic Chapter 1 submission check is:
- `.github/workflows/corrected-submission-geography.yml`

This keeps old model/figure reproduction available without allowing retired publication surfaces to repromote themselves.
