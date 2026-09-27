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


The manual-only boundary also covers the retired effect-fingerprint / P1 paired-support audits and the pre-v13 Figure 3 renderers. These remain reproducible but cannot reactivate from historical branch pushes.

## Superseded scientific-design documents still retained on main

These files remain in place because historical runners or provenance documents reference their paths, but they are **not active scientific contracts**:

- `docs/chapter1_progressive_analysis_contract.md` — old universal-vs-branching H1–H5 ladder;
- `docs/chapter1_explanation_gap_validation_plan.md` — old prospective falsification plan for that ladder;
- `docs/chapter1_progressive_pollination_concordance.md` — old sampled guild-architecture concordance layer;
- `config/chapter1_global_branching.yml` — old branching configuration;
- `.github/workflows/validate-chapter1-explanation-gap-contract.yml` — historical replay only, now manual-only.

The current scientific hierarchy is the corrected H1–H4 surface selected by `config/chapter1_submission_current.json`.

## Retired automatic Chapter 1 analysis workflows

The pre-corrected all-data / H3 / H5 analysis workflows that previously reacted to pushes on `ch1-all-data-primary`, together with the old latest-trait / robustness checkpoint workflows, are now **manual-only historical replay**. The current automatic Chapter 1 check remains `.github/workflows/corrected-submission-geography.yml`.

The manual-only invariant is enforced by `tests/test_chapter1_historical_workflows_manual_only.py`, which now covers these retired analysis workflows in addition to the earlier v8–v14 publication/render workflows.

The following active-looking paths are also explicitly historical and must not be used to select the present paper analysis:

- `analysis/v2/README.md`;
- `config/chapter1_progressive_analysis.yml`;
- `config/chapter1_all_data_progressive_analysis.yml`;
- `docs/CHAPTER1_DATABASE_RELEASE.md` when it refers to the progressive H1–H5 dispatcher;
- `config/chapter1_database_versions/README.md` as an execution selector.

Trait-database manifests remain valid provenance objects; they simply no longer select the scientific analysis contract. The current selector is `config/chapter1_submission_current.json`.

## Retained pre-corrected working manuscripts and freezes

These files remain on their original paths because historical runners, figure locks and provenance notes reference them. They are **not current submission documents** and now carry an explicit historical banner:

- `docs/chapter1_manuscript_full_v7_submission_order_20260909.md`;
- `docs/chapter1_manuscript_v8_reframing_20260913.md`;
- `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`;
- `docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md`;
- `docs/chapter1_v9_submission_figure_sync_20260915.md`;
- `docs/chapter1_submission_freeze_20260915_p1_defended.md`;
- `docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`;
- `docs/chapter1_submission_freeze_20260915_p1_p2_defended.md`;
- `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`;
- `docs/chapter1_latest_trait_reanalysis_20260908.md`.

Historical numerical values in those documents, including the old 8,265-unit geography, are retained intentionally as provenance. The present paper uses the corrected 8,264-unit surface selected by `config/chapter1_submission_current.json`.

Additional retained submission-freeze / figure-sync documents are also historical and bannered:

- `docs/chapter1_submission_freeze_20260909.md`;
- `docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md`;
- `docs/chapter1_v8_submission_figure_sync_20260915.md`;
- `docs/chapter1_v10_submission_figure_sync_20260915.md`;
- `docs/chapter1_v11_submission_figure_sync_20260915.md`;
- `docs/chapter1_v13_submission_figure_sync_20260917.md`.

## Final retired analysis runners

The following pre-corrected scientific runners are also manual-only historical replay. They no longer react to branch pushes or pull requests:

- response-geometry power / calibrated / observed runners;
- support-limited H4 prospective temporal and PolLimCrop runners, including the legacy v1 PolLimCrop preflight;
- H5c/H5d specificity and identifiability runners;
- NEE channel qualification, GloBI supported-catalog, N1/N2 discovery and source-positive-scan runners.

Their implementation and frozen artifacts remain available for provenance. Validator-only workflows may still run on relevant code changes; they test retained historical code but do not select or publish the current Chapter 1 result surface.

## Stale canonical-label cleanup

Historical execution documentation no longer advertises removed or superseded workflows as current/canonical routes. In particular:

- `analysis/v2/README.md` is explicitly a historical WHEN/WHERE surface and no longer points to a canonical workflow;
- `docs/CHAPTER1_DATABASE_RELEASE.md` treats database manifests as provenance identities, not scientific selectors;
- `config/chapter1_database_versions/README.md` distinguishes the legacy database pointer from the paper-analysis selector;
- `data/v2/README.md` states that retained Bombus/pollinator assets do not define the current H1–H4 mechanism.

The only current Chapter 1 scientific selector remains `config/chapter1_submission_current.json`.

## Retired NEE acquisition/result runners

The remaining NEE **result-generating** runners are also manual-only historical replay:

- GBIF background campaign runner;
- source-proxy-from-full-scan runner;
- island Search observation runner and sharded variant;
- island Search transport audit runner;
- sharded source-positive-scan runner.

NEE validator workflows remain branch-triggered where useful because they only validate retained historical contracts/implementations. They do not execute or promote a current Chapter 1 scientific result.
