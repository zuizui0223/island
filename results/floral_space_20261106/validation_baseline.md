# Validation of pre-existing test failures

Date: 2026-11-06 (JST).

**All five failures from the full project test run reproduce on the unchanged baseline commit.** No existing source, tests, or submission documents were modified during this verification. The full suite remains red; these failures are not reported as passing or silently excluded.

## Verification scope

- Baseline commit: `3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5`.
- Isolation: a detached Git worktree at `/workspace/scratch/625333e4fa64/island-baseline`, containing none of the new floral-space modules, tests, or results.
- Runtime: the same Python environment as the full project run; the project's pytest configuration selects its own `src` directory.
- Full project run: **1,460 passed, 1 skipped, 5 failed**, with 8 warnings, in 49.73 seconds.
- Independent baseline run: **the same 5 tests failed**, with 1 warning, in 0.57 seconds. Only the five named tests were run in the baseline worktree.
- `README.md`, `docs/PAPER_PIPELINE.md`, `src/island_v2/v1_category_search.py`, and all three affected test files were byte-identical between the working repository and the baseline worktree when checked. The baseline worktree had no tracked or untracked changes after the test run.

## Reproduced failures and immediate causes

| Test | Failure reproduced on baseline | Cause and classification |
| --- | --- | --- |
| `tests/test_chapter1_v13_legacy_boundary.py::test_active_handoff_points_to_archive_boundary` | Missing `legacy/chapter1-pre-v13/` in `docs/PAPER_PIPELINE.md`. | Existing submission-document contract mismatch. The current pipeline describes the corrected-geography submission surface. |
| `tests/test_chapter1_v13_submission_surface.py::test_readme_is_v13_only_and_points_pre_v13_to_legacy` | Missing `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md` in `README.md`. | Existing v13-specific test expectation differs from the current corrected submission handoff. |
| `tests/test_chapter1_v13_submission_surface.py::test_pipeline_is_v13_only_and_points_pre_v13_to_legacy` | The same v13 manuscript path is missing from `docs/PAPER_PIPELINE.md`. | Existing v13-specific test expectation differs from the current corrected submission pipeline. |
| `tests/test_chapter1_v13_submission_surface.py::test_current_handoff_is_global_only` | Missing exact phrase `floral accessibility` in `README.md`. | Existing text-contract mismatch; the failed assertion does not involve floral-space calculations. |
| `tests/test_v1_category_search.py::test_disabled_web_provider_does_not_parse_sitemap_snapshot` | `ImportError`: HTTPX detects a SOCKS proxy but the optional `socksio` package is absent. | Existing environment-dependent initialization failure. `collect_free_sources()` constructs `httpx.Client` even when every provider is disabled. The failure occurs at client construction, before any provider request or sitemap parsing. |

The first four failures establish disagreement between existing tests and existing handoff text; this check does not decide which submission contract should replace the other. The fifth reproduces in the supplied runtime and is not a regression introduced by the new analysis.

## Reproduction command

Run from a clean checkout of the baseline commit:

```bash
python -m pytest -q \
  tests/test_chapter1_v13_legacy_boundary.py::test_active_handoff_points_to_archive_boundary \
  tests/test_chapter1_v13_submission_surface.py::test_readme_is_v13_only_and_points_pre_v13_to_legacy \
  tests/test_chapter1_v13_submission_surface.py::test_pipeline_is_v13_only_and_points_pre_v13_to_legacy \
  tests/test_chapter1_v13_submission_surface.py::test_current_handoff_is_global_only \
  tests/test_v1_category_search.py::test_disabled_web_provider_does_not_parse_sitemap_snapshot
```

Captured local evidence:

- `work/floral_space_inputs/full_pytest.log`; SHA-256 `0e6589dc9fe80f5fe93b6dbe11660cb5c0d295820bf0409c5152e1be20b2feb2`.
- `work/floral_space_inputs/baseline_selected_pytest.log`; SHA-256 `901b24a408dd8a70d3d973af18c9d9e3451504e95ca2ab11527e0d2977d25976`.

These unrelated existing failures were documented without changing the submission handoff or network acquisition code.
