# Floral-profile pair sharing implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Test conditional floral-profile sharing by nonshared taxa on matched island pairs and deliver reproducible results in the existing draft PR.

**Architecture:** Reuse the tested regional permutation design and add a pair-scoring layer that removes shared taxon indices before profile lookup. A separate runner selects matched pairs, executes the fixed family, and saves aggregate outputs. A renderer and report communicate the measured result without changing the earlier richness analysis.

**Tech Stack:** Existing Python, NumPy, pandas, SciPy and Matplotlib.

**Spec:** `docs/superpowers/specs/2026-11-06-floral-profile-pair-sharing-design.md`

## Global Constraints

- Work on `analysis/floral-space-20261106`, based on `506ea1f6792d9bf63cfdd3b6ec0daafddbff2f04`; no edits to existing Chapter 1 or prior analytical outputs.
- Four scopes, four exact regime keys, 50 complete names per island, 20 nonshared names per pair side, at most 2,000 pairs, 999 permutations, seed 20261106.
- All and Direct share pair IDs within each encoding and regime. Pair selection is blind to profile overlap.
- Shared names are removed before profile lookup; the full regional pool is retained for genus permutations.
- New results go to `results/floral_space_pairs_20261106/`; public outputs omit species-to-profile rows.
- Reuse the existing finite-tail and BH helpers; report BY with the same fixed 16-comparison family.
- No new dependencies. All implementation uses test-first verification of meaningful mathematical and integration behavior.
- The controller owns all git commits and publication; workers edit only their assigned files and do not dispatch subagents.

## Review Focus

- A shared taxon removed from one pair but present uniquely in another retains one global regional assignment: exact overlapping-pair fixture.
- Taxa found only on small islands remain in the permutation pool: reuse and exercise the existing regional preparation in an end-to-end fixture.
- Nondefault CLI thresholds and permutation counts propagate into outputs and figure labels: parameterized integration fixture and metadata-driven rendering.
- Empty or wholly overlapping pairs never become Jaccard=1 or fictitious p-values: rejection/empty-family fixtures.
- Input row order does not change sampled pairs or seeded inference: shuffled-input fixture and deterministic sorting.

### Task 1: Pair preparation, scoring and permutation kernel

**Files:**
- Create: `src/island_v2/floral_space_pairs.py`
- Test: `tests/test_floral_space_pairs.py`

**Interfaces:**
- Consumes `RegionPermutationDesign` from `floral_space_null`.
- Produces `prepare_pair_design(region_design, pairs: pandas.DataFrame, *, min_unique: int = 20) -> PairPermutationDesign`. `pairs` has `island_a`, `island_b`; preserve the given row order after canonicalizing endpoint order; reject repeated unordered/self/unknown pairs and insufficient support. Expose `.region`, `.pairs`, `.member_side_indices`, `.member_taxon_indices` and `.n_pairs` or equivalent documented properties.
- Produces `pair_overlap_metrics(design, assignment: numpy.ndarray) -> dict[str, numpy.ndarray]`, with `n_profiles_a`, `n_profiles_b`, `profile_intersection`, `profile_union`, and `profile_jaccard`.
- Produces `run_pair_comparison(design, *, permutations: int = 999, seed=20261106) -> tuple[dict, pandas.DataFrame, pandas.DataFrame]`. The second item has one row per draw and columns `permutation`, `median_profile_jaccard`, `median_profile_intersection`, `median_profile_union`, `median_n_profiles_a`, `median_n_profiles_b`. The third item includes pair IDs, taxon counts, observed metrics, `null_mean_jaccard`, and `null_sd_jaccard`. Summary fields include `observed_median_jaccard`, `observed_minus_null_median`, `observed_to_null_median`, `n_pairs`, `n_islands`, `n_taxa_pool`, `n_genera`, `n_permutable_genera`, `n_permutations`, `status`, and all finite-tail helper fields. Add `null_min`, `null_max`, `null_n_distinct_sampled_values` and observed/null medians for the four descriptive component diagnostics. If the null median is zero, the ratio is NaN with an explicit status/flag; keep the absolute difference.

- [x] Write tests before implementation: exact enumeration with repeated labels; shared-before-profile removal; one assignment across overlapping pairs; genus multiset preservation; invalid/empty pairs; deterministic draw streams; sufficient-support rejection.
- [x] Use this exact fixture: genus G taxa a,b,c,d carry profile IDs 0,0,0,1; genus H taxon s carries 2. Islands I={s,a,b}, J={s,c}, K={d}; eligibility=2, min_unique=1. Metrics I: n_species=3,n_profiles=2,expected=2.5; J:2,2,2; K:1,1,1. Pair I,J has exact labelled-permutation Jaccards 0 six times, 0.5 twelve times, 1 six times. K is ineligible but d remains in the pool. For coherence replace K with L={s,d}, metrics2,2,2, and use pairs(I,J),(I,L): every joint draw has median0.5, although individual pair values vary. Include a separate exclusive carrier of the shared taxon's profile to catch removal of profiles instead of taxa.
- [x] Run the new focused tests and record the expected failing output.
- [x] Implement vectorized occupancy across pair sides using the original regional profile IDs; do not sample each pair independently. Stream pair-level first and second moments rather than saving a large pair-by-draw matrix.
- [x] Run the focused tests and Ruff; provide exact commands and passing counts in the report.
- [x] Controller reviews and commits the accepted kernel before integrating the runner.

### Task 2: Matched pair selection and reproducible suite

**Files:**
- Create: `src/island_v2/floral_space_pair_analysis.py`
- Test: `tests/test_floral_space_pair_analysis.py`

**Interfaces:**
- Consumes Task 1's `prepare_pair_design` and `run_pair_comparison`.
- Produces `select_matched_pairs(all_design, direct_design, *, min_unique=20, max_pairs=2000, seed=20261106) -> tuple[pandas.DataFrame, dict]`. Use the intersection of eligible island IDs and uniform reservoir sampling over sorted qualifying unordered pairs. Check that the designs refer to the same regime. Audit counts before support filtering, after filtering, and after sampling. Use both scopes' membership; do not assume nesting.
- Produces `run_suite(results: Path, flora_path: Path, geography_path: Path, output: Path, *, permutations=999, min_species=50, min_unique=20, max_pairs=2000, seed=20261106) -> dict` and a module CLI exposing every argument.
- Fixed encoding pairs are (`all`, `direct`) and (`all_single_state`, `direct_single_state`); fixed summary order follows `SCOPES` then `REGIMES`. Pair selection and permutation seeds are those in the spec.

- [x] Write failing tests for matched support filtering, zero supported pairs, row-order invariance, seeded selection, fixed 16 family slots, same All/Direct pairs, and CLI output on small fixtures with nondefault thresholds/draw counts.
- [x] Implement preparation with existing `build_membership` and `prepare_region`; retain their saved metric checks and full pools.
- [x] Save the five tabular/manifest outputs named in the spec; deterministic gzip mtime=0. Each source file and both new modules are hashed. The manifest records software, CLI arguments, seed streams, parent commit and source result provenance.
- [x] Run focused tests plus the earlier 57 floral-space tests; validate BH/BY against SciPy on fixture outputs.
- [x] Controller reviews, runs the production suite, reconstructs inference from saved draws, and commits verified outputs.

### Task 3: Figure, scientific report and release verification

**Files:**
- Create: `scripts/render_floral_space_pairs.py`
- Create: `results/floral_space_pairs_20261106/README.md`
- Create: `results/floral_space_pairs_20261106/validation.json`
- Generate: `results/floral_space_pairs_20261106/pair_sharing_results.png` and `.svg`

**Interfaces:**
- Consumes the actual outputs of Task 2. All counts, parameters and labels come from the saved manifest and tables.

- [x] Render observed median Jaccard against the 95% null range for all 16 comparisons, with separate panels for the two encodings, consistent scales and clear sample counts.
- [x] Inspect the final raster figure for readability, clipping, correct interval meaning and appropriate scale; use no invented data or confidence intervals.
- [x] Write a concise report with the exact tested population, matched pairs, new results of either direction, relation to prior within-island richness results, limitations and replay commands.
- [x] Run the full project test suite once; distinguish existing baseline failures from regressions. Record focused tests, lint, numerical reconstruction, hashes and public release checks.
- [x] Obtain an independent whole-change review, fix material issues, publish a new commit to the existing draft PR using an expected-head guard, and verify remote/local tree equality.

## Completion record

All three tasks are complete. The measured result, full test status, independent numerical checks, visual checks and whole-change review are recorded in [the pair-sharing report](../../../results/floral_space_pairs_20261106/README.md) and [validation.json](../../../results/floral_space_pairs_20261106/validation.json). The output is delivered through the existing draft PR #321; the earlier result directory and Chapter 1 submission files remain unchanged.
