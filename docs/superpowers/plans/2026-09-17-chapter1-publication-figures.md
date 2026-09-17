# Chapter 1 Publication Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a deterministic publication-figure and source-data pipeline for the current Chapter 1 paper, rendering three main figures plus Extended Data figures/tables directly from frozen locks and pinned GitHub Actions artifacts.

**Architecture:** Separate data materialization from rendering. `chapter1_publication_data.py` validates and converts frozen artifacts/locks into tidy source-data tables; `chapter1_publication_figures.py` renders publication artwork from those tables only. A dedicated GitHub Actions workflow downloads pinned artifacts, verifies counts/digests, renders PDF/SVG/PNG output, runs tests, and uploads one publication bundle.

**Tech Stack:** Python 3.11, pandas, numpy, matplotlib, geopandas/shapely for map boundaries, pytest, Ruff, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-17-chapter1-publication-figures-design.md`

## Global Constraints

- Publication-facing artwork must not contain `v13`, workflow IDs, artifact IDs, branch names, or repository implementation jargon.
- Main figures use 180 mm double-column width, <=170 mm body height, Arial/Helvetica, 5–7 pt final text, RGB, >=0.25 pt lines.
- Save each main/Extended Data figure as PDF, SVG, and >=300 dpi PNG.
- Every plotted numeric value must come from frozen locks, pinned artifacts, or deterministic summaries of them.
- Every figure/table gets companion source data and provenance metadata.
- Do not refit scientific models in the figure workflow.
- Keep 4,453 trait-informed island inputs distinct from primary H1 complete-model support (4,334 context-specific islands) and Direct-only support (4,256).

---

### Task 1: Freeze publication-data contracts with failing tests

**Files:**
- Create: `tests/test_chapter1_publication_figures.py`

**Interfaces:**
- Consumes: none.
- Produces: behavioral contract for `materialize_publication_data(...)`, `render_publication_bundle(...)`, and publication labels.

- [ ] **Step 1: Write failing tests for source-data invariants**

Tests must assert:

```python
assert geography["island_id"].nunique() == 8265
assert trait_inputs["island_id"].nunique() == 4453
assert primary_omnibus["n_unique_islands"].sum() == 4334
assert direct_omnibus["n_unique_islands"].sum() == 4256
assert glopl_sites["site_key"].nunique() == 1248
assert island_glopl_sites["site_key"].nunique() == 197
assert glopl_tested_islands["island_id"].nunique() == 37
assert set(glopl_tested_islands["island_id"]) <= set(trait_inputs["island_id"])
```

Also assert Figure 2 source data has exactly `6 * 4 * 2 = 48` all-observed atomic rows and Figure 3 primary trait rows reproduce the frozen functional-bridge estimates.

- [ ] **Step 2: Write failing tests for publication labels**

```python
for text in all_publication_labels:
    assert "v13" not in text.lower()
    assert "workflow" not in text.lower()
    assert "artifact" not in text.lower()
```

- [ ] **Step 3: Run the focused tests on the branch**

Run in CI: `pytest -q tests/test_chapter1_publication_figures.py`

Expected RED: import/module errors because publication modules do not yet exist.

- [ ] **Step 4: Commit the RED tests**

Commit message: `test: define Chapter 1 publication figure contract`

---

### Task 2: Materialize tidy publication source data

**Files:**
- Create: `src/island_v2/chapter1_publication_data.py`
- Modify: `tests/test_chapter1_publication_figures.py`

**Interfaces:**
- Consumes directories containing the five pinned input artifacts plus the two committed result locks.
- Produces `PublicationData` paths and CSV/JSON files under `figures/chapter1/source_data/`.

- [ ] **Step 1: Implement input validation only**

Define constants for pinned artifact identities:

```python
PURPOSE_RUN = 29228212586
PURPOSE_ARTIFACT = "purpose-shortest-distance-regime-29228212586"
PURPOSE_DIGEST = "sha256:695f35b97bae07e81b05deab537dd73fa687b2d99e9efdb1cb3babd2fa12dfb6"
ATOMIC_RUN = 34961775336
ATOMIC_ARTIFACT = "chapter1-all-data-probability-34961775336"
ATOMIC_DIGEST = "sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b"
GLOPL_OVERLAP_RUN = 35085378171
GLOPL_OVERLAP_ARTIFACT = "chapter1-h5-glopl-island-overlap-preflight-35085378171"
GLOPL_OVERLAP_DIGEST = "sha256:5ab7aecfbeb4fac2db26bcdc13edd4b47e953b1d1e6f35255330a0b2f70fe041"
GLOPL_GLOBAL_RUN = 35090599662
GLOPL_GLOBAL_ARTIFACT = "chapter1-h5-glopl-global-distance-35090599662"
GLOPL_GLOBAL_DIGEST = "sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623"
GLOPL_SHAPE_RUN = 35091385147
GLOPL_SHAPE_ARTIFACT = "chapter1-h5-glopl-global-shape-audit-35091385147"
GLOPL_SHAPE_DIGEST = "sha256:ee7700a03dcbe9b0d6184195c482e3ec7915a5e7a9e8f1c6a85d0e4fe2de8e44"
FUNCTIONAL_RUN = 35141624253
FUNCTIONAL_ARTIFACT = "chapter1-v13-functional-bridge-35141624253"
FUNCTIONAL_DIGEST = "sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e"
```

Validate required files and committed lock values before reading plot data.

- [ ] **Step 2: Build Figure 1 source data**

From `purpose_shortest_island_data.csv`, output all 8,265 island coordinates with `trait_informed = n_trait_species > 0`. From `GLOPL_METADATA_ISLAND_MATCH.csv.gz`, output all 1,248 unique GloPL coordinates and exact frozen-island matches. Assert 197 island sites and 37 unique matched islands.

- [ ] **Step 3: Build Figure 2 source data**

Read `all/beta_binomial_within_slopes.csv`, `direct/beta_binomial_within_slopes.csv`, and corresponding omnibus tables. Keep `stratum == "all_observed"`. Add evidence scope, 95% CI as `estimate ± 1.96 * cluster_robust_se`, family labels, and model-specific omnibus support. Output 48 atomic rows.

Build deterministic descriptive summaries from the committed unified lock for classic-orientation means and two family means; label these rows `descriptive` rather than inferential.

- [ ] **Step 4: Build Figure 3 source data**

Read the committed unified and functional locks plus GloPL result/shape artifacts. Create one tidy table with rows for global distance, supplemental-only, no-zero-constant, mainland-to-offshore step, within-offshore gradient, offshore-only gradient, four primary trait-state effects, and their available sensitivities. Compute 95% CI only when a frozen SE is available. Preserve `not_evaluable` rows explicitly.

- [ ] **Step 5: Build Extended Data tables**

Write `extended_data_table1_data_layers.csv` through `extended_data_table6_glopl_tested_islands.csv` as defined in the spec.

- [ ] **Step 6: Write provenance manifest**

Record source paths, expected counts, lock SHAs supplied by the workflow, pinned artifact metadata, and SHA-256 digests of generated source-data files.

- [ ] **Step 7: Run focused tests**

`pytest -q tests/test_chapter1_publication_figures.py -k data`

Expected GREEN.

- [ ] **Step 8: Commit**

Commit message: `feat: materialize Chapter 1 publication source data`

---

### Task 3: Render the three main figures

**Files:**
- Create: `src/island_v2/chapter1_publication_figures.py`
- Modify: `tests/test_chapter1_publication_figures.py`

**Interfaces:**
- Consumes source-data CSVs from Task 2 only.
- Produces main figures in `figures/chapter1/main/`.

- [ ] **Step 1: Add failing rendering tests**

Assert the renderer creates PDF/SVG/PNG for all three main figures; inspect SVG text for forbidden internal terminology and required biological labels; assert physical figure width is 180 mm within tolerance.

- [ ] **Step 2: Implement shared publication style**

Use Matplotlib rcParams with Arial/Helvetica fallback, 5–7 pt text at final size, >=0.25 pt lines, accessible Okabe-Ito-derived palette, lower-case bold panel labels, and no in-artwork figure title.

- [ ] **Step 3: Implement Figure 1**

Panel a real global map: 8,265 light-grey geography points, 4,453 trait-informed inputs, 1,248 GloPL sites, 37 GloPL-tested frozen islands. Panel b compact biological evidence design. Keep map dominant and remove repository terminology.

- [ ] **Step 4: Implement Figure 2**

Panels a/b coefficient forests for primary and Direct-only scopes with identical x limits. Panels c/d descriptive classic orientation and two-family summaries. Display context model n and joint P values once per context; no significance stars.

- [ ] **Step 5: Implement Figure 3**

Panel a global GloPL estimate/CI; panel b sensitivity + offshore diagnostic forest; panel c functional trait forest with support columns; panel d robustness/negative-result matrix. Distinguish primary, sensitivity, post-hoc diagnostic, and non-supported evidence using shape/line style and labels.

- [ ] **Step 6: Verify focused rendering tests GREEN**

`pytest -q tests/test_chapter1_publication_figures.py -k render`

- [ ] **Step 7: Commit**

Commit message: `feat: render publication-ready Chapter 1 main figures`

---

### Task 4: Render Extended Data figures and tables

**Files:**
- Modify: `src/island_v2/chapter1_publication_figures.py`
- Modify: `tests/test_chapter1_publication_figures.py`

- [ ] **Step 1: Add failing tests for four Extended Data figures and six tables**
- [ ] **Step 2: Implement atomic support/sample-size audit figure**
- [ ] **Step 3: Implement geographic sampling/GloPL overlap audit figure**
- [ ] **Step 4: Implement detailed GloPL sensitivity figure**
- [ ] **Step 5: Implement functional-bridge sensitivity matrix**
- [ ] **Step 6: Verify all Extended Data outputs and source tables**

Run: `pytest -q tests/test_chapter1_publication_figures.py`

- [ ] **Step 7: Commit**

Commit message: `feat: add Chapter 1 Extended Data display items`

---

### Task 5: Add deterministic GitHub Actions rendering workflow

**Files:**
- Create: `.github/workflows/render-chapter1-publication-figures.yml`

**Interfaces:**
- Downloads the six pinned artifacts listed above.
- Calls `python -m island_v2.chapter1_publication_data` then `python -m island_v2.chapter1_publication_figures`.
- Uploads `chapter1-publication-figures-${{ github.run_id }}`.

- [ ] **Step 1: Add workflow trigger for the feature branch and pull requests to main**
- [ ] **Step 2: Download pinned artifacts with `gh run download`**
- [ ] **Step 3: Verify required artifact files and manifest digests**
- [ ] **Step 4: Run focused tests before rendering**
- [ ] **Step 5: Render bundle and run structural verification**
- [ ] **Step 6: Upload complete `figures/chapter1/` directory**
- [ ] **Step 7: Commit**

Commit message: `ci: render Chapter 1 publication figures from frozen inputs`

---

### Task 6: Synchronize manuscript legends and submission documentation

**Files:**
- Modify: `docs/chapter1_v13_submission_figure_sync_20260917.md`
- Modify: `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
- Modify: `README.md` only if a short publication-output pointer is needed.

- [ ] **Step 1: Rewrite publication-facing figure legends without internal version labels**
- [ ] **Step 2: Describe uncertainty, n, and evidence role panel-by-panel**
- [ ] **Step 3: Keep internal provenance in the reproducibility section only**
- [ ] **Step 4: Run submission-surface tests plus publication tests**

Run:
`pytest -q tests/test_chapter1_v13_submission_surface.py tests/test_chapter1_publication_figures.py`

- [ ] **Step 5: Commit**

Commit message: `docs: sync Chapter 1 publication figures and legends`

---

### Task 7: Final branch verification and PR

- [ ] **Step 1: Run the publication workflow on branch HEAD**
- [ ] **Step 2: Verify workflow artifact contains all expected PDF/SVG/PNG/CSV/JSON outputs**
- [ ] **Step 3: Run Ruff for new modules/tests**
- [ ] **Step 4: Verify no generated publication label contains internal `v13`/workflow/artifact jargon**
- [ ] **Step 5: Compare branch against main and ensure scientific locks are unchanged**
- [ ] **Step 6: Open a PR to main with exact CI and artifact evidence; do not merge without a separate integration decision**
