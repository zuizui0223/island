# Chapter 1 v13 Unified Island Syndrome Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze and promote a reproducible v13 Chapter 1 paper in which geographic isolation is associated with a global floral/reproductive island syndrome, increasing experimental pollen limitation, two partially independent plant response pathways, a post-hoc functional bridge, and flora-layer-specific taxonomic realization.

**Architecture:** Preserve v11/v12 as immutable provenance. Add one new post-hoc functional-bridge analysis that consumes the already frozen Route A/B exact-species GloPL overlap artifacts, verifies their pinned provenance, and estimates trait-state differences in current pollen limitation with the same publication weighting, measurement fixed effects, context adjustment, and publication-cluster covariance. Then assemble a v13 result lock from existing canonical locks plus the newly verified artifact, generate the v13 manuscript/figure map, validate the claim ceiling, and only then repoint submission-facing documentation.

**Tech Stack:** Python 3, pandas, numpy, scipy/statistical utilities already in `island_v2`, pytest, Ruff, GitHub Actions, Markdown/JSON/YAML.

**Spec:** `docs/superpowers/specs/2026-09-17-chapter1-v13-unified-island-syndrome-design.md`

## Global Constraints

- Never edit or delete v11/v12 manuscript or result-lock files.
- Functional-bridge evidence is labelled `posthoc_functional_triangulation`, never confirmatory.
- Reuse the exact frozen Route A/B trait recodes and exact-species overlaps; no synonym/genus fallback and no recoding after outcome inspection.
- Publication total analysis weight remains exactly 1.0 within each trait analysis.
- Primary covariance is publication-cluster robust.
- Primary model: `PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects`.
- Sensitivities: supplemental-only, no-zero-constant, within-publication, within-publication×site.
- `selfing_mating_system` and `shallow_open_tube` remain non-evaluable unless the existing frozen support artifacts already admit them; no threshold relaxation.
- GloBI remains supplementary and cannot carry a causal mechanism claim.
- No statement may claim historical causal mediation from pollen limitation to trait evolution.
- Do not merge automatically; finish with a PR against `main`.

---

### Task 1: Freeze the v13 functional-bridge contract and write RED tests

**Files:**
- Create: `config/chapter1_v13_functional_bridge_v1.yml`
- Create: `tests/test_chapter1_v13_functional_bridge.py`
- Later create: `src/island_v2/chapter1_v13_functional_bridge.py`

**Interfaces:**
- Consumes two canonical parent artifacts containing `out/MATCHED_EFFECT_ROWS.csv.gz` from Route A and Route B.
- Produces `load_parent_rows(path) -> pandas.DataFrame`, `aggregate_species_measurement_cells(rows) -> pandas.DataFrame`, `fit_global_trait_level(cells) -> dict`, `fit_within_group_trait_level(cells, group_columns) -> dict`, and `run_functional_bridge(...) -> dict`.

- [ ] **Step 1: Write the frozen config**

Pin:

```yaml
contract: chapter1_v13_functional_bridge_v1
inferential_role: posthoc_functional_triangulation
parents:
  reproductive_assurance:
    workflow_run_id: 35093274622
    artifact_id: 10444749163
    required_contract: chapter1_h5_glopl_reproductive_assurance_moderation_v1
  floral_architecture:
    workflow_run_id: 35094588521
    artifact_id: 10445189257
    artifact_digest: sha256:ef81551d41cc747b07e927b9a5c4f13eaf0fc013a5cf8ac3fb1a5491aadec8c4
    required_contract: chapter1_h5_glopl_floral_architecture_moderation_v1
primary_model: PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects
publication_total_weight: 1.0
cluster: study_key
sensitivities: [supplemental_only, no_zero_constant, within_publication, within_publication_site]
```

The config must list the four evaluable frozen traits explicitly: `self_compatibility`, `autonomous_selfing`, `generalized_form`, `actinomorphic_symmetry`, while retaining the two support-failed traits as non-evaluable provenance.

- [ ] **Step 2: Write failing unit tests for publication weighting and aggregation**

```python
def test_aggregate_cells_gives_each_publication_total_weight_one():
    cells = aggregate_species_measurement_cells(_toy_rows())
    totals = cells.groupby("study_key")["analysis_weight"].sum()
    assert np.allclose(totals.to_numpy(), 1.0)
```

- [ ] **Step 3: Write failing test for global trait coefficient direction**

Use synthetic data with two contexts and measurement levels where protected/accessible state lowers PL by 0.5 after adjustment:

```python
def test_global_trait_level_recovers_negative_trait_effect():
    result = fit_global_trait_level(_synthetic_cells(effect=-0.5))
    assert result["evaluable"] is True
    assert result["trait_state_estimate"] < 0
    assert result["trait_state_one_sided_negative_p"] < 0.05
```

- [ ] **Step 4: Write failing within-publication and within-site tests**

```python
def test_within_publication_requires_both_trait_states():
    result = fit_within_group_trait_level(_one_state_per_publication(), ["study_key"])
    assert result["evaluable"] is False
    assert result["reason"] == "no_groups_with_both_trait_states"


def test_within_site_recovers_paired_negative_difference():
    result = fit_within_group_trait_level(
        _paired_site_cells(effect=-0.4), ["study_key", "site_key"]
    )
    assert result["evaluable"] is True
    assert result["trait_state_estimate"] < 0
```

- [ ] **Step 5: Run RED test**

Run:

```bash
pytest tests/test_chapter1_v13_functional_bridge.py -q
```

Expected: FAIL because `chapter1_v13_functional_bridge` is not implemented.

- [ ] **Step 6: Commit RED contract/tests**

```bash
git add config/chapter1_v13_functional_bridge_v1.yml tests/test_chapter1_v13_functional_bridge.py
git commit -m "test: freeze v13 functional bridge contract"
```

---

### Task 2: Implement the functional-bridge estimator and GREEN unit tests

**Files:**
- Create: `src/island_v2/chapter1_v13_functional_bridge.py`
- Modify only if necessary: `tests/test_chapter1_v13_functional_bridge.py`

**Interfaces:**
- Reads canonical `MATCHED_EFFECT_ROWS.csv.gz` files only after parent provenance is verified by the workflow.
- Writes `functional_bridge_trait_results.csv`, `functional_bridge_manifest.json`, and `functional_bridge_summary.md`.

- [ ] **Step 1: Implement cell aggregation exactly as Route A/B**

```python
def aggregate_species_measurement_cells(rows: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "study_key", "site_key", "species_key", "analysis_regime", "z_distance",
        "trait", "trait_state", *MEASUREMENT_COLUMNS,
    ]
    out = rows.groupby(group_cols, as_index=False, dropna=False).agg(
        PL_Effect_Size=("PL_Effect_Size", "mean"),
        n_effect_rows=("PL_Effect_Size", "size"),
    )
    out["analysis_weight"] = 1.0 / out.groupby("study_key")["site_key"].transform("size")
    return out
```

- [ ] **Step 2: Implement global adjusted trait-level model**

Build the design in this exact order:

```python
[intercept, context_dummies..., z_distance, trait_state, measurement_dummies...]
```

Use the existing publication-cluster robust WLS semantics from `chapter1_h5_glopl_global_distance._clustered_wls` and return two-sided plus one-sided-negative p-values for `trait_state`.

- [ ] **Step 3: Implement within-group fixed-effect contrasts**

Restrict to publications or publication×site groups containing both trait states. Weighted-demean `PL_Effect_Size` and all non-intercept design columns within the selected fixed-effect group, drop zero-variance columns, and fit the remaining design with publication clustering. Return the number of qualifying fixed-effect groups.

- [ ] **Step 4: Implement frozen sensitivities**

For each admitted trait compute:

```text
primary
supplemental_only: PL_Effect_Size_Type2 == "Sup"
no_zero_constant: Constant_added_bool == False
within_publication
within_publication_site
```

No sensitivity may change trait coding, matching, support thresholds, or weighting.

- [ ] **Step 5: Run unit tests**

```bash
pytest tests/test_chapter1_v13_functional_bridge.py -q
ruff check src/island_v2/chapter1_v13_functional_bridge.py tests/test_chapter1_v13_functional_bridge.py
```

Expected: all pass.

- [ ] **Step 6: Commit implementation**

```bash
git add src/island_v2/chapter1_v13_functional_bridge.py tests/test_chapter1_v13_functional_bridge.py
git commit -m "feat: add v13 GloPL functional triangulation"
```

---

### Task 3: Add reproducible GitHub Actions execution and lock the observed result

**Files:**
- Create: `.github/workflows/run-chapter1-v13-functional-bridge.yml`
- Create after successful run: `config/chapter1_v13_functional_bridge_result_lock.json`
- Create: `docs/chapter1_v13_functional_bridge_result_20260917.md`

**Interfaces:**
- Workflow downloads exact historical artifacts `10444749163` and `10445189257`, verifies the configured artifact provenance/digest where GitHub exposes it, runs tests before analysis, runs the estimator, and uploads one v13 artifact.
- Result lock records workflow run ID, artifact ID, artifact digest, branch head SHA, and every reported estimate/p-value.

- [ ] **Step 1: Write workflow**

Use `gh api` with `${{ github.token }}` to download each historical artifact ZIP by numeric ID, unzip into separate parent directories, verify required `out/MATCHED_EFFECT_ROWS.csv.gz`, run the v13 module, then upload its output.

- [ ] **Step 2: Push workflow and observe first run**

The workflow is branch-scoped to `ch1-v13-unified-island-syndrome` and triggered on changes to the v13 config/module/tests/workflow.

- [ ] **Step 3: Diagnose any failure using job logs; do not loosen the scientific contract**

Allowed fixes are implementation/provenance bugs only. Do not change trait recodes, source artifacts, or statistical model because of observed outcomes.

- [ ] **Step 4: On GREEN, freeze reproduced values**

The lock must distinguish:

```json
{
  "inferential_role": "posthoc_functional_triangulation",
  "historical_trait_evolution_causally_identified": false,
  "primary": {},
  "sensitivities": {},
  "within_group": {}
}
```

Use the reproduced artifact values even if they differ from preliminary chat values.

- [ ] **Step 5: Commit result lock and result note**

```bash
git add config/chapter1_v13_functional_bridge_result_lock.json docs/chapter1_v13_functional_bridge_result_20260917.md
git commit -m "results: lock v13 functional bridge"
```

---

### Task 4: Assemble and validate the v13 paper-level result lock

**Files:**
- Create: `config/chapter1_v13_unified_island_syndrome_result_lock.json`
- Create: `src/island_v2/chapter1_v13_submission_lock.py`
- Create: `tests/test_chapter1_v13_submission_lock.py`
- Create: `docs/chapter1_unified_hypothesis_20260917.md`

**Interfaces:**
- Consumes existing canonical all-data/v12/GloPL/H3 locks plus the new functional-bridge lock.
- Produces one fail-closed paper-level summary with H1-H5 evidence classes and prohibited claims.

- [ ] **Step 1: Write RED validator tests**

Tests must fail if:
- v11/v12 parent contract names change;
- H4 functional bridge is labelled confirmatory;
- a causal historical-selection flag is true;
- GloBI is promoted to mechanism;
- the new functional-bridge artifact provenance is absent;
- H5 taxonomic-realization layers are collapsed into one causal mediation claim.

- [ ] **Step 2: Implement validator**

```python
def validate_v13_lock(lock: dict) -> dict:
    assert lock["architecture"]["H1"] == "global_island_syndrome"
    assert lock["H4"]["inferential_role"] == "posthoc_functional_triangulation"
    assert lock["claim_ceiling"]["historical_pollen_limitation_selected_traits"] is False
    return {"verified": True, ...}
```

- [ ] **Step 3: Build the lock from canonical parent locks**

H1/H3 must point to the all-data/two-panel canonical locks; H2 to the full-global GloPL lock; H4 to the newly generated functional-bridge lock; H5 to H3A/H3B and native-status boundary locks.

- [ ] **Step 4: Run validator tests**

```bash
pytest tests/test_chapter1_v13_submission_lock.py -q
ruff check src/island_v2/chapter1_v13_submission_lock.py tests/test_chapter1_v13_submission_lock.py
```

- [ ] **Step 5: Commit lock, validator, and unified hypothesis note**

```bash
git add config/chapter1_v13_unified_island_syndrome_result_lock.json src/island_v2/chapter1_v13_submission_lock.py tests/test_chapter1_v13_submission_lock.py docs/chapter1_unified_hypothesis_20260917.md
git commit -m "results: freeze v13 unified island syndrome"
```

---

### Task 5: Write the v13 manuscript and figure contract

**Files:**
- Create: `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
- Create: `docs/chapter1_v13_submission_figure_sync_20260917.md`

**Interfaces:**
- Manuscript may only state values present in the v13 paper lock or explicitly labelled parent robustness notes.
- Figure contract maps the four figures in the approved design to exact result-lock fields/artifacts.

- [ ] **Step 1: Rewrite Abstract/Introduction around the v13 H1-H5 architecture**

The opening claim is global recurrence under a shared pollination-service constraint, not rejection of a universal syndrome.

- [ ] **Step 2: Rewrite Methods in the new inferential order**

Order: global plant response → GloPL → dual pathways → post-hoc functional triangulation → taxonomic realization → robustness boundaries.

- [ ] **Step 3: Rewrite Results using only locked values**

Explicitly label H4 as post-hoc. Report failed Route A/B slope-moderation tests as a distinct negative result, not as contradiction of the trait-level functional association.

- [ ] **Step 4: Rewrite Discussion and claim ceiling**

State triangulation rather than historical mediation. GloBI and named pollinator identities remain supplementary/non-promoted.

- [ ] **Step 5: Add four-figure sync document**

Map each panel to the exact lock/source artifact and mark confirmatory/post-hoc/descriptive evidence visibly.

- [ ] **Step 6: Validate manuscript claims against lock**

Add a lightweight test or script assertion that every headline numeric value appears in the v13 lock and that prohibited phrases such as `caused the evolution`, `pollinator abundance declined globally`, and `GloBI proves` do not occur in the manuscript.

- [ ] **Step 7: Commit manuscript and figure contract**

```bash
git add docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md docs/chapter1_v13_submission_figure_sync_20260917.md
git commit -m "docs: write v13 global pollination constraint manuscript"
```

---

### Task 6: Promote the v13 submission surface and verify provenance preservation

**Files:**
- Modify: `README.md`
- Modify: `docs/PAPER_PIPELINE.md`
- Create: `.github/workflows/audit-chapter1-v13-submission-surface.yml`
- Create: `tests/test_chapter1_v13_submission_surface.py`

**Interfaces:**
- Promotion occurs only after Tasks 1-5 are GREEN.
- v11/v12 files must remain byte-identical to `main`.

- [ ] **Step 1: Write RED surface test**

Assert README/PAPER_PIPELINE list the v13 result lock, manuscript, and figure contract first, while retaining v11/v12 as historical provenance.

- [ ] **Step 2: Update README and PAPER_PIPELINE**

Do not delete old links; relabel them historical.

- [ ] **Step 3: Add provenance guard**

Compare the known immutable v11/v12 files on the v13 branch against `main` via hashes or git diff and fail if any changed.

- [ ] **Step 4: Run full v13 verification**

```bash
pytest tests/test_chapter1_v13_functional_bridge.py tests/test_chapter1_v13_submission_lock.py tests/test_chapter1_v13_submission_surface.py -q
ruff check src/island_v2/chapter1_v13_functional_bridge.py src/island_v2/chapter1_v13_submission_lock.py tests/test_chapter1_v13_*.py
```

Run the GitHub Actions submission-surface audit and require GREEN.

- [ ] **Step 5: Open PR**

Open a PR from `ch1-v13-unified-island-syndrome` to `main`, include exact workflow/artifact/digest provenance, and leave it unmerged for explicit integration review.
