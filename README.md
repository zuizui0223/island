# Island — Chapter 1 v14 candidate reanalysis

> **PR #237 reorganizes the publication logic without rewriting v13 provenance.**
> The v14 analysis order is now:
>
> 1. **H1 — global floral/reproductive island syndrome:** reproductive assurance + plain colour + accessible/generalized floral structure.
> 2. **H2 — floral-shift decomposition:** selfing-syndrome/reproductive-assurance route versus selfing-adjusted pollinator-facing colour/accessibility route.
> 3. **H3 — independent global pollen limitation:** the frozen GloPL isolation analysis, formerly v13 H2.
> 4. **H4 — functional bridge:** exact-species trait × GloPL triangulation, unchanged in role.
>
> H1 is refit as a seven-response beta-binomial multivariate analysis. H2 explicitly
> conditions colour and accessibility responses on selfing_core. The detailed design
> is in [docs/chapter1_v14_hypothesis_reorder_20260918.md](docs/chapter1_v14_hypothesis_reorder_20260918.md)
> and the frozen-input preflight results are in
> [docs/chapter1_v14_preflight_results_20260918.md](docs/chapter1_v14_preflight_results_20260918.md).
>
> **v13 remains immutable historical/submission provenance until the v14 CI artifact is
> reproduced and locked.**

---

# Island — Chapter 1 v13 paper repository

This repository's **current publication-facing surface is Chapter 1 v13**. The active paper remains a global-only **H1–H4** synthesis. Earlier Chapter 1 hypothesis ladders, defended submission states and submission-surface audits are historical provenance, not the current scientific spine.

## 1. Current scientific claim

The v13 paper now tests one upgraded integrated idea:

> **Geographic isolation is associated globally with stronger experimental pollen limitation and with a recurrent functional floral/reproductive island-syndrome core. Plants express that core through partially separable reproductive-assurance and floral-accessibility responses, while pollinator-facing floral display can be reorganized in context-dependent ways.**

The four geographic strata remain replication strata for recurrence of the functional core. The raw flower-colour and colour–architecture analyses add a reproduced extended H3 layer showing that recurrence at the functional level does not require identical detailed phenotype vectors or one universal pollinator mechanism.

## 2. Current H1–H4 architecture

1. **H1 — recurrent global functional island syndrome.** Six identically oriented floral/reproductive responses show supported multivariate isolation responses in all four predeclared geographic replication strata, with a positive descriptive classic-island direction in both evidence scopes.
2. **H2 — global pollination constraint.** Experimental pollen limitation increases with geographic separation in the full-global GloPL analysis (`beta=+0.07937`; one-sided `p=0.01772`).
3. **H3 — partially separable plant-response modules.** The frozen recurrent core contains reproductive assurance and floral accessibility/generalization. A reproduced raw-colour extension adds a third pollinator-facing display module: colour composition and, in some contexts, colour–architecture coupling change with isolation after adjustment for measured reproductive assurance.
4. **H4 — functional bridge.** Exact-species GloPL analyses provide explicitly post-hoc functional triangulation. Autonomous selfing is the strongest reproducible association with lower current pollen limitation; architecture evidence is weaker but directionally concordant.

The raw-colour extension does **not** retroactively become a seventh confirmatory H1 atomic trait. It upgrades the biological interpretation of H3 while preserving the frozen H1–H4 result lock.

## 3. Canonical v13 paper surface

Read these first, in this order:

1. [`docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`](docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md) — current manuscript;
2. [`config/chapter1_v13_unified_island_syndrome_result_lock.json`](config/chapter1_v13_unified_island_syndrome_result_lock.json) — frozen H1–H4 results and claim ceiling;
3. [`docs/chapter1_v13_submission_figure_sync_20260917.md`](docs/chapter1_v13_submission_figure_sync_20260917.md) — three-main-figure evidence-role contract;
4. [`config/chapter1_v13_functional_bridge_result_lock.json`](config/chapter1_v13_functional_bridge_result_lock.json) — post-hoc exact-species functional triangulation;
5. [`docs/chapter1_unified_hypothesis_20260917.md`](docs/chapter1_unified_hypothesis_20260917.md) — upgraded biological hypothesis and evidence map;
6. [`docs/chapter1_v13_raw_colour_audit_20260917.md`](docs/chapter1_v13_raw_colour_audit_20260917.md) — reproduced raw flower-colour extension;
7. [`docs/chapter1_v13_raw_colour_coupling_audit_20260917.md`](docs/chapter1_v13_raw_colour_coupling_audit_20260917.md) — reproduced colour-conditioned architecture extension;
8. [`docs/chapter1_submission_freeze_v13_20260917.md`](docs/chapter1_submission_freeze_v13_20260917.md) — current submission freeze.

The compact reproducibility map is [`docs/PAPER_PIPELINE.md`](docs/PAPER_PIPELINE.md).

## 4. Frozen database and parent evidence

The paper uses a fixed universe of **8,265 islands** and **106,295 analysis-applicable plant species**. The three raw trait-coverage axes are flower colour, floral structural complexity and reproductive assurance.

Final trait snapshot:

- source run: **34191508045**;
- artifact: `source-scale-batch-integration-34191508045`;
- resolved cells: **222,688 / 318,885 = 69.83%**;
- reproductive-assurance cells: **48,497 / 106,295 = 45.63%**.

Database construction is documented in [`docs/DATABASE_BUILD.md`](docs/DATABASE_BUILD.md), and the versioned release route is in [`DATABASE_RELEASE.md`](DATABASE_RELEASE.md).

The frozen parent **result locks required by v13 remain active** so the submission contract can validate their exact contracts. Retired pre-v13 publication and mechanism implementations are no longer importable active code or workflows; their historical provenance is indexed under [`legacy/chapter1-pre-v13/`](legacy/chapter1-pre-v13/) and is recoverable from the pre-cleanup Git commit recorded there.

## 5. Evidence roles and claim ceiling

v13 now combines four evidence layers:

- **plant-side functional recurrence:** isolation-associated floral/reproductive response vectors across four geographic replication strata;
- **independent ecological pressure:** global pollen-supplementation experiments measuring pollen limitation;
- **pollinator-facing display extension:** raw flower-colour composition and colour-conditioned architecture, analyzed without one-colour-one-pollinator scoring;
- **functional compatibility:** frozen trait states compared with current experimental pollen limitation.

The synthesis is **triangulation, not historical causal mediation**. It does not establish that pollen limitation caused the observed trait evolution, that pollinator abundance globally declines with isolation, that a named pollinator was lost, or that the observed assemblage pattern arose through within-lineage evolution. Raw colour does not uniquely identify realized pollinator identity or directly measure attraction investment.

## 6. Active validation

The current submission contract is enforced by:

- [`tests/test_chapter1_v13_submission_lock.py`](tests/test_chapter1_v13_submission_lock.py);
- [`tests/test_chapter1_v13_submission_surface.py`](tests/test_chapter1_v13_submission_surface.py);
- [`tests/test_chapter1_v13_functional_bridge.py`](tests/test_chapter1_v13_functional_bridge.py);
- [`tests/test_chapter1_v13_raw_colour_audit.py`](tests/test_chapter1_v13_raw_colour_audit.py);
- [`tests/test_chapter1_v13_colour_architecture_audit.py`](tests/test_chapter1_v13_colour_architecture_audit.py);
- [`tests/test_chapter1_v13_raw_colour_coupling_audit.py`](tests/test_chapter1_v13_raw_colour_coupling_audit.py);
- [`tests/test_chapter1_v13_legacy_boundary.py`](tests/test_chapter1_v13_legacy_boundary.py);
- [`tests/test_v13_active_surface_legacy_boundary.py`](tests/test_v13_active_surface_legacy_boundary.py);
- [`.github/workflows/audit-chapter1-v13-submission.yml`](.github/workflows/audit-chapter1-v13-submission.yml);
- [`.github/workflows/run-chapter1-v13-raw-colour-audit.yml`](.github/workflows/run-chapter1-v13-raw-colour-audit.yml).

The audits run for changes targeting `main` and check the repository suite, v13 linting, the v13-only paper surface, the legacy boundary and preservation of frozen historical provenance.

## 7. Chapter 1 / Chapter 2 handoff

Chapter 1 now establishes **recurrent but non-uniform response**: isolation is associated with a global functional core and stronger pollen limitation, yet the detailed pollinator-facing display response is context dependent.

Chapter 2 (`izu-core`) asks the next-level question: why can a broad ecological constraint recur while detailed response trajectories and determinant rankings differ among systems? Its current transportability framework treats receiving functional state, interaction regime and nonlinear response geometry as reasons that one determinant ranking need not transport unchanged across contexts.

The programme-level handoff is:

`Ch1: what recurs? -> functional core + context-dependent display`

`Ch2: why need responses not be identical? -> state/context/response geometry -> different response branches`

This is a conceptual bridge, not cross-repository causal validation.

## 8. Repository map

```text
README.md
  <- current v13 entry point
docs/PAPER_PIPELINE.md
  <- compact v13 reproducibility map
docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md
config/chapter1_v13_unified_island_syndrome_result_lock.json
config/chapter1_v13_functional_bridge_result_lock.json
docs/chapter1_v13_submission_figure_sync_20260917.md
docs/chapter1_unified_hypothesis_20260917.md
docs/chapter1_v13_raw_colour_audit_20260917.md
docs/chapter1_v13_raw_colour_coupling_audit_20260917.md
docs/chapter1_submission_freeze_v13_20260917.md
docs/DATABASE_BUILD.md
DATABASE_RELEASE.md
src/island_v2/
data/v2/
legacy/chapter1-pre-v13/
legacy/v1/
```

## 9. Legacy boundary

The pre-v13 publication and mechanism branches are archived/indexed under [`legacy/chapter1-pre-v13/`](legacy/chapter1-pre-v13/). Representative historical files are physically retained there with their former relative paths; the complete retired implementation is exactly recoverable from the pre-cleanup commit recorded in [`legacy/chapter1-pre-v13/ARCHIVE_MANIFEST.md`](legacy/chapter1-pre-v13/ARCHIVE_MANIFEST.md).

The active root intentionally retains only the v13 publication surface, shared data/trait/geography infrastructure, and frozen parent result locks required by the v13 contract. Retired pre-v13 implementations must not re-enter the active package without explicitly reopening the publication architecture.

`legacy/v1/` remains the separately frozen historical v1 analysis.