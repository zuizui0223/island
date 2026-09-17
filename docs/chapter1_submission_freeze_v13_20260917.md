# Chapter 1 v13 submission freeze — 2026-09-17

Status: review candidate; not merged.

## Canonical publication surface

1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
3. `docs/chapter1_v13_submission_figure_sync_20260917.md`
4. `config/chapter1_v13_functional_bridge_result_lock.json`
5. `docs/chapter1_unified_hypothesis_20260917.md`

Historical v11/v12 manuscripts and result locks remain immutable provenance and are not replaced or rewritten by v13.

## Frozen paper architecture

- **H1 — recurrent global island syndrome:** isolation-associated plant responses recur in the classic island-syndrome direction across all four broad contemporary-flora contexts; regional response vectors may differ and are not assumed equal.
- **H2 — global pollination constraint:** experimental pollen limitation increases with geographic separation in the full-global GloPL analysis.
- **H3 — dual plant response pathways:** reproductive assurance and generalized/accessibility floral architecture are partially separable plant-side response families.
- **H4 — functional bridge:** exact-species GloPL trait-state comparisons are explicitly `posthoc_functional_triangulation`, not confirmatory mediation. Autonomous selfing is the strongest reproducible functional bridge to lower current pollen limitation; architecture evidence is weaker and heterogeneous.
- **H5 — taxonomic realization:** broad contemporary all-observed and defended native Palearctic responses have different taxonomic representation and must not be collapsed into one assembly mechanism.

## Canonical v13 functional bridge

- workflow run: `35141624253`
- artifact: `10465048981`
- artifact digest: `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`
- inferential role: `posthoc_functional_triangulation`
- the previously frozen Route A/B distance-by-trait moderation failures remain negative results and are not reclassified.

## Reused frozen evidence

- all-data primary probability run: `34961775336`, artifact `10394245237`, digest `sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b`
- full-global GloPL run: `35090599662`, artifact `10444156159`, digest `sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623`
- broad source-free H3A run: `35043562731`, artifact `10426327347`
- defended native Palearctic P1/P3 provenance remains frozen in the historical v11/v12 surfaces.

## Claim ceiling

v13 may state association, recurrence, partial pathway separation, functional compatibility/triangulation, and flora-layer-specific taxonomic realization.

v13 must not state that pollen limitation historically caused trait evolution, that pollen limitation statistically mediates the global syndrome, that pollinator abundance or visitation globally declines with isolation, that GloBI proves the causal mechanism, that Palearctic genus structuring explains the broad all-observed response, or that post-hoc functional triangulation is confirmatory.

## Submission gate

Promotion to review-ready requires a fresh branch-head audit after this freeze commit that passes:

- all `tests/test_chapter1_v13_*.py`;
- Ruff on v13 source/tests;
- historical v11/v12 immutability guard;
- v13 submission-surface guards for README, `docs/PAPER_PIPELINE.md`, manuscript headlines, evidence-role labels, and prohibited causal wording.

No merge to `main` is authorized by this freeze document.
