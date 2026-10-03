# Chapter 1 submission package — corrected geography draft

Status: **Ecology Letters first-shot draft conditional on data-policy clearance; Global Ecology and Biogeography fallback if that gate cannot be closed**.

This directory is the clean submission-facing package for Chapter 1. The scientific baseline is now the corrected geography merged in PR #242 and selected by `config/chapter1_submission_current.json`. The older v14 result lock is retained only as superseded provenance.

## Current package

- `MANUSCRIPT.md` — clean submission draft updated to the corrected geography baseline.
- `FIGURE_CAPTIONS.md` — captions for Main Figures 1–6 using corrected H1–H4 estimates.
- `COVER_LETTER_DRAFT.md` — Ecology Letters first-shot cover letter.
- `ECOLOGY_LETTERS_TARGET.md` — current journal fit, limits and editorial positioning.
- `NOVELTY_STATEMENT.md` — concise conceptual novelty statement.
- `TITLE_PAGE_TEMPLATE.md` — Ecology Letters title-page requirements and current counts.
- `DATA_ACCESSIBILITY_DRAFT.md` — rights-aware data/code availability wording.
- `ECOLOGY_LETTERS_DATA_GATE.md` — hard submission gate created by third-party redistribution restrictions.
- `ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md` — pre-submission editor inquiry for that gate.
- `figures/` — complete Main Figures 1–6 in SVG/PDF plus `FIGURE_MANIFEST.json` and deterministic renderer.
- `SUPPLEMENT_PLAN.md` — Supplementary Information assembly map.
- `SUPPLEMENTARY_INFORMATION_DRAFT.md` — internal source map for SI assembly.
- `SUPPLEMENTARY_INFORMATION.md` — submission-facing S1–S7 text bound to deterministic tables and corrected outputs.
- `supplement/` — deterministic supplementary tables and hash manifest.
- `supplement/figures/` — Supplementary Figures S1–S7 in SVG/PDF plus figure manifest.
- `GRAPHICAL_ABSTRACT_BRIEF.md` — constraint–response triangle concept.
- `GRAPHICAL_ABSTRACT_SHORT_TEXT.md` — <=500-character graphical-abstract text draft.
- `GEB_FALLBACK.md` — ready-to-convert GEB structured-abstract/double-anonymous fallback.
- `SUBMISSION_CHECKLIST.md` — remaining journal-specific work and claim boundary.
- `SUBMISSION_READINESS.md` — current complete/blocked submission state and final external-input gates.
- `PACKAGE_MANIFEST.md` — package contents and source-of-truth pointers.

## Current baseline

- corrected analysis universe: **8,264 island units**;
- broad H1 union: **4,379 islands**;
- corrected distance: minimum minor-great-circle arc separation to source-matched GSHHG 2.3.7 continental coastlines on a mean-radius sphere;
- 1,113 formerly spurious island zero distances are now positive;
- 996 true continental GloPL site zeros remain zero;
- corrected H2/H3/H4 outputs are selected from `results/geography_20260924/`;
- final finite-cluster H1 outputs are frozen under `results/h1_final_directional_20261003/`; the older high-dimensional H1 Wald tables remain provenance/reference only.

## Scientific spine

1. **H1a — global-average tendency:** the frozen classic-island directional score is positive on average across the four predeclared regions under Paule-Mandel random effects with modified Hartung-Knapp inference (all-analysis mean = 0.0691, one-sided p = 0.0256; Direct-only mean = 0.0635, p = 0.0238).
2. **H1b — heterogeneity:** regional magnitude varies substantially. Heterogeneity is supported in the primary sandwich analysis (I² = 0.819 all-analysis; 0.685 Direct-only), while an exact delete-cluster stress test retains it in all-analysis (p = 0.00789) but not Direct-only (p = 0.131). The stronger claim that all four regions independently support the same direction fails because northern mid-latitudes are weak (intersection-union p = 0.134 / 0.106).
3. **H2 — decomposition:** measured reproductive assurance does not statistically absorb all floral accessibility change; finite-cluster FDR support is concentrated in northern high latitudes and tropical all-analysis, while colour is more region dependent. This is conditional decomposition, not mediation.
4. **H3 — pressure:** experimental pollen limitation increases with corrected geographic isolation (beta = 0.0919, finite-publication p = 0.01594); the offshore-only robustness gradient also remains positive (p = 0.02459).
5. **H4 — function:** the literal H2 reproductive-assurance and generalized-accessibility scores are associated with lower current pollen limitation in exact-species post-hoc functional triangulation (finite-publication p = 0.00417 and 0.02334).

## Submission boundary

The manuscript may claim a positive but leverage-sensitive global-average island-syndrome direction, strong regional heterogeneity, conditional H2 decomposition, an independent pollen-limitation gradient and post-hoc functional compatibility.

It must not claim:
- causal mediation from pollen limitation to trait evolution;
- a globally observed decline in pollinator abundance or visitation;
- a universal named bee/butterfly/bird mechanism;
- a universal four-region floral island syndrome or independent support in every region;
- uniform positive change in all seven H1 traits;
- that the tropical Direct-only H2 accessibility result is FDR-supported;
- that colour is a global H4 functional bridge;
- within-lineage evolution rather than assemblage composition;
- that the corrected geography is a prospective confirmation.

Prospective H4 validation audits remain outside the H4 result.

## Figure status

Main Figures 1–6 are committed under `figures/` in SVG/PDF and are reproducible from the corrected result surface with `scripts/submission/render_chapter1_main_figures.py`. Supplementary Figures S1–S7 are committed under `supplement/figures/` and are reproducible with `scripts/submission/render_chapter1_supplement_figures.py`. CI validates both figure manifests and file integrity.

## Remaining work

Ecology Letters remains the first-shot target **only after the data-policy gate is closed**. The scientific narrative, references, main figures, SI text/tables and SI figures are assembled. Remaining blockers are: EL data-policy clearance or rights closure; final author/corresponding-author/ORCID metadata; funding, acknowledgements, conflicts and author contributions; freezing the final submission commit/tag and minting the software/results DOI; and compiling the final journal-ready manuscript/SI files.
