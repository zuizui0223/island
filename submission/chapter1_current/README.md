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
- H1/H2/raw colour/architecture/H3/exact H4/atomic H4 outputs are selected from `results/geography_20260924/`.

## Scientific spine

1. **H1 — pattern:** isolation is associated with a recurrent multivariate floral/reproductive response across four regions, but individual traits are not uniformly positive.
2. **H2 — decomposition:** reproductive assurance increases, while floral accessibility remains associated with isolation after reproductive-assurance adjustment in the primary analysis; detailed colour × architecture responses are context dependent.
3. **H3 — pressure:** experimental pollen limitation increases with corrected geographic isolation (beta = 0.0919, p = 0.0157).
4. **H4 — function:** the literal H2 reproductive-assurance and generalized-accessibility scores are associated with lower current pollen limitation in exact-species post-hoc functional triangulation.

## Submission boundary

The manuscript may claim pattern, conditional decomposition, an independent pollen-limitation gradient and post-hoc functional compatibility.

It must not claim:
- causal mediation from pollen limitation to trait evolution;
- a globally observed decline in pollinator abundance or visitation;
- a universal named bee/butterfly/bird mechanism;
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
