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
- `figures/Figure3_constraint_response_triangle.svg` — publication-facing conceptual Figure 3 source.
- `SUPPLEMENT_PLAN.md` — Supplementary Information assembly map.
- `GRAPHICAL_ABSTRACT_BRIEF.md` — constraint–response triangle concept.
- `GRAPHICAL_ABSTRACT_SHORT_TEXT.md` — <=500-character graphical-abstract text draft.
- `GEB_FALLBACK.md` — ready-to-convert GEB structured-abstract/double-anonymous fallback.
- `SUBMISSION_CHECKLIST.md` — remaining journal-specific work and claim boundary.
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

Main Figures 1–6 and the poster workflow were regenerated on 25 September 2026 from the corrected PR #242 result surface. The repository keeps the source tables and captions as the scientific source of truth; exported PDF/PNG/PPTX files are submission artifacts rather than inferential inputs.

## Remaining work

Ecology Letters remains the first-shot target **only after the data-policy gate is closed**. The scientific narrative and format fit a Letter. Remaining blockers are: EL data-policy clearance or rights closure, final author metadata, permanent archive DOI(s), funding/acknowledgements/conflicts/contributions, final reference styling and assembly of the planned Supplementary Information.
