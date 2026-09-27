# Chapter 1 submission package — Ecology Letters first shot

Current package:
- `MANUSCRIPT.md` — clean submission-facing article draft updated to PR #242 corrected geography.
- `FIGURE_CAPTIONS.md` — captions for Main Figures 1–6 using corrected estimates.
- `COVER_LETTER_DRAFT.md` — Ecology Letters first-shot cover letter.
- `ECOLOGY_LETTERS_TARGET.md` — target-journal scope, constraints and fallback logic.
- `NOVELTY_STATEMENT.md` — editorial novelty statement.
- `TITLE_PAGE_TEMPLATE.md` — required title-page fields and current counts.
- `DATA_ACCESSIBILITY_DRAFT.md` — rights-aware data-accessibility language.
- `SUPPLEMENT_PLAN.md` — Supplementary Information structure.
- `SUPPLEMENTARY_INFORMATION.md` — assembled S1–S7 text with Tables S1–S6 and Supplementary Data manifest.
- `GRAPHICAL_ABSTRACT_BRIEF.md` — graphical-abstract concept.
- `SUBMISSION_CHECKLIST.md` — completion checklist and claim boundary.

Primary scientific source of truth:
- `config/chapter1_submission_current.json`
- `docs/chapter1_corrected_submission_20260924.md`
- `results/geography_20260924/`

Superseded provenance:
- `config/chapter1_v14_canonical_result_lock.json`
- the uncorrected v14 manuscript/results are retained for audit only and must not be submitted unchanged.

Main figure sequence:
1. Global geographic and data scope
2. Database construction and analytical workflow
3. Inferential structure / working hypothesis
4. H1 recurrent multivariate response
5. H2 conditional pathway decomposition
6. H3 pollen-limitation pressure and H4 functional bridge

Main Figures 1–6 were regenerated from the corrected tables under `results/geography_20260924/` on 25 September 2026. Exported PNG/PDF/PPTX files are delivery artifacts; the committed corrected tables remain the scientific source of truth.

The first-shot target is Ecology Letters. The current manuscript is within the Letter limits (main text <5,000 words; 6 display items; abstract <150 words). Global Ecology and Biogeography is the fallback if editorial rejection is based on generality/novelty rather than scientific validity.
