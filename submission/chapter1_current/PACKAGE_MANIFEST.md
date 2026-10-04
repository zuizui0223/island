# Chapter 1 submission package — Ecology Letters first shot

Current package:
- `MANUSCRIPT.md` — clean submission-facing article draft synchronized to the final finite-cluster H1 and corrected geography.
- `FIGURE_CAPTIONS.md` — captions for Main Figures 1–6 using corrected estimates.
- `COVER_LETTER_DRAFT.md` — Ecology Letters first-shot cover letter.
- `ECOLOGY_LETTERS_TARGET.md` — target-journal scope, constraints and fallback logic.
- `NOVELTY_STATEMENT.md` — editorial novelty statement.
- `TITLE_PAGE_TEMPLATE.md` — required title-page fields and current counts.
- `DATA_ACCESSIBILITY_DRAFT.md` — rights-aware data-accessibility language.
- `ECOLOGY_LETTERS_DATA_GATE.md` — EL submission-policy gate for incomplete trait-ledger redistribution rights.
- `ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md` — draft pre-submission editorial inquiry.
- `figures/FIGURE_MANIFEST.json` — complete Main Figures 1–6 package manifest.
- `figures/Figure1_global_scope_corrected.*` through `Figure6_H3_H4_functional_bridge.*` — Main Figures 1–6 in SVG/PDF.
- `SUPPLEMENT_PLAN.md` — Supplementary Information structure.
- `SUPPLEMENTARY_INFORMATION_DRAFT.md` — internal SI source map.
- `SUPPLEMENTARY_INFORMATION.md` — submission-facing S1–S7 text.
- `supplement/SUPPLEMENT_TABLES_MANIFEST.json` — deterministic table/source hash manifest.
- `supplement/Table_S1_data_summary.csv` through `Table_S6b_H4_atomic.csv` — generated SI tables.
- `supplement/figures/FIGURE_MANIFEST.json` — Supplementary Figure S1–S7 package manifest.
- `supplement/figures/Figure_S1_*` through `Figure_S7_*` — Supplementary Figures S1–S7 in SVG/PDF.
- `GRAPHICAL_ABSTRACT_BRIEF.md` — graphical-abstract concept.
- `GRAPHICAL_ABSTRACT_SHORT_TEXT.md` — graphical-abstract short text.
- `GEB_FALLBACK.md` — structured-abstract and double-anonymous fallback conversion.
- `SUBMISSION_CHECKLIST.md` — completion checklist and claim boundary.
- `SUBMISSION_READINESS.md` — final current-state gate for submission.

Primary scientific source of truth:
- `config/chapter1_h1_final_directional_result_lock.json`
- `config/chapter1_submission_current.json`
- `submission/chapter1_current/MANUSCRIPT.md`
- `results/h1_final_directional_20261003/`
- `results/geography_20260924/`

Superseded provenance:
- `config/chapter1_v14_canonical_result_lock.json`
- the uncorrected v14 manuscript/results are retained for audit only and must not be submitted unchanged.

Main figure sequence:
1. Global geographic and data scope
2. Database construction and analytical workflow
3. Constraint–response triangle / inferential boundary
4. H1 global-average directional tendency and regional modulation
5. H2 conditional pathway decomposition
6. H3 pollen-limitation pressure and H4 functional bridge

Main Figures 1–6 are committed in SVG/PDF and bound to the final H1 result lock plus corrected H2–H4 tables by deterministic renderer + CI. Supplementary Figures S1–S7 are likewise committed and guarded. The committed corrected tables remain the inferential source of truth.

The first-shot target is Ecology Letters **conditional on closing the data-policy gate**. SI text, deterministic tables and Supplementary Figures S1–S7 are assembled; only final journal-file compilation remains after metadata freeze. The manuscript is within the Letter limits (main text <5,000 words; 6 display items; abstract <150 words). If third-party redistribution restrictions cannot be cleared or accepted by the EL editors, switch before submission to Global Ecology and Biogeography, whose policy explicitly permits editorial exceptions for legal data-sharing restrictions.
