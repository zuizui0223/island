# Chapter 1 submission readiness — updated 2026-10-15

## Current decision

**Q1-only empirical paper.** The scientific submission surface remains corrected geography plus final separate-trait H1; the *editorial narrative* now tests the pollen-supply-constraint premise and uses constraint (source H3) → plant reproductive/floral response (sources H1 + H2) → contemporary functional compatibility (source H4). `izu-core` Ch2/Q2 simulations are a separate paper. Composite/domain H1 scores are historical provenance only.

First-shot journal: **Ecology Letters**, conditional on closing the journal's data-policy gate.  
Fallback: **Global Ecology and Biogeography** if the data-sharing restriction cannot be accepted.

## Complete

### Science

- corrected 8,264-island geographic baseline;
- final H1: seven separate traits × four regions, broad All primary;
- all 112 H1 fits converged;
- self-compatibility increases with isolation in all four regions in broad All and WCVP-compatible All analyses;
- other reproductive, colour and structural responses are region dependent;
- Direct-only evidence-quality sensitivity and WCVP regional-native-compatible origin sensitivity completed;
- H2 reproductive-assurance / accessibility conditional decomposition completed;
- H3 independent GloPL pollen-limitation gradient completed;
- H4 exact-species functional triangulation completed;
- causal and claim boundaries fixed.

### Manuscript package

- working title and manuscript reorganized around a three-stage pollen-constraint question;
- working abstract: **133 words**; working main body: **4,675 words** by whitespace count including headings;
- five Q1 Main PDF/SVG pairs are rendered and machine-verified; Main captions and references are now numbered 1–5, with original six-source figures archived as replay provenance; 120-dpi PDF previews and the 178-mm, five-page A4 proof have been visually inspected; 178-mm minimum text audit passes (7.159–8.144 pt), but final composed Main/SI still needs QA;
- editing-only Main DOCX and 23-page text PDF, 28-page Main-plus-five-figures review PDF, SI DOCX and 14-page text PDF, 21-page SI-plus-seven-figures review PDF are committed under `submission/chapter1_current/editorial_proofs/`. All **37 independently rendered text-PDF pages** passed the canvas/text audit; in addition, LibreOffice rendered the actual DOCX files to **19 Main portrait pages and 12 SI landscape pages**, with nine SI tables reflowed at readable widths and all 31 Word-derived pages independently checked. The files are not submission-authorized;
- cover letter, novelty statement, comparative study-scale audit and prior literature placement synchronized to traitwise H1 (PR #315 merged);
- discussion now identifies reduced *functional pollinator service* as a literature-supported, **unmeasured** possible explanation for H3; no pollinator declines, guild turnover or historical mediation are inferred from the Chapter 1 results;
- current interpretation no longer relies on a pooled H1 score.

### Current scientific claim

> **Independent GloPL experiments reveal an isolation-associated pollen-limitation gradient; global island floras show selective reproductive and nonuniform floral responses, some persisting after selfing adjustment. Exact-species trait–PL associations provide post-hoc functional compatibility but do not identify historical causal mediation.**

## Status of prior post-PR #315 verification (historical 8 October check)

- Main branch contains the 8 October literature/scale clarification (squash commit `33ceede43a551a5e89df35905047fa405bf066a5`), with H1–H4 methods/results and frozen numerical tables unchanged.
- Final H1 and corrected-geography workflows triggered by the final PR head (`37726796591`, `37726796582`) were **queued, not passed** when checked on 8 October 2026. Do not describe these checks as green unless their completed status is independently confirmed.
- Official Ecology Letters data/code policy was rechecked on 8 October 2026 and still requires complete analysis-used data access at submission. The 46,274-cell redistribution-approved derivative does not replace the full 222,688 resolved-cell analysis ledger.
- Pre-submission data-policy inquiry is prepared; editorial permission is not documented in the repository.

## Blocking external or author inputs

1. **Ecology Letters data-policy clearance.**  
   The complete analysis ledger is not wholly redistributable. Do not submit until rights are closed or the editors explicitly approve the proposed restricted-data arrangement.

2. **Author metadata.**  
   Final author order, affiliations, emails, ORCIDs, corresponding-author details, funding, acknowledgements, conflicts, contributions and final approval remain author-confirmed inputs. The cover letter now holds the approval statement as an explicit confirmation placeholder.

3. **Final software/results DOI and tag.**  
   Freeze the final submission commit, mint the release DOI, then insert it into the final data-accessibility text.

4. **Journal-composed Main/SI QA.**  
   Five PDFs/SVGs, PNG previews, and a five-page A4 proof with all figures at 178 mm are committed and visually checked. Actual print-width text is at least 7.159 pt and the original figures remain replayable. The remaining gate is their placement and reference consistency in the final compiled journal Main/SI, not a standalone figure size issue.

5. **Final journal files.**  
   Editorial review DOCX/PDF files and checked figure-inclusive PDFs now exist. After author metadata, rights, and DOI freeze, perform actual journal-template production, typography/word-count checks, approval and final submission.

## Claim ceiling

May claim:
- recurrent self-compatibility enrichment with isolation;
- regionally contingent selfing, colour and structural responses;
- selected selfing-adjusted accessibility responses;
- increasing experimental pollen limitation with isolation;
- exact-species functional compatibility of reproductive assurance/accessibility with lower current pollen limitation.

Must not claim:
- family-wise proof of a universal syndrome;
- historical causal mediation;
- global pollinator decline;
- a universal named pollinator mechanism;
- within-lineage evolution rather than assemblage filtering.

## Working branch editorial status (15 October)

- Q1 manuscript and SI crosswalk have been rewritten on PR #319's branch. The six-edge PR #319 proposal remains an archival relationship index, **not** the Main Results outline.
- No frozen H1–H4 tables, source scripts, p-value rules or result locks were modified by this editorial restructuring.
- [Main/Supp and figure migration contract](CH1_MAIN_SUPP_FIGURE_PLAN_20261015.md). Completion of this PR is not equivalent to readiness for Ecology Letters submission.
