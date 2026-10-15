# Chapter 1 Q1 five-figure visual QA — 2026-10-15

## Assets inspected

- Five PDF figures in `figures/q1_main5/`, rendered from frozen result tables.
- Five SVG equivalents validated against `FIGURE_MANIFEST.json`.
- Five 120-dpi PNG previews in `figures/q1_main5/qa_png/`, rasterized from PDFs using `pdftoppm` in the [figure CI job](https://github.com/zuizui0223/island/actions/runs/38710540900).
- Original six source figures remain unchanged in `figures/`.

## Direct preview checks

| Figure | PDF-to-raster check | Remaining journal-layout risk |
| --- | --- | --- |
| 1: independent data streams | Two different geographic samples, numerators/denominators and untested dashed historical causal path legible; no clipping observed | Text may be compact if scaled to a narrow column |
| 2: pollen limitation | Coefficients, CIs and publication denominators readable; last-row coefficient label moved above point and rechecked | Confirm exact dimensions and label spacing at print scale |
| 3: seven traits × four regions | Original six-source display remains byte-for-byte archived. Q1 figure has been refaceted vertically from the same 112 frozen fit rows; values were not refit | Full-width 178-mm print-layout legibility check pending refreshed assets |
| 4: selfing-adjusted floral effects | Three scales explicitly separated and q values legible; southern q label repositioned and rechecked | Small q-value text should be checked at final page scale |
| 5: post-hoc GloPL link | Separately matched sample sizes, estimates and intervals legible; no visible overlap | Ensure legends and notes remain legible at final width |

## Scientific and submission boundary

Redrawing figure assets and label positions did not modify source H1–H4 result tables, inference families, coefficients or sample denominators. Neither pollinator disappearance nor historical mediation is inferred.

**Checked:** vector-file integrity, five PDFs rasterized, five 120-dpi visual previews directly inspected, visible axis/label clashes corrected.

**Still pending:** compile the complete journal-ready Main/SI files, then inspect the figures in that final composition; close authorship, trait-data rights and DOI gates. The A4 proof verifies actual print width but is NOT a completed journal submission.

## New print-width audit (follow-up)

The earlier widescreen raster previews were larger than the actual journal display width. At a standardized 178-mm width, minimum effective text in the old print-facing versions was only approximately 4.4 pt (Fig. 3), 4.0 pt (Fig. 4) and 4.4 pt (Fig. 5). The renderer was refactored to 7-inch-wide, vertically stacked panels; the updated workflow now generates `PRINT_WIDTH_AUDIT.json` with a fail-closed 6.5-pt minimum. The numbers and manuscript inference limits were not changed. Revised PDF-to-PNG inspection is still required after CI passes.


## True paper-width A4 proof — completed 2026-10-15

- [Five-page vector proof](figures/q1_main5/qa_print/Q1_Main5_178mm_A4_print_proof.pdf) created using `pypdf` without rasterizing the figure content. Each Q1 figure is centered at **178 mm** on an A4 page (16-mm horizontal margins).
- [Proof manifest](figures/q1_main5/qa_print/PRINT_PROOF_MANIFEST.json) independently verifies five pages, source PDF names, placements, and A4 fit; Figure 3 is 240.87 mm high and fits A4 with 28.06-mm vertical margins.
- [Text-size audit](figures/q1_main5/PRINT_WIDTH_AUDIT.json) passed all five figures: minimum effective glyph sizes at 178 mm are **7.159, 8.141, 8.071, 8.126, 8.144 pt** (Figures 1–5). Threshold: 6.5 pt.
- Pages 1–5 were visually inspected at 120 dpi in `figures/q1_main5/qa_print/page-1.png` … `page-5.png`. No figure or label clipping was apparent. The four-color geographic legend was restored to Figure 1; the seven-outcome × four-region Figure 3 is presented as four readable stacked regional panels using the unchanged 112 fitted values.
- [CI figure-proofs run](https://github.com/zuizui0223/island/actions/runs/38712444070) succeeded. Historical six-source and frozen numerical artifacts remain unchanged.

**Scope:** This proof establishes readable figures and page fitting on A4 at the intended 178-mm double-column width. It is not the journal's final typeset PDF, and it cannot address data-sharing rights or author approval.
