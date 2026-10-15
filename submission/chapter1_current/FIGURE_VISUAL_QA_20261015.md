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

**Still pending:** compile journal-ready Main/Supp PDFs; check at actual print scale; close authorship, data-rights and DOI gates. Raster preview QA is not final journal acceptance.

## New print-width audit (follow-up)

The earlier widescreen raster previews were larger than the actual journal display width. At a standardized 178-mm width, minimum effective text in the old print-facing versions was only approximately 4.4 pt (Fig. 3), 4.0 pt (Fig. 4) and 4.4 pt (Fig. 5). The renderer was refactored to 7-inch-wide, vertically stacked panels; the updated workflow now generates `PRINT_WIDTH_AUDIT.json` with a fail-closed 6.5-pt minimum. The numbers and manuscript inference limits were not changed. Revised PDF-to-PNG inspection is still required after CI passes.
