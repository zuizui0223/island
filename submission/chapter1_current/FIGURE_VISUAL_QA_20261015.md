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
| 3: seven traits × four regions | Promoted byte-for-byte from frozen H1; trait rows, regional estimates, intervals and evidence scopes visible | Four-region panel is information-dense and may need full-width placement |
| 4: selfing-adjusted floral effects | Three scales explicitly separated and q values legible; southern q label repositioned and rechecked | Small q-value text should be checked at final page scale |
| 5: post-hoc GloPL link | Separately matched sample sizes, estimates and intervals legible; no visible overlap | Ensure legends and notes remain legible at final width |

## Scientific and submission boundary

Redrawing figure assets and label positions did not modify source H1–H4 result tables, inference families, coefficients or sample denominators. Neither pollinator disappearance nor historical mediation is inferred.

**Checked:** vector-file integrity, five PDFs rasterized, five 120-dpi visual previews directly inspected, visible axis/label clashes corrected.

**Still pending:** compile journal-ready Main/Supp PDFs; check at actual print scale; close authorship, data-rights and DOI gates. Raster preview QA is not final journal acceptance.
