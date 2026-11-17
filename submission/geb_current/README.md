# Global Ecology and Biogeography — Chapter 1 first-shot

**Active editorial target, 17 November 2026. Not submitted.** The underlying Chapter 1 H1–H4 scientific estimates and data-version lock remain `config/chapter1_submission_current.json` (4 October 2026).

## Start here

- [Double-anonymous research article](MANUSCRIPT_BLINDED.md) — seven-item structured abstract, alphabetical keywords, source-locked H1–H4 Results and stated inferential limits
- [GEB figure legends](FIGURE_CAPTIONS_GEB.md) — six figures, lowercase panel convention and equal-area projection
- [Supplementary Information](SUPPLEMENTARY_INFORMATION.md) — corrected H1–H4 tables and methodological safeguards
- [Identifying title-page template](TITLE_PAGE_IDENTIFYING_TEMPLATE.md) — requires verified corresponding/coauthor information
- [GEB cover-letter draft](COVER_LETTER_GEB.md) — editorial significance under 250 words, approval pending
- [Rights-aware Data and Code Availability](DATA_CODE_AVAILABILITY.md) — source permissions and editor decision are not presumed
- [GEB submission checklist](SUBMISSION_CHECKLIST_GEB.md) — completion versus external blockers
- [GEB review exporter](../../scripts/submission/build_chapter1_geb_review.py) and [Mollweide figure renderer](../../scripts/submission/render_chapter1_geb_figures.py)

## Editorial changes, not new analyses

*Ecology Letters* was the former target and its files under `submission/chapter1_current/` are retained for scientific provenance. GEB was selected as first shot on 17 November 2026. The current GEB text keeps the original body (4,958 words by repository Markdown count) and H1–H4 estimates, changing only the target-specific abstract, keywords, anonymity, data statement, figure format and publication package.

The GEB Research Article instructions request around 5,000 main-text words, a seven-field structured abstract of at most 300 words, 6–10 alphabetic keywords, a separate identifying title page and double-anonymous manuscript, line numbers, and a cover letter uploaded separately as PDF. Worldwide distribution maps must use an equal-area projection. Official instructions: https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html

## Not closed

- **Data permissions:** 222,688 analysed trait-axis cells; only 46,274 currently authorized for public redistribution. The other 176,414 cells include both third-party legal restrictions and review-required records. Wiley's data-sharing policy permits discretionary exceptions for legal restrictions, not blanket release of all data. This must be resolved with the editor.
- **Authorship:** final names/order, affiliations, ORCID, funding, COI and approval pending.
- **Anonymous repository access:** no reviewer link minted; do not expose a personal GitHub profile or author metadata in blinded documents.
- **Software/results archive:** no final source DOI or release tag yet; explicit user/coauthor approval is necessary.
- **File QA:** complete review rendering, anonymization, maps and figure legends must pass CI before GEB package is promoted.

The package builder may export **review copies**, but must not submit to Wiley, email editors, release a dataset, or identify author metadata without approval.
