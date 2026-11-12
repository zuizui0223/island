# Chapter 1 editorial review bundle — 12 November 2026

**Status: editorial review artifact, not a journal submission or archival release.**

The GitHub Actions workflow [Chapter 1 — rights-safe editorial review bundle](../../.github/workflows/build-chapter1-editorial-review.yml) runs from the branch `submission/chapter1-reviewpack-20261112`.

It checks the selected current submission contract, abstract/body word counts, and main-figure manifest. It then compiles:

- DOCX and review PDF for the main manuscript, including six main figures and captions.
- DOCX and review PDF for Supplementary Information, including seven supplementary figures.
- Incomplete title-page and unapproved cover-letter DOCX/PDF copies.
- Original Markdown sources, main and supplementary figure PDFs, six aggregate supplemental CSV tables, and a SHA-256 manifest.

**What it does not contain:** the restricted 222,688-cell scientific trait ledger, source excerpts, raw GloPL/GBIF files, or source-material records lacking redistribution permission. The separately released 46,274-cell derivative is not the complete analysed ledger.

**Open submission gates:**
1. Ecology Letters written clarification/rights closure concerning the 176,414 redistribution-restricted or review-required resolved cells.
2. Complete author metadata, approval, funding and conflict-of-interest declarations.
3. Explicit approval for the final software/results release/tag and permanent DOI.
4. Final journal-specific typesetting and full visual author review.

This review packet is generated from already frozen results. It does not change any H1–H4 estimates, release data publicly, send email, or submit the paper to a journal.

Current source: `config/chapter1_submission_current.json`. Historical H1–H5 composite hypotheses remain excluded from the active submission.
