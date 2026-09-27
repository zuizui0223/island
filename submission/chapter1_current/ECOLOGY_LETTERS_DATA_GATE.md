# Ecology Letters data-policy gate

Checked against the Ecology Letters author guidelines on 27 September 2026.

Official source:
- https://onlinelibrary.wiley.com/page/journal/14610248/homepage/forauthors.html

## Why this is a real submission gate

Ecology Letters requires manuscripts using data/code to make the analysed raw data (or the exact subset of existing data used), metadata, analysis code and derived products accessible to editors/reviewers at submission and permanently archived in an external repository before publication. The journal explicitly states that acceptance is contingent on data-editor reproducibility.

The current Chapter 1 scientific database contains:

- 222,688 resolved species × axis cells;
- 46,274 cells currently classified as redistributable;
- 176,414 resolved cells currently blocked or review-required for redistribution;
- full-database SHA-256:
  `a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`.

The 46,274-cell rights-filtered public subset is therefore **not** the complete analysis input and cannot by itself satisfy a reproducibility claim for the paper.

## Ecology Letters decision rule

Ecology Letters remains the first-shot target only if **one** of the following is closed before submission:

1. **Rights closure:** the analysis-used trait cells can be archived under valid redistribution permissions; or
2. **Editorial exception:** Ecology Letters explicitly confirms that the third-party legal/licensing restrictions can be handled through an approved exception and that the proposed reviewer-access/reconstruction package is sufficient.

Do not submit to Ecology Letters while this gate is unresolved.

## What can be archived now

Publicly archivable now:

- exact code commit / software release;
- corrected geography outputs and receipts;
- model-ready derived result tables that do not reproduce restricted source rows;
- 46,274-cell rights-filtered trait subset;
- machine-readable source/provenance/rights inventory;
- scripts that reconstruct analysis inputs from openly retrievable sources where licensing permits;
- GloPL/GBIF/GSHHG version/source receipts.

## GEB fallback

Global Ecology and Biogeography also mandates data/code sharing, but its current policy explicitly permits editorial exceptions when sharing conflicts with legal requirements and asks authors to justify restrictions in the Data and Code Availability Statement.

Therefore:

- **EL first shot:** conditional on closing the data-policy gate;
- **GEB fallback:** currently the safer policy fit if third-party redistribution restrictions cannot be cleared.

Official GEB source:
- https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html

## Immediate next action

Send the pre-submission data-policy inquiry in `ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md` before creating the final EL submission archive.
