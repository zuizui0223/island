# Chapter 1 submission software/results archive candidate

Status: **prepared but not published**.

This record is intended to supply the permanent code/derived-results DOI required by the submission package without redistributing the complete third-party trait ledger.

## Proposed record type

Zenodo **Software** record linked to the exact GitHub submission release/tag.

Draft metadata:
- `zenodo/chapter1_submission_zenodo_metadata.json`

## Inclusion boundary

Include:

- exact repository source code at the frozen submission commit;
- `config/chapter1_submission_current.json`;
- corrected geography code and replay instructions under `scripts/geography_correction/`;
- submission builders under `scripts/submission/`;
- frozen final H1 directional result lock and replay code/results, including the WCVP regional-native-compatible provenance sensitivity;
- aggregate corrected H2/H3/H4 result tables under `results/geography_20260924/`;
- geometry validation/audit receipts;
- H1 optimizer/convergence audit;
- deterministic supplementary tables and manifests;
- submission-facing manuscript, captions, SI and data-accessibility documents;
- figure source files and rendered figures once their rendering PRs are merged.

Exclude from this software/results archive:

- the full 222,688-cell Chapter 1 trait ledger;
- raw third-party source excerpts or source files whose redistribution rights are unresolved;
- any file that would silently turn the software archive into redistribution of the restricted scientific trait database.

The rights-filtered public trait derivative is already archived separately:

- DOI `10.5281/zenodo.22704973`
- 46,274 redistributable cells

## Preferred citation architecture

The paper should cite three distinct objects:

1. **software + corrected derived results:** this release DOI, once minted;
2. **rights-filtered trait derivative:** DOI `10.5281/zenodo.22704973`;
3. **GloPL source data:** Bennett et al. 2018a,b / Dryad DOI `10.5061/dryad.dt437`.

This separation prevents a software/results DOI from being misdescribed as public release of the complete trait ledger.

## Freeze gate

Do not create the final GitHub release / Zenodo deposit until all of the following are true:

- final manuscript title is frozen;
- main Figures 1–6 are committed and guarded;
- Supplementary Figures S1–S7 are committed and guarded;
- final SI text/tables are on main;
- author metadata intended for the software record is confirmed;
- the submission commit SHA is known;
- the user explicitly approves the external release.

## Candidate tag

Suggested tag after freeze:

`chapter1-submission-v1`

The tag is only a proposed name. It has not been created.

## Post-publication update

After a DOI is minted:

- insert it into `submission/chapter1_current/DATA_ACCESSIBILITY_DRAFT.md`;
- insert it into `TITLE_PAGE_TEMPLATE.md`;
- update `SUBMISSION_CHECKLIST.md`;
- cross-link the software DOI with trait-subset DOI `10.5281/zenodo.22704973`;
- rerun current submission CI.
