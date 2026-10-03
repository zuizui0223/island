# Chapter 1 submission readiness — 2026-09-28

## Current decision

The corrected Chapter 1 scientific package is **analysis-complete; main-figure regeneration is in progress after the final inference audit**.

First-shot journal: **Ecology Letters**, conditional on closing the journal's data-policy gate.

Fallback: **Global Ecology and Biogeography** if the third-party trait-data restriction cannot be accepted by Ecology Letters.

## Complete

### Science

- corrected 8,264-island geographic baseline;
- H1a positive global-average predeclared island-syndrome direction supported in all-analysis and Direct-only random-effects synthesis;
- H1b regional heterogeneity supported; strict four-region recurrence explicitly fails because northern mid-latitudes are weak;
- raw three-axis/state models retained as descriptive reorganization and floristic-origin diagnostics, not confirmatory H1 tests;
- H2 reproductive-assurance / accessibility conditional decomposition re-audited with finite-cluster inference;
- H3 independent GloPL pollen-limitation gradient retained under finite-publication inference, including offshore-only robustness;
- H4 exact-species functional triangulation retained under finite-publication inference;
- causal/claim boundaries fixed in manuscript and SI.

### Manuscript package

- title and 140-word abstract;
- main text 4,977 words;
- 18 references in Ecology Letters style;
- novelty positioned against recent island-colonization and range-edge literature;
- Ecology Letters cover letter;
- novelty statement;
- running title and keywords;
- GEB fallback structured abstract.

### Figures and Supplementary Information

- deterministic Main Figure 1–6 renderer updated to the final directional H1 and finite-cluster H2–H4 result surface;
- committed SVG/PDF figures must be regenerated once from that updated renderer before submission;
- Supplementary Information S1–S7 assembled;
- deterministic Tables S1, S2a–S2h, S3, S5, S6a and S6b;
- Supplementary Figures S1–S7 committed in SVG/PDF;
- SI figure manifest and deterministic renderer;
- CI verifies figure files, manifests, current corrected values and inferential boundaries.

### Data/reproducibility preparation

- rights-filtered Chapter 1 trait derivative published at Zenodo:
  DOI `10.5281/zenodo.22704973`;
- GloPL source data cited with Dryad DOI `10.5061/dryad.dt437`;
- unpublished software/results Zenodo metadata prepared;
- explicit inclusion/exclusion boundary prevents redistribution of the restricted full trait ledger;
- current submission CI runs manuscript/package, SI, figure and archive-candidate guards.

## Blocking items that require external or author input

### 1. Ecology Letters data-policy clearance

Current full scientific trait ledger:
- 222,688 resolved cells;
- 46,274 currently redistributable;
- 176,414 blocked or review-required for redistribution.

Before Ecology Letters submission, either:

1. close rights for the analysis-used trait cells; or
2. obtain explicit written editorial clearance confirming that the third-party licensing restriction and proposed reviewer-access/reconstruction plan are acceptable.

Prepared inquiry:
- `ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md`
- To: `ecolets2@cefe.cnrs.fr`
- Cc: `ecolets@cefe.cnrs.fr`

Do not submit to Ecology Letters before this gate is closed.

### 2. Author metadata

Required inputs still absent from the repository:

- final author order;
- affiliations;
- author emails;
- corresponding-author postal/contact details;
- ORCID identifiers;
- funding;
- acknowledgements;
- conflict-of-interest confirmation;
- author-contribution statement.

These should be filled only from author-confirmed information.

### 3. Final permanent software/results DOI

Prepared but not published:

- `zenodo/CHAPTER1_SUBMISSION_RELEASE.md`
- `zenodo/chapter1_submission_zenodo_metadata.json`

After the target-journal/data-policy decision and author metadata freeze:

1. freeze final submission commit;
2. create proposed tag `chapter1-submission-v1`;
3. create GitHub/Zenodo software-results release;
4. insert the minted DOI into Data Accessibility/title page;
5. rerun submission CI.

External publication requires explicit approval.

### 4. Final journal files

After metadata and DOI freeze:

- compile final manuscript file;
- compile final Supplementary Information PDF/file;
- check journal upload naming;
- perform final visual inspection of all figures at submission scale.

## Current scientific claim in one sentence

**Geographic isolation is associated with a positive global-average floral/reproductive island-syndrome direction whose strength differs strongly among regions; reproductive assurance and accessibility are the most consistent functional components, and the same isolation axis is independently associated with increasing pollen limitation.**
