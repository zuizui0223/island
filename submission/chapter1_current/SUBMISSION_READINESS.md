# Chapter 1 submission readiness — 2026-09-28

## Current decision

The corrected Chapter 1 analysis is **complete under the final finite-cluster H1 inference**; submission prose and Figure 4 are being synchronized to that result.

First-shot journal: **Ecology Letters**, conditional on closing the journal's data-policy gate.

Fallback: **Global Ecology and Biogeography** if the third-party trait-data restriction cannot be accepted by Ecology Letters.

## Complete

### Science

- corrected 8,264-island geographic baseline;
- H1a: positive global-average classic-island direction under Paule-Mandel + modified Hartung-Knapp inference;
- H1b: strong regional heterogeneity; strict four-region recurrence is unsupported;
- raw three-axis result retained as descriptive phenotype/provenance audit, not confirmatory H1;
- H1 Direct-only northern-high optimizer warning independently closed;
- H2 reproductive-assurance / accessibility conditional decomposition;
- H3 independent GloPL pollen-limitation gradient;
- H4 exact-species functional triangulation;
- causal/claim boundaries fixed in manuscript and SI.

### Manuscript package

- title and abstract revised for the final H1a/H1b claim; recount pending;
- main-text recount pending after final inferential revision;
- 18 references in Ecology Letters style;
- novelty positioned against recent island-colonization and range-edge literature;
- Ecology Letters cover letter;
- novelty statement;
- running title and keywords;
- GEB fallback structured abstract.

### Figures and Supplementary Information

- Main Figures 1–6 are committed; Figure 4 requires final H1a/H1b visual synchronization;
- main-figure manifest and deterministic renderer;
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

**Geographic isolation is associated with a positive average shift toward reproductive assurance and accessible floral architecture, but the magnitude and phenotype are strongly region dependent; the same gradient is independently associated with increasing pollen limitation.**
