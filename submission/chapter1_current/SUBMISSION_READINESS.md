# Chapter 1 submission readiness — 2026-10-04

## Current decision

The Chapter 1 scientific analysis and submission-facing inference are **frozen and synchronized on main**.

First-shot journal: **Ecology Letters**, conditional on closing the journal's data-policy gate.

Fallback: **Global Ecology and Biogeography** if the third-party trait-data restriction cannot be accepted by Ecology Letters.

No further exploratory analysis should alter the frozen H1–H4 claims unless a new analysis phase is explicitly opened.

## Complete

### Science

- corrected 8,264-island geographic baseline;
- H1a: positive global-average classic-island direction under Paule-Mandel + modified Hartung-Knapp inference;
- strict four-region recurrence is unsupported because northern mid-latitudes are weak;
- H1b: regional heterogeneity is supported in the primary all-analysis; the Direct-only exact cluster-jackknife sensitivity is weaker and does not support formal heterogeneity;
- exact delete-one-spatial-cluster jackknife sensitivity retains H1a in both evidence scopes;
- WCVP regional-native-compatible replay applies the same final one-dimensional seven-indicator finite-cluster H1 to 513,320 island×species rows on 2,372 islands and reproduces the positive global-average tendency plus weak northern mid-latitudes;
- northern-midlatitude isolation-support diagnostic is complete: the region is much more mainland-proximate, but adjusted residual isolation variation is not depleted, so the weak slope is not explained by a simple lack of distance range or island-count power;
- raw three-axis and floristic-origin results are retained as descriptive phenotype/provenance audits, not confirmatory substitutes for H1;
- H1 Direct-only northern-high optimizer warning independently closed;
- H2 reproductive-assurance / accessibility conditional decomposition;
- H3 independent GloPL pollen-limitation gradient;
- H4 exact-species functional triangulation;
- causal/claim boundaries fixed in manuscript, SI and result locks.

### Manuscript package

- final title and **137-word abstract** synchronized to H1a/H1b;
- main text **4,996 words** by repository Markdown count; final formatted journal recount still required;
- 18 references in Ecology Letters style;
- novelty positioned against recent island-colonization and range-edge literature;
- Ecology Letters cover letter;
- novelty statement;
- running title and keywords;
- GEB fallback structured abstract.

### Figures and Supplementary Information

- Main Figures 1–6 are committed and synchronized;
- Figure 4 is generated directly from the final H1 result lock and shows regional directional scores, H1a/H1b, weighting/leverage sensitivities and provenance divergence;
- main-figure manifest and deterministic renderer are current;
- Supplementary Information S1–S7 assembled;
- deterministic Tables S1, S2a–S2h, S3, S5, S6a and S6b;
- Supplementary Figures S1–S7 committed in SVG/PDF;
- SI figure manifest and deterministic renderer;
- main submission CI passes at commit `4e8ce8f40210f00c7e64d5aba53dc91122591df0`.

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
