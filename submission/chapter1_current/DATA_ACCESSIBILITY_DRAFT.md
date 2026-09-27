# Data accessibility draft — Ecology Letters conditional version

## Status

**Not submission-ready for Ecology Letters until `ECOLOGY_LETTERS_DATA_GATE.md` is closed.** The wording below records the truthful current state; it is not a substitute for the journal's requirement that the analysed data/subset be accessible to reviewers and archived appropriately.

## Submission-facing draft

The rights-filtered Chapter 1 Database 1.0 public derivative is already archived at Zenodo (DOI: `10.5281/zenodo.22704973`). It contains 46,274 redistribution-authorized cells and is explicitly not the complete scientific analysis ledger. The exact analysis code, corrected geography outputs and model-ready derived result tables still require a submission-specific permanent archive before submission. The full frozen Chapter 1 trait ledger cannot be redistributed wholesale because it integrates third-party source material with heterogeneous reuse terms. Its exact scientific identity is preserved by an immutable SHA-256 digest and machine-readable source/provenance records. A rights-filtered public subset contains only cells explicitly authorized for redistribution; omission from that subset reflects redistribution rights rather than biological missingness or scientific exclusion.

The GloPL pollen-limitation data are from Bennett et al. (2018) and are version-pinned in the analysis provenance. GBIF occurrence inputs and GSHHG geography are likewise referenced through their source/version receipts.

## Before submission

Before submission, add permanent archive DOI(s) for:

1. exact submission code commit / software archive;
2. corrected derived H1–H4 result tables and geography receipts;
3. machine-readable submission provenance/rights manifest identifying restricted source lineages.

Already resolved:

- rights-filtered Chapter 1 trait subset: Zenodo DOI `10.5281/zenodo.22704973`.

Do **not** state that the rights-filtered public trait subset is the exact analysis database.

## Scientific database identity

- frozen full trait ledger SHA-256: `a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`;
- resolved scientific cells: 222,688;
- currently redistribution-authorized public cells: 46,274;
- current paper selector: `config/chapter1_submission_current.json`.


## GEB fallback wording if EL clearance fails

Global Ecology and Biogeography explicitly permits editorial exceptions where sharing conflicts with legal requirements. For a GEB submission, retain the same public code/derived products/provenance package, describe the third-party restrictions explicitly, and request the legal-restriction exception rather than implying full public redistribution.


## Existing permanent dataset citation

Chapter 1 Database 1.0 rights-filtered public derivative. Zenodo. DOI: `10.5281/zenodo.22704973`. This deposit contains 46,274 of the 222,688 resolved scientific cells and must be described as a rights-filtered derivative rather than the complete analysis ledger.
