> **Current H1 surface:** final separate-trait finite-cluster inference in `results/h1_final_traitwise_t_20261004/`; no pooled/domain H1 score is active.

> **Implementation status:** Submission-facing S1–S7 text is assembled in `SUPPLEMENTARY_INFORMATION.md`. Active deterministic tables are S1, traitwise H1 S2, finite-cluster H2 S3a/S3b, H3 S5 and H4 score S6a. Table S4 remains the complete raw-pattern file family.

# Supplementary Information plan — Ecology Letters first shot

The main paper should remain conceptual and compact. Technical audit detail, exhaustive coefficients and provenance belong here rather than in the main narrative.

## Appendix S1 — Geographic exposure and correction

- GSHHG 2.3.7 coastline construction;
- continental sibling reconstruction;
- exclusion of one continental split component;
- 1,113 repaired spurious island zero distances;
- 996 retained true continental GloPL zeros;
- spherical minor-great-circle arc implementation;
- numerical validation and geometry receipts;
- statement that the correction is a post-hoc measurement repair.

## Appendix S2 — Island flora and trait evidence

- GBIF occurrence assignment;
- taxonomic normalization;
- trait ontology;
- evidence precedence;
- missingness handling;
- all-analysis versus Direct-only scope;
- 106,295 species × 3 axes;
- 222,688 / 318,885 resolved cells;
- provenance and redistribution-rights boundary.

## Appendix S3 — H1 traitwise recurrence and provenance sensitivity

- seven individual reproductive, colour and structural traits;
- broad contemporary flora as primary scope;
- WCVP regional-native-compatible origin sensitivity;
- Direct-only evidence-quality sensitivity;
- spatial-cluster finite-sample t inference;
- pointwise intervals and two-sided P values;
- explicit statement that the 28 individual tests do not provide a family-wise syndrome test;
- historical composite/joint/domain analyses retained only as provenance.

## Appendix S4 — H2 full decomposition

- selfing_core construction;
- generalized_accessible construction;
- adjusted plain-colour model;
- five raw colour states;
- form/symmetry/tube components;
- colour × architecture joint prevalence;
- architecture conditional on colour;
- complete FDR families;
- tropical Direct-only accessibility remains positive but is not FDR-supported (q = 0.1267).

## Appendix S5 — H3 pollen-limitation analysis

- GloPL version receipt;
- aggregation to measurement cells;
- publication-total weight one;
- measurement-condition controls;
- publication-cluster-robust covariance;
- primary, no-zero and supplemental-only sensitivities.

## Appendix S6 — H4 exact-species bridge

- exact-match normalization;
- no synonym rescue / no genus fallback;
- literal H2 score reconstruction;
- matched species/publication counts;
- score-model coefficients;
- atomic-trait sensitivities;
- explicit post-hoc status.

## Appendix S7 — Inferential and validation boundary

- conditional decomposition ≠ mediation;
- pollen limitation ≠ pollinator abundance;
- trait architecture ≠ realized pollinator identity;
- sorting versus within-lineage evolution unresolved;
- support-limited prospective H4 audit stopped before outcome unblinding;
- historical pre-corrected workflows excluded from the current submission surface.

## Supplementary tables

- Table S1: deterministic data/geography summary;
- Table S2: final 112-row traitwise H1 table;
- Table S3a/S3b: All and Direct-only H2 finite-cluster audits;
- Table S4: complete corrected raw-colour/architecture CSV family;
- Table S5: finite-publication H3 primary and sensitivity models;
- Table S6a: finite-publication H4 exact H2-score bridge;
- atomic H4 sensitivities remain reported in SI text from the corrected atomic source.

## Supplementary figures

Prioritize diagnostics and robustness, not additional narrative figures:

- S1 corrected versus old distance distribution;
- S2 regional H1 Direct-only comparison;
- S3 H2 raw colour-state response matrix;
- S4 colour-conditioned architecture heatmap;
- S5 H3 sensitivity estimates;
- S6 H4 atomic-trait sensitivity forest plot;
- S7 support/missingness and evidence-scope diagnostics.
