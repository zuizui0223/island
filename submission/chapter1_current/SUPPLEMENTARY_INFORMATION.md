# Supplementary Information

## A recurrent global floral island syndrome extends beyond the selfing syndrome

This Supplementary Information accompanies the corrected Chapter 1 submission selected by config/chapter1_submission_current.json. All numerical results below are bound to results/geography_20260924/ and to the deterministic supplementary tables under submission/chapter1_current/supplement/.

The geography correction is a post-hoc repair of a known exposure mismatch. H4 is post-hoc functional triangulation. Neither is reclassified here as prospective confirmation.

---

# Appendix S1. Corrected geography and analysis scope

## S1.1 Corrected island universe

The plant analysis contains 8,264 island units. One GSHHG split component corresponding to continental Eurasia was excluded after the geography audit. The broad H1 analytic union contains 4,379 islands.

The original exposure calculation compared GSHHG island geometry with coarser Natural Earth continental polygons. Within the broad H1 union, 1,113 islands therefore received spurious zero distance. Recalculation against continental coastlines reconstructed from the same GSHHG 2.3.7 high-resolution archive assigned positive distance to all 1,113.

For GloPL, exposure is calculated from study-site coordinates. Of 1,248 unique GloPL sites, 996 lie on seeded continental land and legitimately retain zero distance.

### Table S1. Data and geography summary

The deterministic machine-readable version is:

- submission/chapter1_current/supplement/Table_S1_data_summary.csv

Key values are:

| Quantity | Corrected value |
| --- | ---: |
| Global island units | 8,264 |
| Broad H1 union | 4,379 |
| Accepted angiosperm species | 106,295 |
| Possible species × axis cells | 318,885 |
| Resolved species × axis cells | 222,688 |
| GloPL experimental rows | 2,969 |
| GloPL measurement cells | 1,408 |
| GloPL unique sites | 1,248 |
| GloPL publications | 919 |
| True continental GloPL zero-distance sites | 996 |

## S1.2 Distance algorithm and validation

Isolation is the minimum coastline-to-coastline separation between minor great-circle arcs on a mean-radius sphere (R = 6371.0088 km). Continental siblings are reconstructed before distance calculation and artificial dateline split edges are omitted. Segment indexing uses triangle-inequality pruning followed by exact arc minima.

Machine-readable receipts:

- results/geography_20260924/spherical_geometry_validation.json
- results/geography_20260924/distance_same_coastline_audit.json
- results/geography_20260924/corrected_universe_exclusion.csv
- results/geography_20260924/formerly_zero_islands_recalculated.csv
- results/geography_20260924/glopl_corrected_site_distances.csv

The numerical validation checks the implementation, not metre-scale accuracy of the coastline source.

---

# Appendix S2. Island floras, trait evidence and redistribution boundary

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was normalized to a fixed accepted-species axis and retained source provenance.

The scientific trait database contains 106,295 accepted angiosperm species and three raw evidence axes:

1. flower colour;
2. floral structural complexity;
3. reproductive assurance.

Of 318,885 possible species × axis cells, 222,688 are resolved (69.83%). Species-direct High/Medium evidence forms the Direct-only sensitivity. The primary all-analysis scope additionally admits validated lower-confidence evidence where direct evidence is unavailable. Missing trait information remains missing rather than becoming trait absence.

The full scientific ledger is identified by SHA-256:

a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a

The rights-filtered public derivative contains 46,274 redistribution-authorized cells and is archived at Zenodo:

DOI 10.5281/zenodo.22704973

This public derivative is not the complete scientific analysis ledger. Omitted cells are omitted because of redistribution-rights status, not biological missingness or exclusion from the analysis.

Supporting documentation:

- docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md
- docs/CHAPTER1_DATABASE_PUBLIC_SUBSET.md
- config/chapter1_database_versions/v1.0.0.yml
- submission/chapter1_current/ECOLOGY_LETTERS_DATA_GATE.md

---

# Appendix S3. H1 recurrent multivariate island response

H1 models seven atomic responses, each coded so that a positive isolation coefficient is in the predicted island-syndrome direction:

1. self-compatibility;
2. predominantly or obligately selfing mating system;
3. autonomous or delayed autonomous selfing;
4. plain colour;
5. generalized floral form;
6. actinomorphic symmetry;
7. shallow or open floral tube.

Within each geographic stratum, the seven-dimensional isolation vector is tested by a joint Wald test after beta-binomial models with standardized log isolation, island area and climate PC1–PC4. Spatial-block cluster-robust covariance is used.

### Table S2a. Complete H1 atomic coefficients

Machine-readable table:

- submission/chapter1_current/supplement/Table_S2a_H1_atomic.csv

This table contains 112 rows spanning all-analysis and Direct-only evidence scopes, all retained strata and all seven atomic responses. It preserves the frozen Direct-only northern-high shallow/open-tube optimizer flag and annotates it rather than silently rewriting the frozen corrected output.

### Table S2b. H1 joint response-vector tests

Machine-readable table:

- submission/chapter1_current/supplement/Table_S2b_H1_joint.csv

Broad all-observed results:

| Region | All-analysis q | Direct-only q |
| --- | ---: | ---: |
| Northern mid-latitude | 3.216 × 10^-10 | 1.370 × 10^-4 |
| Northern high latitude | 2.433 × 10^-5 | 3.794 × 10^-8 |
| Tropical | 3.498 × 10^-7 | 3.528 × 10^-5 |
| Southern extratropical | 2.504 × 10^-18 | 3.155 × 10^-31 |

All four regions support the seven-response vector in both evidence scopes. This is a multivariate recurrence claim, not a claim that every atomic coefficient is positive or individually supported.

In the primary all-analysis scope, 26 of 28 regional atomic coefficients are positive. The southern shallow/open-tube coefficient is negative (β = -0.21660, SE = 0.06083, nominal P = 0.000370). Northern-high generalized form is positive but weak (P = 0.1018), and southern selfing mating system is positive but weak (P = 0.0540).

## S3.1 Direct-only northern-high optimizer audit

The frozen corrected Direct-only northern-high table marked shallow/open tube as optimizer_success=false even though the seven-response vector was FDR-supported. Because the historical vector_supported flag was based on q-value and did not itself require every optimizer flag to be true, we performed a dedicated numerical audit.

Audit sources:

- results/geography_20260924/h1_direct_northern_high_convergence_audit.json
- scripts/geography_correction/audit_h1_direct_convergence.py

Results:

- frozen shallow/open-tube estimate: 0.3544814;
- enhanced re-fit estimate: 0.3544784;
- absolute coefficient change: 2.93 × 10^-6;
- enhanced re-fit: converged;
- fully converged seven-response replay: q = 3.793 × 10^-8;
- fully converged six-response sensitivity excluding shallow/open tube: q = 1.430 × 10^-8.

An independent Python 3.11 multistart confirmation also passed. All multistart fits succeeded, the seven-response test gave P = 1.898 × 10^-8, the six-response sensitivity gave P = 7.154 × 10^-9, and the maximum slope deviation from the frozen solution was 2.85 × 10^-6.

The six-response result is a sensitivity only and does not replace the predeclared seven-response H1 estimand.

---

# Appendix S4. H2 conditional decomposition and raw floral patterns

H2 separates measured reproductive assurance from additional floral response.

The species-level selfing_core score uses compatibility, mating system and autonomous-selfing capacity only. The generalized_accessible score summarizes generalized floral form, actinomorphy and shallow/open tube.

The fitted model family includes:

selfing_core ~ isolation + area + climate

generalized_accessible ~ isolation + selfing_core + area + climate

plain_colour ~ isolation + selfing_core + area + climate

Persistence of an isolation coefficient after adjustment for selfing_core is conditional decomposition, not mediation.

### Table S3. Complete H2 conditional decomposition

Machine-readable table:

- submission/chapter1_current/supplement/Table_S3_H2_decomposition.csv

Selfing-adjusted accessibility estimates:

| Region | All-analysis β | SE | P | q | Direct-only β | SE | P | q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Northern mid-latitude | 0.02068 | 0.01281 | 0.1065 | 0.1703 | 0.01715 | 0.01144 | 0.1340 | 0.2143 |
| Northern high latitude | 0.12557 | 0.03390 | 0.000212 | 0.000847 | 0.11657 | 0.03716 | 0.001709 | 0.006835 |
| Tropical | 0.06846 | 0.02494 | 0.006061 | 0.01616 | 0.05249 | 0.02616 | 0.04483 | 0.1196 |
| Southern extratropical | 0.04929 | 0.03979 | 0.2154 | 0.2873 | 0.05779 | 0.03473 | 0.09611 | 0.1922 |

All eight estimates are positive. FDR support is region dependent. In the primary scope it is strongest in northern high latitudes and the tropics. Tropical Direct-only remains nominally positive but is not FDR-supported.

## S4.1 Raw colour and colour × architecture results

Raw colour, joint colour × architecture and architecture conditional on colour are retained in full rather than compressed into a very large static table.

All-analysis files:

- results/geography_20260924/all/raw_patterns/raw_colour_model_results.csv
- results/geography_20260924/all/raw_patterns/raw_colour_joint_omnibus.csv
- results/geography_20260924/all/raw_patterns/raw_colour_architecture_model_results.csv
- results/geography_20260924/all/raw_patterns/raw_colour_conditioned_architecture_model_results.csv

Direct-only files:

- results/geography_20260924/direct/raw_patterns/raw_colour_model_results.csv
- results/geography_20260924/direct/raw_patterns/raw_colour_joint_omnibus.csv
- results/geography_20260924/direct/raw_patterns/raw_colour_architecture_model_results.csv
- results/geography_20260924/direct/raw_patterns/raw_colour_conditioned_architecture_model_results.csv

These descriptive trait combinations do not identify realized pollinator identity.

---

# Appendix S5. H3 experimental pollen limitation

GloPL pollen limitation is the log response ratio of reproduction after supplemental pollen addition versus natural pollen receipt. Positive values indicate increased reproduction after added pollen.

The corrected analysis contains 2,969 experimental rows, 1,408 publication × coordinate × measurement cells, 1,248 unique sites and 919 publications. Each publication contributes total analysis weight one. Uncertainty is publication-cluster robust.

### Table S5. H3 pollen-limitation models

Machine-readable table:

- submission/chapter1_current/supplement/Table_S5_H3_pollen_limitation.csv

| Analysis | β distance | SE | Two-sided P | Cells | Publications | Sites |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Primary | 0.09191 | 0.03806 | 0.01575 | 1,408 | 919 | 1,248 |
| Supplemental-only | 0.04410 | 0.04351 | 0.31082 | 828 | 470 | 736 |
| No-zero-constant | 0.09089 | 0.03865 | 0.01869 | 1,375 | 912 | 1,238 |

The primary and no-zero-constant models support a positive isolation association. Supplemental-only remains positive but unsupported. H3 is an isolation-associated pollen-delivery constraint, not direct evidence of declining pollinator abundance or visitation.

Primary source files:

- results/geography_20260924/h3_original_corrected_comparison.json
- results/geography_20260924/h3_corrected_measurement_cells.csv.gz
- results/geography_20260924/h3_corrected_effect_rows.csv.gz

---

# Appendix S6. H4 exact-species functional triangulation

H4 exact-matches the literal H2 species scores to GloPL. Matching uses only case normalization and underscore/space normalization; synonym rescue and genus fallback are not used.

The model is:

pollen limitation ~ H2 score + corrected isolation + geographic context + measurement conditions

Each publication again contributes total weight one and uncertainty is publication-cluster robust.

### Table S6a. Exact H2-score bridge

Machine-readable table:

- submission/chapter1_current/supplement/Table_S6a_H4_scores.csv

| Trait family / score | Analysis | β trait score | SE | Two-sided P | Species | Publications |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Reproductive assurance / selfing_core | Primary | -0.29830 | 0.10352 | 0.00396 | 455 | 409 |
| Reproductive assurance / selfing_core | No-zero-constant | -0.30078 | 0.10419 | 0.00389 | 453 | 408 |
| Reproductive assurance / selfing_core | Supplemental-only | -0.09953 | 0.09398 | 0.28962 | 283 | 241 |
| Accessibility / generalized_accessible | Primary | -0.29566 | 0.12896 | 0.02187 | 143 | 143 |
| Accessibility / generalized_accessible | No-zero-constant | -0.30366 | 0.12991 | 0.01942 | 142 | 142 |
| Accessibility / generalized_accessible | Supplemental-only | -0.29950 | 0.14562 | 0.03971 | 101 | 98 |

Negative coefficients indicate that stronger expression of the island-associated H2 score is associated with lower current experimental pollen limitation.

### Table S6b. Atomic-trait H4 sensitivities

Machine-readable table:

- submission/chapter1_current/supplement/Table_S6b_H4_atomic.csv

Primary atomic estimates:

| Trait | β trait state | SE | Two-sided P | Species | Publications |
| --- | ---: | ---: | ---: | ---: | ---: |
| Self-compatibility | -0.11885 | 0.08246 | 0.1495 | 499 | 429 |
| Autonomous selfing | -0.44492 | 0.08020 | 2.89 × 10^-8 | 558 | 469 |
| Generalized form | -0.18448 | 0.09174 | 0.04434 | 246 | 246 |
| Actinomorphic symmetry | -0.38076 | 0.08697 | 1.20 × 10^-5 | 582 | 479 |

Selfing mating system and shallow/open tube did not pass their parent support gates for atomic H4 promotion.

H4 is explicitly post-hoc functional triangulation. It does not identify historical mediation.

---

# Appendix S7. Inferential boundaries, reproducibility and data availability

The submission distinguishes four inferential levels.

1. Assemblage pattern (H1): contemporary island-flora composition changes with isolation.
2. Conditional decomposition (H2): measured reproductive assurance does not statistically absorb all floral-accessibility response. This is not mediation.
3. Independent ecological pressure (H3): experimental pollen limitation increases with isolation. This is not a direct measure of pollinator abundance or visitation.
4. Functional compatibility (H4): reproductive-assurance and accessibility states associated with isolation are associated with lower current pollen limitation in exact-species overlap. This is post-hoc functional triangulation.

The unobserved historical edge is:

past isolation-associated pollen limitation → selection / sorting / persistence → present trait composition

The current analysis does not distinguish species sorting, differential colonization/persistence and within-lineage evolutionary change.

A separate prospective H4 validation effort stopped at the support gate before outcome unblinding because the prespecified minimum sample/publication requirements were not met. It is a design/support result, not a biological null.

## S7.1 Deterministic supplementary tables

The generated table manifest is:

- submission/chapter1_current/supplement/SUPPLEMENT_TABLES_MANIFEST.json

It records source and output SHA-256 hashes and supports byte-for-byte table rebuilds through:

- scripts/submission/build_chapter1_supplement_tables.py

Generated tables:

- Table_S1_data_summary.csv
- Table_S2a_H1_atomic.csv
- Table_S2b_H1_joint.csv
- Table_S3_H2_decomposition.csv
- Table_S5_H3_pollen_limitation.csv
- Table_S6a_H4_scores.csv
- Table_S6b_H4_atomic.csv

Table S4 is intentionally represented by the complete corrected raw-pattern file family rather than by one oversized duplicated table.

## S7.2 Submission data-policy boundary

The rights-filtered trait derivative is publicly archived at DOI 10.5281/zenodo.22704973. The complete 222,688-cell scientific analysis ledger is not yet wholly redistributable.

Ecology Letters submission remains blocked until either:

1. redistribution rights for the analysis-used trait cells are closed; or
2. the editors explicitly approve a legal/licensing exception and reviewer-access/reconstruction plan.

See:

- submission/chapter1_current/ECOLOGY_LETTERS_DATA_GATE.md
- submission/chapter1_current/DATA_ACCESSIBILITY_DRAFT.md
- submission/chapter1_current/ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md

## S7.3 Reproducibility receipts

- results/geography_20260924/repository_replay_verification.json
- results/geography_20260924/refit_completion.json
- results/geography_20260924/raw_pattern_refit_completion.json
- results/geography_20260924/h1_direct_northern_high_convergence_audit.json

This Supplementary Information is descriptive of the current corrected scientific surface; it does not replace permanent data/code archiving.
