# Supplementary Information — GEB research article (anonymous review)

## Island isolation is associated with recurrent reproductive assurance but regionally contingent floral change

This Supplementary Information accompanies the GEB-first-shot research article and the unchanged corrected scientific result selector config/chapter1_submission_current.json. Appendices S1–S7 retain the frozen primary results. Appendix S8 documents the new, explicitly post-hoc paired-trait extension.

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

### S2.1 Interpretation of GBIF-observed island floras

GBIF occurrence records document reported presences and are not standardized equal-effort species inventories, measures of plant abundance or reliable evidence that an unreported species is absent. The island-by-species layer therefore describes **recorded assemblages**, not complete island-wide species lists. Islands with no usable records were not converted to trait zeros.

Coverage can vary with access to islands and collecting institutions, research history, digitization, survey season and taxonomic interest. Non-random collecting practices can also be related to the traits of plants that observers collect or identify (Meyer et al. 2016; Daru et al. 2018). Thus the probability that a species enters the observed flora may depend on island isolation and the traits being analysed. On small islands, species-identification or coordinate uncertainty may affect assignment to a particular island polygon; this is a distinct problem from the GSHHG coastline-distance repair.

The present H1 regressions include island area, climate and spatial-block robust uncertainty, but **do not estimate occurrence-detection probabilities, equalize survey effort or establish comparable sampling completeness**. These controls do not remove effort-related observation bias. Neither the size nor sign of residual bias is estimable from the occurrence compilation alone. Independent floristic checklists or standardized repeat surveys would be required to assess coverage and detection directly.

### S2.2 Trait evidence, missingness and floristic origin


The scientific trait database contains 106,295 accepted angiosperm species and three raw evidence axes:

1. flower colour;
2. floral structural complexity;
3. reproductive assurance.

Of 318,885 possible species × axis cells, 222,688 are resolved (69.83%). Species-direct High/Medium evidence forms the Direct-only sensitivity. The primary all-analysis scope additionally admits validated lower-confidence evidence where direct evidence is unavailable. Missing trait information remains missing rather than becoming trait absence.

The full scientific ledger is identified by SHA-256:

a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a

The rights-filtered public derivative contains 46,274 redistribution-authorized cells and is archived at Zenodo:

[Public derivative DOI withheld pending anonymous-access decision]

This public derivative is not the complete scientific analysis ledger. Omitted cells are omitted because of redistribution-rights status, not biological missingness or exclusion from the analysis.

Supporting documentation:

- docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md
- docs/CHAPTER1_DATABASE_PUBLIC_SUBSET.md
- config/chapter1_database_versions/v1.0.0.yml
- Data-rights plan (author-only; not supplied with blinded materials)

---

# Appendix S3. Final traitwise t inference without multiplicity correction

### H1: seven separate trait responses

The primary analysis is broad contemporary island flora with All evidence. Each of seven binary traits is fitted separately in each of four regions by beta-binomial logit regression, adjusting for standardized corrected log isolation, log island area and climate PC1–4. Reproduction, colour and structure are labels, not aggregate scores. Spatial-block sandwich standard errors use a finite-cluster t reference with G−1 degrees of freedom and pointwise 95% intervals. All individual p values are two-sided and unadjusted; no Holm correction is applied. The final rule was user-approved after comparison of alternatives on 4 October 2026 and is retrospective, not preregistered. Individual significance does not establish family-wise significance or formal differences between regional slopes.

WCVP regional-native-compatible flora is the sole active origin sensitivity; Direct-only is a separate evidence-quality sensitivity. Both use the identical model and inference rule. WCVP compatibility retains existing source-native records and upgrades unresolved records only under accepted TDWG-L3 native-range compatibility; introduced records are not overwritten. This establishes regional compatibility rather than exact focal-island nativity. H2–H4 and their previously defined functional covariates are unchanged.

### Table S2. Final traitwise H1 results

The active machine-readable H1 table is:

- submission/chapter1_current/supplement/Table_S2_H1_traitwise.csv

It contains all 112 final fits (All/Direct-only × broad/WCVP × four regions × seven traits), their finite-cluster t intervals and pointwise two-sided P values. All 112 fits converged. Earlier joint/vector H1 tables are retained only as historical provenance and are not active submission inference.

### H1: recurrent reproductive responses with regional floral differences

All 112 fits converged. Primary All evidence supports 17 of 28 individual associations at nominal two-sided P < 0.05, the same set of significant traits as the original poster's normal approximation. Self-compatibility increases with isolation in all four regions (P = 0.00569, 0.00671, 0.01812 and 0.01123 in northern mid-latitude, northern high-latitude, tropical and southern extratropical floras, respectively). Selfing mating system increases in the first three regions. Autonomous selfing increases in northern high latitudes and the tropics. Generalized form increases in northern mid-latitudes, the tropics and southern extratropics; actinomorphy increases in northern high latitudes and the tropics. Shallow/open tubes increase in northern high latitudes but decrease in southern extratropics. Plain colour increases in southern extratropics.

In WCVP All, 10 of 28 individual associations meet the same nominal threshold. Increasing self-compatibility remains supported in all four regions. Northern-high autonomous selfing, actinomorphy and shallow/open tubes remain positive and supported; tropical actinomorphy remains supported. Plain colour increases in both tropical and southern extratropical WCVP floras. Several broad-flora associations, including the southern shallow/open-tube decrease, no longer meet the threshold. This can reflect changed composition, precision and coverage, rather than proving an introduced-species mechanism. Direct-only yields 15/28 supported associations in broad flora and 11/28 in WCVP flora. Complete coefficients, intervals and unadjusted p values are reported in results/h1_final_traitwise_t_20261004/traitwise_results.csv.

# Appendix S4. H2 conditional decomposition and raw floral patterns

H2 separates measured reproductive assurance from additional floral response.

The species-level selfing_core score uses compatibility, mating system and autonomous-selfing capacity only. The generalized_accessible score summarizes generalized floral form, actinomorphy and shallow/open tube.

The fitted model family includes:

selfing_core ~ isolation + area + climate

generalized_accessible ~ isolation + selfing_core + area + climate

plain_colour ~ isolation + selfing_core + area + climate

Persistence of an isolation coefficient after adjustment for selfing_core is conditional decomposition, not mediation.

### Table S3. Complete H2 conditional decomposition

Active machine-readable finite-cluster tables:

- submission/chapter1_current/supplement/Table_S3a_H2_all_finite_cluster.csv
- submission/chapter1_current/supplement/Table_S3b_H2_direct_finite_cluster.csv

Selfing-adjusted accessibility estimates use t(G−1) P values and BH q values from the final finite-cluster audit:

| Region | All-analysis β | SE | P | q | Direct-only β | SE | P | q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Northern mid-latitude | 0.02068 | 0.01281 | 0.1103 | 0.1764 | 0.01715 | 0.01144 | 0.1378 | 0.2205 |
| Northern high latitude | 0.12557 | 0.03390 | 0.000474 | 0.00243 | 0.11657 | 0.03716 | 0.002744 | 0.01097 |
| Tropical | 0.06846 | 0.02494 | 0.007174 | 0.01913 | 0.05249 | 0.02616 | 0.04751 | 0.1267 |
| Southern extratropical | 0.04929 | 0.03979 | 0.2232 | 0.2977 | 0.05779 | 0.03473 | 0.1046 | 0.2091 |

All eight estimates are positive. FDR support is region dependent. In the primary scope it is strongest in northern high latitudes and the tropics. Tropical Direct-only remains nominally positive but is not FDR-supported.

## S4.1 Complete-three-component selfing sensitivity

To test whether incomplete measurement of `selfing_core` creates the residual accessibility association, a stricter Direct-only mediator was rebuilt using only species with all three reproductive components observed and informative: self-incompatibility/compatibility, mating system and autonomous-selfing capacity. The complete score was available for 564 species and yielded island-level scores on 2,610 islands before covariate filtering.

Machine-readable results:

- results/h2_complete_selfing_corrected_20261003/h2_complete_selfing_sensitivity.csv
- results/h2_complete_selfing_corrected_20261003/h2_complete_selfing_sensitivity_summary.json
- results/h2_complete_selfing_corrected_20261003/baseline_reconstruction_gate.json

The baseline Direct-only score was independently reconstructed to machine precision before applying the complete-three-component restriction.

| Region | Complete-score islands | β isolation | SE | P | q |
| --- | ---: | ---: | ---: | ---: | ---: |
| Northern mid-latitude | 1,725 | 0.02296 | 0.01107 | 0.03812 | 0.05082 |
| Northern high latitude | 221 | 0.08508 | 0.03433 | 0.01321 | 0.03759 |
| Tropical | 438 | 0.05058 | 0.02153 | 0.01880 | 0.03759 |
| Southern extratropical | 157 | 0.03487 | 0.03165 | 0.27060 | 0.27060 |

All four coefficients remain positive. FDR support persists in northern high latitudes and the tropics, and northern mid-latitudes lie immediately above the 0.05 FDR threshold. This weakens a missing-component explanation for H2 but is not a formal errors-in-variables correction.

## S4.2 Raw colour and colour × architecture results

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
| Primary | 0.09191 | 0.03806 | 0.01594 | 1,408 | 919 | 1,248 |
| Supplemental-only | 0.04410 | 0.04351 | 0.31134 | 828 | 470 | 736 |
| No-zero-constant | 0.09089 | 0.03865 | 0.01891 | 1,375 | 912 | 1,238 |

The primary and no-zero-constant models support a positive isolation association. Supplemental-only remains positive but unsupported. H3 tests current pollen limitation, not a measured decline in pollinator abundance or visitation.

### Table S5b. Offshore-gradient robustness

Because 996 of 1,248 GloPL sites are true continental zero-distance sites, we tested whether H3 reduces to a mainland/offshore step. The P values below use finite-publication t inference; older asymptotic P values (0.02319 and 0.14364) are retained in the machine-readable audit. Offshore means positive corrected coast distance, not independently confirmed membership in an island polygon. Offshore distance was re-standardized within the restricted sample; its β magnitude is not directly comparable with the full-sample β.

| Analysis | β | SE | Two-sided P | Cells | Publications | Sites |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Offshore continuous gradient | 0.22031 | 0.09704 | 0.02459 | 276 | 153 | 252 |
| Mainland vs offshore indicator | 0.12743 | 0.08714 | 0.14398 | 1,408 | 919 | 1,248 |

A leave-one-publication-out jackknife of the offshore gradient retained positive estimates in all 153 deletions. The minimum coefficient was 0.17850, the maximum was 0.26930 and the weakest asymptotic two-sided P value was 0.04892.

Machine-readable results:

- results/h3_offshore_gradient_20261003/h3_offshore_gradient_summary.json
- results/h3_offshore_gradient_20261003/h3_offshore_leave_one_publication.csv

These are post-hoc robustness analyses. They show that the positive H3 association persists within offshore sites and is not captured by a binary mainland/offshore contrast, but they do not establish causation.

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
| Reproductive assurance / selfing_core | Primary | -0.29830 | 0.10352 | 0.00417 | 455 | 409 |
| Reproductive assurance / selfing_core | No-zero-constant | -0.30078 | 0.10419 | 0.00410 | 453 | 408 |
| Reproductive assurance / selfing_core | Supplemental-only | -0.09953 | 0.09398 | 0.29068 | 283 | 241 |
| Accessibility / generalized_accessible | Primary | -0.29566 | 0.12896 | 0.02334 | 143 | 143 |
| Accessibility / generalized_accessible | No-zero-constant | -0.30366 | 0.12991 | 0.02082 | 142 | 142 |
| Accessibility / generalized_accessible | Supplemental-only | -0.29950 | 0.14562 | 0.04239 | 101 | 98 |

Negative coefficients indicate that stronger expression of the island-associated H2 score is associated with lower current experimental pollen limitation.

### Table S6b. Atomic-trait H4 sensitivities

Primary atomic estimates are reported with finite-publication t inference:

| Trait | β trait state | SE | Two-sided P | Species | Publications |
| --- | ---: | ---: | ---: | ---: | ---: |
| Self-compatibility | -0.11885 | 0.08246 | 0.1502 | 499 | 429 |
| Autonomous selfing | -0.44492 | 0.08020 | 4.85 × 10^-8 | 558 | 469 |
| Generalized form | -0.18448 | 0.09174 | 0.04543 | 246 | 246 |
| Actinomorphic symmetry | -0.38076 | 0.08697 | 1.47 × 10^-5 | 582 | 479 |

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

The two earlier predeclared GloPL distance-by-trait moderation families were also replayed after the 24 September geography correction, using the exact archived matched effect rows, trait states, support decisions, models and sensitivities. Before distance replacement, the replay reproduced the frozen estimates to a maximum absolute difference of 3.93 × 10^-14. Correcting the exposure did not rescue either buffering family.

| Frozen family / trait | Corrected distance × trait interaction | SE | Two-sided P | One-sided buffering P | Frozen-rule support |
| --- | ---: | ---: | ---: | ---: | --- |
| Reproductive assurance: self-compatibility | +0.04444 | 0.07414 | 0.5489 | 0.7256 | No |
| Reproductive assurance: autonomous selfing | -0.08329 | 0.07237 | 0.2498 | 0.1249 | No |
| Floral architecture: generalized form | -0.05204 | 0.09589 | 0.5873 | 0.2937 | No |
| Floral architecture: actinomorphy | +0.03413 | 0.07961 | 0.6682 | 0.6659 | No |

Selfing mating system and shallow/open tube remain support-limited under their frozen preflight gates. Autonomous selfing retains the predicted negative interaction in both frozen measurement sensitivities, but it remains statistically unsupported. Generalized form is negative in the primary and no-zero-constant fits but reverses sign in the supplemental-only sensitivity. Therefore H4 should be interpreted as an association between island-enriched states and lower average current pollen limitation, not evidence that those states flatten the isolation-associated pollen-limitation gradient.

Reproducibility surface:

- results/geography_20260924/corrected_trait_moderation_replay_20261003/
- validation workflow run 37094533233; artifact 11263028661

## S7.1 Deterministic supplementary tables

The generated table manifest is:

- submission/chapter1_current/supplement/SUPPLEMENT_TABLES_MANIFEST.json

It records source and output SHA-256 hashes and supports byte-for-byte table rebuilds through:

- scripts/submission/build_chapter1_supplement_tables.py

Generated active tables:

- Table_S1_data_summary.csv
- Table_S2_H1_traitwise.csv
- Table_S3a_H2_all_finite_cluster.csv
- Table_S3b_H2_direct_finite_cluster.csv
- Table_S5_H3_pollen_limitation.csv
- Table_S6a_H4_scores.csv

Atomic H4 sensitivities are reported in the SI text from the corrected atomic source table; earlier H1 joint/vector and domain tables remain historical provenance only.

Table S4 is intentionally represented by the complete corrected raw-pattern file family rather than by one oversized duplicated table.

## S7.2 Submission data-policy boundary

The rights-filtered trait derivative is publicly archived at [Public derivative DOI withheld pending anonymous-access decision]. The complete 222,688-cell scientific analysis ledger is not yet wholly redistributable.

Data and code access must be reviewed with the GEB editor. The third-party legal/reuse permissions are not automatically waived; permissible paths require either:

1. redistribution rights for the analysis-used trait cells are closed; or
2. the editors explicitly approve a legal/licensing exception and reviewer-access/reconstruction plan that is lawful for the affected source records.

See:

- submission/geb_current/DATA_CODE_AVAILABILITY.md
- config/chapter1_geb_submission.json

## S7.3 Reproducibility receipts

- results/geography_20260924/repository_replay_verification.json
- results/geography_20260924/refit_completion.json
- results/geography_20260924/raw_pattern_refit_completion.json
- results/geography_20260924/h1_direct_northern_high_convergence_audit.json

This Supplementary Information is descriptive of the current corrected scientific surface; it does not replace permanent data/code archiving.


---

# Appendix S8. Post-hoc depth-dependent floral symmetry response

This appendix accompanies the post-hoc H1–H2 extension; the original analyses in Appendices S1–S7 are unchanged. All outputs below derive from the own GBIF occurrence-based island flora and own species-Direct trait evidence, with the pinned corrected geography. No official GIFT island checklist is used. The separate island–candidate-mainland descriptive comparison uses official GIFT mainland donor lists only. Original data rights restrictions remain in force.

## S8.1 Estimand and support

Among species with both unambiguous recorded radial/bilateral symmetry and shallow/open or deep tube states, define for island i:

D_i = p(radial | deep, i) − p(radial | shallow/open, i).

Exclude intermediate depth, asymmetric, unresolved and multistate records; require at least five species in each depth group and at least ten paired carriers. A companion marginal contrast is p(radial) − p(shallow/open) on a common carrier set. Regress these proportions against the raw H1 log-distance transform with area, four climate PCs, log carrier count and documented coverage, using equal island weights and finite spatial-cluster t inference. Require at least 50 islands, 10 blocks and full rank. The eight regional contrast slots include failures in Holm adjustment. The separate tropical-minus-northern-midlatitude contrast family uses matched 5–95% exposure support; it is not a robust demonstrated region interaction. The H1 primary traitwise beta-binomial results are unchanged. High-only, literal-native and n≥30 variants are nested sensitivities.

### Table S8a. All regional and cross-region pair contrasts

Machine-readable full 60 slots: `Table_S8a_flower_pair_isolation.csv`; detailed individual cell slopes in original anonymous receipt [review provenance available via approved anonymous access]. The following is the broad Direct conditional n≥10 subset, not all fit slots.

| Region | Islands | Blocks | Slope | 95% spatial interval | Holm p | Joint sensitivity |
|---|---:|---:|---:|---|---:|---|
| northern_midlatitude | 604 | 53 | 0.00198 | [-0.01007, 0.01404] | 1.000 | [-0.01772, 0.02087] |
| northern_high_latitude | 6 | 2 | — | — | not fitted | — |
| tropical | 257 | 66 | 0.02762 | [0.01539, 0.03985] | 0.0002206 | [0.00510, 0.04943] |
| southern_extratropical | 42 | 17 | — | — | not fitted | — |
| tropical_minus_north_mid | 570 | 79 | 0.03603 | [0.00926, 0.06281] | 0.01798 | [-0.00846, 0.08188] |

Unsupported slots retain p=1 for the planned adjustment but are not reported as tested null effects. Joint ranges use 499 coherent positive genus weights and geographical block resampling, not calibrated confidence intervals. Results are outcome-motivated, not independent confirmations. The original source artifact is SHA-pinned, not newly redistributed.

## S8.2 Evidence-grade and matched-island diagnosis

The broad tropical conditional slope was 0.02762 (257 islands); on the 62 High-only eligible islands the broad estimate was 0.00970. Changing to High-only on the same islands and broad covariates gave 0.00147. The paired High-minus-broad slope difference was -0.00822, interval [-0.02564, 0.00920], P=0.3423. Original High-only covariates give -0.00072. These are distinct species rosters, not repeat measurements or evidence of grade equivalence. The tropical literal-native conditional n=44 fails support. Original WCVP-compatible replay artifacts were unavailable during the extension, so no substitute or claim of native validation is made.

Table S8d, all grade diagnostic estimates: `Table_S8d_grade_matched_diagnostic.csv`. Signed change identity and source receipt: grade result [review provenance available via approved anonymous access].

## S8.3 Measured reproductive assurance

Use the original species-Direct H2 state audit and equal-weight selfing_core definition (at least two informative components); complete-three-component and separate SC and autonomous/delayed capacity sensitivities. The reproductive mean uses either all scored flora or the reproductive-known subset of the morphology carriers. Require at least ten scored reproductive species and the original fit support. Baseline and adjusted slopes use identical islands within each row; the rows are not mutually independent. Joint sensitivity reweights both response and reproductive mean with the same genus multipliers; nuisance covariates and cohort remain fixed. Holm over eight score × roster slots per scope.

### Table S8b. Reproductive conditioning

| Roster / score | Islands | Matched baseline | Adjusted slope | Holm p | Joint sensitivity (adjusted) |
|---|---:|---:|---:|---:|---|
| whole_flora / min2 | 244 | 0.02796 | 0.01967 | 0.003613 | [-0.00123, 0.04565] |
| whole_flora / complete3 | 72 | 0.01193 | 0.01213 | 0.1220 | [-0.01357, 0.03883] |
| whole_flora / SC | 257 | 0.02762 | 0.02356 | 0.001115 | [0.00192, 0.04727] |
| whole_flora / autonomy | 257 | 0.02762 | 0.02043 | 0.003613 | [0.00019, 0.04545] |
| pair_carriers / min2 | 62 | 0.00689 | 0.00561 | 0.1979 | [-0.01725, 0.02453] |
| pair_carriers / complete3 | 2 | — | — | not fitted | — |
| pair_carriers / SC | 188 | 0.02889 | 0.02463 | 0.00004177 | [0.00273, 0.04588] |
| pair_carriers / autonomy | 106 | 0.02226 | 0.02089 | 0.00005184 | [0.00087, 0.04294] |

Full 24 slots including unsupported and High-only/native sensitivities: `Table_S8b_reproductive_conditioning.csv`; original H2 [review provenance available via approved anonymous access] and joint weights [review provenance available via approved anonymous access]. The four atomic adjustment ranges stay positive; the whole-flora min2 joint range crosses zero. No equivalence, causal independence, mediation fraction or independent replication inferred.

## S8.4 Depth-specific decomposition

On the exact same mask and covariates, the difference of the two component slopes equals the slope of D_i to numerical precision. The two additional unconditioned H1 component tests form a Holm2 family. H2 component estimates are descriptive. All four fitted difference slopes reproduce their parent outputs.

### Table S8c. Depth-specific slopes

| Roster / score | Islands | Radial among deep | Radial among shallow/open | Difference |
|---|---:|---:|---:|---:|
| none / none | 257 | 0.01988 | -0.00774 | 0.02762 |
| whole_flora / min2 | 244 | 0.01516 | -0.00451 | 0.01967 |
| pair_carriers / SC | 188 | 0.01837 | -0.00626 | 0.02463 |
| pair_carriers / autonomy | 106 | 0.01382 | -0.00707 | 0.02089 |

All 12 component/contrast rows with intervals: `Table_S8c_depth_specific_components.csv`. The deep component is the better-supported contributor. The shallow/open negative unconditioned slope is nominally and Holm2 supported in the spatial model but its joint-weight range crosses zero. The whole-flora min2 and pair-carrier SC deep component remain positive in the joint sensitivity; autonomy-adjusted deep component does not, although the contrast stays positive.

**Figure S8. Depth-dependent symmetry response.** Figure S8. A: radial prevalence slopes in the two tube-depth groups. B: paired difference. Solid intervals denote 95% spatial-cluster intervals; dashed ranges denote 2.5–97.5% joint genus/geography sensitivity, not confidence intervals. Each row has its own fixed cohort, shared by both components. Slopes are absolute proportions per unit raw H1 log-isolation transform, not the standardized original H1 beta-binomial slopes.

## S8.5 H4 joint architecture support

The exact-species GloPL morphology bridge was assessed before outcome fitting. Require 50 total species, 30 publications, and at least ten distinct species and five publications per each of the four morphology cells. All four planned panels fail, with 32 global and eight tropical matched species before selfing adjustment. No effect estimated.

Table S8e: all 16 cell-support rows: `Table_S8e_GloPL_joint_support.csv`; verified support receipt [review provenance available via approved anonymous access]. Original H4 remains the exact-score functional bridge, not evidence of a depth-dependent pollen-limitation mechanism. H3/H4 share GloPL experiments.

## S8.6 Reproducibility and claim ceiling

Each receipt records the original successful run/head, artifact ID, reported ZIP SHA256 and anonymous stdout results. Aggregate table assembly checks original receipt identity and numerical schema; no new biological model is fitted. The supplemental results are exploratory, scope-dependent and do not identify evolutionary change, a particular pollinator, approach mechanics, nectar accessibility or historical mediation. All H1–H4 original submission files and selectors are unchanged.
