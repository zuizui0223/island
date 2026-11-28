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

The rights-filtered public derivative contains 46,274 redistribution-authorized cells and is archived separately. The public derivative DOI is withheld pending the anonymous-access decision.

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

The rights-filtered trait derivative is publicly archived; the DOI is withheld pending the anonymous-access decision. The complete 222,688-cell scientific analysis ledger is not yet wholly redistributable.

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

Machine-readable full 60 slots: `Table_S8a_flower_pair_isolation.csv`; detailed individual cell slopes in original anonymous receipt [anonymous reviewer access pending editorial approval]. The following is the broad Direct conditional n≥10 subset, not all fit slots.

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

Table S8d, all grade diagnostic estimates: `Table_S8d_grade_matched_diagnostic.csv`. Signed change identity and source receipt: grade result [anonymous reviewer access pending editorial approval].

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

Full 24 slots including unsupported and High-only/native sensitivities: `Table_S8b_reproductive_conditioning.csv`; original H2 [anonymous reviewer access pending editorial approval] and joint weights [anonymous reviewer access pending editorial approval]. The four atomic adjustment ranges stay positive; the whole-flora min2 joint range crosses zero. No equivalence, causal independence, mediation fraction or independent replication inferred.

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

**Figure S8. Depth-dependent symmetry response.** (a) radial prevalence slopes in the two tube-depth groups. (b) Paired difference. Solid intervals denote 95% spatial-cluster intervals; dashed ranges denote 2.5–97.5% joint genus/geography sensitivity, not confidence intervals. Each row has its own fixed cohort, shared by both components. Slopes are absolute proportions per unit raw H1 log-isolation transform, not the standardized original H1 beta-binomial slopes.

## S8.5 H4 joint architecture support

The exact-species GloPL morphology bridge was assessed before outcome fitting. Require 50 total species, 30 publications, and at least ten distinct species and five publications per each of the four morphology cells. All four planned panels fail, with 32 global and eight tropical matched species before selfing adjustment. No effect estimated.

Table S8e: all 16 cell-support rows: `Table_S8e_GloPL_joint_support.csv`; verified support receipt [anonymous reviewer access pending editorial approval]. Original H4 remains the exact-score functional bridge, not evidence of a depth-dependent pollen-limitation mechanism. H3/H4 share GloPL experiments.

## S8.6 Reproducibility and claim ceiling

Each receipt records the original successful run/head, artifact ID, reported ZIP SHA256 and anonymous stdout results. Aggregate table assembly checks original receipt identity and numerical schema; no new biological model is fitted. The supplemental results are exploratory, scope-dependent and do not identify evolutionary change, a particular pollinator, approach mechanics, nectar accessibility or historical mediation. All H1–H4 original submission files and selectors are unchanged.

## S8.7 Local genus-margin diagnostic

On the original tropical paired-carrier masks, write R for the number of radial species, D and S for the numbers of deep and shallow/open species, and X for the radial-deep count. The observed response is X/D − (R−X)/S. Within each island and genus, the fixed-margin expectation of X is the product of the radial and deep counts divided by the number of paired carriers in that genus. Sum over genera and substitute this expectation into the conditional contrast. Repeating within genus and joint evidence grades provides a second expectation. The observed slope on the original covariate design is the exact sum of the genus expectation, incremental grade term and remaining association. These expectations use observed local trait margins and therefore are not held-out predictions.

A planned species-coherent null would permute symmetry once per species within genus and joint grades, with depth and the island occurrence matrix fixed, recomputing the local expectation for each of 1,999 draws. No stratum in the primary cohort contained variation in both symmetry and depth, so the test was not estimable. The local expected slopes equalled the observed slopes in every fitted cohort. This is an identifiability limitation, not evidence of absence of a within-lineage effect or proof that genus turnover caused the observed pattern. Zero incremental grade component does not exclude documentation bias.

### Table S8f. Local genus-margin accounting

| Cohort | Islands | Blocks | Observed slope | Genus expectation | Within-stratum residual | Exchangeable strata |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Broad n≥10 | 257 | 66 | 0.02762 | 0.02762 | 0 | 0 |
| Broad n≥30 | 122 | 49 | 0.02345 | 0.02345 | 0 | 0 |
| High-only | 62 | 30 | -0.00072 | -0.00072 | 0 | 0 |
| Literal native | 44 | 17 | — | — | — | — |

All four slots, including grade increment and null support, are retained in Table_S8f_local_genus_margins.csv. There is no permutation P value for any fitted slot. The primary broad cohort contains 718 species in 371 genera; 14 islands contain a genus–grade stratum with both depth classes, but there is no joint-variable exchangeable stratum. The literal-native cohort fails the original island-count threshold.

## S8.8 Geographic transport of paired floral composition

In each held-out spatial block, exclude all accepted paired-trait species recorded anywhere in that block from training, even if the same species occurs in other blocks. Primary donors are documented non-tropical floras; other tropical blocks form a nested sensitivity. Use unique donor species, rather than occurrence counts. Require at least 100 donor species, ten donor blocks and both depth categories. For each genus, estimate four joint category frequencies with smoothing of five pseudo-observations toward the donor-only Jeffreys-smoothed prior (0.5 per cell). Unseen genera use the donor prior. This analysis supplies observed paired-carrier genera at prediction time and uses the original trait-conditioned eligibility mask. It is not a geography-only forecast, a prediction of genus entry or establishment, or a new sample independent of the original database. No first-source-region embargo is available in this pinned geography.

Compare three mixtures: donor-only blind prior, equally weighted represented genera, and represented genera weighted by their observed paired species counts. The observed category vector q and predicted vector p give squared category error ||q−p||²; differences in this error equal differences in multiclass Brier loss on the same carriers. Also report log loss and squared error in the conditional contrast. Average paired gains within geographic blocks, then equally across blocks. The 499 fixed-fold resampling ranges describe sensitivity to block weighting, not calibrated confidence intervals, because folds and taxa overlap. Species- and genus-equal gains use per-species mean loss across held-out occurrences. Fit predicted conditional contrasts to the same H1 covariate design; the slope differences telescope exactly. No new significance family is assigned.

### Table S8g. Held-out composition and isolation response

| Cohort | Donors | Scored / planned islands | Brier gain weighted vs blind | Observed slope | Weighted prediction slope |
| --- | --- | ---: | ---: | ---: | ---: |
| Broad n≥10 | Non-tropical | 257 / 257 | 0.05087 | 0.02762 | 0.00397 |
| Broad n≥10 | Other tropical | 257 / 257 | 0.04232 | 0.02762 | 0.00371 |
| Broad n≥30 | Non-tropical | 122 / 122 | 0.04696 | 0.02345 | 0.00262 |
| Broad n≥30 | Other tropical | 122 / 122 | 0.03875 | 0.02345 | 0.00105 |
| High-only | Non-tropical | 62 / 62 | 0.01828 | -0.00072 | -0.00234 |
| High-only | Other tropical | 31 / 62 | — | — | — |
| Literal native | Non-tropical | — / 44 | — | — | — |
| Literal native | Other tropical | — / 44 | — | — | — |

The primary non-tropical donor panel yields conditional-contrast RMSE of 0.25145, 0.21317 and 0.15992 for the blind, equal-genus and species-count-weighted mixtures respectively. Block-equal Brier gains are 0.01905 for equal-genus versus blind and 0.03182 for weighted versus equal-genus; corresponding fixed-fold sensitivity ranges are [0.01587, 0.02256] and [0.02675, 0.03752]. Both increments are positive under species- and genus-equal weighting for non-tropical donors. Under other tropical donors, genus-equal weighted-versus-equal increments are negative in both broad cohorts. The High-only other-tropical comparison is not completed because six folds fail donor support; no reduced-support result is substituted. Both native slots fail the original parent gate. Table_S8g_genus_transport.csv retains all eight panels, score comparisons, resampling ranges, taxon-equal sensitivities and donor support. Composition skill is not successful transfer of the isolation gradient.

## S8.9 Evaluation-only donor support decomposition

Split each island's observed and predicted four-cell proportions into contributions from carriers with at least one eligible donor congener and those with none. Denote predicted contributions Pk and Pu, and observed contributions Qk and Qu. For the nonlinear conditional contrast f, the supported contribution to residual f(Qk+Qu)−f(Pk+Pu) is the average of f(Qk+Pu)−f(Pk+Pu) and f(Qk+Qu)−f(Pk+Qu). The unsupported contribution is the analogous average over the two orders. These two contributions sum exactly to the residual; the same identity holds for their regression slopes on the fixed design. Observed target categories are used solely after prediction in this two-factor Shapley accounting. For the category error, project each signed vector difference onto their sum; the two signed projections sum to squared error, but are not independent or causal fractions. All archived original forecast fields were reproduced within absolute tolerance 10^−10 before the decomposition.

### Table S8h. Donor support and residual isolation slope

| Cohort | Donors | Residual slope | Supported genera | Unsupported genera | Unsupported carrier mass |
| --- | --- | ---: | ---: | ---: | ---: |
| Broad n≥10 | Non-tropical | 0.02365 | 0.01298 | 0.01066 | 0.53476 |
| Broad n≥10 | Other tropical | 0.02391 | 0.01383 | 0.01008 | 0.47052 |
| Broad n≥30 | Non-tropical | 0.02083 | 0.00874 | 0.01209 | 0.54862 |
| Broad n≥30 | Other tropical | 0.02241 | 0.00858 | 0.01383 | 0.47821 |
| High-only | Non-tropical | 0.00162 | 0.00187 | -0.00025 | 0.70116 |

Full eight-slot accounting, including unsupported panels, is in Table_S8h_genus_transport_coverage.csv. In the primary non-tropical panel, mean carrier proportions with zero, one, two to four, or at least five donor species are 53.48%, 14.55%, 13.26% and 18.71%. Both residual components are positive in the broad cohorts, but their relative sizes depend on carrier threshold. The unsupported component is larger at n≥30; no universally dominant source of slope error is established. Supported genera can be sparsely sampled; shrinkage and documentation remain explanations of prediction error. Neither component measures evolutionary change or the effect of collecting new data. No new statistical test is assigned to the decomposition. All results remain exploratory same-database diagnostics, and no revised H4 joint morphology model was fitted.

---

# Appendix S9. Reproductive strategies and source-relative floral composition

The following archived analyses use the project-assembled GBIF occurrence-derived target flora (4,505 parent islands) and integrated species-level Direct evidence. The historical identifier GIFT4505 refers to this target roster, not official GIFT island checklists. Official GIFT native mainland floras (748 candidate regions in the pinned donor input) serve as comparison pools (Weigelt et al. 2020). They are not demonstrated ancestral sources. Geographic nearest-five (k5) and geography-filtered climate-nearest-ten definitions, pooled versus first-ranked pools, inclusive versus focal-name-disjoint rosters and origin/grade sensitivities are retained only where present in the original design. Different cohorts and weights do not constitute independent replications. All comparisons are exploratory unless already designated part of the frozen H1–H4. The exact source and software provenance inventory is maintained in the author-only audit; a lawful anonymous access route is pending.

## S9.1 Reproductive strategies on common carriers

Use high/medium species-Direct unambiguous compatibility, mating system and autonomous/delayed selfing states. For the common-carrier SC versus mating-system test, retain the same species with both known reproductive states and colour, form and symmetry recorded. Each island and its nearest-five potential mainland pool is sampled at 20 species, without replacement, in 20 replicates. Each draw supplies both reproductive fractions; subtract source from island and then compare the two differences. Geographic blocks and first-rank source signatures are alternative equal-weight sensitivity groupings.

### Table S9a. Northern-midlatitude same-carrier reproductive contrasts

| Origin scope | Islands | SC | Mating-system selfing | Selfing minus SC |
| --- | ---: | ---: | ---: | ---: |
| All recorded | 606 | −0.05023 | +0.01958 | +0.06981 |
| Literal native | 45 | +0.04821 | +0.10879 | +0.06057 |

The corresponding block intervals for the differences of contrasts are [0.04827, 0.09228] and [0.01080, 0.10930]; first-source-equal values are +0.07268 [0.05527, 0.09003] and +0.06617 [0.02350, 0.10596]. Broad SC interval [−0.07488, −0.02838]; broad selfing [−0.01065, 0.05061]. Native SC [0.01592, 0.07516] and selfing [0.05830, 0.15193]. Northern-high same-carrier n=31, five blocks and three sources, does not support a further regional inference; tropical and southern intersections fail the original report threshold. Trait-specific (not common-carrier) broad northern-midlatitude contrasts are SC −0.03926 (721 islands), mating-system +0.03577 (638) and autonomy +0.03051 (57; interval crosses zero). On the same 67 native-SC eligible islands, native-minus-all SC share +0.06213 [0.03719, 0.09331]; for 49 mating-system islands +0.08499 [0.05790, 0.11391]. Absolute source contrast and roster selection must not be confused with the recurrent H1 isolation slope. A positive difference does not imply self-incompatible plants self or that one reproductive state causes the other.

## S9.2 Joint architecture relative to candidate sources

With truly paired, single-state symmetry and depth, define the joint frequency J = p(symmetry) p(depth) + covariance, computing the product and covariance within each roster. Island-minus-source covariance is decomposed into between- and within-genus components. The between term depends on genus weights and genus-level marginals; it is not a causal turnover estimate.

### Table S9b. Joint morphology contrasts (percentage points)

| Scope / pair | Islands | Joint difference | Covariance difference | Between genera | Within genera |
| --- | ---: | ---: | ---: | ---: | ---: |
| Tropical radial × deep | 666 | +10.989 | +4.018 | +4.057 | −0.039 |
| Tropical bilateral × deep | 666 | −1.368 | −3.425 | −3.639 | +0.214 |
| Tropical tubular/bell × deep | 598 | +5.172 | −1.097 | −1.634 | +0.537 |
| Native tropical radial × deep | 96 | +12.63 | +7.08 | — | — |
| Native tropical bilateral × deep | 96 | −6.58 | −6.60 | — | — |

Broad radial-deep block interval [8.958,13.313] pp; bilateral-deep [−2.68,−0.05]; tubular/bell-deep [3.82,6.61]. Native joint intervals [8.64,17.11] and [−8.82,−4.17]. Covariance intervals for broad radial, bilateral and tubular/bell pairs are [3.1737,4.9234], [−4.3149,−2.5685], [−1.9404,−0.2753] pp. The form category overlaps symmetry categories; do not sum all three joint frequencies or infer a conditional radial proportion. A later k5-disjoint cohort (1,579 northern, 670 tropical islands) gave raw tropical-minus-northern radial-deep minus bilateral-deep difference +0.108934 [0.077743,0.141316]. The four donor variants gave +0.100653 to +0.124492; all 96 support-qualified residual tests failed Holm correction (minimum 0.384). For High-only tropical, residual +0.046200, raw permutation p=0.004, Holm p=0.384; five genera account for 74.97% of absolute contributions, and subtracting their fixed-design contributions leaves +0.003122. Matched-grade results did not survive correction. Native regional gates failed. Historical oceanic–continental and LGM-unconnected–connected raw radial/deep contrasts (+6.69/+6.90 pp) are retained as exploratory source-history sensitivity; overlapping geological histories and incomplete ancestry audit preclude independent geological inference.

## S9.3 Depth conditional composition and regional factorization

For deep and nondeep categories require at least ten carriers in both groups on target and donor sides. Full 3 × 4 categorical tables retain asymmetric symmetry and intermediate depth where recorded; do not equate these cohorts to S8's two-depth restriction. Source comparisons use disjoint and inclusive pools, geographic and climate filtering and common eligible islands. Block and source group bootstraps and shared genus weights assess sensitivity. For J = d q, with d = P(deep) and q = P(radial | deep), the exact per-island difference is ΔJ = Δd(q_target+q_source)/2 + Δq(d_target+d_source)/2.

### Table S9c. Climate-disjoint radial-deep factorization

| Region | Common islands | Joint difference | Depth prevalence term | Conditional symmetry term |
| --- | ---: | ---: | ---: | ---: |
| Northern midlatitude | 264 | +0.0580128 | +0.0263188 | +0.0316940 |
| Tropical | 201 | +0.1087307 | +0.0730687 | +0.0356620 |
| Tropical minus northern | — | +0.0507179 | +0.0467499 | +0.0039681 |

In the original mode-specific conditional cohort, conditional radial proportion was higher than the source among deep carriers in northern (+0.1703343 [0.1415206,0.1978840]) and tropical (+0.1165609 [0.0850246,0.1485838]) floras. The northern direction survives joint genus/group weighting for all four nested variants, tropical only for k5-disjoint. The common-cohort regional differences in depth prevalence (+0.0783053) and conditional radial proportion (−0.0472502) do not both survive joint genus/group sensitivity across both groupings. Among conditional-eligible tropical k5-disjoint islands, bilateral-deep joint difference is +0.0126495 (interval crosses zero), unlike the broader cohort's negative joint difference, despite falling conditional bilateral prevalence (−0.0951204). No robust between-region conditional difference or within-lineage change is established. Native and High-only conditional cohorts fail the original regional gate.

---

# Appendix S10. Categorical functional assembly, rare profiles and isolation

## S10.1 Source-relative floral diversity and profile representation

The original triple-complete analysis retains 2,234 species with unambiguous colour, form and symmetry, 928 island floras with at least 30 carriers and 926 with an eligible candidate source. Equal-30 samples without replacement (30 rarefactions) are compared. The original plug-in marginal metric is the mean of the three 1−sum(p²) values, and joint Q is 1−sum(joint profile p²); joint Hill N2 is the reciprocal concentration. These nominal measures are distinct from later finite-sample distinct-species-pair Q and hierarchical with-replacement Q. Mutual information measures association of colour with the form–symmetry pair; excess MI subtracts its exact fixed-margin expectation. These are recorded categorical states, not continuous hypervolume or direct functional performance.

### Table S10a. Island minus candidate source: block-equal mean

| Region | Islands / blocks | Marginal Q | Joint Q | Open-radial share |
| --- | ---: | ---: | ---: | ---: |
| North midlatitude | 743 / 53 | −0.02571 | −0.01420 | +0.04582 |
| North high latitude | 53 / 7 | −0.08249 | −0.04614 | +0.13718 |
| Tropical | 106 / 42 | −0.00417 | +0.00268 | +0.03920 |
| Southern extratropical | 26 / 13 | +0.04309 | +0.03028 | +0.17101 |

Northern-midlatitude block intervals are [−0.03707,−0.01508], [−0.02189,−0.00687], [0.02357,0.06756]; first-source-equal counterparts −0.01746, −0.00985 and +0.02904 retain direction. Northern-high has only three distinct nearest-source signatures; tropical Q intervals cross zero and open-form source-group direction fails; southern is small. On the same 67 native-eligible northern islands, median marginal Q is +0.01304 all and +0.01546 native; native-minus-all +0.00136 [−0.00553,0.00786]. Native-minus-all joint Q and open-radial intervals also cross zero. The broad contraction is not a native-confirmed general result.

Across 2,000 sampled island pairs median taxon Jaccard 0.0350 and profile Jaccard 0.2845 (0.2124 after shared taxa removed); across 2,324 colour–form floras mean frequency overlap is 0.386901, species-disjoint 0.247183 and genus-disjoint 0.214879. Recurrence despite taxonomic dissimilarity is descriptive, not above-null adaptive convergence. In 669 northern Direct single-state islands profile concentration result −11.4% relative to genus-fixed permutation median (q=0.016); multistate sensitivity lacks adjusted support. Profile differentiation passes 12/16 overlapping fixed-binary-margin BY comparisons but 0/16 genus-fixed, shared-name-exclusive comparisons. In 928 triple-complete floras 909 show positive excess MI; northern-midlatitude median 0.182 nats, tropical 0.066. These medians are not a formal regional contrast.

## S10.2 Discrete profile-space and matched genus count

Reference set: 100 actual three-trait profiles, with 29 rare (1–2 species), 25 common (≥26) and 46 intermediate profiles. A reference profile is supported within one Hamming mismatch. Exact-k reference subsets account for number of observed profile centres. This is categorical neighbourhood coverage, not a continuous hypervolume. Geographic k5 candidate donors are excluded if sharing focal accepted species in the disjoint analyses. The original tropical disjoint-source k-conditioned contrast −0.02469 [−0.03207,−0.01766] (source-equal −0.02435 [−0.03257,−0.01638]) is sensitive to species-weighted references and native status.

Fixing genus richness: sample 15 distinct genera and one species per genus on each side, repeat 40 independent draws, then compare rare-minus-common coverage relative to exact-k references. Ten-genus and identical-genus identity designs retain their own cohorts. The conditional block and first-source-signature resampling groups overlap and are not independent evolutionary histories.

### Table S10b. Rare-minus-common profile gaps

| Scope | Islands | Block-equal difference | Conditional interval |
| --- | ---: | ---: | --- |
| Tropical, 15 genera | 106 | −0.03669 | [−0.05308,−0.01896] |
| Northern midlatitude, 15 genera | 742 | −0.02694 | [−0.04257,−0.01153] |
| Southern extratropical, 15 genera | 26 | −0.00879 | [−0.03587,0.01868] |
| Native northern midlatitude, 15 genera | 67 | +0.01276 | [−0.01511,0.03674] |
| Native tropical, 15 genera | 17 | −0.00865 | [−0.04995,0.02809] |
| Tropical, same 5 genus identities | 81 | +0.00261 | [−0.01579,0.01982] |
| Tropical, same 8 genus identities | 50 | +0.00914 | [−0.01475,0.03532] |

Source-equal n15 tropical/north gaps −0.03972/−0.02227; both negative intervals. At n10 distinct genera tropical/north −0.03160/−0.01913, also negative under source weighting. Exact rare-support genus accounting on 106 tropical islands: total −0.02970, shared genera +0.00945, one-sided genera −0.03914 (rounding); these are descriptive Shapley terms, not mediation. The same-genus designs select different species, islands and fewer genera. WCVP regional-native sensitivity: northern all-recorded rare gap −0.02627 versus WCVP +0.01851; same-island WCVP-minus-all +0.02256 [0.01186,0.03462]. These results establish dependence on the recorded genus identity and origin scope, not universal native contraction.

## S10.3 Corrected-isolation functional-space associations

The corrected replay uses the same frozen observed-flora parent and Direct traits, 928 triple-eligible islands, 927 matched geography, 87 strict-native before final support gates. Historical source-relative metrics and older geographic gradients are not interchangeable with these regressions. The measured-coverage and same-morphology-roster reproductive adjustments retain their original multiplicity families. A 399-draw coherent regional genus whole-profile null retains common species labels across occurrences. It tests exchangeability within genus, not phylogenetic independence.

### Table S10c. Coverage and reproductive conditioning

| Region / response | Islands | Slope | BY q |
| --- | ---: | ---: | ---: |
| Northern MI, coverage | 742 | +0.02003 | 0.02353 |
| Tropical marginal Q, coverage | 105 | −0.01108 | 0.00850 |
| Northern MI, same-roster selfing | 739 | +0.01770 | 0.03109 |
| Tropical marginal Q, same-roster selfing | 68 | −0.01087 | 0.01846 |

Coverage-adjusted intervals: [0.00825,0.03180], [−0.01661,−0.00555]. Same-roster intervals [0.00699,0.02840] and [−0.01662,−0.00512]. Neither baseline slope is BY-supported; all-flora reproductive conditioning yields BY 0.05051 for both. Joint Q is unsupported. The coherent genus null expected northern MI +0.012683 and tropical marginal Q −0.008676; residuals +0.007343 and −0.002408, focal Holm 0.10 for both; residual BY 0.05378. Same-roster selfing-conditioned residuals +0.006585/−0.002353 have pointwise-model BY 0.04837/0.04782 but both coherent p=0.055 and BY 0.68271; the null does not preserve all morphology–reproduction coupling. The native subset (67 north,17 tropical,3 south) provides no adjusted support. These are conditional associations, not within-lineage selection or causal independence.

Earlier, differently specified tropical slopes for zygomorphic prevalence (−0.559 pp per log2(1+km), Holm 0.0295) and SC (+0.593 pp, conditional genus residual +0.297 pp, Holm 0.0476) are historical sensitivity only; autonomous +1.053 pp has residual +0.058 with interval crossing zero. The SC and autonomy cohorts differ, and native robustness is absent. The corrected H1 estimates supersede these historical slopes. Missing specialist profile nulls did not support unusual depletion.

## S10.4 Carrier completeness and marginal versus joint structure

Four-axis colour/form/symmetry/depth complete carriers yield 1,326 primary islands, 370 species and 116 profiles; tropical n264 and northern n1,002, southern below gate. Across five outcomes and three regional slots all 15 primary Holm tests equal 1. Tropical raw additive Q +0.020244, joint Q −0.005984, fixed-margin departure −0.006267 and rarefied repertoire −0.010752 have pointwise conditional intervals excluding zero, but are not supported residual discoveries. Requiring all four axes changes the tropical broad radial-minus-bilateral deep source contrast +0.122974 to −0.086662; ordered island eligibility and symmetric target/donor carrier terms −0.048012, −0.103593, −0.058031 sum to the change −0.209636. Requiring colour alone reverses direction; requiring form alone does not. Colour availability among radial-deep versus bilateral-deep is 0.148988 versus 0.390079 in tropical target carriers and 0.341036 versus 0.333777 in donors; the difference of differences −0.248349 is nearly the genus/pair-grade expectation −0.248889 (p=0.988). All seven primary documentation-null Holm values equal 1. No missing-data correction is established.

Six pair-specific spaces, 24 panels: all 18 primary regional residual Holm values equal 1. On common three-axis carriers (north 1,460/tropical 492), the regional difference between colour–symmetry and form–symmetry joint-Q source contrasts is −0.043328 [−0.072109,−0.013644]. Exact marginal part −0.041820 and coupling part −0.001508 (interval crosses zero); coherent genus-grade p=0.364. The follow-up within-block source-relative isolation-gradient regional contrast −0.0011705 [−0.0132426,0.0130399], coherent p=0.726; native/high subsets fail support. Local fixed-margin partition (43 of 96 panels supported for fitting, not 43 significant) found mean joint Q departures largely reproduced by genus-grade strata, with only one joint-Q isolation slope having a pointwise interval excluding zero. Form–symmetry slope −0.001757 is pointwise. The northern open-radial and tubular-zygomorphic cell directions largely match genus-grade expectation; one genus dominates the lip residual. No universal within-genus dissolution of trait coupling.

---

# Appendix S11. Repertoire mismatch and scale dependence

## S11.1 Equal-effort repertoire mismatch

For each profile c, let h_T(c) and h_S(c) be exact probabilities of observing that profile in independent samples of k species without replacement from recorded island and candidate-source rosters. Expected exclusive counts divided by k are sum h_T(1−h_S)/k and sum h_S(1−h_T)/k. Balanced mismatch is twice the smaller expected exclusive count divided by k (equivalently twice the smaller already-normalized exclusive count). Signed deficit is source-only minus target-only. These are expected-count statistics, not probabilities, samplewise turnover or historical replacement.

### Table S11a. Tropical-minus-northern balanced mismatch, k=5 disjoint

| Pair | Difference | Pointwise block interval | Primary genus-grade Holm |
| --- | ---: | --- | ---: |
| Colour × form | +0.134580 | [0.099908,0.169181] | 0.102 |
| Symmetry × depth | +0.129073 | [0.106249,0.150335] | 0.088 |

Signed deficit contrasts −0.001914 and −0.014204 have intervals spanning zero; both Holm=1. All 24 nested regional mismatch contrasts are positive; 20/20 support-qualified pointwise intervals above zero. First-source-weighted differences +0.064493/+0.070853. The four primary genus-grade tests fail; the nested panels are not independent replications.

## S11.2 Sampling-baseline and systematic-distance decomposition

The exact sum of expected between-roster exclusive counts divided by k equals the mean of within-roster independent-draw mismatch plus sum((h_T−h_S)²)/k. This identity separates finite-sample mismatch from systematic occupancy differences, not an evolutionary process. Subsequent analyses were selected after the original mismatch result; the two-endpoint residual Holm family is separate.

### Table S11b. Systematic symmetry–depth regional contrast

| Component | Geographic-block weighting | First-source weighting |
| --- | ---: | ---: |
| Observed tropical minus north | +0.038840 | +0.046410 |
| Genus-grade expectation | +0.025879 | +0.035874 |
| Residual | +0.012962 | +0.010536 |

Primary k5 disjoint n1,579 north/102 blocks/98 source groups and 670 tropical/91 blocks/76 source groups. North systematic distance 0.077988, tropical 0.116828. Observed block interval [0.021211,0.058194], residual [0.007163,0.018836], coherent p=0.014, Holm p=0.028 over the two post-result endpoints; leave-block signs positive. First-source residual interval [0.004859,0.016865], p=0.014. Colour–form residual +0.005649, p=0.414. Climate-source weighting, High-only, n10 and native sensitivity limit generalization; native fails support and residual reverses direction. No new independently confirmed effect.

Radial-state contributions sum +0.012104 (93.4% of net block residual): deep +0.004373, shallow +0.002875, intermediate +0.002663, absent/open +0.002192. Other states can offset the total; these are exact arithmetic contributions, not individual-cell significance tests. The observed systematic-distance regional contrast +0.038840 comprises directional mean difference +0.016069 and between-island heterogeneity +0.022771; first-source +0.046410 = +0.006184 +0.040227. These follow E(δ²)=E(δ)²+Var(δ); separate nonlinear component uncertainty is not available. The balanced mismatch contrast +0.129073 includes within-roster baseline +0.077736, systematic distance +0.038840 and reduced absolute deficit +0.012498 (rounding).

## S11.3 Rao Q, hierarchical diversity and source aggregation

For distinct species pairs, joint symmetry–depth Q equals P(depth differs) plus the sum across depth classes of P(both in class) times P(symmetry differs | both in class). Separate exposure weights from conditional symmetry diversity, retaining one-sided support separately. On the conditional-common climate-disjoint cohort, joint Q source deltas north +0.02483285 = +0.02671208 depth Q −0.00187923 same-depth contribution; tropical +0.01738004 = +0.00731118 +0.01006886. Northern positive persists across both grouping genus sensitivities; tropical does not. Broader n≥5 cohort north −0.03111390, tropical +0.00651876. Ordered northern eligibility change +0.05595 partitions +0.04635 group support, +0.00355 richness and +0.00604 depth balance. This is order-dependent accounting, not a unique causal partition.

A separate with-replacement profile measure decomposes regional gamma into equal-group/equal-island alpha plus within-group beta plus between-group beta. Expected equal-size samples use m=min(target, source carriers) for each island, with exact finite-sample mixture corrections. The pooled broad cohort is 2,244 islands; first-ranked comparison retains 2,190, with common m within island across source definitions.

### Table S11c. Joint hierarchical diversity differences (percentage points)

| Cohort / reference | Alpha | Within-group beta | Between-group beta | Gamma |
| --- | ---: | ---: | ---: | ---: |
| North, broad observed | −8.133 | +2.244 | +4.025 | −1.864 |
| Tropical, broad observed | −5.302 | +2.892 | +2.177 | −0.233 |
| North, broad matched size | −2.630 | −0.869 | +1.659 | −1.841 |
| Tropical, broad matched size | +0.647 | −0.525 | −0.327 | −0.205 |
| North, common-size pooled | −2.073 | −0.750 | +1.441 | −1.382 |
| North, common-size first-ranked | −0.827 | −0.265 | +1.133 | +0.042 |
| Tropical, common-size pooled | +0.725 | −0.533 | −0.323 | −0.131 |
| Tropical, common-size first-ranked | +0.491 | −0.419 | −0.812 | −0.739 |

Climate-disjoint block-weighted examples. Raw alpha negative and both beta terms positive across the four broad donor variants and two groupings, including observed-profile genus/group weighting. However observed gamma sensitivity crosses zero, and after size matching within-group beta reverses. First-ranked-minus-pooled common-size gamma differences +1.424 [0.939,1.900] north and −0.608 [−0.924,−0.266] tropical demonstrate comparison-reference dependence on the same islands. Tropical High-only first-ranked matched gamma ranges −3.023 to −2.381 pp and between-group beta −5.532 to −2.815 pp across eight nested comparisons; observed-profile genus-weighted gamma sensitivity usually spans zero. There is no jointly rarefied/genus-reweighted inference. Native regional gate fails. No stable ecological beta compensation or equivalence claim.

---

# Appendix S12. Prediction, syndrome proxies and independent functional constraints

## S12.1 Predictability and support

Endpoint-specific Direct coverage in the 4,505-island parent: colour 3,132, form 2,806, symmetry 3,145 and SC 2,958 islands; true colour–form 2,324, colour–symmetry 2,699 and form–symmetry 2,611 with ≥5 carriers. Observed genus labels are not predicted colonization histories. Earlier local prediction (606 islands/44 archipelagos) reduced exact error 0.007351→0.005473 and kernel error 0.011137→0.006574 (BY approximately 0.00093), but geographical nearest reference performed comparably and taxa can overlap. Matched donor support remains narrow. Symmetry was variable in only 81/2,200 multi-species genera in the recorded 32,039-name Direct roster; held-out Brier 0.2381→0.01768 with known genus, with post-hoc smoothing. Nine/nine matched Direct symmetry component wins and five/nine SC wins do not establish general transport.

The older geographically excluded known-genus forecast reported seven of 13 endpoints improving versus external blind. The oracle-free geography→genus→trait analysis is distinct. Candidate-source increment in genus Brier +0.000390, +0.000122, +0.000348 for north/tropical/south, respectively; geographical nearness performs better in north. In 39 all-recorded endpoint-region panels, 13 improved corrected nominal Q error relative to local blind, nine proper category scores, six both; tropical symmetry n925 Q-error gain +0.008968 and southern primary colour n184 +0.016851 are strict benchmark candidates. Yet the taxon-equal analysis yields zero of 39 genus-equal wins over all three reference models, and no native Southern replication. In 59 qualified panels congener-profile mismatch has the largest signed projection onto total category prediction error, not a causal variance share. Evidence-grade category association beyond genus appears in 43/55 supported panels; different species by grade are not repeated measurements.

### Table S12a. Held-out profile prediction on corrected-isolation cohorts

| Outcome | Support | Observed slope | Predicted slopes (two fixed seeds) | MSE decrease |
| --- | --- | ---: | --- | --- |
| North MI | 742 islands | +0.02003 | +0.00979/+0.00862 | 54.6%/51.2% |
| Tropical marginal Q | 105 islands | −0.01108 | −0.00444/−0.00399 | 53.0%/49.5% |
| Native northern MI | 67 islands | +0.00484 | — | 39.3%/25.5% |

Held-out species proper Brier north 0.9380→0.6479, tropical 0.9597→0.8356; log loss 3.5984→2.5026, 3.8737→3.2662. These older taxon holdouts are not geographic/phylogenetic holdouts; mean donor coverage 72–74% north and 45–47% tropical. Forecast draw ranges include zero for tropical Q. They complement S8.8 but are not that four-cell geographic holdout. Native tropical below fit gate. The own-database reproductive-capacity support audit (576 slots,335 fits) found one BY-supported documentation signal, SC-known colour support in northern high latitudes +0.02635 [0.01491,0.03778], BY 0.01738 (171 islands); autonomy-known triple profile coverage median 0.88505/0.70036 north/tropical is descriptive with sparse 132 global names and no supported autonomous slope. This is potential capacity and documentation, not actual reproduction without pollinators.

## S12.2 Plant-defined syndromes

Five nonexclusive plant morphology proxy templates were explored, not assigned animal guilds. The archived pair atlas contains 80 panels and 64 scope/pair/template groups. Retrospective consistency across four donor/overlap variants and both grouping intervals yields 17/64 covariance screens, 10/64 also meeting phi support (six broad, four literal native). No endpoint passes both origin scopes. The following screen-qualified phi ranges are descriptive, not multiplicity-adjusted tests.

| Scope | Proxy / trait pair | Phi difference range |
| --- | --- | --- |
| Broad | Bird: colour × depth | +0.07629 to +0.09605 |
| Broad | Bird: colour × size | +0.05027 to +0.05481 |
| Broad | Bat: colour × size | −0.07006 to −0.06064 |
| Broad | Bird: form × depth | +0.07190 to +0.08216 |
| Broad | Specialist bee: symmetry × depth | +0.03731 to +0.05120 |
| Broad | Open generalist: symmetry × depth | −0.03949 to −0.03620 |
| Literal native | Specialist bee: colour × symmetry | −0.05746 to −0.05395 |
| Literal native | Specialist bee: form × symmetry | −0.11535 to −0.08899 |
| Literal native | Open generalist: form × symmetry | −0.14960 to −0.12919 |
| Literal native | Specialist bee: symmetry × size | +0.07350 to +0.08544 |

Both specialist and open-generalist native form–symmetry values are negative, not evidence of substitution. Colour–size coherent-genus joint sensitivity passes 0/72 contrasts under both grouping definitions. All are proxy patterns, not realized effective visitors.

## S12.3 Animal channel and functional bridge limitations

The local five-guild by four-region positivity audit qualified 0/20 comparisons; inadequate detection cannot be equated with absence. A hypothetical EuPPollNet all-bee removal comparison (76 networks, 11 islands, four blocks) shows +9.36 pp median excess disconnection relative to equal-species random removal; simulation, not observed decline or reproductive effect. The 21-island Bombus-nondetected/bird-detected candidate set yields no floral Q, joint-profile coordination or open-radial direction robust to block resampling. Birds tubular association is positive in 11/11 nested public comparisons; bird-family permutation supports 6/24 exploratory BH contrasts but none under genus restriction. In a 30-island common plant cohort, bird diet specialist-minus-partial source-relative deep residuals +0.195580/+0.187680, form −0.093692/−0.074534 and joint form/depth −0.035551/−0.030217 for geo/climate references; no 160-test Holm pass or native robustness. Original broader form/depth family had 0/500 corrected tests. Origin-uncertainty bounds: 460 feasible endpoints, zero sign-identified; no effect inferred. Bombus lip/papilionaceous +0.0148127 [−0.015426,0.037922]; bat tubular +0.06251 across 49 islands/four blocks Holm 0.273; original bat source interaction −0.0263522, p=0.205, Holm 1. Four-guild common support at most seven islands. These data cannot distinguish pollinator replacement or actual floral selection. The direct GloPL–FD island bridge (30 islands/26 blocks) had all three Q/open-form intervals crossing zero and no same-island same-species match. The original H4 component-score bridge and S8.5 insufficient joint-morphology fit are separate endpoints.

## S12.4 Evidence coverage and archive boundary

Electronic Table_S12b_migrated_evidence_inventory.csv records all eligible and excluded entries in the historical source-first registry, with source support, reported result, inference ceiling and appendix disposition; Table_S12c_recent_evidence_inventory.csv records the later analysis families and source-aggregation/eligibility extensions. These inventories summarize archived reports and are not full statistical panel tables. Full archived numerical panels, code and provenance require a lawful editor-approved reviewer access mechanism. Two incomplete-ancestry records and one independent external-only readiness entry are excluded from empirical synthesis, not counted as negative biological effects. The historical directional results and limited support assessments do not alter the corrected H1–H4 selector. All new appendices are outcome-motivated; no new models were fitted to compile this version.

Reference for the additional candidate floras: Weigelt, P., König, C. & Kreft, H. (2020). GIFT – A Global Inventory of Floras and Traits for macroecology and biogeography. J. Biogeogr., 47, 16–43. https://doi.org/10.1111/jbi.13623
