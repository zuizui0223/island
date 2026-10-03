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

## S3.2 Current seven-response taxonomic representation depth

We asked whether the current H1 vector is represented primarily by family/genus composition or whether a finer-grained residual remains. For each atomic outcome, scored species were retained only when both their family and genus contained at least one other scored species. Family and genus expectations were leave-one-species-out means from the fixed scored species pool. Taxonomy therefore never filled missing traits. The same island-species observations were used at observed, family-residual and genus-residual stages.

A beta-binomial common-support gate was fitted before interpreting residualization. H1 remained FDR-supported on this exact taxonomically eligible species set in all four regions in each of four surfaces: all-observed and WCVP regional-native-compatible floras, each under all-analysis and Direct-only trait evidence.

Across the four broad and WCVP regional-native-compatible surfaces, the tropical seven-response vector remained supported after genus residualization:

| Flora / evidence | Post-genus q |
| --- | ---: |
| All-observed, all-analysis | 0.00733 |
| All-observed, Direct-only | 0.00559 |
| Regional-native-compatible, all-analysis | 0.00144 |
| Regional-native-compatible, Direct-only | 0.000311 |

Self-compatibility remained positive and individually supported after genus residualization in all four surfaces (P = 0.00293, 0.00355, 9.93 × 10^-5 and 7.09 × 10^-5). Plain colour and generalized form were also supported post-genus in both regional-native-compatible surfaces.

Paired 300-draw spatial-block bootstrap quantified attenuation uncertainty. In the tropical all-analysis scopes, the additional family-to-genus attenuation was positive: median 0.458, 95% interval 0.212–0.662 in all-observed flora and median 0.341, interval 0.150–0.552 in regional-native-compatible flora. Direct-only intervals crossed zero. Thus genus composition contributes to the tropical vector in these broad surfaces, but does not exhaust it. Section S3.3 tests whether that conclusion is stable to exact floristic-status provenance.

Northern high latitudes showed substantial family-level attenuation: the family attenuation interval excluded zero in all-observed all-analysis (0.300–0.702), all-observed Direct-only (0.181–0.573) and regional-native-compatible all-analysis (0.156–0.670); the regional-native-compatible Direct-only lower bound was -0.0019. The incremental family-to-genus attenuation crossed zero in every northern-high surface, so a sharp genus breakpoint is not identified.

Northern mid-latitude depth was model/evidence sensitive because the equal-island observed stage did not reproduce the H1 vector in two of four surfaces. Southern extratropical depth was floristic-status sensitive: genus residuals remained supported in all-observed flora but not in regional-native-compatible flora, while the paired attenuation intervals did not establish a precise native family/genus breakpoint.

This audit localizes taxonomic representation depth only. Persistence after genus residualization does not prove within-lineage evolution; it can reflect within-genus species sorting, unmeasured source composition, persistence filters, introductions in the broad observed flora, or genuine evolutionary change. Conversely, attenuation after family/genus adjustment does not prove dispersal or colonization filtering.

Reproducibility surface:

- results/h1_taxonomic_depth_current_20261003/
- src/island_v2/chapter1_h1_taxonomic_depth_current.py
- config/chapter1_h1_taxonomic_depth_current.yml
- validation workflow run 37095612710; artifact 11264117016; digest sha256:d2b46479d6eb40457157490b667660bf9728361faae0ab7321ec1cc8d6d1db5b

## S3.3 Tropical floristic-status provenance diagnostic

The broad/WCVP tropical below-genus result was further partitioned by floristic-status provenance. This diagnostic uses the same current H1 taxonomic-depth implementation but disables the already-completed bootstrap because the purpose is support localization, not a second attenuation-uncertainty analysis.

Four non-overlapping provenance groups were examined:

- source-backed native non-endemics: 118,429 island-by-species rows, 425 islands, 26,200 species;
- source-backed endemics: 17,444 rows, 143 islands, 17,429 species;
- source-backed native records with unresolved endemism: 14,658 rows, 380 islands, 8,469 species;
- records originally unresolved for origin but upgraded by WCVP regional-native compatibility: 362,789 rows, 2,281 islands, 60,117 species.

The tropical results were:

| Status provenance / evidence | Parent common-support H1 q | Post-genus q | Interpretation |
| --- | ---: | ---: | --- |
| Native non-endemic / all-analysis | 0.003758 | 0.1177 | compatible with genus structuring |
| Native non-endemic / Direct-only | 1.85 × 10^-5 | 0.0754 | compatible with genus structuring |
| Endemic / all-analysis | 0.7246 | 0.2114 | parent H1 not reproduced |
| Endemic / Direct-only | not testable | not testable | one retained outcome |
| Source-native, endemism unresolved / all-analysis | 0.01266 | 0.000648 | below-genus residual retained |
| Source-native, endemism unresolved / Direct-only | 0.3600 | 0.00202 | parent beta-binomial H1 gate fails; no depth inference |
| WCVP-upgraded origin-unresolved / all-analysis | 9.88 × 10^-6 | 0.00215 | below-genus residual retained |
| WCVP-upgraded origin-unresolved / Direct-only | 0.002415 | 0.000603 | below-genus residual retained |

Thus the strict source-backed native-nonendemic tropical H1 is strongly supported before genus adjustment but is not supported after genus residualization. The source-backed endemic subset cannot adjudicate a finer-grained response because its parent H1 is unsupported or not testable. The most stable below-genus residual instead occurs in the large WCVP-upgraded group whose exact focal-island origin status remains unresolved.

This changes the interpretation of S3.2. The tropical residual cannot be promoted as evidence of repeated native within-lineage evolution. It is a status-provenance-sensitive finer-than-genus component of the broad contemporary/regional-native-compatible flora. Within-genus species sorting, regional source-pool structure, imperfect focal-island status assignment, persistence filtering and genuine evolutionary change remain unresolved alternatives.

Reproducibility surface:

- results/tropical_status_depth_20261003/
- validation workflow run 37096771148; artifact 11264951254

## S3.4 Frozen source-matched genus entry versus within-genus loading

A separately frozen historical lineage-representation bridge was replayed after the
24 September geography correction to refine the strict tropical native-nonendemic
assembly interpretation. This bridge is not the current seven-response H1 estimand. It
uses the previously frozen floral functional position
`(-large_bee_like + generalized_accessible) / 2` and asks whether source-matched change
is expressed through representation of genera or through extra species weighting within
genera that are already represented.

All historical design choices were retained: tropical context, source-backed
`native_nonendemic` flora, broad and Direct evidence scopes, four source modes,
`prevalence_richness` source matching, a minimum of five represented genera, equal-island
OLS, island area and climate PC1–PC4 controls, and spatial-block cluster-robust
covariance. Only geographic distance was replaced.

Before the replacement, all 24 frozen slopes were reproduced to numerical precision
(maximum absolute slope and SE differences 1.15 × 10^-16; maximum P-value difference
2.17 × 10^-15).

On corrected geography, genus-entry enrichment was negative and FDR-supported in all
four source modes in both evidence scopes:

| Evidence | Entry β range | Maximum entry q | Loading β range | Maximum loading q |
| --- | ---: | ---: | ---: | ---: |
| Broad | -0.06861 to -0.05552 | 0.02925 | -0.02046 to -0.01657 | 0.15221 |
| Direct | -0.06713 to -0.06269 | 0.02059 | -0.02122 to -0.01617 | 0.17060 |

Species-weighted enrichment was also FDR-supported in all four source modes
(maximum q = 0.02860 broad; 0.02606 Direct), whereas the loading increment was not
FDR-supported in any source mode. The supported species-weighted shift is therefore
already expressed at genus entry, with no supported additional signal from weighting
species within represented genera.

This result is concordant with S3.3, where the current seven-response strict
native-nonendemic tropical H1 becomes unsupported after genus residualization. However,
the estimands differ: the entry/loading bridge uses an older frozen floral functional
position. It therefore provides supplementary assembly concordance rather than a direct
mechanistic decomposition of current H1. Genus entry can reflect arrival, establishment,
persistence, habitat filtering or biotic interactions and does not identify any one of
those processes.

Reproducibility surface:

- results/corrected_lineage_entry_loading_20261003/
- scripts/geography_correction/replay_corrected_lineage_entry_loading.py
- validation workflow run 37097462474; artifact 11264653941; digest sha256:afc980240819a0ce11fe8269fb7ec71d89e001184481dd1f0970e55de37625fb

## S3.5 Literal additive decomposition of the current H1

The historical scalar bridge in S3.4 is not the current seven-response H1 estimand. We
therefore rebuilt the source-matched decomposition using the literal current H1 binary
state separately for each atomic outcome. An initial v1 current-H1 bridge still inherited
the historical all-island-species weighting inside genera; because current H1 is defined
only on trait-resolved species for each outcome, that version is not used for scientific
interpretation.

The final v2 analysis uses only trait-resolved, source-candidate, source-genus-scored
native non-endemic species and enforces the exact identity

`raw H1 mean = source expectation + genus-entry enrichment + within-genus loading + within-genus trait residual`.

The row-level identity passed to machine precision in both evidence scopes (maximum
absolute error 1.11 × 10^-16). Regressing each additive component on corrected
isolation, island area and climate PC1–PC4 with spatial-block cluster-robust covariance
also closed exactly: maximum slope identity error was 6.59 × 10^-17 across 28
all-analysis atomic identities and 1.04 × 10^-16 across 24 Direct-only identities.

The result is an antagonistic hierarchy rather than a serial assembly route:

| Component | All-analysis mean slope across source modes | All-analysis max q | Direct-only mean slope across source modes | Direct-only max q |
| --- | ---: | ---: | ---: | ---: |
| Raw H1 mean | -0.00358 to -0.00180 | 0.000219 | -0.00724 to -0.00559 | 0.000798 |
| Source-availability-matched expectation | +0.00637 to +0.01019 | 0.799 | +0.00702 to +0.00798 | 0.495 |
| Genus-entry enrichment | -0.01392 to -0.01013 | 0.01287 | -0.01649 to -0.01467 | 0.07951 |
| Within-genus loading | -0.000245 to +0.00170 | 0.01470 | -0.000663 to -0.000028 | 0.340 |
| Within-genus trait residual | -0.000025 to +0.00116 | 0.05726 | +0.00132 to +0.00286 | 0.02048 |

The restricted raw H1 vector is jointly nonzero in every source mode but its oriented
mean is negative, so this source-evaluable strict native-nonendemic subset cannot be
treated as a miniature version of the globally recurrent H1 direction. Source
expectations are directionally positive but are not 4/4 FDR-supported. Genus entry is
negative in every source mode and is 4/4 FDR-supported only in the all-analysis scope.
Additional within-genus loading has no stable positive direction.

The only 4/4 positive Direct-only component is the within-genus trait residual. Its
joint vector is FDR-supported in all four source modes (maximum q = 0.02048). Autonomous
selfing is the strongest positive atomic residual (beta = +0.0267 to +0.0310 across
source modes; P = 0.0035–0.0049). Actinomorphy and generalized form also point positive,
whereas plain colour, self-compatibility and shallow/open tube point negative in this
restricted frame.

Thus current H1 is not generated here by preferential entry of more H1-like source
genera. Instead, source expectations, realized genus entry and finer-grained species
trait composition oppose one another. The positive Direct-only within-genus residual
does not demonstrate within-lineage evolution: it can still arise from within-genus
species sorting, persistence, unmeasured source structure or genuine evolutionary
change. Because the restricted raw vector itself differs in orientation from the broad
global H1, this decomposition is retained as a mechanistic boundary rather than promoted
as the generator of global recurrence.

Reproducibility surface:

- results/current_h1_source_lineage_decomposition_v2_20261003/
- config/chapter1_current_h1_source_lineage_decomposition_v2.yml
- src/island_v2/chapter1_current_h1_source_lineage_decomposition_v2.py
- validation workflow run 37101224737; artifact 11266406362; digest sha256:57e5c4319a47eb0c2f55c01ee17c9afc2085f62142e7cc91bc014dafbcfdae70

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
| Primary | 0.09191 | 0.03806 | 0.01575 | 1,408 | 919 | 1,248 |
| Supplemental-only | 0.04410 | 0.04351 | 0.31082 | 828 | 470 | 736 |
| No-zero-constant | 0.09089 | 0.03865 | 0.01869 | 1,375 | 912 | 1,238 |

The primary and no-zero-constant models support a positive isolation association. Supplemental-only remains positive but unsupported. H3 is an isolation-associated pollen-delivery constraint, not direct evidence of declining pollinator abundance or visitation.

### Table S5b. Offshore-gradient robustness

Because 996 of 1,248 GloPL sites are true continental zero-distance sites, we tested whether H3 reduces to a mainland/offshore step.

| Analysis | β | SE | Two-sided P | Cells | Publications | Sites |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Offshore continuous gradient | 0.22031 | 0.09704 | 0.02319 | 276 | 153 | 252 |
| Mainland vs offshore indicator | 0.12743 | 0.08714 | 0.14364 | 1,408 | 919 | 1,248 |

A leave-one-publication-out jackknife of the offshore gradient retained positive estimates in all 153 deletions. The minimum coefficient was 0.17850, the maximum was 0.26930 and the weakest two-sided P value was 0.04892.

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
