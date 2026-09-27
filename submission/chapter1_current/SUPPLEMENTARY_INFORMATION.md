# Supplementary Information

## A recurrent global floral island syndrome extends beyond the selfing syndrome

This Supplementary Information accompanies the Ecology Letters first-shot manuscript. It expands the corrected Chapter 1 submission baseline selected by `config/chapter1_submission_current.json`. The source of truth for corrected numerical results is `results/geography_20260924/`.

The geography correction is a post-hoc repair of a known exposure mismatch. H4 is post-hoc functional triangulation. Neither is reclassified here as prospective confirmation.

---

# Appendix S1. Corrected geographic exposure

## S1.1 Island universe

The corrected plant analysis contains **8,264 island units**. The original locked universe contained 8,265 units; one unit was the western dateline-split component (`0-W`) of the Eurasian continental sibling rather than an island and was excluded. The broad H1 analytic union contains **4,379 islands**.

The original exposure calculation compared GSHHG island geometry with coarser Natural Earth continental polygons. Within the broad H1 union, **1,113 islands** had therefore received spurious zero distance. Recalculation against continental coastlines reconstructed from the same GSHHG 2.3.7 high-resolution archive assigned positive distance to all 1,113.

For GloPL, exposure is calculated from study-site coordinates rather than from plant island identities. Of 1,248 unique GloPL sites, **996 lie on seeded continental land** and legitimately retain zero distance.

### Table S1. Corrected geography summary

| Quantity | Corrected value |
|---|---:|
| Global island units | 8,264 |
| Broad H1 union | 4,379 |
| Spurious island zero distances repaired | 1,113 |
| Excluded continental split components | 1 |
| GloPL unique sites | 1,248 |
| True continental GloPL zero-distance sites | 996 |
| Sphere radius | 6371.0088 km |
| Independent numerical geometry checks | 24 |
| Maximum analytic–numerical discrepancy | 1.13 × 10^-8 km |

## S1.2 Distance algorithm

Distance is minimum coastline-to-coastline separation between minor great-circle arcs on a mean-radius sphere. Continental siblings are reconstructed before distance calculation and artificial dateline split edges are omitted. Segment indexing uses triangle-inequality pruning followed by exact arc minima.

The implementation passed analytic crossing/endpoint/dateline tests, subdivision checks and 24 independent numerical minimizations. The numerical tolerance validates the algorithm only; it does not imply metre-scale accuracy of the source coastline data.

Machine-readable receipts:

- `results/geography_20260924/spherical_geometry_validation.json`
- `results/geography_20260924/distance_same_coastline_audit.json`
- `results/geography_20260924/corrected_universe_exclusion.csv`
- `results/geography_20260924/formerly_zero_islands_recalculated.csv`
- `results/geography_20260924/glopl_corrected_site_distances.csv`

---

# Appendix S2. Island floras and trait evidence

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was retained with source provenance and normalized to a fixed accepted-species axis. The final scientific trait database contains **106,295 species × 3 evidence axes = 318,885 possible cells**, of which **222,688 (69.83%)** are resolved.

The three raw evidence axes are:

1. flower colour;
2. floral structural complexity;
3. reproductive assurance.

Species-direct High/Medium evidence forms the Direct-only sensitivity. The primary all-analysis scope additionally admits validated lower-confidence evidence where direct evidence is unavailable. Missing trait information remains missing rather than becoming trait absence.

The full immutable trait ledger contains third-party material with heterogeneous redistribution rights and cannot be redistributed wholesale. The current rights-filtered public derivative contains **46,274** redistribution-authorized resolved cells. Omission from the public derivative indicates rights status rather than biological missingness or scientific exclusion.

Scientific database SHA-256:

`a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`

---

# Appendix S3. H1 full multivariate response

H1 models seven atomic responses, coded so that a positive isolation coefficient is in the predicted island-syndrome direction:

1. self-compatibility;
2. predominantly/obligately selfing mating system;
3. autonomous/delayed autonomous selfing;
4. plain colour;
5. generalized floral form;
6. actinomorphic symmetry;
7. shallow/open floral tube.

Within each geographic stratum, the seven-dimensional isolation vector is tested by a joint Wald test after fitting beta-binomial models with standardized log isolation, island area and climate PC1–PC4. Spatial-block cluster-robust covariance is used.

### Table S2. H1 joint response-vector tests in the broad all-observed analysis

| Region | All-analysis islands | All-analysis q | Direct-only islands | Direct-only q |
|---|---:|---:|---:|---:|
| Northern mid-latitude | 2,173 | 3.216 × 10^-10 | 2,161 | 1.370 × 10^-4 |
| Northern high latitude | 411 | 2.433 × 10^-5 | 408 | 3.794 × 10^-8 |
| Tropical | 1,493 | 3.498 × 10^-7 | 1,467 | 3.528 × 10^-5 |
| Southern extratropical | 302 | 2.504 × 10^-18 | 292 | 3.155 × 10^-31 |

All four regions support the seven-response vector in both evidence scopes. This is a multivariate recurrence claim, not a claim that every atomic coefficient is positive or individually supported.

In the primary all-analysis scope, 26 of 28 regional atomic coefficients are positive. The clearest opposite-sign response is southern shallow/open tube (β = -0.21660, SE = 0.06083, nominal P = 0.000370). Northern-high-latitude generalized form is positive but weak (P = 0.1018), and southern selfing mating system is positive but weak (P = 0.0540).

Complete H1 tables:

- `results/geography_20260924/all/beta_binomial_within_omnibus.csv`
- `results/geography_20260924/all/beta_binomial_within_slopes.csv`
- `results/geography_20260924/direct/beta_binomial_within_omnibus.csv`
- `results/geography_20260924/direct/beta_binomial_within_slopes.csv`

Native and native-nonendemic rows are retained in the same files. Where a stratum does not pass support requirements it is labelled not testable rather than biological null.

---

# Appendix S4. H2 conditional decomposition and raw floral patterns

H2 separates measured reproductive assurance from additional floral response. The species-level `selfing_core` score uses compatibility, mating system and autonomous-selfing capacity only. The `generalized_accessible` score summarizes generalized floral form, actinomorphy and shallow/open tube.

Primary models are:

`selfing_core ~ isolation + area + climate`

`generalized_accessible ~ isolation + selfing_core + area + climate`

`plain_colour ~ isolation + selfing_core + area + climate`

Persistence of an isolation coefficient after adjustment for `selfing_core` is conditional decomposition, not mediation.

### Table S3. Selfing-adjusted accessibility response

| Region | All-analysis β | SE | P | q | Direct-only β | SE | P | q |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Northern mid-latitude | 0.02068 | 0.01281 | 0.1065 | 0.1703 | 0.01715 | 0.01144 | 0.1340 | 0.2143 |
| Northern high latitude | 0.12557 | 0.03390 | 0.000212 | 0.000847 | 0.11657 | 0.03716 | 0.001709 | 0.006835 |
| Tropical | 0.06846 | 0.02494 | 0.006061 | 0.01616 | 0.05249 | 0.02616 | 0.04483 | 0.1196 |
| Southern extratropical | 0.04929 | 0.03979 | 0.2154 | 0.2873 | 0.05779 | 0.03473 | 0.09611 | 0.1922 |

All eight estimates are positive. FDR support is region dependent. In the primary scope it is strongest in northern high latitudes and the tropics. Direct-only northern high latitude is also FDR-supported; tropical Direct-only remains nominally positive but not FDR-supported.

The manuscript therefore uses the bounded statement that accessibility is not statistically exhausted by reproductive assurance, while the strength of the additional floral response is geographically contingent.

H2 additionally evaluates:

- five raw colour states after reproductive-assurance adjustment;
- raw floral form;
- symmetry;
- tube depth;
- colour × form joint prevalence;
- colour × tube-depth joint prevalence;
- architecture conditional on colour.

These analyses do not assign realized pollinator identity. Complete machine-readable outputs are under:

- `results/geography_20260924/all/raw_patterns/`
- `results/geography_20260924/direct/raw_patterns/`
- `results/geography_20260924/all/h2_decomposition_models.csv`
- `results/geography_20260924/direct/h2_decomposition_models.csv`

---

# Appendix S5. H3 experimental pollen limitation

GloPL pollen limitation is the log response ratio of reproduction after supplemental pollen addition versus natural pollen receipt. Positive values indicate improved reproduction after added pollen.

The corrected analysis contains:

- 2,969 experiment-level rows;
- 1,408 publication × coordinate × measurement cells;
- 1,248 unique sites;
- 919 publications.

Each publication contributes total analysis weight one. Uncertainty is publication-cluster robust.

### Table S4. Corrected H3 geographic-isolation models

| Analysis | β distance | SE | Two-sided P | Cells | Publications | Sites |
|---|---:|---:|---:|---:|---:|---:|
| Primary | 0.09191 | 0.03806 | 0.01575 | 1,408 | 919 | 1,248 |
| No-zero-constant sensitivity | 0.09089 | 0.03865 | 0.01869 | 1,375 | 912 | 1,238 |
| Supplemental-only sensitivity | 0.04410 | 0.04351 | 0.31082 | 828 | 470 | 736 |

The primary and no-zero-constant models support a positive isolation association. Supplemental-only remains positive but unsupported. H3 is interpreted as an isolation-associated pollen-delivery constraint, not as direct evidence of declining pollinator abundance or visitation.

Source:

- `results/geography_20260924/h3_original_corrected_comparison.json`

---

# Appendix S6. H4 exact-species functional triangulation

H4 exact-matches species-level H2 scores to GloPL. Matching uses only case normalization and underscore/space normalization; synonym rescue and genus fallback are not used.

The model is:

`pollen limitation ~ H2 score + corrected isolation + geographic context + measurement conditions`

Each publication again has total weight one and uncertainty is publication-cluster robust.

### Table S5. Exact H2-score bridge

| Trait family / score | Analysis | β trait score | SE | Two-sided P | Species | Publications |
|---|---|---:|---:|---:|---:|---:|
| Reproductive assurance / selfing_core | Primary | -0.29830 | 0.10352 | 0.00396 | 455 | 409 |
| Reproductive assurance / selfing_core | No-zero-constant | -0.30078 | 0.10419 | 0.00389 | 453 | 408 |
| Reproductive assurance / selfing_core | Supplemental-only | -0.09953 | 0.09398 | 0.28962 | 283 | 241 |
| Accessibility / generalized_accessible | Primary | -0.29566 | 0.12896 | 0.02187 | 143 | 143 |
| Accessibility / generalized_accessible | No-zero-constant | -0.30366 | 0.12991 | 0.01942 | 142 | 142 |
| Accessibility / generalized_accessible | Supplemental-only | -0.29950 | 0.14562 | 0.03971 | 101 | 98 |

Negative coefficients indicate that stronger expression of the island-associated H2 score is associated with lower current experimental pollen limitation.

### Table S6. Primary atomic-trait H4 sensitivities

| Trait | β trait state | SE | Two-sided P | Species | Publications |
|---|---:|---:|---:|---:|---:|
| Self-compatibility | -0.11885 | 0.08246 | 0.1495 | 499 | 429 |
| Autonomous selfing | -0.44492 | 0.08020 | 2.89 × 10^-8 | 558 | 469 |
| Generalized form | -0.18448 | 0.09174 | 0.04434 | 246 | 246 |
| Actinomorphic symmetry | -0.38076 | 0.08697 | 1.20 × 10^-5 | 582 | 479 |

Selfing mating system and shallow/open tube did not pass their parent support gates for atomic H4 promotion.

Complete source tables:

- `results/geography_20260924/h4_exact_corrected.csv`
- `results/geography_20260924/h4_atomic_corrected.csv`

H4 is explicitly **post-hoc functional triangulation**. It does not identify historical mediation.

---

# Appendix S7. Inferential boundaries and stopped validation

The submission distinguishes four levels of inference.

1. **Assemblage pattern (H1).**  
   Contemporary island-flora composition changes with geographic isolation.

2. **Conditional decomposition (H2).**  
   Measured reproductive assurance does not statistically absorb all floral-accessibility response. This is not causal mediation.

3. **Independent ecological pressure (H3).**  
   Experimental pollen limitation increases with geographic isolation. This is not a direct measure of pollinator abundance or visitation.

4. **Functional compatibility (H4).**  
   Reproductive-assurance and accessibility states associated with isolation are associated with lower current pollen limitation in exact-species overlap. This is post-hoc and does not identify the historical causal sequence.

The unobserved historical edge is:

`past isolation-associated pollen limitation -> selection / sorting / persistence -> present trait composition`

The current analysis does not distinguish species sorting, differential colonization/persistence and within-lineage evolutionary change.

A separate prospective H4 validation effort stopped at the support gate before outcome unblinding because the prespecified minimum sample/publication requirements were not met. It is a design/support result, not a biological null and is not used to promote or refute H4.

---

# Supplementary Data manifest

The following committed files constitute the machine-readable Supplementary Data surface for the current corrected analysis.

## Geography

- `results/geography_20260924/corrected_geography_covariates.csv`
- `results/geography_20260924/gshhg_spherical_distances_all.csv`
- `results/geography_20260924/glopl_corrected_site_distances.csv`
- `results/geography_20260924/spherical_geometry_validation.json`
- `results/geography_20260924/distance_same_coastline_audit.json`

## H1

- `results/geography_20260924/all/beta_binomial_within_omnibus.csv`
- `results/geography_20260924/all/beta_binomial_within_slopes.csv`
- `results/geography_20260924/direct/beta_binomial_within_omnibus.csv`
- `results/geography_20260924/direct/beta_binomial_within_slopes.csv`

## H2

- `results/geography_20260924/all/h2_decomposition_models.csv`
- `results/geography_20260924/direct/h2_decomposition_models.csv`
- `results/geography_20260924/all/raw_patterns/`
- `results/geography_20260924/direct/raw_patterns/`

## H3

- `results/geography_20260924/h3_original_corrected_comparison.json`
- `results/geography_20260924/h3_corrected_measurement_cells.csv.gz`
- `results/geography_20260924/h3_corrected_effect_rows.csv.gz`

## H4

- `results/geography_20260924/h4_exact_corrected.csv`
- `results/geography_20260924/h4_atomic_corrected.csv`
- `results/geography_20260924/h4_ra_corrected_effect_rows.csv.gz`
- `results/geography_20260924/h4_arch_corrected_effect_rows.csv.gz`

---

# Supplementary references

References are shared with the main manuscript. No supplementary-only biological claim depends on an uncited external source.
