# Floral-space geometry and hypervolume feasibility

**Decision, 6 November 2026:** use the common square-root-distance PCoA for visualization and retain discrete recorded-profile diversity and species support as the first analysis. No Gaussian hypervolume estimates or tests of total floral-space volume were produced.

This is an exploratory diagnostic. The choice of distance transformation followed inspection of the observed distance geometry; it was not a prospective model selection rule.

## 1. Inputs and counting units

The common space uses flower colour and floral form. A profile is the pair of recorded state sets. A taxon with several reported colours or forms remains one profile; the Cartesian product is never counted as several observed flowers. Counts below are retained taxon-name units, including unchanged hybrid formula keys, rather than a new claim of fully reconciled biological species identities. The analysis-specific taxonomic gate and exact input hashes are retained with the results.

| Scope | Complete name units | Complete units in corrected flora | Profiles in scope | Islands with complete profiles | Islands with ≥50 complete units | Maximum profiles on one island |
|---|---:|---:|---:|---:|---:|---:|
| All | 16,739 | 16,725 | 496 | 3,858 | 1,532 | 301 |
| Direct | 8,455 | 8,452 | 467 | 3,779 | 1,442 | 284 |

All and Direct share the same 496-profile reference coordinates. Fitting occurs once on unique All profiles, then name units are mapped to coordinates by profile ID. Identical profiles therefore have identical coordinates regardless of their number of carriers. The All table has 7,222 complete units with exactly one colour and one form; multistate records are retained as recorded sets.

## 2. What the ordination preserves

For each trait, dissimilarity is one minus intersection size divided by union size. The colour and form dissimilarities are then averaged with equal weight. Single-state inputs reduce to nominal Gower mismatch. The primary PCoA embeds the square root of this distance. This is an explicit transformation of the geometry, not additional trait evidence. Ordered depth/size traits are not part of this distance definition.

The original distance matrix has 23.7935% negative inertia, defined as the absolute sum of negative eigenvalues divided by the absolute sum of all eigenvalues. The square-root matrix has no material negative inertia and 75 positive axes. Lingoes also removes negative eigenvalues, but its eigenvalue shift c=1.590764 adds 2c=3.181528 to squared off-diagonal distances and leaves 494 positive axes. That definition follows the [official ape PCoA documentation](https://search.r-project.org/CRAN/refmans/ape/html/pcoa.html); our stored c is half the squared-distance addition.

| Retained common axes | Cumulative positive inertia | RMSE against original distance | RMSE against square-root distance |
|---|---:|---:|---:|
| 2 | 15.35% | 0.50154 | 0.58363 |
| 3 | 22.57% | 0.42233 | 0.50355 |
| 4 | 29.25% | 0.35954 | 0.44037 |
| 5 | 35.73% | 0.30741 | 0.38663 |
| 6 | 42.11% | 0.26045 | 0.33861 |
| 17 | 80.49% | 0.06652 | 0.10226 |
| 34 | 95.38% | 0.10121 | 0.02992 |
| 75 | 100.00% | 0.11810 | 0.00000 |

Two displayed axes retain 15.35% of positive inertia. Reaching 80% requires 17 axes, and reaching 95% requires 34 axes. Two-dimensional plots are therefore projections, not estimates of the entire profile space. The minimum original-distance RMSE occurs at 17 axes (0.06652); this is a diagnostic and is not by itself a dimension-selection rule. Both original-distance and transformed-distance error are retained because the transformation changes the target.

Functional-space dimensionality should be assessed from distance fidelity, rather than fixed at the number of axes convenient for plotting; poor representations can alter diversity conclusions. This application follows that principle from [Maire et al. (2015)](https://doi.org/10.1111/geb.12299). Eigenvalue percentages alone do not replace the reported distance errors.

## 3. Dimensional fidelity versus island support

[Blonder et al. (2018), Table 2](https://doi.org/10.1111/2041-210X.12865) proposes the KDE screening heuristic ln(m)>d, with m input points and d dimensions. Here m counts unique profiles because the proposed endpoint is profile breadth, with one point per type. This conservative project-specific use of the heuristic is not a theorem proving that any fitted density would be impossible or invalid.

| Dimensions | Minimum unique profiles | Common inertia retained | All islands meeting rule | All meeting rule and ≥50 complete units | Direct meeting rule and ≥50 complete units |
|---|---:|---:|---:|---:|---:|
| 2 | 8 | 15.35% | 2,601 | 1,532 | 1,442 |
| 3 | 21 | 22.57% | 1,977 | 1,532 | 1,441 |
| 4 | 55 | 29.25% | 1,136 | 1,136 | 1,085 |
| 5 | 149 | 35.73% | 79 | 79 | 74 |
| 6 | 404 | 42.11% | 0 | 0 | 0 |

The dimensions supported by many islands omit most of the common geometry. At six dimensions no island meets the unique-profile screening rule, while even that projection retains less than half of the common inertia. The present analysis therefore does not interpret a low-dimensional Gaussian hypervolume as total floral-profile breadth.

## 4. Meaning of a future limited projection analysis

A deliberately limited analysis could estimate a probability-content region of an equally weighted Gaussian mixture over unique profiles in a fixed common projection. Its volume would describe a smoothed distribution in those named axes. It would not count observed combinations or establish the extent of the full floral space. PCoA does not make the spaces between categorical states into observed or biologically feasible flowers.

Such a comparison would need fixed axes and scaling, comparable bandwidths, a declared probability threshold, sampling controls, and explicit projection-loss and coverage reports. Giving repeated carrier species additional KDE weight would instead measure how species are allocated among profiles. Blonder et al. discuss these bandwidth and sample-size comparability requirements and caution against applying hypervolume algorithms directly to unordered categorical inputs; their discussion of ordination does not remove the biological discreteness of the original observations.

## 5. Inspectable evidence

- [Machine-readable feasibility record](hypervolume_feasibility.json) contains input hashes, both transformations, scope counts, the sampling screen and the exact decision.
- [Primary ordination diagnostics](all/ordination_diagnostics.csv) and [Lingoes comparison](all/lingoes_ordination_diagnostics.csv) preserve all tested dimensions.
- [Common reference coordinates](all/profile_coordinates.csv), [All island metrics](all/island_metrics.csv) and [Direct island metrics](direct/island_metrics.csv) retain the values behind the tables.
- [Taxonomic exclusions](all/excluded_taxon_names.csv) and [All input manifest](all/input_manifest.json) identify the analysis population and frozen inputs.

The feasibility record was rebuilt from the final taxonomy-gated outputs. Every saved primary diagnostic value was reproduced from the final profiles with `island_v2.floral_space_geometry`. The full coordinate-distance reconstruction was additionally checked against the square-root dissimilarity matrix using SciPy's pairwise Euclidean-distance calculation.
