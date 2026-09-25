# A recurrent global floral island syndrome separates reproductive assurance from pollinator-facing floral reorganization under increasing pollen limitation

## Abstract

Island floras are expected to experience unreliable reproduction as geographic isolation increases, yet it remains unclear whether reproductive assurance and floral simplification form one recurrent global syndrome, and whether floral change is merely a consequence of increased selfing. We assembled a global island plant database containing 106,295 angiosperm species and 222,688 resolved species-by-trait cells across flower colour, floral structure and reproductive assurance. After a geography audit, the primary submission baseline uses 8,264 island units and source-matched GSHHG 2.3.7 continental coastlines. We tested four linked hypotheses: whether isolation is associated with a recurrent floral/reproductive island syndrome (H1); whether floral change remains after conditioning on reproductive assurance (H2); whether experimental pollen limitation increases with isolation in the independent GloPL database (H3); and whether the same species-level trait states enriched with isolation are associated with lower current pollen limitation (H4).

A seven-response multivariate isolation response was supported in all four predeclared geographic regions and in both evidence scopes, although individual trait responses were not uniformly positive. Reproductive-assurance components generally increased with isolation, and selfing-adjusted floral accessibility remained positive in all four regions, with FDR support in northern high latitudes and the tropics in the primary analysis. Selfing-adjusted flower-colour composition and colour–architecture coupling also changed with isolation, but detailed directions differed among regions. Independent pollen-supplementation experiments showed that pollen limitation increased with geographic isolation (standardized slope = 0.0919 ± 0.0381, two-sided p = 0.0157). In exact-species post-hoc functional comparisons, the H2 reproductive-assurance score was associated with lower current pollen limitation (beta = -0.2983, p = 0.0040), as was the H2 generalized-accessibility score (beta = -0.2957, p = 0.0219).

These results identify recurring functional components of an island floral syndrome while rejecting a compulsory pathway in which floral change is only a downstream consequence of selfing. Island isolation is associated with stronger pollen limitation, whereas the reproductive-assurance and accessibility states enriched with isolation are themselves associated with lower current pollen limitation. The detailed floral trajectory remains region dependent, and the historical causal sequence is not identified by these cross-sectional analyses.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; flower colour; pollination

## Introduction

Oceanic and continental islands impose recurrent challenges to plant reproduction. Geographic separation can reduce connectivity with source regions, constrain repeated immigration and expose colonists to low mate availability and uncertain mutualistic interactions. These conditions motivate a classic expectation of island biology: traits that permit reproduction with fewer external partners should become increasingly advantageous as isolation increases.

Baker's law formalized this logic for mating systems. A colonist capable of uniparental reproduction can establish when mates or compatible pollen donors are rare, creating a route from geographic isolation to reproductive assurance. Consistent with this expectation, self-compatible species are frequently over-represented on islands. Yet reproductive assurance is only one way in which plants may respond to unreliable pollination.

Floral phenotype can also change. Specialized floral architectures restrict access to particular visitors, whereas open, radial and shallow flowers can be used by a broader set of potential pollen vectors. Similarly, flower colour and floral architecture jointly structure interactions with functional pollinator groups. If reliable pollination becomes more difficult under isolation, then two non-exclusive pathways could generate an island floral syndrome. First, greater dependence on self-fertilization could produce a selfing syndrome in which floral investment and specialization change as a downstream consequence of reproductive assurance. Second, pollinator-facing floral traits could reorganize independently of selfing if the benefits of attraction or specialization change when the realized pollinator community becomes less predictable.

Distinguishing these routes is important because superficially similar island phenotypes can arise for different reasons. A more open or generalized flower could accompany increased selfing, but it could also remain advantageous in outcrossing plants if broader visitor access buffers pollen delivery. Likewise, flower-colour composition can shift through changes in reproductive investment, pollinator-facing selection, species sorting or colonization filtering. A global island syndrome therefore should not be assumed to represent a single serial mechanism.

A second problem is that floral-trait patterns alone cannot demonstrate that pollination actually becomes more limiting with isolation. Direct evidence requires experiments that compare natural reproduction with reproduction after supplemental pollen addition. Pollen limitation provides such an outcome: if supplemental pollen increases seed or fruit production, then natural pollen receipt constrains reproduction. A global increase in pollen limitation with island isolation would therefore provide an independent ecological counterpart to a plant-side island syndrome, without requiring a claim that pollinator abundance itself has been measured.

Here we combine a global island plant-trait database with the independent GloPL pollen-supplementation database to test four linked hypotheses. **H1** asks whether geographic isolation is associated with a recurrent floral/reproductive island syndrome spanning reproductive assurance, flower colour and floral accessibility. **H2** asks whether floral change is reducible to a selfing-syndrome pathway or whether floral colour and accessibility responses remain after measured reproductive assurance is conditioned on. **H3** tests whether experimental pollen limitation independently increases with geographic isolation. **H4** asks whether the same species-level reproductive-assurance and accessibility scores used in H2 are associated with lower current pollen limitation. Together, these analyses distinguish global recurrence, conditional pathway decomposition, ecological pressure and functional compatibility.

## Materials and Methods

### Global island plant database

We used a corrected global analysis universe of **8,264 island units**. Candidate units originated from the frozen GSHHG 2.3.7 high-resolution island geometry used by the project. A geography audit showed that one locked unit was exactly the western split component (`0-W`) of the Eurasian continental sibling and was therefore excluded. The corrected broad H1 union contains 4,379 islands. The trait database itself retained its frozen taxonomic and provenance structure, containing 106,295 accepted angiosperm species.

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was assembled source by source from floras, monographs, public trait resources and primary literature. Every accepted record retained source provenance. Species-direct high- and medium-confidence evidence had priority, while validated lower-confidence evidence was used only where direct evidence was unavailable. Missing trait information was retained as missing rather than converted to trait absence, and family-level inference was not used as a general fill rule.

The final frozen trait snapshot contained 222,688 resolved cells of 318,885 possible species-by-axis cells (69.83%) across three raw evidence axes: flower colour, floral structural complexity and reproductive assurance. The broad all-analysis-eligible evidence scope was used as the primary plant analysis, with species-direct high/medium evidence retained as a Direct-only sensitivity.

### Geographic isolation and covariates

The primary geographic exposure is the minimum separation between island and continental coastlines measured on a mean-radius sphere (R = 6371.0088 km). Island and continental geometries were taken from the same GSHHG 2.3.7 high-resolution archive. Continental siblings were reconstructed before distance calculation, artificial dateline split edges were omitted, and distance was calculated as the minimum minor-great-circle arc separation between coastline segments. This is a spherical coastline distance, not a WGS84 ellipsoidal distance.

The geography audit identified 1,113 broad-H1 islands that had been assigned spurious zero distance when GSHHG island polygons were compared with a coarser Natural Earth continental geometry. Under the source-matched coastline metric, all 1,113 have positive distance. This correction is a post-hoc measurement repair selected as the primary submission baseline, not a new prospective confirmation.

Plant models used standardized log distance, log island area and climate PC1-PC4. Spatial dependence was addressed with spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of source separation, connectivity and colonization opportunity rather than as a randomized treatment.

Four predeclared geographic replication strata were analysed separately: northern mid-latitude, northern high-latitude, tropical and southern extratropical islands. These strata were used to assess recurrence of response structure, not to require identical coefficients among regions.

For GloPL, geographic exposure was recomputed from each site's coordinates against the same continental geometry. Sites lying on seeded continental land were assigned zero distance; 996 of 1,248 sites are true continental zeros and remain zero in the corrected baseline.

### H1: recurrent floral/reproductive island syndrome

H1 contained seven atomic responses, coded so that positive isolation coefficients represented the predicted island-syndrome direction:

1. self-compatibility;
2. predominantly or obligately selfing mating system;
3. autonomous or delayed autonomous selfing;
4. plain colour;
5. generalized floral form;
6. actinomorphic symmetry; and
7. shallow or open floral tube.

The responses were grouped into three biological domains: reproductive assurance (responses 1-3), colour (response 4) and accessibility/generalization (responses 5-7). Plain colour contrasted white and green/brown/inconspicuous flowers against yellow/orange, red/pink and blue/purple flowers. It is therefore a colour-composition contrast, not a direct measure of pigment investment or animal-perceived conspicuousness.

Island-level trait counts were modelled using beta-binomial logit models with standardized log isolation distance, island area and climate PC1-PC4. Within each geographic stratum, a seven-dimensional Wald test evaluated the joint isolation-response vector. Individual response coefficients were retained so that joint support could not be misread as uniform support for every trait.

### H2: separating reproductive assurance from additional floral change

We first defined a species-level reproductive-assurance score, `selfing_core`, using self-incompatibility/compatibility, mating system and autonomous-selfing capacity only. Flower colour, form, symmetry, tube depth and flower size were excluded from this score.

The reproductive-assurance response was evaluated as:

`selfing_core ~ isolation + island area + climate`.

We then asked whether floral responses persisted after reproductive assurance was conditioned on:

`plain_colour ~ isolation + selfing_core + island area + climate`

and

`generalized_accessible ~ isolation + selfing_core + island area + climate`.

The `generalized_accessible` score summarizes open/generalized floral form, actinomorphy and shallow/open tube depth. Persistence of an isolation coefficient after adjustment for `selfing_core` was interpreted as evidence that the floral response is not reducible to measured reproductive assurance. This is a conditional decomposition rather than a formal causal mediation analysis.

To avoid assigning plants to named pollinator syndromes from a single composite score, we additionally analysed raw reported colour states and raw architecture. Five colour states were analysed after adjustment for `selfing_core`: white, red/pink, yellow/orange, blue/purple and green/brown/inconspicuous. We then tested raw colour × floral-form and colour × tube-depth combinations both as joint prevalence and as conditional architecture given colour. These analyses evaluate display-architecture reorganization but do not identify realized pollinator identity.

### H3: experimental pollen limitation

We used the version-pinned GloPL database of pollen-supplementation experiments. Pollen limitation was represented by the log response ratio of reproduction after supplemental pollen addition to reproduction under natural pollen receipt; positive values therefore indicate that reproduction increases when pollen is added.

The frozen analysis contained 2,969 experimental rows, 1,408 measurement cells, 1,248 sites and 919 publications. Duplicate measurements were aggregated at the publication-by-coordinate-by-measurement level. Each publication contributed total analysis weight one. The primary model related pollen limitation to standardized corrected log geographic distance while controlling for broad geographic context and experimental measurement conditions. Uncertainty was estimated with publication-cluster-robust covariance.

### H4: functional bridge from H2 traits to current pollen limitation

H4 was an explicitly post-hoc functional triangulation analysis. Species in the H2 trait dataset were exact-matched to GloPL species after only case normalization and underscore/space normalization; synonym rescue and genus fallback were not used.

The primary H4 predictors were the literal Direct-only species-level H2 scores used in the plant analysis: `selfing_core` and `generalized_accessible`. Both are soft-membership scores scaled from 0 to 1, with higher values representing stronger concordance with reproductive assurance or generalized accessibility, respectively. Missing component traits were not coded as zero.

For each H2 score, we fitted:

`pollen limitation ~ H2 trait score + standardized distance + geographic context + measurement conditions`.

Each publication again received total weight one and uncertainty was publication-cluster robust. A negative H2-score coefficient indicates that species expressing the island-associated trait state more strongly tend to experience lower current experimental pollen limitation, conditional on geographic distance and measured study structure.

We also retained atomic-trait sensitivities. Because the global colour response was geographically heterogeneous and lacked a single directional functional prediction, colour was not promoted into the global H4 bridge.

## Results

### H1: a recurrent multivariate island response with region-dependent expression

The seven-response isolation vector remained supported in all four broad geographic strata after the geography correction. Joint FDR-adjusted q-values were **3.216 × 10^-10** in northern mid-latitudes (2,173 islands), **2.433 × 10^-5** in northern high latitudes (411 islands), **3.498 × 10^-7** in the tropics (1,493 islands) and **2.504 × 10^-18** in southern extratropical islands (302 islands). Direct-only joint support also remained in all four regions.

Joint support did not imply a uniformly positive seven-trait syndrome. Twenty-six of 28 broad all-analysis coefficients were positive, but the southern shallow/open-tube coefficient was negative (beta = -0.2166, nominal p = 0.000370). Northern high-latitude generalized form weakened to p = 0.1018, and southern selfing mating system weakened to p = 0.0540. The corrected baseline therefore supports recurring functional components with region- and trait-dependent expression rather than one universal response of every floral trait (Figure 4).

### H2: floral accessibility is not reducible to reproductive assurance

The `selfing_core` isolation slope was positive in all four regions in the primary analysis, although nominal support remained concentrated in the tropical and southern strata.

After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. FDR-adjusted q-values were 0.1703 in northern mid-latitudes, 0.000847 in northern high latitudes, 0.01616 in the tropics and 0.2873 in southern extratropical islands. Thus the primary evidence for a selfing-adjusted accessibility response is strongest in northern high latitudes and the tropics. In the Direct-only sensitivity, the tropical accessibility estimate remained nominally positive (p = 0.0448) but no longer survived FDR correction (q = 0.1196); it is therefore not treated as a replicated FDR-supported result.

Colour responses remained more heterogeneous. Southern adjusted plain colour remained supported (q = 0.000847). Raw colour, architecture and colour-conditioned architecture analyses showed region-specific responses. Northern-high-latitude blue/purple architecture associations remained negative, while tropical Direct evidence retained a positive yellow/orange × butterfly/deep-tube coupling (q = 0.03949). These descriptive trait combinations do not establish realized pollinator identity.

The corrected H2 results therefore support a partially separable floral response: accessibility cannot be reduced to measured reproductive assurance in all regions, but the strength and detailed form of the residual floral response are geographically contingent (Figure 5).

### H3: pollen limitation increases with isolation

Across 2,969 pollen-supplementation experiments, the corrected standardized distance coefficient was positive (**beta = 0.09191, SE = 0.03806; two-sided p = 0.01575; one-sided positive p = 0.00787**). The no-zero-constant sensitivity also retained a supported positive coefficient (beta = 0.09089, p = 0.01869).

The supplemental-only sensitivity was positive but unsupported (beta = 0.04410, p = 0.31082). The primary result therefore supports an isolation-associated increase in experimental pollen limitation, while the evidence is not equally strong under every measurement definition (Figure 6A).

### H4: island-associated response traits are linked to lower current pollen limitation

At the atomic-trait level, autonomous selfing provided the strongest functional association. Species with autonomous or delayed selfing had lower current pollen limitation after adjustment for distance, broad context and measurement structure (beta = -0.44492, p = 2.89 × 10^-8). Actinomorphy was also strongly negative (beta = -0.38076, p = 1.20 × 10^-5), generalized form was weaker but supported (beta = -0.18448, p = 0.04434), and self-compatibility was negative but imprecise (beta = -0.11885, p = 0.1495).

The literal H2-to-H4 bridge gave the same overall interpretation. The Direct-only `selfing_core` score matched 455 GloPL species across 409 publications and was negatively associated with current pollen limitation (**beta = -0.29830, SE = 0.10352, p = 0.00396**). The association remained supported in the no-zero-constant sensitivity, whereas the supplemental-only subset remained negative but unsupported (p = 0.28962).

The Direct-only `generalized_accessible` score matched 143 GloPL species across 143 publications. Higher accessibility was associated with lower current pollen limitation (**beta = -0.29566, SE = 0.12896, p = 0.02187**). The supplemental-only sensitivity was also negative (p = 0.03971), and the no-zero-constant sensitivity remained supported.

Thus the same two response families enriched along the island-isolation gradient are associated, in exact-species overlap, with lower current pollen limitation (Figure 6B-C). This functional alignment does not establish that historical pollen limitation caused the observed trait distributions.

## Discussion

### A recurrent functional core with non-uniform floral realization

The corrected baseline identifies a recurrent multivariate floral/reproductive response across four geographically distinct island regions while making the non-uniformity of individual traits explicit. Reproductive assurance and floral accessibility/generalization provide the clearest recurring functional components, whereas colour and specific architecture traits vary more strongly among regions. The negative southern shallow/open-tube response is especially important: the island syndrome is a multivariate tendency, not a rule that every floral component must move in the same direction everywhere.

This distinction separates recurrence at the functional level from uniformity at the phenotype level. Islands can repeatedly favour, sort or retain plants that are less dependent on precise external pollen delivery without converging on one globally identical flower phenotype. Similar broad ecological constraints can therefore generate a common functional core alongside geographically contingent display and architecture.

### Floral change is not explained by reproductive assurance alone

The selfing-syndrome route remains an important explanation. Isolation was associated with greater reproductive assurance, consistent with Baker-type establishment and reproductive-assurance processes. Selfing can reduce dependence on pollen vectors and can be accompanied by changes in floral investment.

Yet the H2 decomposition shows that this is not the whole story. In the primary analysis, floral accessibility remained positively associated with isolation after measured reproductive assurance was controlled, with strongest support in northern high latitudes and the tropics. Therefore the observed accessibility response cannot be treated simply as a statistical by-product of greater selfing.

This residual association is compatible with an additional pollinator-facing route, but it does not prove that route causally. When reliable interaction with a restricted set of visitors becomes less dependable, the relative value of accessible floral architecture could change even without a shift in mating system. The region-specific colour × architecture results are consistent with such context dependence while also showing why one universal named pollination syndrome is too restrictive.

### Pollen limitation provides an independent ecological counterpart

The corrected GloPL analysis supplies a separate empirical layer: experimental pollen limitation increases with geographic isolation. This matters because plant trait patterns alone cannot distinguish adaptation, species sorting, colonization filtering or other historical processes from changes in pollination service.

Pollen limitation should not be equated with pollinator abundance. It integrates the reproductive consequences of pollen quantity, pollen quality, visitation, mate availability and other processes affecting natural pollen receipt. Our result therefore supports an isolation-associated pollination constraint rather than a global decline of one pollinator group.

### Functional compatibility links the syndrome to current pollen limitation

The H4 bridge asks whether the traits enriched with isolation are functionally aligned with lower dependence on natural pollen receipt. Both literal H2 species scores showed the predicted negative association with current pollen limitation. Reproductive assurance provides the most direct interpretation: autonomous selfing can permit reproduction when external pollen delivery is insufficient. Floral accessibility is less mechanistically specific, but open and generalized floral architectures may broaden the set of visitors capable of contacting reproductive organs and transferring pollen.

These associations connect the global trait pattern to a plausible functional consequence. They do not identify historical mediation. The data do not observe an ancestral state, an increase in pollen limitation, subsequent selection and later trait change within the same lineages. H4 should therefore be interpreted as post-hoc functional triangulation rather than proof of the evolutionary sequence.

### Why detailed floral responses differ among regions

The recurrent functional core coexists with strongly context-dependent colour and colour–architecture responses. This non-uniformity is biologically informative rather than a failure of the island-syndrome hypothesis. The selective or filtering value of a floral phenotype depends on the plant's starting state and on the realized interaction community. Similar broad constraints can therefore generate different detailed phenotypic trajectories across regions.

This result creates a natural next question: why does the same broad island-like change produce different pollinator-facing response branches? Addressing that question requires models or local systems that explicitly represent starting plant state, realized pollinator composition and their interaction, rather than treating a coarse island gradient as sufficient to determine one phenotype.

### Limitations

The primary plant responses describe contemporary island-flora composition. They do not by themselves distinguish species sorting, colonization filtering and within-lineage evolutionary change. Geographic distance is a composite exposure that covaries with connectivity and source supply. Although plant models adjust for island area and climate, residual confounding remains possible.

The corrected distance metric is itself a model of geography. It uses source-matched GSHHG shorelines and exact spherical arc minimization, but it assumes a mean-radius sphere and inherits positional limits from the underlying coastline data. Numerical precision should not be interpreted as metre-scale geographic accuracy. The geography correction was selected post hoc after an audit identified the original source mismatch; it is the primary measurement baseline because it repairs a known exposure error, not because it constitutes new prospective confirmation.

The H2 conditional analyses do not constitute causal mediation. A floral association that remains after adjustment for `selfing_core` shows that the pattern is not statistically reducible to measured reproductive assurance, but it does not prove direct selection by pollinators. Likewise, raw colour and architecture do not identify realized visitor identity.

Finally, H4 is post-hoc functional triangulation. The exact H2 species scores were tested against an already inspected GloPL outcome base. A separate prospective validation audit was support-limited before outcome unblinding and is not part of the H4 result. Stronger causal inference will require lineage-resolved or longitudinal observations linking pollination service, reproductive success and trait change through time.

## Conclusion

Geographic isolation is associated with a recurrent multivariate floral/reproductive island response characterized by recurring reproductive-assurance and floral-accessibility components, but not by uniform change in every floral trait. Floral accessibility remains associated with isolation after reproductive assurance is accounted for in the primary analysis, showing that the floral response is not reducible to a selfing syndrome alone. Independent pollen-supplementation experiments show that pollen limitation increases with isolation, while the same reproductive-assurance and accessibility scores enriched with isolation are associated with lower current pollen limitation. Together, these results support a recurrent functional response to island isolation while leaving the detailed floral trajectory contingent on ecological context and the historical causal pathway unresolved.

## References

Baker, H. G. 1955. Self-compatibility and establishment after long-distance dispersal. *Evolution* 9:347-349.

Bennett, J. M. et al. 2018. GloPL, a global data base on pollen limitation of plant reproduction. *Scientific Data* 5:180249.

Fenster, C. B. et al. 2004. Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics* 35:375-403.

Grossenbacher, D. L. et al. 2017. Self-compatibility is over-represented on islands. *New Phytologist* 215:469-478.

Rosas-Guerrero, V. et al. 2014. A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters* 17:388-400.

Sicard, A. & Lenhard, M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany* 107:1433-1443.
