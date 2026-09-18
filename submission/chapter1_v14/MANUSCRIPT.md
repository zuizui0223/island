# A recurrent global floral island syndrome separates reproductive assurance from pollinator-facing floral reorganization under increasing pollen limitation

## Abstract

Island floras are expected to experience unreliable reproduction as geographic isolation increases, yet it remains unclear whether reproductive assurance and floral simplification form one recurrent global syndrome, and whether floral change is merely a consequence of increased selfing. We assembled a global island plant database spanning 8,265 islands and 106,295 angiosperm species, with 222,688 resolved species-by-trait cells across flower colour, floral structure and reproductive assurance. We tested four linked hypotheses. First, we asked whether isolation is associated with a recurrent floral/reproductive island syndrome. Second, we separated a reproductive-assurance route from floral changes that remain after conditioning on a strict selfing score. Third, we tested whether experimental pollen limitation increases with isolation using the independent GloPL database. Finally, we asked whether the same species-level trait scores that increase with isolation are associated with lower current pollen limitation.

A seven-response multivariate island-syndrome vector was supported in all four predeclared geographic regions and in both evidence scopes. Reproductive assurance increased with isolation, and generalized floral accessibility remained positively associated with isolation after controlling for selfing. Selfing-adjusted flower-colour composition and colour–architecture coupling also changed with isolation, but the detailed directions differed among regions. Independent pollen-supplementation experiments showed that pollen limitation increased with geographic isolation (standardized slope = 0.0794 ± 0.0377, two-sided p = 0.0354). In exact-species post-hoc functional comparisons, the H2 selfing score was associated with lower current pollen limitation (beta = -0.2971, p = 0.0041), as was the H2 generalized-accessibility score (beta = -0.2960, p = 0.0222).

These results identify a recurrent global floral island-syndrome direction while rejecting a compulsory pathway in which floral change is only a downstream consequence of selfing. Island isolation is associated with stronger pollen limitation, whereas the reproductive-assurance and accessibility states enriched with isolation are themselves associated with lower current pollen limitation. The evidence supports functional convergence on reproductive assurance and floral accessibility while leaving the historical causal sequence and realized pollinator mechanisms unresolved.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; flower colour; pollination

## Introduction

Oceanic and continental islands impose recurrent challenges to plant reproduction. Geographic separation can reduce connectivity with source regions, constrain repeated immigration and expose colonists to low mate availability and uncertain mutualistic interactions. These conditions motivated one of the classic expectations of island biology: traits that permit reproduction with fewer external partners should become increasingly advantageous as isolation increases.

Baker's law formalized this logic for mating systems. A colonist capable of uniparental reproduction can establish when mates or compatible pollen donors are rare, creating a route from geographic isolation to reproductive assurance. Consistent with this expectation, self-compatible species are frequently over-represented on islands. Yet reproductive assurance is only one way in which plants may respond to unreliable pollination.

Floral phenotype can also change. Specialized floral architectures restrict access to particular visitors, whereas open, radial and shallow flowers can be used by a broader set of potential pollen vectors. Similarly, flower colour and floral architecture jointly structure interactions with functional pollinator groups. If reliable pollination becomes more difficult under isolation, then two non-exclusive pathways could generate an island floral syndrome. First, greater dependence on self-fertilization could produce a selfing syndrome in which floral investment and specialization decline as a downstream consequence of reproductive assurance. Second, pollinator-facing floral traits could reorganize independently of selfing if the benefits of attraction or specialization change when the realized pollinator community becomes less predictable.

Distinguishing these routes is important because superficially similar island phenotypes can arise for different reasons. A more open or generalized flower could accompany increased selfing, but it could also remain advantageous in outcrossing plants if broader visitor access buffers pollen delivery. Likewise, a shift toward less conspicuous flower colours could accompany reduced investment under selfing, but colour can also respond directly to changes in the composition of effective pollinators. Therefore a global island syndrome should not be assumed to represent a single serial mechanism.

A second problem is that floral-trait patterns alone cannot demonstrate that pollination actually becomes more limiting with isolation. Direct evidence requires experiments that compare natural reproduction with reproduction after supplemental pollen addition. Pollen limitation provides such an outcome: if supplemental pollen increases seed or fruit production, then natural pollen receipt constrains reproduction. A global increase in pollen limitation with island isolation would therefore provide an independent ecological counterpart to a plant-side island syndrome, without requiring a claim that pollinator abundance itself has been measured.

Here we combine a global island plant-trait database with the independent GloPL pollen-supplementation database to test four linked hypotheses. **H1** asks whether geographic isolation is associated with a recurrent floral/reproductive island syndrome spanning reproductive assurance, flower colour and floral accessibility. **H2** asks whether floral change is reducible to a selfing-syndrome pathway or whether pollinator-facing colour and accessibility responses remain after measured reproductive assurance is conditioned on. **H3** tests whether experimental pollen limitation independently increases with geographic isolation. **H4** asks whether the same species-level reproductive-assurance and accessibility scores used in H2 are associated with lower current pollen limitation. Together, these analyses distinguish global recurrence, pathway decomposition, ecological pressure and functional compatibility.

## Materials and Methods

### Global island plant database

We constructed a fixed global island universe of 8,265 islands using GSHHG-based island geometry and associated geographic covariates. Observed island floras were assembled from GBIF occurrence records assigned to islands. The final analysis universe contained 106,295 accepted angiosperm species.

Trait evidence was assembled source by source from floras, monographs, public trait resources and primary literature. Every accepted record retained source provenance. Species-direct high- and medium-confidence evidence had priority, while validated lower-confidence evidence was used only where direct evidence was unavailable. Missing trait information was retained as missing rather than converted to trait absence, and family-level inference was not used as a general fill rule.

The final frozen snapshot contained 222,688 resolved cells of 318,885 possible species-by-axis cells (69.83%) across three raw evidence axes: flower colour, floral structural complexity and reproductive assurance. The broad all-analysis-eligible evidence scope was used as the primary plant analysis, with species-direct high/medium evidence retained as a Direct-only sensitivity.

### Geographic isolation and covariates

The primary plant exposure was log-transformed distance to the nearest continental source boundary. Models included log island area and four climate principal components. Spatial dependence was addressed with spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of source separation, connectivity and colonization opportunity rather than as a randomized treatment.

Four predeclared geographic replication strata were analysed separately: northern mid-latitude, northern high-latitude, tropical and southern extratropical islands. These strata were used to assess recurrence of the direction of response, not to require identical coefficients among regions.

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

Island-level trait counts were modelled using beta-binomial logit models with standardized log isolation distance, island area and climate PC1-PC4. Within each geographic stratum, a seven-dimensional Wald test evaluated the joint isolation-response vector. For descriptive comparison among domains, standardized coefficients were averaged within each biological domain and then the three domain means were given equal weight.

### H2: separating selfing from pollinator-facing floral change

We first defined a species-level reproductive mechanism score, `selfing_core`, using self-incompatibility/compatibility, mating system and autonomous-selfing capacity only. Flower colour, form, symmetry, tube depth and flower size were excluded from this score.

The reproductive-assurance route was evaluated as:

`selfing_core ~ isolation + island area + climate`.

We then asked whether floral responses persisted after reproductive assurance was conditioned on:

`plain_colour ~ isolation + selfing_core + island area + climate`

and

`generalized_accessible ~ isolation + selfing_core + island area + climate`.

The `generalized_accessible` score summarizes open/generalized floral form, actinomorphy and shallow/open tube depth. Persistence of an isolation coefficient after adjustment for `selfing_core` was interpreted as evidence that the floral response is not reducible to measured reproductive assurance. This is a conditional decomposition rather than a formal causal mediation analysis.

To avoid assigning plants to named pollinator syndromes from a single composite score, we additionally analysed raw reported colour states and raw architecture. Five colour states were analysed after adjustment for `selfing_core`: white, red/pink, yellow/orange, blue/purple and green/brown/inconspicuous. We then tested raw colour × floral-form and colour × tube-depth combinations both as joint prevalence and as conditional architecture given colour. These analyses evaluate display-architecture reorganization but do not identify realized pollinator identity.

### H3: experimental pollen limitation

We used the version-pinned GloPL database of pollen-supplementation experiments. Pollen limitation was represented by the log response ratio of reproduction after supplemental pollen addition to reproduction under natural pollen receipt; positive values therefore indicate that reproduction increases when pollen is added.

The frozen analysis contained 2,969 experimental rows from 1,248 sites and 919 publications. Duplicate measurements were aggregated at the publication-by-coordinate-by-measurement level. Each publication contributed total analysis weight one. The primary model related pollen limitation to standardized log geographic distance while controlling for broad geographic context and experimental measurement conditions. Uncertainty was estimated with publication-cluster-robust covariance.

### H4: functional bridge from H2 traits to current pollen limitation

H4 was an explicitly post-hoc functional triangulation analysis. Species in the H2 trait dataset were exact-matched to GloPL species after only case normalization and underscore/space normalization; synonym rescue and genus fallback were not used.

The primary H4 predictors were the literal Direct-only species-level H2 scores used in the plant analysis: `selfing_core` and `generalized_accessible`. Both are soft-membership scores scaled from 0 to 1, with higher values representing stronger concordance with reproductive assurance or generalized accessibility, respectively. Missing component traits were not coded as zero.

For each H2 score, we fitted:

`pollen limitation ~ H2 trait score + standardized distance + geographic context + measurement conditions`.

Each publication again received total weight one and uncertainty was publication-cluster robust. A negative H2-score coefficient indicates that species expressing the island-associated trait state more strongly tend to experience lower current experimental pollen limitation, conditional on geographic distance and measured study structure.

We also retained atomic-trait and reconstructed-family sensitivities. Because the global colour response was geographically heterogeneous and lacked a single directional functional prediction, colour was not promoted into the global H4 bridge.

## Results

### H1: a recurrent multivariate island syndrome

The seven-response isolation vector was supported in all four geographic strata in the primary all-analysis scope and in the Direct-only sensitivity. Joint false-discovery-rate adjusted q-values were 2.36 × 10^-6 in northern mid-latitudes, 2.61 × 10^-6 in northern high latitudes, 3.78 × 10^-7 in the tropics and 1.60 × 10^-12 in southern extratropical islands. Direct-only evidence retained joint support in all four strata.

The equal-weight three-domain orientation was positive in every region. Reproductive assurance showed positive mean isolation responses in all four strata, as did floral accessibility/generalization. The colour component was less uniform: the plain-colour coefficient was approximately zero in northern mid-latitudes but positive in northern high-latitude, tropical and southern extratropical strata. Thus the global signal is a recurrent multivariate syndrome direction rather than an identical response of every component in every region.

### H2: floral accessibility persists beyond measured selfing

The `selfing_core` isolation coefficient was positive in all four geographic strata under both evidence scopes, consistent with increasing reproductive assurance along the isolation gradient.

However, the floral response was not absorbed by reproductive assurance. After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. Support was strongest in northern high latitudes (beta = 0.1248, q = 0.00144) and the tropics (beta = 0.0733, q = 0.00494). Direct-only estimates retained the same positive direction in all four strata, with northern high latitudes supported after correction and the tropical estimate close to the correction boundary.

Selfing-adjusted colour responses were more heterogeneous. The raw five-colour vector changed with isolation in northern mid-latitude, tropical and southern extratropical floras, but not in northern high latitudes. Northern mid-latitudes showed a replicated decline in red/pink flowers. In northern high latitudes, where the overall five-colour vector was not supported, blue/purple flowers nevertheless showed a replicated decline in specialized/deep architecture. Tropical Direct evidence instead showed increased yellow/orange coupling to a deep-tube component, while southern extratropical flowers showed mixed yellow/orange restructuring.

These results reject a compulsory model in which all floral change is only a downstream consequence of selfing. Floral accessibility forms an additional response dimension, while colour and architecture reorganize in context-dependent ways.

### H3: pollen limitation increases with isolation

Across 2,969 pollen-supplementation experiments, the standardized distance coefficient was positive (beta = 0.07937, SE = 0.03773; two-sided p = 0.03543; one-sided positive p = 0.01772). The no-zero-constant sensitivity retained a supported positive coefficient, whereas the supplemental-only estimate was smaller but remained positive.

A post-hoc shape audit indicated that the pattern was not attributable solely to a mainland-versus-offshore contrast. The offshore-only distance coefficient was also positive (beta = 0.23545, p = 0.00587), supporting a continuing gradient among sampled offshore sites.

### H4: island-associated response traits are linked to lower current pollen limitation

At the atomic-trait level, autonomous selfing provided the strongest functional bridge. Species with autonomous or delayed selfing had lower current pollen limitation after adjustment for distance, broad context and measurement structure (beta = -0.44672, p = 2.60 × 10^-8). Actinomorphy and generalized floral form were directionally concordant, although support and robustness were weaker.

The direct H2-to-H4 bridge gave the same overall interpretation. The literal Direct-only `selfing_core` score matched 455 GloPL species across 409 publications and was negatively associated with current pollen limitation (beta = -0.29706, SE = 0.10353, p = 0.00411). The association remained supported in the no-zero-constant sensitivity, whereas the supplemental-only subset remained negative but was not supported.

The literal Direct-only `generalized_accessible` score matched 143 GloPL species across 143 publications. Higher accessibility was associated with lower current pollen limitation (beta = -0.29601, SE = 0.12940, p = 0.02216). The association remained negative and supported in both the supplemental-only and no-zero-constant sensitivities.

Thus the same two response families that increase with island isolation are associated, in exact-species overlap, with lower current pollen limitation. This functional alignment does not establish that historical pollen limitation caused the observed trait distributions.

## Discussion

### A recurrent functional core with non-uniform floral realization

Our results identify a recurrent floral/reproductive island-syndrome direction across four geographically distinct island regions. The strongest recurring components are reproductive assurance and floral accessibility/generalization. The colour component contributes to the multivariate syndrome but does not respond uniformly across regions.

This combination matters because it separates recurrence at the functional level from uniformity at the phenotype level. Islands can repeatedly favour, sort or retain plants that are less dependent on precise external pollen delivery without converging on one globally identical flower phenotype. The same broad ecological problem can therefore produce a common functional core alongside geographically contingent display and architecture.

### Floral simplification is not explained by selfing alone

The selfing-syndrome route remains an important explanation. Isolation was associated with greater reproductive assurance, consistent with Baker-type establishment and reproductive-assurance processes. Selfing can reduce dependence on pollen vectors and can be accompanied by reductions in floral investment.

Yet the H2 decomposition shows that this is not the whole story. Floral accessibility remained positively associated with isolation after measured reproductive assurance was controlled. Therefore the tendency toward open, radial or shallow floral architecture cannot be treated simply as a statistical by-product of greater selfing in the observed data.

A second pollinator-facing route is therefore required at the level of explanation. When reliable interaction with a restricted set of visitors becomes less dependable, the relative value of accessible floral architecture may increase even without a shift in mating system. The region-specific colour × architecture results reinforce this interpretation while also showing why one universal named pollination syndrome is too restrictive.

### Pollen limitation provides an independent ecological counterpart

The GloPL analysis supplies a separate empirical layer: experimental pollen limitation increases with geographic isolation. This is important because the plant trait analysis alone cannot distinguish adaptation, sorting, colonization filtering or other historical processes from a direct change in pollination service.

Pollen limitation should not be equated with pollinator abundance. It integrates the reproductive consequences of pollen quantity, pollen quality, visitation, mate availability and other processes affecting natural pollen receipt. Our result therefore supports an isolation-associated pollination constraint rather than a global decline of one pollinator group.

### Functional compatibility links the syndrome to current pollen limitation

The H4 bridge asks whether the traits enriched with isolation are at least functionally aligned with reduced dependence on natural pollen receipt. Both literal H2 species scores showed the predicted negative association with current pollen limitation. Reproductive assurance provides the most direct interpretation: autonomous selfing can permit reproduction when external pollen delivery is insufficient. Floral accessibility is less mechanistically specific, but open and generalized floral architectures may broaden the set of visitors capable of contacting reproductive organs and transferring pollen.

These associations help connect the global trait pattern to a plausible functional consequence. They do not, however, identify historical mediation. The data do not observe an ancestral state, an increase in pollen limitation, subsequent selection and later trait change within the same lineages. H4 should therefore be interpreted as functional triangulation rather than proof of the evolutionary sequence.

### Why detailed floral responses differ among regions

The recurrent functional core coexists with strongly context-dependent colour and colour–architecture responses. This non-uniformity is biologically informative rather than a failure of the island-syndrome hypothesis. The selective or filtering value of a floral phenotype depends on the plant's starting state and on the realized interaction community. Similar broad constraints can therefore generate different detailed phenotypic trajectories across regions.

This result creates a natural next question: why does the same broad island-like change produce different pollinator-facing response branches? Addressing that question requires models or local systems that explicitly represent starting plant state, realized pollinator composition and their interaction, rather than treating a coarse island gradient as sufficient to determine one phenotype.

### Limitations

The primary plant responses describe contemporary island-flora composition. They do not by themselves distinguish species sorting, colonization filtering and within-lineage evolutionary change. Geographic distance is also a composite exposure that covaries with connectivity and source supply. Although plant models adjust for island area and climate, residual confounding remains possible.

The H2 conditional analyses do not constitute causal mediation. A floral association that remains after adjustment for `selfing_core` shows that the pattern is not statistically reducible to measured reproductive assurance, but it does not prove direct selection by pollinators. Likewise, raw colour and architecture do not identify realized visitor identity.

Finally, H4 is post-hoc functional triangulation. The exact H2 species scores were tested against an already inspected GloPL outcome base. A separate prospective validation audit was support-limited before outcome unblinding and is not part of the H4 result. Stronger causal inference will require lineage-resolved or longitudinal observations linking pollination service, reproductive success and trait change through time.

## Conclusion

Geographic isolation is associated with a recurrent global floral/reproductive island syndrome characterized by increasing reproductive assurance and floral accessibility, with a more heterogeneous colour component. Floral accessibility remains associated with isolation after controlling for reproductive assurance, showing that island floral change is not reducible to a selfing syndrome alone. Independent pollen-supplementation experiments show that pollen limitation increases with isolation, while the same reproductive-assurance and accessibility scores enriched with isolation are associated with lower current pollen limitation. Together, these results identify a recurrent functional response to island isolation while leaving the detailed floral trajectory contingent on ecological context and the historical causal pathway unresolved.

## References

Baker, H. G. 1955. Self-compatibility and establishment after long-distance dispersal. *Evolution* 9:347-349.

Bennett, J. M. et al. 2018. GloPL, a global data base on pollen limitation of plant reproduction. *Scientific Data* 5:180249.

Fenster, C. B. et al. 2004. Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics* 35:375-403.

Grossenbacher, D. L. et al. 2017. Self-compatibility is over-represented on islands. *New Phytologist* 215:469-478.

Rosas-Guerrero, V. et al. 2014. A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters* 17:388-400.

Sicard, A. & Lenhard, M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany* 107:1433-1443.
