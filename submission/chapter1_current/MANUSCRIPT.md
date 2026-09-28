# A recurrent global floral island syndrome extends beyond the selfing syndrome

## Abstract

Geographic isolation can make plant reproduction unreliable, but whether island floral change is simply a by-product of selfing is unresolved. Across 8,264 island units and 106,295 angiosperm species, a seven-trait floral–reproductive response to isolation recurred in four regions despite non-uniform individual traits. After conditioning on reproductive assurance, floral accessibility remained positively associated with isolation in all regions and was FDR-supported in northern high latitudes and the tropics. Independently, 2,969 pollen-supplementation experiments showed increasing pollen limitation with isolation (β = 0.0919, P = 0.0157). In exact-species post-hoc comparisons, reproductive-assurance and accessibility scores enriched with isolation were each associated with lower current pollen limitation. Island isolation therefore aligns stronger pollen limitation with a recurrent functional response that is not reducible to a serial selfing syndrome; reproductive assurance and pollinator-facing accessibility are partially separable while detailed floral trajectories remain geographically contingent.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; flower colour; pollination

## Introduction

Islands are natural tests of how dispersal, establishment, biotic interactions and evolution reshape ecological strategies. Functional island biogeography predicts that insularity can repeatedly filter traits linked to colonization and persistence, while the detailed outcome depends on the regional species pool and interaction environment (Schrader et al. 2021). Plant reproduction is especially exposed to these filters because successful establishment may require both compatible mates and effective pollen transfer.

Baker's law formalized one solution: colonists capable of uniparental reproduction can establish when mates or compatible pollen donors are scarce (Baker 1955). Later theory and synthesis have emphasized that this prediction concerns enrichment of reproductive assurance rather than universally high selfing rates, and that mate limitation and pollinator limitation are distinct processes (Cheptou 2012; Pannell et al. 2015). Consistent with this colonization-filter expectation, self-compatibility is over-represented on islands in broad comparative data (Grossenbacher et al. 2017). A recent global analysis of 3,222 flowering-plant species further showed that self-compatibility, lifespan and floral symmetry/generalization jointly predict the probability of island colonization, with arrival opportunity modifying that filter (Zell et al. 2025). Reproductive assurance therefore provides a well-supported starting point for a floral island syndrome, but not necessarily its complete explanation. Regional island studies have also described subdued colours and accessible floral forms as components of an insular pollination syndrome (Abe 2006), whereas sister-taxon comparisons across Pacific islands found no universal evolutionary reduction in flower size and instead emphasized island- and lineage-dependent outcomes (Hetherington-Rauth & Johnson 2020). What remains unresolved is whether, among already assembled island floras, increasing source isolation is associated with a recurrent multivariate reproductive response and whether that same geographic gradient carries an independently measured reproductive constraint.

Floral phenotype can respond along a second axis. Greater selfing is often accompanied by a selfing syndrome involving reduced or reorganized floral investment (Sicard & Lenhard 2011). Yet floral form, symmetry, tube depth and colour also mediate access and effectiveness of animal pollen vectors, and combinations of these traits often predict functional pollinator groups better than any single character (Fenster et al. 2004; Rosas-Guerrero et al. 2014). An accessible flower could therefore become common because selfing relaxes the value of specialized attraction, or because broader visitor access remains advantageous even in plants that do not shift mating system. These alternatives make a key distinction: a floral island syndrome could be a single serial pathway from isolation to selfing to floral change, or a multicomponent response in which reproductive assurance and pollinator-facing architecture are only partly coupled.

Trait distributions alone cannot reveal whether natural pollen delivery actually becomes more constraining with isolation. Pollen-supplementation experiments provide an independent outcome by comparing natural reproduction with reproduction after added pollen. Pollen limitation integrates pollen quantity, quality, mate availability and pollinator service, and can influence both demography and selection on floral traits (Ashman et al. 2004; Harder & Aizen 2010; Knight et al. 2005). The GloPL synthesis makes this outcome available globally (Bennett et al. 2018a,b), and global analyses already show that pollen limitation varies systematically with ecological context and pollinator dependence (Bennett et al. 2020). A previous GloPL–GBIF synthesis found little evidence that pollen limitation generally increases towards species' range edges and excluded remote islands more than 200 km from mainland from its range-edge analysis (Dawson-Glass & Hargreaves 2022). The island source-isolation gradient tested here is therefore a distinct biogeographic axis. Linking island trait patterns to this experimental outcome can test whether an inferred island syndrome is functionally aligned with reduced dependence on external pollen delivery, without claiming that pollinator abundance itself has been measured.

Here we combine a global island plant-trait database with independent pollen-supplementation data to ask whether island reproductive change forms one serial syndrome or a recurrent but internally separable response. **H1** tests whether geographic isolation is associated with a recurring multivariate floral/reproductive response across four geographic regions. **H2** asks whether floral accessibility and colour responses remain after measured reproductive assurance is conditioned on. **H3** tests whether experimental pollen limitation independently increases with isolation. **H4** asks whether the same species-level reproductive-assurance and accessibility states enriched with isolation are associated with lower current pollen limitation. Together, these analyses test a constraint–response triangle linking geographic isolation, contemporary pollen limitation and plant reproductive strategy while keeping historical causation explicitly unresolved.

## Materials and Methods

### Global island plant database

We used a corrected global analysis universe of **8,264 island units**. Candidate units originated from GSHHG 2.3.7 high-resolution island geometry. A geography audit identified one continental split component that had been retained as an island; excluding it yielded the corrected 8,264-unit universe and a broad H1 union of 4,379 islands. The trait database retained the same taxonomic and provenance structure, containing 106,295 accepted angiosperm species.

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was assembled source by source from floras, monographs, public trait resources and primary literature. Every accepted record retained source provenance. Species-direct high- and medium-confidence evidence had priority, while validated lower-confidence evidence was used only where direct evidence was unavailable. Missing trait information was retained as missing rather than converted to trait absence, and family-level inference was not used as a general fill rule.

The final trait dataset contained 222,688 resolved cells of 318,885 possible species-by-axis cells (69.83%) across three raw evidence axes: flower colour, floral structural complexity and reproductive assurance. The broad all-analysis-eligible evidence scope was used as the primary plant analysis, with species-direct high/medium evidence retained as a Direct-only sensitivity.

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

Island-level trait counts were modelled using beta-binomial logit models with standardized log isolation distance, island area and climate PC1-PC4. Within each geographic stratum, a seven-dimensional Wald test evaluated the joint isolation-response vector. Individual response coefficients were retained so that joint support could not be misread as uniform support for every trait. Optimizer convergence was audited separately from statistical support; a single Direct-only northern-high-latitude shallow/open-tube fit with a frozen termination warning was re-fitted under expanded line-search settings and checked with a six-response leave-one-component sensitivity.

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

We used the version-pinned GloPL database of pollen-supplementation experiments (Bennett et al. 2018a,b). Pollen limitation was represented by the log response ratio of reproduction after supplemental pollen addition to reproduction under natural pollen receipt; positive values therefore indicate that reproduction increases when pollen is added.

The analysis contained 2,969 experimental rows, 1,408 measurement cells, 1,248 sites and 919 publications. Duplicate measurements were aggregated at the publication-by-coordinate-by-measurement level. Each publication contributed total analysis weight one. The primary model related pollen limitation to standardized corrected log geographic distance while controlling for broad geographic context and experimental measurement conditions. Uncertainty was estimated with publication-cluster-robust covariance.

### H4: functional bridge from H2 traits to current pollen limitation

H4 was an explicitly post-hoc functional triangulation analysis. Species in the H2 trait dataset were exact-matched to GloPL species after only case normalization and underscore/space normalization; synonym rescue and genus fallback were not used.

The primary H4 predictors were the literal Direct-only species-level H2 scores used in the plant analysis: `selfing_core` and `generalized_accessible`. Both are soft-membership scores scaled from 0 to 1, with higher values representing stronger concordance with reproductive assurance or generalized accessibility, respectively. Missing component traits were not coded as zero.

For each H2 score, we fitted:

`pollen limitation ~ H2 trait score + standardized distance + geographic context + measurement conditions`.

Each publication again received total weight one and uncertainty was publication-cluster robust. A negative H2-score coefficient indicates that species expressing the island-associated trait state more strongly tend to experience lower current experimental pollen limitation, conditional on geographic distance and measured study structure.

We also retained atomic-trait sensitivities. Because the global colour response was geographically heterogeneous and lacked a single directional functional prediction, colour was not promoted into the global H4 bridge.

## Results

### H1: a recurrent multivariate island response with region-dependent expression

The seven-response isolation vector remained supported in all four broad geographic strata after the geography correction. Joint FDR-adjusted q-values were **3.216 × 10^-10** in northern mid-latitudes (2,173 islands), **2.433 × 10^-5** in northern high latitudes (411 islands), **3.498 × 10^-7** in the tropics (1,493 islands) and **2.504 × 10^-18** in southern extratropical islands (302 islands). Direct-only joint support also remained in all four regions. The frozen northern-high-latitude Direct-only table carried an optimizer warning for shallow/open tube, but an independent re-fit converged to essentially the same coefficient (absolute change = 2.93 × 10^-6) and the seven-response vector remained supported (q = 3.79 × 10^-8). Excluding that component entirely also retained a fully converged six-response joint result (q = 1.43 × 10^-8; Appendix S3).

Joint support did not imply a uniformly positive seven-trait syndrome. Twenty-six of 28 broad all-analysis coefficients were positive, but the southern shallow/open-tube coefficient was negative (beta = -0.2166, nominal p = 0.000370). Northern high-latitude generalized form weakened to p = 0.1018, and southern selfing mating system weakened to p = 0.0540. The corrected baseline therefore supports recurring functional components with region- and trait-dependent expression rather than one universal response of every floral trait (Figure 4).

### H2: floral accessibility is not reducible to reproductive assurance

The `selfing_core` isolation slope was positive in all four regions in the primary analysis, although nominal support remained concentrated in the tropical and southern strata.

After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. FDR-adjusted q-values were 0.1703 in northern mid-latitudes, 0.000847 in northern high latitudes, 0.01616 in the tropics and 0.2873 in southern extratropical islands. Thus the primary evidence for a selfing-adjusted accessibility response is strongest in northern high latitudes and the tropics. In the Direct-only sensitivity, the tropical accessibility estimate remained nominally positive (p = 0.0448) but no longer survived FDR correction (q = 0.1196); it is therefore not treated as a replicated FDR-supported result.

Colour responses remained more heterogeneous. Southern adjusted plain colour remained supported (q = 0.000847). Raw colour, architecture and colour-conditioned architecture analyses showed region-specific responses. Northern-high-latitude blue/purple architecture associations remained negative, while tropical Direct evidence retained a positive yellow/orange × butterfly/deep-tube coupling (q = 0.03949). These descriptive trait combinations do not establish realized pollinator identity.

The corrected H2 results therefore support a partially separable floral response: accessibility estimates remained positive after reproductive-assurance adjustment in all four regions, but FDR support was concentrated in northern high latitudes and the tropics. The strength and detailed form of the additional floral response are geographically contingent (Figure 5).

### H3: pollen limitation increases with isolation

Across 2,969 pollen-supplementation experiments, the corrected standardized distance coefficient was positive (**beta = 0.09191, SE = 0.03806; two-sided p = 0.01575; one-sided positive p = 0.00787**). The no-zero-constant sensitivity also retained a supported positive coefficient (beta = 0.09089, p = 0.01869).

The supplemental-only sensitivity was positive but unsupported (beta = 0.04410, p = 0.31082). The primary result therefore supports an isolation-associated increase in experimental pollen limitation, while the evidence is not equally strong under every measurement definition (Figure 6A).

### H4: island-associated response traits are linked to lower current pollen limitation

At the atomic-trait level, autonomous selfing provided the strongest functional association. Species with autonomous or delayed selfing had lower current pollen limitation after adjustment for distance, broad context and measurement structure (beta = -0.44492, p = 2.89 × 10^-8). Actinomorphy was also strongly negative (beta = -0.38076, p = 1.20 × 10^-5), generalized form was weaker but supported (beta = -0.18448, p = 0.04434), and self-compatibility was negative but imprecise (beta = -0.11885, p = 0.1495).

The literal H2-to-H4 bridge gave the same overall interpretation. The Direct-only `selfing_core` score matched 455 GloPL species across 409 publications and was negatively associated with current pollen limitation (**beta = -0.29830, SE = 0.10352, p = 0.00396**). The association remained supported in the no-zero-constant sensitivity, whereas the supplemental-only subset remained negative but unsupported (p = 0.28962).

The Direct-only `generalized_accessible` score matched 143 GloPL species across 143 publications. Higher accessibility was associated with lower current pollen limitation (**beta = -0.29566, SE = 0.12896, p = 0.02187**). The supplemental-only sensitivity was also negative (p = 0.03971), and the no-zero-constant sensitivity remained supported.

Thus the same two response families enriched along the island-isolation gradient are associated, in exact-species overlap, with lower current pollen limitation (Figure 6B-C). This functional alignment does not establish that historical pollen limitation caused the observed trait distributions.

## Discussion

### A recurrent functional response, not a uniform phenotype

The central result is not that every island flower changes in the same direction. Instead, geographic isolation is associated with a recurrent multivariate shift in reproductive strategy whose detailed floral realization varies among regions. This reconciles earlier evidence for regional island pollination syndromes (Abe 2006) with comparative evidence that individual floral traits can evolve idiosyncratically among islands and lineages (Hetherington-Rauth & Johnson 2020). This distinction is consistent with functional island biogeography: common filters can repeatedly favour similar ecological functions while regional species pools, histories and interaction networks produce different trait combinations (Schrader et al. 2021). The negative southern shallow/open-tube coefficient is therefore informative rather than anomalous. It shows that a floral island syndrome is best understood as multivariate recurrence, not an all-or-nothing checklist of trait changes.

The strongest recurring components are reproductive assurance and floral accessibility/generalization. That combination extends classic island reproductive theory. Baker's law predicts enrichment of uniparental reproductive capacity during colonization, but neither Baker's law nor the selfing syndrome requires all pollinator-facing traits to collapse into a single downstream response (Cheptou 2012; Pannell et al. 2015; Sicard & Lenhard 2011). Our global result places that distinction at assemblage scale.

### Floral reorganization extends beyond the selfing syndrome

Reproductive assurance remains one plausible route through which isolated floras reduce dependence on uncertain pollen delivery. The over-representation of self-compatible species on islands provides independent comparative support for this colonization-filter logic (Grossenbacher et al. 2017), and Zell et al. (2025) showed that breeding system and floral symmetry jointly predict island colonization probability. Our analysis asks a different question: among established island floras, does increasing geographic isolation continue to organize reproductive function? H2 shows that floral accessibility is not statistically exhausted by measured reproductive assurance. After adjustment for selfing_core, accessibility estimates remain positive in all four regions and are FDR-supported in northern high latitudes and the tropics.

This pattern rejects a compulsory serial interpretation in which isolation first increases selfing and floral reorganization is only its by-product. It is instead compatible with partially separable responses: reproductive assurance can reduce the need for external pollen, while accessible floral architectures may alter the range of visitors capable of contacting reproductive organs. Pollination-syndrome research shows that floral traits work in combinations and that functional groups can exert different selective pressures, while also warning against identifying realized pollinators from phenotype alone (Fenster et al. 2004; Rosas-Guerrero et al. 2014). The region-specific colour × architecture results fit that view: the broad function recurs, but the detailed display does not.

### Independent experiments identify an isolation-associated pollen constraint

The GloPL analysis supplies an ecological layer that is independent of the island trait database. Experimental pollen limitation increases with geographic isolation even after the corrected geography is used. This result matters because pollen limitation is a reproductive outcome rather than a floral proxy. It captures shortfalls in successful pollen receipt arising from pollen quantity, quality, mate availability and pollinator service (Ashman et al. 2004; Knight et al. 2005). It therefore supports an isolation-associated constraint on pollen delivery without requiring a claim of global pollinator decline.

The effect is modest and not equally strong under every measurement sensitivity, but its direction aligns with a large literature showing that pollen limitation responds to ecological context and can shape floral adaptation (Bennett et al. 2020; Harder & Aizen 2010). The contrast with prior GloPL geography is informative: pollen limitation showed little general increase towards species' range edges (Dawson-Glass & Hargreaves 2022), whereas it increases along the island-to-continent isolation axis tested here. H3 therefore identifies a geographic reproductive constraint specific to island source isolation rather than a generic tendency for pollen limitation to rise at all distributional margins. The independent H3 result changes the interpretation of H1–H2: the island trait gradient occurs along a geographic axis that is also associated with experimentally measured reproductive constraint.

### A constraint–response triangle, with one causal edge still missing

The strongest synthesis comes from combining three associations. First, isolation is associated with recurrent reproductive-assurance and accessibility responses. Second, isolation is independently associated with greater pollen limitation. Third, in exact-species overlap, stronger reproductive-assurance and accessibility scores are associated with lower current pollen limitation. Autonomous selfing provides the clearest atomic example, consistent with its direct capacity to reproduce when external pollen is insufficient.

Together these results form a constraint–response triangle. They show that the traits enriched along the island-isolation gradient are functionally aligned with reduced contemporary pollen limitation, while the same geographic gradient is associated with stronger pollen limitation. What they do **not** observe is the historical edge from past pollen limitation through selection, sorting or persistence to present-day trait composition. H4 was explicitly post-hoc, and these data contain neither ancestral states nor temporal changes in pollen service. The triangle therefore strengthens functional interpretation without converting correlation into mediation.

### Why regional floral trajectories diverge

The recurrent functional core coexists with strongly context-dependent colour and colour–architecture responses. This is expected if the payoff of a floral phenotype depends on the starting plant strategy and on the realized interaction community. Pollinator groups differ in their effectiveness and in the floral traits they favour, but the same phenotype can interact with multiple groups and the same functional outcome can be achieved by different trait combinations (Fenster et al. 2004; Rosas-Guerrero et al. 2014). Regional divergence is therefore a prediction of a multicomponent response rather than evidence against recurrence.

This result also defines the next mechanistic question. Global cross-sectional data can identify recurrent functions and their ecological correlates, but explaining why one region changes colour, tube depth or architecture differently from another requires direct information on realized visitors, pollen transfer, mating outcomes and lineage history.

### Limitations

The primary plant responses describe contemporary island-flora composition. They do not distinguish species sorting, colonization filtering and within-lineage evolutionary change. Geographic distance is a composite exposure that covaries with connectivity and source supply, and residual confounding remains possible despite adjustment for island area and climate.

The corrected distance metric repairs a known geometry mismatch using source-matched GSHHG shorelines, but the repair was selected post hoc after audit. It is therefore the preferred measurement baseline, not a prospective confirmation. Numerical distance precision also exceeds the positional accuracy implied by the underlying shoreline data.

H2 is a conditional decomposition, not causal mediation. A floral association that remains after adjustment for selfing_core shows that measured reproductive assurance does not statistically absorb the pattern; it does not prove direct selection by pollinators. H4 is likewise post-hoc functional triangulation. A separate prospective validation audit was support-limited before outcome unblinding and is not part of the H4 result. Stronger causal inference will require lineage-resolved or longitudinal data linking pollination service, reproductive success and trait change through time.

## Conclusion

Global island floras show a recurrent reproductive response to geographic isolation, but that response is not a single serial selfing syndrome. Reproductive assurance and floral accessibility form partially separable components: accessibility estimates remain positive after reproductive-assurance adjustment across regions, with strongest statistical support in northern high latitudes and the tropics. Independently, pollen-supplementation experiments show increasing pollen limitation with isolation, and the same reproductive-assurance and accessibility states enriched along the island gradient are associated with lower current pollen limitation. Island isolation therefore aligns ecological constraint with functionally compatible plant responses while leaving the historical causal route—and the reason detailed floral trajectories diverge among regions—unresolved.

## References

Abe, T. (2006). Threatened pollination systems in native flora of the Ogasawara (Bonin) Islands. *Ann. Bot.*, 98, 317–334.

Ashman, T.-L., Knight, T.M., Steets, J.A., Amarasekare, P., Burd, M., Campbell, D.R. et al. (2004). Pollen limitation of plant reproduction: ecological and evolutionary causes and consequences. *Ecology*, 85, 2408–2421.

Baker, H.G. (1955). Self-compatibility and establishment after ‘long-distance’ dispersal. *Evolution*, 9, 347–349.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018a). Data from: GloPL, a global data base on pollen limitation of plant reproduction. Dryad Digital Repository. Available at: https://doi.org/10.5061/dryad.dt437.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018b). GloPL, a global data base on pollen limitation of plant reproduction. *Sci. Data*, 5, 180249.

Bennett, J.M., Steets, J.A., Burns, J.H., Burkle, L.A., Vamosi, J.C., Wolowski, M. et al. (2020). Land use and pollinator dependency drives global patterns of pollen limitation in the Anthropocene. *Nat. Commun.*, 11, 3999.

Cheptou, P.-O. (2012). Clarifying Baker's Law. *Ann. Bot.*, 109, 633–641.

Dawson-Glass, E. & Hargreaves, A.L. (2022). Does pollen limitation limit plant ranges? Evidence and implications. *Philos. Trans. R. Soc. B*, 377, 20210014.

Fenster, C.B., Armbruster, W.S., Wilson, P., Dudash, M.R. & Thomson, J.D. (2004). Pollination syndromes and floral specialization. *Annu. Rev. Ecol. Evol. Syst.*, 35, 375–403.

Grossenbacher, D.L., Brandvain, Y., Auld, J.R., Burd, M., Cheptou, P.-O., Conner, J.K. et al. (2017). Self-compatibility is over-represented on islands. *New Phytol.*, 215, 469–478.

Harder, L.D. & Aizen, M.A. (2010). Floral adaptation and diversification under pollen limitation. *Philos. Trans. R. Soc. B*, 365, 529–543.

Hetherington-Rauth, M.C. & Johnson, M.T.J. (2020). Floral trait evolution of angiosperms on Pacific islands. *Am. Nat.*, 196, 87–100.

Knight, T.M., Steets, J.A., Vamosi, J.C., Mazer, S.J., Burd, M., Campbell, D.R. et al. (2005). Pollen limitation of plant reproduction: pattern and process. *Annu. Rev. Ecol. Evol. Syst.*, 36, 467–497.

Pannell, J.R., Auld, J.R., Brandvain, Y., Burd, M., Busch, J.W., Cheptou, P.-O. et al. (2015). The scope of Baker's law. *New Phytol.*, 208, 656–667.

Rosas-Guerrero, V., Aguilar, R., Martén-Rodríguez, S., Ashworth, L., Lopezaraiza-Mikel, M., Bastida, J.M. et al. (2014). A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecol. Lett.*, 17, 388–400.

Schrader, J., Wright, I.J., Kreft, H. & Westoby, M. (2021). A roadmap to plant functional island biogeography. *Biol. Rev.*, 96, 2851–2870.

Sicard, A. & Lenhard, M. (2011). The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Ann. Bot.*, 107, 1433–1443.
Zell, A.N., Miranda, C.H., Grady, E.L., Grossenbacher, D.L. & Igić, B. (2025). Island colonization in flowering plants is determined by the interplay of breeding system, lifespan, floral symmetry, and arrival opportunity. *New Phytol.*, 245, 420–432.
