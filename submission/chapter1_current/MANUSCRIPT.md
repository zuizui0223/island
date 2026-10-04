> **Superseded inference draft:** the t/Holm interpretation below is no longer primary. The user restored original seven-trait poster inference on 4 October 2026: see `results/h1_poster_original_wcvp_20261004/README.md` and the current selector. This draft must not be submitted or used to replace the original poster conclusions.

# Traitwise floral responses to island isolation and their functional alignment with pollen limitation

## Abstract

Across a corrected universe of 8,264 island units and a database of 106,295 angiosperms, we examine seven individual reproductive and floral responses to isolation. Broad-flora traitwise inference supports increasing selfing in northern mid-latitudes, increasing plain colour in southern extratropical floras and decreasing shallow/open tubes in the latter region after correction for 28 comparisons. Only the southern colour response is supported in every evidence and WCVP sensitivity scope. Independently, pollen limitation increases with isolation (P = 0.0159), while reproductive assurance and accessibility are associated with lower current pollen limitation (P = 0.00417 and 0.0233). These assemblage and functional associations reveal components of an island syndrome without establishing universal floral simplification or historical causation.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; flower colour; pollination

## Introduction

Islands are natural tests of how dispersal, establishment, biotic interactions and evolution reshape ecological strategies. Functional island biogeography predicts that insularity can repeatedly filter traits linked to colonization and persistence, while the detailed outcome depends on the regional species pool and interaction environment (Schrader et al. 2021). Plant reproduction is especially exposed to these filters because successful establishment may require both compatible mates and effective pollen transfer.

Baker's law formalized one solution: colonists capable of uniparental reproduction can establish when mates or compatible pollen donors are scarce (Baker 1955). Later theory and synthesis have emphasized that this prediction concerns enrichment of reproductive assurance rather than universally high selfing rates, and that mate limitation and pollinator limitation are distinct processes (Cheptou 2012; Pannell et al. 2015). Consistent with this colonization-filter expectation, self-compatibility is over-represented on islands in broad comparative data (Grossenbacher et al. 2017). A recent global analysis of 3,222 flowering-plant species further showed that self-compatibility, lifespan and floral symmetry/generalization jointly predict the probability of island colonization, with arrival opportunity modifying that filter (Zell et al. 2025). Reproductive assurance therefore provides a well-supported starting point for a floral island syndrome, but not necessarily its complete explanation. Regional island studies have also described subdued colours and accessible floral forms as components of an insular pollination syndrome (Abe 2006), whereas sister-taxon comparisons across Pacific islands found no universal evolutionary reduction in flower size and instead emphasized island- and lineage-dependent outcomes (Hetherington-Rauth & Johnson 2020). What remains unresolved is whether, among already assembled island floras, increasing source isolation is associated with a recurrent multivariate reproductive response and whether that same geographic gradient carries an independently measured reproductive constraint.

Floral phenotype can respond along a second axis. Greater selfing is often accompanied by a selfing syndrome involving reduced or reorganized floral investment (Sicard & Lenhard 2011). Yet floral form, symmetry, tube depth and colour also mediate access and effectiveness of animal pollen vectors, and combinations of these traits often predict functional pollinator groups better than any single character (Fenster et al. 2004; Rosas-Guerrero et al. 2014). An accessible flower could therefore become common because selfing relaxes the value of specialized attraction, or because broader visitor access remains advantageous even in plants that do not shift mating system. These alternatives make a key distinction: a floral island syndrome could be a single serial pathway from isolation to selfing to floral change, or a multicomponent response in which reproductive assurance and pollinator-facing architecture are only partly coupled.

Trait distributions alone cannot reveal whether natural pollen delivery actually becomes more constraining with isolation. Pollen-supplementation experiments provide an independent outcome by comparing natural reproduction with reproduction after added pollen. Pollen limitation integrates pollen quantity, quality, mate availability and pollinator service, and can influence both demography and selection on floral traits (Harder & Aizen 2010; Knight et al. 2005). The GloPL synthesis makes this outcome available globally (Bennett et al. 2018a,b), and global analyses already show that pollen limitation varies systematically with ecological context and pollinator dependence (Bennett et al. 2020). A previous GloPL–GBIF synthesis found little evidence that pollen limitation generally increases towards species' range edges and excluded remote islands more than 200 km from mainland from its range-edge analysis (Dawson-Glass & Hargreaves 2022). The island source-isolation gradient tested here is therefore a distinct biogeographic axis. Linking island trait patterns to this experimental outcome can test whether an inferred island syndrome is functionally aligned with reduced dependence on external pollen delivery, without claiming that pollinator abundance itself has been measured.

Here we combine global island floras and independent pollen-supplementation experiments. H1 tests individual reproductive, colour and structural responses to isolation across four regions. H2 asks whether floral responses persist after adjustment for reproductive assurance. H3 tests isolation-associated pollen limitation. H4 tests functional associations between the H2 states and pollen limitation in exact-species overlap. The design separates observed patterns from their potential ecological mechanisms.

## Materials and Methods

### Global island plant database

We used a corrected global analysis universe of **8,264 island units**. Candidate units originated from GSHHG 2.3.7 high-resolution island geometry. A geography audit identified one continental split component that had been retained as an island; excluding it yielded the corrected 8,264-unit universe and a broad H1 union of 4,379 islands. The trait database retained the same taxonomic and provenance structure, containing 106,295 accepted angiosperm species.

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was assembled source by source from floras, monographs, public trait resources and primary literature. Every accepted record retained source provenance. Species-direct high- and medium-confidence evidence had priority, while validated lower-confidence evidence was used only where direct evidence was unavailable. Missing trait information was retained as missing rather than converted to trait absence, and family-level inference was not used as a general fill rule.

The final trait dataset contained 222,688 resolved cells of 318,885 possible species-by-axis cells (69.83%) across three raw evidence axes: flower colour, floral structural complexity and reproductive assurance. The broad all-analysis-eligible evidence scope was used as the primary plant analysis, with species-direct high/medium evidence retained as a Direct-only sensitivity.

The primary H1 sampling frame deliberately retained all observed island-flora records, irrespective of floristic origin, to avoid discarding the majority of records whose origin status was unresolved. Among 1,039,757 island-by-species records in the frozen flora ledger, 150,531 were source-backed native, 19,442 were source-backed introduced and 869,784 (83.65%) had unresolved origin status. The all-observed analysis therefore estimates trait composition in the contemporary observed flora; it is not, by itself, an analysis of native colonization or in-situ evolution. WCVP regional-native compatibility is the sole active floristic-origin sensitivity.

Because those strict subsets were severely support-limited, we additionally performed a regional-native compatibility sensitivity using the World Checklist of Vascular Plants (WCVP; Govaerts et al. 2021). Accepted species were exact-matched to the WCVP bulk names table, and native distributions were read at TDWG level 3. An unresolved island-by-species record was upgraded to regional-native-compatible only when the focal island mapped unambiguously to one TDWG level-3 unit and WCVP listed that accepted species as native in the same unit. Source-backed introduced records were never overwritten. The WCVP archive used for this replay had SHA-256 `d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa`. This procedure establishes regional native compatibility rather than exact island-level nativeness. The earlier complementary-origin and TDWG-area analyses are archived provenance diagnostics and are not included as additional active flora partitions in this replay.

### Geographic isolation and covariates

The primary geographic exposure is the minimum separation between island and continental coastlines measured on a mean-radius sphere (R = 6371.0088 km). Island and continental geometries were taken from the same GSHHG 2.3.7 high-resolution archive. Continental siblings were reconstructed before distance calculation, artificial dateline split edges were omitted, and distance was calculated as the minimum minor-great-circle arc separation between coastline segments. This is a spherical coastline distance, not a WGS84 ellipsoidal distance.

The geography audit identified 1,113 broad-H1 islands that had been assigned spurious zero distance when GSHHG island polygons were compared with a coarser Natural Earth continental geometry. Under the source-matched coastline metric, all 1,113 have positive distance. This correction is a post-hoc measurement repair selected as the primary submission baseline, not a new prospective confirmation.

Plant models used standardized log distance, log island area and climate PC1-PC4. Spatial dependence was addressed with spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of source separation, connectivity and colonization opportunity rather than as a randomized treatment.

Four predeclared geographic replication strata were analysed separately: northern mid-latitude, northern high-latitude, tropical and southern extratropical islands. These strata were used to assess recurrence of response structure, not to require identical coefficients among regions.

For GloPL, geographic exposure was recomputed from each site's coordinates against the same continental geometry. Sites lying on seeded continental land were assigned zero distance; 996 of 1,248 sites are true continental zeros and remain zero in the corrected baseline.

### H1: individual trait responses across four regions

H1 now tests seven previously defined binary outcomes separately: self-compatibility, selfing mating system, autonomous/delayed selfing, plain colour, generalized form, actinomorphy and shallow/open tubes. Reproduction, colour and structure are organizational domains only; neither domain scores nor an omnibus directional score are calculated. This user-requested revision on 4 October 2026 follows inspection of earlier results and is explicitly retrospective.

For each trait and region, counts of species expressing the state are modelled against their informative species denominator using a beta-binomial logit regression. Predictors are standardized corrected log isolation, log area and climate PC1–4. Spatial-block sandwich uncertainty uses a t reference with G−1 degrees of freedom; pointwise 95% intervals accompany every slope. Two-sided Holm adjustment covers all 28 region-by-trait tests within each flora/evidence scope, retaining failed tests in the family. Positive slopes indicate increasing prevalence, not an imposed direction of evolution. Cross-region differences in significance do not constitute a formal interaction test.

Broad contemporary flora is primary. The sole active floristic-origin sensitivity is WCVP regional-native compatibility; Direct-only is separately retained as an evidence-quality sensitivity. Historical pooled-score, strict-native, complementary-origin and high-dimensional omnibus analyses remain archived, not active evidence for H1. H2–H4 retain their distinct estimands and existing functional covariates.

### H2: separating reproductive assurance from additional floral change

We first defined a species-level reproductive-assurance score, `selfing_core`, using self-incompatibility/compatibility, mating system and autonomous-selfing capacity only. Flower colour, form, symmetry, tube depth and flower size were excluded from this score.

The reproductive-assurance response was evaluated as:

`selfing_core ~ isolation + island area + climate`.

We then asked whether floral responses persisted after reproductive assurance was conditioned on:

`plain_colour ~ isolation + selfing_core + island area + climate`

and

`generalized_accessible ~ isolation + selfing_core + island area + climate`.

The `generalized_accessible` score summarizes open/generalized floral form, actinomorphy and shallow/open tube depth. Persistence of an isolation coefficient after adjustment for `selfing_core` was interpreted as evidence that the floral response is not reducible to measured reproductive assurance. This is a conditional decomposition rather than a formal causal mediation analysis.

As a reliability sensitivity to incomplete mediator measurement, we rebuilt the Direct-only `selfing_core` score using only species for which self-incompatibility/compatibility, mating system and autonomous-selfing capacity were all observed and informative. This complete-three-component score was available for 564 Direct-only species and produced island-level scores on 2,610 islands before model covariate filtering. We then refitted the same `generalized_accessible ~ isolation + selfing_core + island area + climate` model. This sensitivity tests whether missing selfing components explain the conditional accessibility signal; it is not a formal errors-in-variables correction.

To avoid assigning plants to named pollinator syndromes from a single composite score, we additionally analysed raw reported colour states and raw architecture. Five colour states were analysed after adjustment for `selfing_core`: white, red/pink, yellow/orange, blue/purple and green/brown/inconspicuous. We then tested raw colour × floral-form and colour × tube-depth combinations both as joint prevalence and as conditional architecture given colour. These analyses evaluate display-architecture reorganization but do not identify realized pollinator identity.

### H3: experimental pollen limitation

We used the version-pinned GloPL database of pollen-supplementation experiments (Bennett et al. 2018a,b). Pollen limitation was represented by the log response ratio of reproduction after supplemental pollen addition to reproduction under natural pollen receipt; positive values therefore indicate that reproduction increases when pollen is added.

The analysis contained 2,969 experimental rows, 1,408 measurement cells, 1,248 sites and 919 publications. Duplicate measurements were aggregated at the publication-by-coordinate-by-measurement level. Each publication contributed total analysis weight one. The primary model related pollen limitation to standardized corrected log geographic distance while controlling for broad geographic context and experimental measurement conditions. Uncertainty was estimated with publication-cluster-robust covariance.

Because 996 of 1,248 GloPL sites are true continental zero-distance sites, we added post-hoc robustness checks to distinguish a continuous offshore gradient from a simple mainland/offshore step. We refitted the model among positive-distance measurement cells only, fitted a separate binary mainland-versus-offshore indicator model, and repeated the offshore-only gradient after deleting each publication in turn. These are robustness analyses rather than new confirmatory tests.

### H4: functional bridge from H2 traits to current pollen limitation

H4 was an explicitly post-hoc functional triangulation analysis. Species in the H2 trait dataset were exact-matched to GloPL species after only case normalization and underscore/space normalization; synonym rescue and genus fallback were not used.

The primary H4 predictors were the literal Direct-only species-level H2 scores used in the plant analysis: `selfing_core` and `generalized_accessible`. Both are soft-membership scores scaled from 0 to 1, with higher values representing stronger concordance with reproductive assurance or generalized accessibility, respectively. Missing component traits were not coded as zero.

For each H2 score, we fitted:

`pollen limitation ~ H2 trait score + standardized distance + geographic context + measurement conditions`.

Each publication again received total weight one and uncertainty was publication-cluster robust. A negative H2-score coefficient indicates that species expressing the island-associated trait state more strongly tend to experience lower current experimental pollen limitation, conditional on geographic distance and measured study structure.

We also retained atomic-trait sensitivities. Because the global colour response was geographically heterogeneous and lacked a single directional functional prediction, colour was not promoted into the global H4 bridge.

## Results

### H1: region- and trait-specific isolation responses

All 112 individual-trait models converged. In broad all-analysis, three of 28 associations survived two-sided Holm correction: selfing mating system increased in northern mid-latitudes (slope 0.04634, adjusted P = 0.04201); plain colour increased in southern extratropical floras (0.07419, P = 0.00559); and shallow/open tubes decreased in the latter region (−0.21660, P = 0.03026). These are associations of assemblage composition with isolation.

The southern plain-colour increase was supported in every evidence/flora scope. Northern-midlatitude selfing was also supported in both Direct-only flora scopes. Other results were scope-specific: broad Direct-only supported tropical autonomous selfing and generalized form; WCVP all-analysis supported tropical plain colour and southern self-compatibility. The southern shallow/open-tube decrease did not survive correction in Direct-only or WCVP analyses. No northern-high-latitude individual contrast passed this 28-test correction. All coefficients and pointwise intervals, including unsupported and opposing estimates, are reported in results/h1_traitwise_20261004/traitwise_results.csv. These tests replace, rather than supplement, the historical omnibus directional score.

### H2: floral accessibility is not reducible to reproductive assurance

The `selfing_core` isolation slope was positive in all four regions in the primary analysis, although nominal support remained concentrated in the tropical and southern strata.

After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. With finite-cluster t references and the same frozen eight-test BH family, q-values were 0.1764 in northern mid-latitudes, 0.00243 in northern high latitudes, 0.01913 in the tropics and 0.2977 in southern extratropical islands. Thus the primary evidence for a selfing-adjusted accessibility response is strongest in northern high latitudes and the tropics. In the Direct-only sensitivity, the tropical accessibility estimate remained nominally positive (finite-cluster p = 0.0475) but did not survive FDR correction (q = 0.1267); it is therefore not treated as a replicated FDR-supported result.

The complete-three-component Direct-only sensitivity reduced the selfing mediator to 564 species but retained a positive isolation coefficient for `generalized_accessible` in all four regions. FDR-adjusted q-values were 0.0508 in northern mid-latitudes, 0.0376 in northern high latitudes, 0.0376 in the tropics and 0.2706 in southern extratropical islands. The corresponding fitted island counts were 1,725, 221, 438 and 157. Thus the conditional accessibility signal in northern high latitudes and the tropics persists when every species contributing to the reproductive-assurance score has all three selfing components directly observed, which weakens a missing-component explanation without eliminating general measurement-error concerns.

Colour responses remained more heterogeneous. Southern adjusted plain colour remained supported after finite-cluster correction (q = 0.00243 in all-analysis and 1.8 × 10^-5 in Direct-only). Raw colour, architecture and colour-conditioned architecture analyses showed region-specific response patterns, including negative northern-high-latitude blue/purple architecture associations and a positive tropical Direct yellow/orange × butterfly/deep-tube pattern. These raw pattern analyses are descriptive and do not establish realized pollinator identity.

The corrected H2 results therefore support a partially separable floral response: accessibility estimates remained positive after reproductive-assurance adjustment in all four regions, but FDR support was concentrated in northern high latitudes and the tropics. The strength and detailed form of the additional floral response are geographically contingent (Figure 5).

### H3: pollen limitation increases with isolation

Across 2,969 pollen-supplementation experiments, the corrected standardized distance coefficient was positive (**beta = 0.09191, SE = 0.03806; finite-publication two-sided p = 0.01594; one-sided positive p = 0.00797**). The no-zero-constant sensitivity also retained a supported positive coefficient (beta = 0.09089, p = 0.01869).

The supplemental-only sensitivity was positive but unsupported (beta = 0.04410, p = 0.31082). The primary result therefore supports an isolation-associated increase in experimental pollen limitation, while the evidence is not equally strong under every measurement definition (Figure 6A).

The post-hoc offshore-only analysis showed that this result was not generated solely by the contrast between 996 mainland zero-distance sites and offshore sites. Among 276 positive-distance measurement cells from 153 publications and 252 sites, pollen limitation still increased with isolation (**beta = 0.22031, SE = 0.09704, finite-publication p = 0.02459**). By contrast, a binary mainland-versus-offshore indicator alone was unsupported (beta = 0.12743, p = 0.14364). Leave-one-publication-out refits kept all 153 offshore-gradient estimates positive; the weakest two-sided result was p = 0.04892. These checks support a continuous offshore gradient while remaining post-hoc robustness analyses.

### H4: island-associated response traits are linked to lower current pollen limitation

At the atomic-trait level, autonomous selfing provided the strongest functional association. Species with autonomous or delayed selfing had lower current pollen limitation after adjustment for distance, broad context and measurement structure (beta = -0.44492, finite-publication p = 4.85 × 10^-8). Actinomorphy was also strongly negative (beta = -0.38076, finite-publication p = 1.47 × 10^-5), generalized form was weaker but supported (beta = -0.18448, finite-publication p = 0.04543), and self-compatibility was negative but imprecise (beta = -0.11885, finite-publication p = 0.1502).

The literal H2-to-H4 bridge gave the same overall interpretation. The Direct-only `selfing_core` score matched 455 GloPL species across 409 publications and was negatively associated with current pollen limitation (**beta = -0.29830, SE = 0.10352, finite-publication p = 0.00417**). The association remained supported in the no-zero-constant sensitivity, whereas the supplemental-only subset remained negative but unsupported (finite-publication p = 0.2907).

The Direct-only `generalized_accessible` score matched 143 GloPL species across 143 publications. Higher accessibility was associated with lower current pollen limitation (**beta = -0.29566, SE = 0.12896, finite-publication p = 0.02334**). The supplemental-only sensitivity was also negative (finite-publication p = 0.04239), and the no-zero-constant sensitivity remained supported (p = 0.02082).

Thus the same two response families enriched along the island-isolation gradient are associated, in exact-species overlap, with lower current pollen limitation (Figure 6B-C). This functional alignment does not establish that historical pollen limitation caused the observed trait distributions.

## Discussion

### Recurrent components without a universal floral checklist

The traitwise analysis makes both recurring and opposing components visible without allowing a positive colour effect to cancel a negative structural effect. Increasing plain colour in southern extratropical floras is the most consistent association across evidence and provenance scopes. The broad-flora decrease in shallow/open tubes in that same region cautions against equating subdued colour with generalized access. Other positive reproductive and structural responses depend on region and evidence scope. These patterns are compatible with regional assembly and interaction environments, but do not identify the visitors responsible or distinguish selection from colonization and persistence.

Removing the aggregate changes the question and its multiplicity burden: a supported average score does not guarantee that individual components survive 28 two-sided tests. We therefore do not carry the old global-average P value forward or claim that all regions independently express one syndrome. The biological domains remain a framework for comparing traits rather than three quantitative indices.

### Floral reorganization extends beyond the selfing syndrome

Reproductive assurance remains one plausible route through which isolated floras reduce dependence on uncertain pollen delivery. The over-representation of self-compatible species on islands provides independent comparative support for this colonization-filter logic (Grossenbacher et al. 2017), and Zell et al. (2025) showed that breeding system and floral symmetry jointly predict island colonization probability. Our analysis asks a different question: among established island floras, does increasing geographic isolation continue to organize reproductive function? H2 shows that floral accessibility is not statistically exhausted by measured reproductive assurance. After adjustment for selfing_core, accessibility estimates remain positive in all four regions and are FDR-supported in northern high latitudes and the tropics.

This pattern is inconsistent with a compulsory serial interpretation in which every floral response must be absorbed by measured `selfing_core`: selfing-adjusted accessibility is supported in selected regions rather than uniformly across regions and evidence scopes. It is therefore compatible with partially separable responses rather than proof of a universal second pathway. Reproductive assurance can reduce the need for external pollen, while accessible floral architectures may alter the range of visitors capable of contacting reproductive organs. Pollination-syndrome research shows that floral traits work in combinations and that functional groups can exert different selective pressures, while also warning against identifying realized pollinators from phenotype alone (Fenster et al. 2004; Rosas-Guerrero et al. 2014). The region-specific colour × architecture results fit that view: the broad function recurs, but the detailed display does not.

### Independent experiments identify an isolation-associated pollen constraint

The GloPL analysis supplies an ecological layer that is independent of the island trait database. Experimental pollen limitation increases with geographic isolation even after the corrected geography is used. This result matters because pollen limitation is a reproductive outcome rather than a floral proxy. It captures shortfalls in successful pollen receipt arising from pollen quantity, quality, mate availability and pollinator service (Knight et al. 2005). It therefore supports an isolation-associated constraint on pollen delivery without requiring a claim of global pollinator decline.

The effect is modest and not equally strong under every measurement sensitivity, but its direction aligns with a large literature showing that pollen limitation responds to ecological context and can shape floral adaptation (Bennett et al. 2020; Harder & Aizen 2010). The contrast with prior GloPL geography is informative: pollen limitation showed little general increase towards species' range edges (Dawson-Glass & Hargreaves 2022), whereas it increases along the island-to-continent isolation axis tested here. H3 therefore identifies a geographic reproductive constraint specific to island source isolation rather than a generic tendency for pollen limitation to rise at all distributional margins. The independent H3 result changes the interpretation of H1–H2: the island trait gradient occurs along a geographic axis that is also associated with experimentally measured reproductive constraint.

### A constraint–response triangle, with one causal edge still missing

The strongest synthesis comes from combining three associations. First, isolation is associated with specific reproductive and floral states, with support and sometimes direction varying among regions. Second, isolation is independently associated with greater pollen limitation. Third, in exact-species overlap, stronger reproductive-assurance and accessibility scores are associated with lower current pollen limitation. Autonomous selfing provides the clearest atomic example, consistent with its direct capacity to reproduce when external pollen is insufficient.

Together these results form a constraint–response triangle. They show that the traits enriched along the island-isolation gradient are functionally aligned with reduced contemporary pollen limitation, while the same geographic gradient is associated with stronger pollen limitation. A corrected-geography replay of the earlier predeclared distance-by-trait tests did not establish that either reproductive-assurance or floral-architecture states buffer this isolation–pollen-limitation slope; both families failed their frozen promotion rules. What the data do **not** observe is the historical edge from past pollen limitation through selection, sorting or persistence to present-day trait composition. H4 was explicitly post-hoc, and these data contain neither ancestral states nor temporal changes in pollen service. The triangle therefore strengthens functional interpretation without converting correlation into mediation.

### Why regional floral trajectories diverge

The supported trait associations coexist with context-dependent colour and colour–architecture responses. This is expected if the payoff of a floral phenotype depends on the starting plant strategy and on the realized interaction community. Pollinator groups differ in their effectiveness and in the floral traits they favour, but the same phenotype can interact with multiple groups and the same functional outcome can be achieved by different trait combinations (Fenster et al. 2004; Rosas-Guerrero et al. 2014). Regional divergence is compatible with a multicomponent response; the present traitwise analysis does not formally test cross-region differences.

This result also defines the next mechanistic question. Global cross-sectional data can identify recurrent functions and their ecological correlates, but explaining why one region changes colour, tube depth or architecture differently from another requires direct information on realized visitors, pollen transfer, mating outcomes and lineage history.

### Limitations

The primary plant responses describe contemporary island-flora composition. They do not distinguish species sorting, colonization filtering and within-lineage evolutionary change. Geographic distance is a composite exposure that covaries with connectivity and source supply, and residual confounding remains possible despite adjustment for island area and climate.

Floristic origin remains an important limitation. The WCVP regional-native-compatible sensitivity retains the southern plain-colour association in both evidence scopes, while other trait associations differ from the broad-flora analysis. WCVP compatibility is regional: a taxon native to a TDWG unit need not be native to every island within it. Consequently, these results cannot identify native-specific historical evolution. Strict-native and complementary-origin analyses are archived rather than retained as additional active sensitivities.

The corrected distance metric repairs a known geometry mismatch using source-matched GSHHG shorelines, but the repair was selected post hoc after audit. It is therefore the preferred measurement baseline, not a prospective confirmation. Numerical distance precision also exceeds the positional accuracy implied by the underlying shoreline data.

H2 is a conditional decomposition, not causal mediation. A floral association that remains after adjustment for selfing_core shows that measured reproductive assurance does not statistically absorb the pattern; it does not prove direct selection by pollinators. The complete-three-component Direct-only sensitivity reduces concern that the result is created only by missing mediator components, but it is not a formal measurement-error model and cannot establish that `selfing_core` is measured without error. H4 is likewise post-hoc functional triangulation. A separate prospective validation audit was support-limited before outcome unblinding and is not part of the H4 result. Stronger causal inference will require lineage-resolved or longitudinal data linking pollination service, reproductive success and trait change through time.

## Conclusion

Traitwise analysis reveals components of the classic island syndrome while preserving regional exceptions: plain colour increases consistently in southern extratropical floras, whereas structural accessibility does not uniformly increase. Reproductive responses depend on region and evidence scope. H2 separates some floral responses from measured reproductive assurance; H3 and H4 independently align isolation-associated pollen limitation with traits associated with lower contemporary limitation. These results make the syndrome empirically specific without constructing an aggregate index or identifying historical selection as its cause.

## References

Abe, T. (2006). Threatened pollination systems in native flora of the Ogasawara (Bonin) Islands. *Ann. Bot.*, 98, 317–334.


Baker, H.G. (1955). Self-compatibility and establishment after ‘long-distance’ dispersal. *Evolution*, 9, 347–349.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018a). Data from: GloPL, a global data base on pollen limitation of plant reproduction. Dryad Digital Repository. Available at: https://doi.org/10.5061/dryad.dt437.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018b). GloPL, a global data base on pollen limitation of plant reproduction. *Sci. Data*, 5, 180249.

Bennett, J.M., Steets, J.A., Burns, J.H., Burkle, L.A., Vamosi, J.C., Wolowski, M. et al. (2020). Land use and pollinator dependency drives global patterns of pollen limitation in the Anthropocene. *Nat. Commun.*, 11, 3999.

Cheptou, P.-O. (2012). Clarifying Baker's Law. *Ann. Bot.*, 109, 633–641.

Dawson-Glass, E. & Hargreaves, A.L. (2022). Does pollen limitation limit plant ranges? Evidence and implications. *Philos. Trans. R. Soc. B*, 377, 20210014.

Fenster, C.B., Armbruster, W.S., Wilson, P., Dudash, M.R. & Thomson, J.D. (2004). Pollination syndromes and floral specialization. *Annu. Rev. Ecol. Evol. Syst.*, 35, 375–403.

Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. (2021). The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Sci. Data*, 8, 215.

Grossenbacher, D.L., Brandvain, Y., Auld, J.R., Burd, M., Cheptou, P.-O., Conner, J.K. et al. (2017). Self-compatibility is over-represented on islands. *New Phytol.*, 215, 469–478.

Harder, L.D. & Aizen, M.A. (2010). Floral adaptation and diversification under pollen limitation. *Philos. Trans. R. Soc. B*, 365, 529–543.

Hetherington-Rauth, M.C. & Johnson, M.T.J. (2020). Floral trait evolution of angiosperms on Pacific islands. *Am. Nat.*, 196, 87–100.

Knight, T.M., Steets, J.A., Vamosi, J.C., Mazer, S.J., Burd, M., Campbell, D.R. et al. (2005). Pollen limitation of plant reproduction: pattern and process. *Annu. Rev. Ecol. Evol. Syst.*, 36, 467–497.

Pannell, J.R., Auld, J.R., Brandvain, Y., Burd, M., Busch, J.W., Cheptou, P.-O. et al. (2015). The scope of Baker's law. *New Phytol.*, 208, 656–667.

Rosas-Guerrero, V., Aguilar, R., Martén-Rodríguez, S., Ashworth, L., Lopezaraiza-Mikel, M., Bastida, J.M. et al. (2014). A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecol. Lett.*, 17, 388–400.

Schrader, J., Wright, I.J., Kreft, H. & Westoby, M. (2021). A roadmap to plant functional island biogeography. *Biol. Rev.*, 96, 2851–2870.

Sicard, A. & Lenhard, M. (2011). The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Ann. Bot.*, 107, 1433–1443.
Zell, A.N., Miranda, C.H., Grady, E.L., Grossenbacher, D.L. & Igić, B. (2025). Island colonization in flowering plants is determined by the interplay of breeding system, lifespan, floral symmetry, and arrival opportunity. *New Phytol.*, 245, 420–432.
