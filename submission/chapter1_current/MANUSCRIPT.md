# Island isolation is associated with recurrent reproductive assurance but regionally contingent floral change

## Abstract

Islands repeatedly challenge plant reproduction, but whether isolation produces one floral syndrome or recurrent functions expressed through different phenotypes remains unresolved. We analysed seven reproductive and floral traits separately across 8,264 island units and 106,295 angiosperms, and paired these data with 2,969 pollen-supplementation experiments. Self-compatibility increased with isolation in all four regions in both broad and WCVP regional-native-compatible analyses. Other selfing, colour and structural responses varied geographically. Floral accessibility remained associated with isolation after conditioning on reproductive assurance in northern high latitudes and tropical all-analysis. Independently, pollen limitation increased with isolation (β = 0.0919, P = 0.0159); in exact-species overlap, reproductive-assurance and accessibility scores were associated with lower current pollen limitation. Thus the most recurrent island response is reproductive function rather than a universal floral phenotype. The results support a constraint–response interpretation while leaving the historical causal path from pollen limitation to contemporary assemblages unresolved.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; functional island biogeography; pollination

## Introduction

Island biogeography is built around repeated ecological problems—arrival, establishment, persistence and interaction under geographic isolation—but repeated problems need not produce identical phenotypes. This distinction is especially important for flowers. Islands are often expected to favour self-compatibility, selfing, reduced floral display or more accessible flowers, and these traits are sometimes grouped as a floral or reproductive island syndrome. Yet a single visible syndrome is not required if different regional floras solve the same reproductive constraint in different ways. The more general question is therefore not whether every isolated flora converges on one flower type, but which reproductive functions recur with isolation, which floral traits remain contingent, and whether the same geographic gradient is associated with an independently measured constraint on pollen delivery.

Baker's law provides the clearest functional prediction. Colonists capable of uniparental reproduction can establish when compatible mates or pollen donors are scarce (Baker 1955). Modern treatments emphasize that this is principally a colonization and establishment filter rather than a claim that island plants must evolve uniformly high selfing rates (Cheptou 2012; Pannell et al. 2015). Global comparative evidence is consistent with this logic: self-compatible species are over-represented on islands (Grossenbacher et al. 2017), and breeding system, lifespan and floral symmetry or generalization jointly predict island colonization probability when arrival opportunity is considered (Zell et al. 2025). These results make reproductive assurance a strong candidate for a recurrent island response. They do not, however, establish that flower colour, shape or accessibility should change in the same direction everywhere. Sister-taxon comparisons across Pacific islands, for example, do not support a universal evolutionary reduction in flower size and instead show substantial island- and lineage-specific variation (Hetherington-Rauth & Johnson 2020).

Floral phenotype also faces a second ecological problem. Selfing can relax selection on traits used to attract or mechanically match animal pollen vectors, generating the familiar selfing syndrome (Sicard & Lenhard 2011). But pollinator-facing traits can also change without being merely downstream consequences of selfing. Floral form, symmetry, tube depth and colour jointly influence which visitors can access rewards and contact reproductive organs, and multivariate floral phenotypes can predict functional pollinator groups even though they do not uniquely identify realized pollinators (Fenster et al. 2004; Rosas-Guerrero et al. 2014). This creates two competing models for island floral change. In a serial model, isolation favours reproductive assurance and floral reorganization follows largely because selfing reduces the value of pollinator-facing investment. In a partially separable model, reproductive assurance changes while floral accessibility or display also responds to the interaction environment. Distinguishing these models requires reproductive and floral traits to be analysed separately rather than collapsed into a single syndrome score.

A second independent evidence layer is needed to ask whether the geographic gradient itself carries a reproductive constraint. Pollen-supplementation experiments quantify the degree to which natural reproduction increases when additional pollen is supplied, integrating pollen quantity, quality, mate availability and pollinator service (Knight et al. 2005; Harder & Aizen 2010). The GloPL synthesis provides this outcome globally (Bennett et al. 2018a,b), and previous analyses show that pollen limitation varies with ecological context and pollinator dependence (Bennett et al. 2020). By contrast, a GloPL–GBIF synthesis found little general increase in pollen limitation towards species' range edges and excluded remote islands beyond 200 km from mainland from that range-edge analysis (Dawson-Glass & Hargreaves 2022). Geographic isolation from continental source therefore represents a distinct biogeographic axis on which to test whether pollen delivery becomes more constraining.

Here we combine global island floras with independent pollen-supplementation experiments to separate recurrent reproductive function from contingent floral phenotype. H1 estimates seven reproductive, colour and structural responses to isolation separately in four geographic regions, with no pooled floral-syndrome score. H2 tests whether floral accessibility and colour responses persist after conditioning on measured reproductive assurance, distinguishing a compulsory selfing-only sequence from partially separable responses. H3 asks whether experimental pollen limitation increases with geographic isolation. H4 then asks whether the reproductive-assurance and accessibility states enriched along the island gradient are associated with lower current pollen limitation in exact-species overlap. Together these analyses test whether island isolation is linked to a repeatable reproductive problem whose functional solutions recur more consistently than their visible floral expression.

## Materials and Methods

### Global island plant database

We used a corrected global analysis universe of **8,264 island units**. Candidate units originated from GSHHG 2.3.7 high-resolution island geometry. A geography audit identified one continental split component that had been retained as an island; excluding it yielded the corrected 8,264-unit universe and a broad H1 union of 4,379 islands. The trait database retained the same taxonomic and provenance structure, containing 106,295 accepted angiosperm species.

Observed island floras were assembled from GBIF occurrence records assigned to island units. Trait evidence was assembled source by source from floras, monographs, public trait resources and primary literature. Every accepted record retained source provenance. Species-direct high- and medium-confidence evidence had priority, while validated lower-confidence evidence was used only where direct evidence was unavailable. Missing trait information was retained as missing rather than converted to trait absence, and family-level inference was not used as a general fill rule.

The final trait dataset contained 222,688 resolved cells of 318,885 possible species-by-axis cells (69.83%) across three raw evidence axes: flower colour, floral structural complexity and reproductive assurance. The broad all-analysis-eligible evidence scope was used as the primary plant analysis, with species-direct high/medium evidence retained as a Direct-only sensitivity.

The primary H1 sampling frame deliberately retained all observed island-flora records, irrespective of floristic origin, to avoid discarding the majority of records whose origin status was unresolved. Among 1,039,757 island-by-species records in the analysis ledger, 150,531 were source-backed native, 19,442 were source-backed introduced and 869,784 (83.65%) had unresolved origin status. The all-observed analysis therefore estimates trait composition in the contemporary observed flora; it is not, by itself, an analysis of native colonization or in-situ evolution. WCVP regional-native compatibility is the sole active floristic-origin sensitivity.

Because those strict subsets were severely support-limited, we additionally performed a regional-native compatibility sensitivity using the World Checklist of Vascular Plants (WCVP; Govaerts et al. 2021). Accepted species were exact-matched to the WCVP bulk names table, and native distributions were read at TDWG level 3. An unresolved island-by-species record was upgraded to regional-native-compatible only when the focal island mapped unambiguously to one TDWG level-3 unit and WCVP listed that accepted species as native in the same unit. Source-backed introduced records were never overwritten. The WCVP archive used for this replay had SHA-256 `d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa`. This procedure establishes regional native compatibility rather than exact island-level nativeness. The earlier complementary-origin and TDWG-area analyses are archived provenance diagnostics and are not included as additional active flora partitions in this replay.

### Geographic isolation and covariates

The primary geographic exposure is the minimum separation between island and continental coastlines measured on a mean-radius sphere (R = 6371.0088 km). Island and continental geometries were taken from the same GSHHG 2.3.7 high-resolution archive. Continental siblings were reconstructed before distance calculation, artificial dateline split edges were omitted, and distance was calculated as the minimum minor-great-circle arc separation between coastline segments. This is a spherical coastline distance, not a WGS84 ellipsoidal distance.

The geography audit identified 1,113 broad-H1 islands that had been assigned spurious zero distance when GSHHG island polygons were compared with a coarser Natural Earth continental geometry. Under the source-matched coastline metric, all 1,113 have positive distance. This correction is a post-hoc measurement repair selected as the primary submission baseline, not a new prospective confirmation.

Plant models used standardized log distance, log island area and climate PC1-PC4. Spatial dependence was addressed with spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of source separation, connectivity and colonization opportunity rather than as a randomized treatment.

Four predeclared geographic replication strata were analysed separately: northern mid-latitude, northern high-latitude, tropical and southern extratropical islands. These strata were used to assess recurrence of response structure, not to require identical coefficients among regions.

For GloPL, geographic exposure was recomputed from each site's coordinates against the same continental geometry. Sites lying on seeded continental land were assigned zero distance; 996 of 1,248 sites are true continental zeros and remain zero in the corrected baseline.

### H1: seven separate trait responses

The primary analysis is broad contemporary island flora with All evidence. Each of seven binary traits is fitted separately in each of four regions by beta-binomial logit regression, adjusting for standardized corrected log isolation, log island area and climate PC1–4. Reproduction, colour and structure are labels, not aggregate scores. Spatial-block sandwich standard errors use a finite-cluster t reference with G−1 degrees of freedom and pointwise 95% intervals. All individual p values are two-sided and unadjusted; no Holm correction is applied. This final inferential specification was selected retrospectively on 4 October 2026 after comparison of alternatives and was not preregistered. Individual significance does not establish family-wise significance or formal differences between regional slopes.

WCVP regional-native-compatible flora is the sole active origin sensitivity; Direct-only is a separate evidence-quality sensitivity. Both use the identical model and inference rule. WCVP compatibility retains existing source-native records and upgrades unresolved records only under accepted TDWG-L3 native-range compatibility; introduced records are not overwritten. This establishes regional compatibility rather than exact focal-island nativity. H2–H4 and their previously defined functional covariates are unchanged.

### Post-baseline source-matched assembly decomposition

We added a post-baseline diagnostic for source-backed native non-endemics using four frozen GIFT mainland-source assignments. For each trait and island, the decomposition separated the source-species expectation, the contribution of represented genera and additional sorting among species within those genera. The identity closed exactly. Components used the H1 covariates and spatial-block covariance; source-mode tests were BH-adjusted and required at least 50 islands per outcome. This localizes assemblage composition but cannot distinguish arrival, establishment, persistence or extinction, or estimate within-species evolution.

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

### H1: recurrent reproductive responses with regional floral differences

All 112 fits converged. Primary All evidence supports 17 of 28 individual associations at nominal two-sided P < 0.05. Self-compatibility increases with isolation in all four regions (P = 0.00569, 0.00671, 0.01812 and 0.01123 in northern mid-latitude, northern high-latitude, tropical and southern extratropical floras, respectively). Selfing mating system increases in the first three regions. Autonomous selfing increases in northern high latitudes and the tropics. Generalized form increases in northern mid-latitudes, the tropics and southern extratropics; actinomorphy increases in northern high latitudes and the tropics. Shallow/open tubes increase in northern high latitudes but decrease in southern extratropics. Plain colour increases in southern extratropics.

In WCVP All, 10 of 28 individual associations meet the same nominal threshold. Increasing self-compatibility remains supported in all four regions. Northern-high autonomous selfing, actinomorphy and shallow/open tubes remain positive and supported; tropical actinomorphy remains supported. Plain colour increases in both tropical and southern extratropical WCVP floras. Several broad-flora associations, including the southern shallow/open-tube decrease, no longer meet the threshold. This can reflect changed composition, precision and coverage, rather than proving an introduced-species mechanism. Direct-only yields 15/28 supported associations in broad flora and 11/28 in WCVP flora. Complete coefficients, intervals and unadjusted P values are reported in Supplementary Table S2.

### Source-matched assembly diagnostic: the testable response is concentrated at genus composition

Source-matched support was adequate in northern mid-latitudes and the tropics only. In northern mid-latitudes, genus-structure enrichment was positive and FDR-supported in all four source modes (All +0.0163 to +0.0221, maximum q = 1.49 × 10^-5; Direct-only +0.0184 to +0.0233, maximum q = 5.54 × 10^-6), whereas additional within-genus species sorting was not robust across modes (maximum q = 0.136 and 0.130). In the tropics, the All source expectation was positive (+0.0185 to +0.0195), but genus structure shifted oppositely (-0.0202 to -0.0147, maximum q = 0.00722); Direct-only genus structure was likewise negative. Thus the testable compositional response was concentrated principally in genus representation, not additional sorting within represented genera.

### H2: floral accessibility is not reducible to reproductive assurance

The `selfing_core` isolation slope was positive in all four regions in the primary analysis, although nominal support remained concentrated in the tropical and southern strata.

After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. With finite-cluster t references and the retained eight-test BH family, q-values were 0.1764 in northern mid-latitudes, 0.00243 in northern high latitudes, 0.01913 in the tropics and 0.2977 in southern extratropical islands. Thus the primary evidence for a selfing-adjusted accessibility response is strongest in northern high latitudes and the tropics. In the Direct-only sensitivity, the tropical accessibility estimate remained nominally positive (finite-cluster p = 0.0475) but did not survive FDR correction (q = 0.1267); it is therefore not treated as a replicated FDR-supported result.

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

### Reproductive function recurs more consistently than floral phenotype

Reproductive traits, especially self-compatibility, recur across regions more consistently than visible floral traits. Self-compatibility increased with isolation in all four primary geographic strata and remained positive in all four WCVP regional-native-compatible analyses, whereas colour and structural responses were more geographically contingent.

This asymmetry shifts the meaning of a floral island syndrome from a universal flower type to a recurrent reproductive function. It is consistent with Baker-type filtering and with evidence that breeding system and floral generalization influence island colonization (Baker 1955; Pannell et al. 2015; Grossenbacher et al. 2017; Zell et al. 2025). The source-matched diagnostic sharpens that interpretation: the northern-midlatitude response is concentrated in genus composition, while in the tropics an increasingly H1-oriented source pool is opposed by realized genus composition. Thus the island pattern is neither a simple reflection of source availability nor evidence that every lineage changes in the same way after colonization.

Functional repeatability therefore need not imply phenotypic uniformity. Isolation may repeatedly favour or filter reproductive independence while regional source pools, interaction communities and colonization histories produce different floral forms.

### Floral reorganization is not exhausted by the selfing syndrome

Reproductive assurance provides one plausible route to reduced dependence on external pollen, but H2 shows that it is not a complete statistical explanation for floral reorganization. After conditioning on selfing_core, the isolation coefficient for generalized accessibility remained positive in all four primary regions and was FDR-supported in northern high latitudes and the tropics. When selfing_core was rebuilt using only species with all three reproductive components directly observed, the accessibility signal persisted in those regions. The Direct-only tropical result weakened after multiplicity correction, so the evidence does not support a universal second pathway. It does, however, reject the stronger claim that measured reproductive assurance necessarily absorbs all isolation-associated floral change.

This is the ecological distinction between a compulsory serial model and partially separable responses. Under the serial model, isolation increases reproductive assurance and floral change is primarily a downstream correlate of that shift. Under the alternative, some floral traits retain an association with isolation after reproductive assurance is held constant. Accessible architectures could broaden the set of visitors capable of contacting reproductive organs, alter mechanical matching, or arise through compositional filtering among lineages that differ in floral access. Our data cannot identify which of these mechanisms operated historically, but they show that the reproductive and pollinator-facing dimensions should not be treated as interchangeable.

The heterogeneous colour and architecture results reinforce this interpretation. Floral syndrome theory predicts that combinations of traits, rather than single colours or shapes, are associated with different functional pollen vectors (Fenster et al. 2004; Rosas-Guerrero et al. 2014). At the same time, those associations are not one-to-one maps from phenotype to realized pollinator. We therefore interpret the raw colour × architecture patterns as evidence that the morphological expression of insularity differs among regions, not as proof that a named pollinator group increased or disappeared. The ecological signal is stronger at the level of reproductive assurance and accessibility than at the level of a universal floral design.

### Isolation is also associated with an independent reproductive constraint

The plant trait data alone cannot show that pollen delivery becomes more difficult with island isolation. H3 addresses that gap with an independent experimental dataset. Across 2,969 pollen-supplementation experiments, pollen limitation increased with corrected geographic isolation. The offshore-only analysis remained positive, while a simple mainland-versus-offshore indicator was unsupported, indicating that the primary association was not generated solely by contrasting continental zero-distance sites with islands. The supplemental-only sensitivity was weaker and unsupported, so the result should be read as a positive primary gradient with measurement sensitivity rather than as an invariant law.

This independent result changes the interpretation of the trait patterns. Geographic isolation is not only associated with reproductive-assurance and accessibility traits; it is also associated with a contemporary experimental shortfall in natural pollen receipt. That does not imply a global decline in pollinator abundance. Pollen limitation can arise from pollen quantity, pollen quality, mate availability, visitor effectiveness or temporal unreliability, and the present data do not separate those processes. The appropriate inference is broader: increasing isolation is associated with a stronger reproductive-service constraint.

This geographic signal is also distinct from a generic range-margin effect. Previous global work found little evidence that pollen limitation rises systematically towards species' range edges (Dawson-Glass & Hargreaves 2022). Here the relevant axis is distance from continental source geography, and the positive offshore gradient suggests that island isolation captures a different ecological process from position within a species' mainland range.

### A constraint–response triangle links pattern, pressure and function

H4 adds a third side to the argument. In exact-species overlap, stronger reproductive-assurance and generalized-accessibility scores were each associated with lower current pollen limitation after accounting for geographic distance and study structure. Autonomous selfing showed the strongest atomic association, as expected for a trait that can directly provide seed production when external pollen is insufficient. Generalized form and actinomorphy were also associated with lower limitation, whereas self-compatibility alone was negative but imprecise.

Taken together, H1–H4 form a constraint–response triangle. First, isolation is associated with recurrent reproductive-assurance traits and selected accessibility responses. Second, isolation is independently associated with greater pollen limitation. Third, species expressing the island-associated reproductive and accessibility states tend to experience lower current pollen limitation. These three associations make the functional interpretation substantially stronger than a trait survey alone: the geographic gradient carries an experimentally measured reproductive constraint, and the traits enriched along that gradient are aligned with reduced contemporary limitation.

One edge of the triangle remains deliberately open. We do not observe historical pollen limitation causing selection, sorting or persistence that generated present-day island assemblages. The exact-species H4 comparison is post-hoc, the island flora data are cross-sectional, and ancestral trait states or temporal changes in pollination service are not available. Earlier prespecified moderation tests also did not show that the H2 trait states significantly buffer the isolation–pollen-limitation slope itself. The present study therefore supports functional compatibility, not causal mediation.

### Island syndromes may be more repeatable in function than in form

A functional interpretation helps reconcile why island studies often recover some elements of a reproductive syndrome without converging on one floral morphology. Isolation can repeatedly increase the value of reproductive independence while the floral consequences depend on regional source pools, climate, interaction communities and the traits already available for colonization. Under that view, regional heterogeneity is not noise around an otherwise universal flower type; it is part of the biological result.

This perspective also clarifies how island macroecology and pollination biology connect. Functional island biogeography emphasizes that insularity filters colonization and persistence through traits rather than acting directly on species richness alone (Schrader et al. 2021). Our results suggest that plant reproduction follows the same logic at two levels. A recurrent functional axis—reduced dependence on uncertain compatible pollen delivery—is visible across broad geography, while the structural and colour traits through which that function is realized vary among regions. The next mechanistic step is therefore not to infer a universal pollinator from floral phenotype, but to measure realized visitor communities, pollen transfer, mating outcomes and lineage histories in systems where the macroecological response is strongest.

### Limitations and inferential boundary

The primary plant analyses describe contemporary island-flora composition. The source-matched diagnostic separates source composition, genus structure and additional within-genus species sorting only in the two contexts with adequate support; it cannot distinguish arrival from establishment, persistence or extinction. Because the active trait state is one value per species, within-species phenotypic change and within-lineage evolution remain unidentified.

Floristic origin is also incompletely resolved: 83.65% of island-by-species records lacked source-backed origin status. WCVP regional-native compatibility preserves the four regional self-compatibility associations but establishes compatibility only at TDWG level 3, not exact focal-island nativeness. Broad–WCVP differences therefore do not identify an introduced-species mechanism.

H1 uses pointwise nominal tests under a retrospective specification, not family-wise proof of a universal syndrome. The geographic exposure is a post-hoc repair of a known geometry mismatch and remains a composite proxy for connectivity and source separation. Finally, H2 is conditional decomposition and H4 is post-hoc functional triangulation; neither identifies historical pollinator selection. Stronger causal inference requires lineage-resolved, temporal or experimental data linking realized pollination, establishment and trait change.

## Conclusion

Island isolation is associated with recurrent reproductive assurance without a universal floral phenotype. Self-compatibility is the most geographically repeatable response, while other reproductive and floral traits vary among regions. Where source support is adequate, the assemblage response is concentrated in genus composition rather than additional sorting within represented genera. Isolation also predicts stronger experimental pollen limitation, and reproductive-assurance and accessibility states are associated with lower current limitation. The result is a recurrent reproductive problem assembled through regionally different solutions; whether island populations additionally change within lineages remains unresolved.

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
