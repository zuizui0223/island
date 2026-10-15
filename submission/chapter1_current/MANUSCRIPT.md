# Do pollen-supply constraints underlie floral island syndromes?

## Abstract

Island floral syndromes are often attributed to reduced pollinator service and greater reproductive assurance, yet the proposed reproductive constraint is rarely tested alongside global floral patterns. We combined a 106,295-species trait database (4,379 sampled islands) with 2,969 previously published GloPL pollen-supplementation experiments. Pollen limitation increased with continental isolation (β = 0.0919, P = 0.0159), although supplemental-only sensitivity was unsupported. In island assemblages, self-compatibility increased with isolation in all four regions; mating-system, colour and structural patterns differed. Floral accessibility retained an isolation association in selected regions after controlling for measured reproductive assurance. Exact-species matching showed post-hoc negative associations of reproductive-assurance and accessibility scores with current pollen limitation. Together these non-paired data are compatible with a pollination-related constraint, but do not establish global pollinator decline, a uniform dulling syndrome, within-lineage floral evolution or historical mediation.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; functional island biogeography; pollination

## Introduction

A familiar explanation for floral island syndromes runs from geographic isolation to fewer or less effective pollinators, then to reproductive assurance and reduced floral attraction. Floral trait patterns alone, however, do not demonstrate insufficient compatible pollen receipt. Nor can pale colour or accessible architecture establish pollinator decline. A test of the classical explanation must distinguish the *ecological constraint* from the reproductive and floral patterns it is invoked to explain.

Pollen-supplementation experiments provide a direct reproductive endpoint: how much does adding pollen improve reproduction relative to natural pollen receipt (Knight et al. 2005; Harder & Aizen 2010)? Such pollen limitation (PL) is closer to the experienced reproductive constraint than pollinator-occurrence proxies, yet also reflects compatible-mate availability, pollen quality and the plant's existing traits; it is not a measure of historical pollinator loss. We use previously published GloPL experiments to examine continental isolation against PL (Bennett et al. 2018a,b). Earlier GloPL analyses related PL to autofertility, floral generalization and pollinator dependence (Burns et al. 2019; Bennett et al. 2020). A global range-edge test addressed a different geographic exposure (Dawson-Glass & Hargreaves 2022).

Baker's law predicts establishment advantages for uniparental reproduction when mates are scarce (Baker 1955; Cheptou 2012; Pannell et al. 2015), and island self-compatibility and mutualism filtering are already documented (Grossenbacher et al. 2017; Razanajatovo et al. 2019; Zell et al. 2025; Delavaux et al. 2024). Yet self-compatibility, autonomous selfing and realized mating system are not interchangeable. In the classical selfing syndrome, increased reproductive independence relaxes selection on floral attraction or mechanical specialization (Sicard & Lenhard 2011). Alternatively, pollinator-facing architecture may shift separately, or assemblage sorting may change observed trait composition without within-lineage evolution. Floral colour, symmetry, form and tube depth do not uniquely identify effective visitors (Fenster et al. 2004; Rosas-Guerrero et al. 2014); island sister-taxon comparisons have also found no universal floral-size reduction (Hetherington-Rauth & Johnson 2020).

We test three components of this explanation using *different, complementary data*: **(1)** whether PL increases with continental isolation; **(2)** which reproductive and floral states covary with isolation across island assemblages and whether selfing adjustment absorbs the floral patterns; and **(3)** whether selected trait states are linked to lower *current* PL. The seven traits, four regions and original H1–H4 estimands are preserved. These analyses assess empirical compatibility with the pollination-constraint hypothesis, not historical causal mediation.

## Materials and Methods

Source-analysis identifiers H1–H4 remain unchanged. The narrative sequence is ecological constraint (source H3), reproductive and floral responses (sources H1 and H2) and post-hoc functional compatibility (source H4); the datasets are not an integrated causal time series.

### Global island plant database

We used a corrected global analysis universe of **8,264 island units** (Figure 1). Candidate units originated from GSHHG 2.3.7 high-resolution island geometry. A geography audit identified one continental split component that had been retained as an island; excluding it yielded the corrected 8,264-unit universe and a broad H1 union of 4,379 islands. The trait database retained the same taxonomic and provenance structure, containing 106,295 accepted angiosperm species.

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

### Published pollen-supplementation experiments and isolation (source H3)

We used the version-pinned GloPL database of pollen-supplementation experiments (Bennett et al. 2018a,b). Pollen limitation was represented by the log response ratio of reproduction after supplemental pollen addition to reproduction under natural pollen receipt; positive values therefore indicate that reproduction increases when pollen is added.

The analysis contained 2,969 experimental rows, 1,408 measurement cells, 1,248 sites and 919 publications. Duplicate measurements were aggregated at the publication-by-coordinate-by-measurement level. Each publication contributed total analysis weight one. The primary model related pollen limitation to standardized corrected log geographic distance while controlling for broad geographic context and experimental measurement conditions. Uncertainty was estimated with publication-cluster-robust covariance.

Because 996 of 1,248 GloPL sites are true continental zero-distance sites, we added post-hoc robustness checks to distinguish a continuous offshore gradient from a simple mainland/offshore step. We refitted the model among positive-distance measurement cells only, fitted a separate binary mainland-versus-offshore indicator model, and repeated the offshore-only gradient after deleting each publication in turn. These are robustness analyses rather than new confirmatory tests.

### Island reproductive and floral responses

#### Seven distinct trait responses (source H1)

The primary analysis is broad contemporary island flora with All evidence. Each of seven binary traits is fitted separately in each of four regions by beta-binomial logit regression, adjusting for standardized corrected log isolation, log island area and climate PC1–4. Reproduction, colour and structure are labels, not aggregate scores: there is no pooled floral-syndrome score. Spatial-block sandwich standard errors use a finite-cluster t reference with G−1 degrees of freedom and pointwise 95% intervals. All individual p values are two-sided and unadjusted; no Holm correction is applied. This final inferential specification was selected retrospectively on 4 October 2026 after comparison of alternatives and was not preregistered. Individual significance does not establish family-wise significance or formal differences between regional slopes.

WCVP regional-native-compatible flora is the sole active origin sensitivity; Direct-only is a separate evidence-quality sensitivity. Both use the identical model and inference rule. WCVP compatibility retains existing source-native records and upgrades unresolved records only under accepted TDWG-L3 native-range compatibility; introduced records are not overwritten. This establishes regional compatibility rather than exact focal-island nativity. H2–H4 and their previously defined functional covariates are unchanged.

#### Floral responses conditional on reproductive assurance (source H2)

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

### Functional associations with contemporary pollen limitation (source H4)

H4 was an explicitly post-hoc functional triangulation analysis. Species in the H2 trait dataset were exact-matched to GloPL species after only case normalization and underscore/space normalization; synonym rescue and genus fallback were not used.

The primary H4 predictors were the literal Direct-only species-level H2 scores used in the plant analysis: `selfing_core` and `generalized_accessible`. Both are soft-membership scores scaled from 0 to 1, with higher values representing stronger concordance with reproductive assurance or generalized accessibility, respectively. Missing component traits were not coded as zero.

For each H2 score, we fitted:

`pollen limitation ~ H2 trait score + standardized distance + geographic context + measurement conditions`.

Each publication again received total weight one and uncertainty was publication-cluster robust. A negative H2-score coefficient indicates that species expressing the island-associated trait state more strongly tend to experience lower current experimental pollen limitation, conditional on geographic distance and measured study structure.

We also retained atomic-trait sensitivities. Because the global colour response was geographically heterogeneous and lacked a single directional functional prediction, colour was not promoted into the global H4 bridge.

## Results

### Constraint: increased pollen limitation with isolation (source H3)

Across 2,969 pollen-supplementation experiments, the corrected standardized distance coefficient was positive (**beta = 0.09191, SE = 0.03806; finite-publication two-sided p = 0.01594; one-sided positive p = 0.00797**). The no-zero-constant sensitivity also retained a supported positive coefficient (beta = 0.09089, p = 0.01869).

The supplemental-only sensitivity was positive but unsupported (beta = 0.04410, p = 0.31082). The primary result therefore supports an isolation-associated increase in experimental pollen limitation, while the evidence is not equally strong under every measurement definition (Figure 2).

The post-hoc offshore-only analysis showed that this result was not generated solely by the contrast between 996 mainland zero-distance sites and offshore sites. Among 276 positive-distance measurement cells from 153 publications and 252 sites, pollen limitation still increased with isolation (**beta = 0.22031, SE = 0.09704, finite-publication p = 0.02459**). By contrast, a binary mainland-versus-offshore indicator alone was unsupported (beta = 0.12743, p = 0.14364). Leave-one-publication-out refits kept all 153 offshore-gradient estimates positive; the weakest two-sided result was p = 0.04892. These checks support a continuous offshore gradient while remaining post-hoc robustness analyses.

### Reproductive and floral patterns and conditional decomposition

#### Reproductive and floral responses of the observed island flora (source H1)

All 112 fits converged. Primary All evidence supports 17 of 28 individual associations at nominal two-sided P < 0.05. Self-compatibility increases with isolation in all four regions (P = 0.00569, 0.00671, 0.01812 and 0.01123 in northern mid-latitude, northern high-latitude, tropical and southern extratropical floras, respectively). Selfing mating system increases in the first three regions. Autonomous selfing increases in northern high latitudes and the tropics. Generalized form increases in northern mid-latitudes, the tropics and southern extratropics; actinomorphy increases in northern high latitudes and the tropics. Shallow/open tubes increase in northern high latitudes but decrease in southern extratropics. Plain colour increases in southern extratropics.

In WCVP All, 10 of 28 individual associations meet the same nominal threshold. Increasing self-compatibility remains supported in all four regions. Northern-high autonomous selfing, actinomorphy and shallow/open tubes remain positive and supported; tropical actinomorphy remains supported. Plain colour increases in both tropical and southern extratropical WCVP floras. Several broad-flora associations, including the southern shallow/open-tube decrease, no longer meet the threshold. This can reflect changed composition, precision and coverage, rather than proving an introduced-species mechanism. Direct-only yields 15/28 supported associations in broad flora and 11/28 in WCVP flora. The traitwise coefficient display is shown in Figure 3; complete estimates, intervals and unadjusted P values are reported in Supplementary Table S2.

#### Does measured reproductive assurance absorb the floral gradient? (source H2)

The `selfing_core` isolation slope was positive in all four regions in the primary analysis, although nominal support remained concentrated in the tropical and southern strata.

After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four primary strata. With finite-cluster t references and the retained eight-test BH family, q-values were 0.1764 in northern mid-latitudes, 0.00243 in northern high latitudes, 0.01913 in the tropics and 0.2977 in southern extratropical islands. Thus the primary evidence for a selfing-adjusted accessibility response is strongest in northern high latitudes and the tropics. In the Direct-only sensitivity, the tropical accessibility estimate remained nominally positive (finite-cluster p = 0.0475) but did not survive FDR correction (q = 0.1267); it is therefore not treated as a replicated FDR-supported result.

The complete-three-component Direct-only sensitivity reduced the selfing mediator to 564 species but retained a positive isolation coefficient for `generalized_accessible` in all four regions. FDR-adjusted q-values were 0.0508 in northern mid-latitudes, 0.0376 in northern high latitudes, 0.0376 in the tropics and 0.2706 in southern extratropical islands. The corresponding fitted island counts were 1,725, 221, 438 and 157. Thus the conditional accessibility signal in northern high latitudes and the tropics persists when every species contributing to the reproductive-assurance score has all three selfing components directly observed, which weakens a missing-component explanation without eliminating general measurement-error concerns.

Colour responses remained more heterogeneous. Southern adjusted plain colour remained supported after finite-cluster correction (q = 0.00243 in all-analysis and 1.8 × 10^-5 in Direct-only). Raw colour, architecture and colour-conditioned architecture analyses showed region-specific response patterns, including negative northern-high-latitude blue/purple architecture associations and a positive tropical Direct yellow/orange × butterfly/deep-tube pattern. These raw pattern analyses are descriptive and do not establish realized pollinator identity.

The corrected H2 results therefore support a partially separable floral response: accessibility estimates remained positive after reproductive-assurance adjustment in all four regions, but FDR support was concentrated in northern high latitudes and the tropics. The strength and detailed form of the additional floral response are geographically contingent (Figure 4).

### Functional consistency: reproductive and floral traits versus current PL (source H4)

At the atomic-trait level, autonomous selfing provided the strongest functional association. Species with autonomous or delayed selfing had lower current pollen limitation after adjustment for distance, broad context and measurement structure (beta = -0.44492, finite-publication p = 4.85 × 10^-8). Actinomorphy was also strongly negative (beta = -0.38076, finite-publication p = 1.47 × 10^-5), generalized form was weaker but supported (beta = -0.18448, finite-publication p = 0.04543), and self-compatibility was negative but imprecise (beta = -0.11885, finite-publication p = 0.1502).

The literal H2-to-H4 bridge gave the same overall interpretation. The Direct-only `selfing_core` score matched 455 GloPL species across 409 publications and was negatively associated with current pollen limitation (**beta = -0.29830, SE = 0.10352, finite-publication p = 0.00417**). The association remained supported in the no-zero-constant sensitivity, whereas the supplemental-only subset remained negative but unsupported (finite-publication p = 0.2907).

The Direct-only `generalized_accessible` score matched 143 GloPL species across 143 publications. Higher accessibility was associated with lower current pollen limitation (**beta = -0.29566, SE = 0.12896, finite-publication p = 0.02334**). The supplemental-only sensitivity was also negative (finite-publication p = 0.04239), and the no-zero-constant sensitivity remained supported (p = 0.02082).

In different exact-species subsets, the selected reproductive-assurance and accessibility scores were each associated with lower current pollen limitation (Figure 5). This functional alignment does not establish that historical pollen limitation caused the observed trait distributions.

## Discussion

### What the pollen-limitation gradient can tell us about isolation

The plant trait data alone cannot show whether successful pollen receipt becomes more limiting with island isolation. The published GloPL experiment synthesis (source H3) addresses this gap through the reproductive response to added pollen rather than pollinator counts. Across 2,969 pollen-supplementation experiments, pollen limitation increased with corrected geographic isolation. The offshore-only analysis remained positive, while a simple mainland-versus-offshore indicator was unsupported, indicating that the primary association was not generated solely by contrasting continental zero-distance sites with islands. The supplemental-only sensitivity was weaker and unsupported, so the result should be read as a positive primary gradient with measurement sensitivity rather than as an invariant law.

An isolation-associated deficit in effective pollination need not reflect fewer pollinator individuals. Oceanic-island networks have fewer interaction types than mainland networks (Traveset et al. 2016). On Japanese islands, reduced long-tongued pollinator biomass coincided with compensatory visits by other guilds but lower fruit set in a focal plant (Hiraiwa & Ushimaru 2017). Across 40 coastal networks, low pollinator functional diversity, rather than species richness, was associated with weaker flower–proboscis matching and pollination success (Hiraiwa & Ushimaru 2024). Thus functional guild replacement could contribute to pollen limitation and contingent floral-access responses. Yet H3 measures neither pollinator identity nor effectiveness, and pollen limitation also reflects compatible-mate availability and pollen quality. Spatial contrasts cannot establish temporal pollinator decline.

This geographic signal is also distinct from a generic range-margin effect. Previous global work found little evidence that pollen limitation rises systematically towards species' range edges (Dawson-Glass & Hargreaves 2022). Here the relevant axis is distance from continental source geography, and the positive offshore gradient suggests that island isolation captures a different ecological process from position within a species' mainland range.

Current PL can itself depend on breeding traits already present in the sampled plant: autonomous selfing may diminish the additional-pollen response despite a history of poor pollinator service. PL therefore cannot be inserted into a contemporaneous mediation model as though it measured pre-evolutionary pollination pressure. The geographic correlations use non-identical plant and experimental-site populations.

### Distinct reproductive and floral responses: a selfing syndrome or parallel pathways?

The clearest result is an asymmetry in repeatability. Reproductive traits, especially self-compatibility, recur across regions more consistently than visible floral traits. In the primary All-evidence analysis, self-compatibility increased with isolation in all four geographic strata, and the same four positive associations were retained in the WCVP regional-native-compatible sensitivity. Selfing mating system and autonomous selfing also increased in several regions, although less uniformly. By contrast, colour and structural responses were more geographically contingent: generalized form and actinomorphy increased in several regions, shallow/open tubes changed in opposite directions between northern high-latitude and southern extratropical floras, and plain colour was concentrated in the south and in tropical WCVP-compatible floras.

This pattern shifts the ecological meaning of a floral island syndrome. The most geographically repeatable individual association is self-compatibility; the other reproductive traits do not constitute a uniformly responding syndrome. That result is consistent with Baker-type filtering, in which reproductive independence can improve establishment under uncertain mate or pollen availability (Baker 1955; Pannell et al. 2015; Grossenbacher et al. 2017). It is also consistent with recent evidence that breeding system and floral generalization jointly influence island colonization probability (Zell et al. 2025). Our analysis extends those studies from the binary question of island occurrence to continuous variation among already assembled island floras: increasing source isolation continues to be associated with the representation of reproductive-assurance traits, whereas the floral phenotype through which insularity is expressed remains region dependent.

The distinction matters because functional repeatability does not require phenotypic uniformity. If isolated floras repeatedly face uncertain compatible pollen delivery, self-compatibility may be favoured or filtered in many regions even when the local interaction community, source flora and colonization history produce different colour or structural outcomes. The data argue against a universal morphological checklist for the island syndrome. A more defensible formulation is a recurrent reproductive function with geographically contingent floral realization.

Reproductive assurance provides one plausible route to reduced dependence on external pollen, but the source H2 conditional models show that measured reproductive assurance does not statistically absorb every floral–isolation association. After conditioning on selfing_core, the isolation coefficient for generalized accessibility remained positive in all four primary regions and was FDR-supported in northern high latitudes and the tropics. When selfing_core was rebuilt using only species with all three reproductive components directly observed, the accessibility signal persisted in those regions. The Direct-only tropical result weakened after multiplicity correction, so the evidence does not support a universal second pathway. It does, however, reject the stronger claim that measured reproductive assurance necessarily absorbs all isolation-associated floral change.

This is the ecological distinction between a compulsory serial model and partially separable responses. Under the serial model, isolation increases reproductive assurance and floral change is primarily a downstream correlate of that shift. Under the alternative, some floral traits retain an association with isolation after reproductive assurance is held constant. Accessible architectures could broaden the set of visitors capable of contacting reproductive organs, alter mechanical matching, or arise through compositional filtering among lineages that differ in floral access. Our data cannot identify which of these mechanisms operated historically, but they show that the reproductive and pollinator-facing dimensions should not be treated as interchangeable.

The regional colour–architecture results show partial rather than uniform concordance with classical pollination syndromes. Corrected northern analyses include declines in specialized or large-bee-associated architecture with isolation. Full-range tropical Direct evidence also includes a positive yellow/orange–deep-tube association, while southern signals are internally mixed. Because regional isolation ranges differ, these patterns should not be read as like-for-like biome contrasts. They are phenotype–syndrome associations, not one-to-one maps to realized pollinators or evidence that named pollinator groups increased, disappeared or replaced one another (Fenster et al. 2004; Rosas-Guerrero et al. 2014).

A more accessible floral architecture is not necessarily less conspicuous or cheaper to produce. Conversely, plain colour is not a direct measure of investment, pollinator preferences or floral reward. We therefore compare serial selfing-syndrome change, partly separable pollinator-facing responses, and colonization or lineage sorting as candidate mechanisms; the cross-sectional analyses do not choose uniquely among them.

### Functional compatibility without a demonstrated historical cause

H4 adds a third side to the argument, although the direction of trait–pollen-limitation association is not itself novel (Burns et al. 2019; Bennett et al. 2020). In exact-species overlap, stronger reproductive-assurance and generalized-accessibility scores were each associated with lower current pollen limitation after accounting for geographic distance and study structure. Autonomous selfing showed the strongest atomic association, as expected for a trait that can directly provide seed production when external pollen is insufficient. Generalized form and actinomorphy were also associated with lower limitation, whereas self-compatibility alone was negative but imprecise.

Taken together, H1–H4 form a constraint–response triangle. First, isolation is associated with recurrent reproductive-assurance traits and selected accessibility responses. Second, isolation is independently associated with greater pollen limitation. Third, species expressing the island-associated reproductive and accessibility states tend to experience lower current pollen limitation. Together these associations extend the plant-trait analysis with an experimentally measured reproductive endpoint across comparable, but not identical, continental-isolation axes. The traits enriched along that gradient also align with reduced contemporary limitation, without requiring identical sampling units.

One edge of the triangle remains deliberately open. We do not observe historical pollen limitation causing selection, sorting or persistence that generated present-day island assemblages. The exact-species H4 comparison is post-hoc and reuses GloPL evidence from H3, rather than constituting a second independent experimental replication. Island flora data are cross-sectional, and ancestral traits or historical pollination service are not observed. Earlier prespecified moderation tests also did not show that the H2 trait states significantly buffer the isolation–pollen-limitation slope itself. The present study therefore supports functional compatibility, not causal mediation.

### Limitations and inferential boundary

Several limits are central to interpretation. First, because the plant-trait states are defined at the species level, their associations with isolation are necessarily compositional and cannot separate three nested processes: between-lineage assemblage sorting, within-species founder-genotype filtering, and post-colonization evolution. Even same-species island–mainland divergence cannot distinguish the latter two, because founders may already carry island-associated states. Resolving this hierarchy requires population-level reproductive measurements with source-matched sampling and, ultimately, genomic tests of whether causal variants predate or postdate colonization.

Second, GBIF occurrences are not complete floristic inventories. Absence from GBIF was never treated as species absence, and H1 therefore estimates trait composition among observed island records rather than complete island floras. Isolation slopes could be biased if observation probability covaries jointly with geographic isolation and reproductive or floral traits; the present data cannot fully exclude geographically structured collecting effort or detectability.

Third, floristic origin is incompletely resolved. The broad analysis intentionally retains contemporary observed flora because 83.65% of island-by-species records lacked source-backed origin status. WCVP regional-native compatibility provides a stricter sensitivity and preserves the four regional self-compatibility associations, but it establishes native compatibility only at TDWG level 3, not exact focal-island nativeness. Differences between broad and WCVP results can reflect composition, coverage and precision and do not by themselves identify an introduced-species mechanism.

Fourth, H1 uses pointwise, nominal two-sided tests for seven traits across four regions under a retrospective final specification. These tests expose trait-specific recurrence and exceptions, but they do not constitute family-wise proof that an entire syndrome is significant. The ecological argument therefore rests on the pattern of effect sizes, recurrence and independent H2–H4 evidence rather than on counting nominally significant coefficients alone.

Fifth, the geographic exposure was repaired after a geometry audit. Using source-matched GSHHG coastlines corrects a known mismatch that had produced 1,113 spurious island zero distances, but the repair is a post-hoc measurement correction rather than a prospective confirmation. Geographic distance also remains a composite proxy for connectivity and source separation rather than a randomized treatment.

Finally, H2 is conditional decomposition rather than mediation, and H4 is post-hoc functional triangulation. Neither analysis identifies historical selection by pollinators. Stronger causal inference will require lineage-resolved comparisons or temporal and experimental data that connect realized pollination service to reproductive success, establishment and trait change.

## Conclusion

The classic pollen-supply-constraint hypothesis has distinguishable empirical components. Published pollen-addition experiments show a positive isolation–PL gradient in the primary analysis, with measurement-dependent uncertainty. Island plant assemblages show a positive self-compatibility association in four regions but no single direction across the other reproductive, colour and architectural states. Selected floral accessibility associations persist after adjustment for measured reproductive assurance. Present-day trait–PL links are post-hoc and cannot establish pollinator decline, selection on duller flowers, historical mediation or within-lineage adaptation. The broad trait database supplies the principal macroecological result; GloPL supplies complementary, independent evidence about the ecological constraint.

## References

Abe, T. (2006). Threatened pollination systems in native flora of the Ogasawara (Bonin) Islands. *Ann. Bot.*, 98, 317–334.

Baker, H.G. (1955). Self-compatibility and establishment after ‘long-distance’ dispersal. *Evolution*, 9, 347–349.

Bennett, J.M., Steets, J.A., Burns, J.H., Burkle, L.A., Vamosi, J.C., Wolowski, M. et al. (2020). Land use and pollinator dependency drives global patterns of pollen limitation in the Anthropocene. *Nat. Commun.*, 11, 3999.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018a). Data from: GloPL, a global data base on pollen limitation of plant reproduction. Dryad Digital Repository. Available at: https://doi.org/10.5061/dryad.dt437.

Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018b). GloPL, a global data base on pollen limitation of plant reproduction. *Sci. Data*, 5, 180249.

Burns, J.H., Bennett, J.M., Li, J., Xia, J., Arceo-Gómez, G., Burd, M. et al. (2019). Plant traits moderate pollen limitation of introduced and native plants: a phylogenetic meta-analysis of global scale. *New Phytol.*, 223, 2063–2075.

Cheptou, P.-O. (2012). Clarifying Baker's Law. *Ann. Bot.*, 109, 633–641.

Dawson-Glass, E. & Hargreaves, A.L. (2022). Does pollen limitation limit plant ranges? Evidence and implications. *Philos. Trans. R. Soc. B*, 377, 20210014.

Delavaux, C.S., Crowther, T.W., Bever, J.D., Weigelt, P. & Gora, E.M. (2024). Mutualisms weaken the latitudinal diversity gradient among oceanic islands. *Nature*, 627, 335–339.

Fenster, C.B., Armbruster, W.S., Wilson, P., Dudash, M.R. & Thomson, J.D. (2004). Pollination syndromes and floral specialization. *Annu. Rev. Ecol. Evol. Syst.*, 35, 375–403.

Govaerts, R., Nic Lughadha, E., Black, N., Turner, R. & Paton, A. (2021). The World Checklist of Vascular Plants, a continuously updated resource for exploring global plant diversity. *Sci. Data*, 8, 215.

Grossenbacher, D.L., Brandvain, Y., Auld, J.R., Burd, M., Cheptou, P.-O., Conner, J.K. et al. (2017). Self-compatibility is over-represented on islands. *New Phytol.*, 215, 469–478.

Harder, L.D. & Aizen, M.A. (2010). Floral adaptation and diversification under pollen limitation. *Philos. Trans. R. Soc. B*, 365, 529–543.

Hetherington-Rauth, M.C. & Johnson, M.T.J. (2020). Floral trait evolution of angiosperms on Pacific islands. *Am. Nat.*, 196, 87–100.

Hiraiwa, M.K. & Ushimaru, A. (2017). Low functional diversity promotes niche changes in natural island pollinator communities. *Proc. R. Soc. B*, 284, 20162218.

Hiraiwa, M.K. & Ushimaru, A. (2024). Loss of functional diversity rather than species diversity of pollinators decreases community-wide trait matching and pollination function. *Funct. Ecol.*, 38, 1296–1308.

Knight, T.M., Steets, J.A., Vamosi, J.C., Mazer, S.J., Burd, M., Campbell, D.R. et al. (2005). Pollen limitation of plant reproduction: pattern and process. *Annu. Rev. Ecol. Evol. Syst.*, 36, 467–497.

Pannell, J.R., Auld, J.R., Brandvain, Y., Burd, M., Busch, J.W., Cheptou, P.-O. et al. (2015). The scope of Baker's law. *New Phytol.*, 208, 656–667.

Razanajatovo, M., van Kleunen, M., Kreft, H., Dawson, W., Essl, F., Pergl, J., Pyšek, P., Winter, M. & Weigelt, P. (2019). Autofertility and self-compatibility moderately benefit island colonization of plants. *Glob. Ecol. Biogeogr.*, 28, 341–352.

Rosas-Guerrero, V., Aguilar, R., Martén-Rodríguez, S., Ashworth, L., Lopezaraiza-Mikel, M., Bastida, J.M. et al. (2014). A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecol. Lett.*, 17, 388–400.

Schrader, J., Wright, I.J., Kreft, H. & Westoby, M. (2021). A roadmap to plant functional island biogeography. *Biol. Rev.*, 96, 2851–2870.

Sicard, A. & Lenhard, M. (2011). The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Ann. Bot.*, 107, 1433–1443.

Traveset, A., Tur, C., Trøjelsgaard, K., Heleno, R., Castro-Urgal, R. & Olesen, J.M. (2016). Global patterns of mainland and insular pollination networks. *Glob. Ecol. Biogeogr.*, 25, 880–890.

Zell, A.N., Miranda, C.H., Grady, E.L., Grossenbacher, D.L. & Igić, B. (2025). Island colonization in flowering plants is determined by the interplay of breeding system, lifespan, floral symmetry, and arrival opportunity. *New Phytol.*, 245, 420–432.
