# A recurrent global floral island syndrome separates selfing and pollinator-facing floral responses under increasing pollen limitation

## Full working manuscript v14 candidate — reordered H1–H4 analysis — 2026-09-18

### Abstract

Island plants are often expected to become more self-reliant and less dependent on specialized pollination as geographic isolation increases, but whether reproductive, colour and structural responses form one recurrent global syndrome remains unresolved. We combined a fixed universe of 8,265 islands and 106,295 accepted angiosperm species with an independent global database of pollen-supplementation experiments. We organized the analysis as a four-step test: **H1**, whether isolation is associated with a recurrent floral/reproductive island syndrome; **H2**, whether the floral component is reducible to a selfing-syndrome route or retains a selfing-independent pollinator-facing component; **H3**, whether experimental pollen limitation independently increases with isolation; and **H4**, whether syndrome trait states are functionally associated with lower current pollen limitation.

In H1, seven identically oriented atomic responses spanning reproductive assurance, plain/inconspicuous colour and floral accessibility showed supported multivariate isolation responses in all four predeclared geographic replication strata in both evidence scopes. Equal-weight three-domain descriptive orientation was positive in all four strata. Colour was not uniform: the plain-colour coefficient was approximately zero in northern mid-latitudes but positive in northern high-latitude, tropical and southern extratropical strata.

H2 separated two non-exclusive plant responses. Reproductive assurance (`selfing_core`) increased with isolation in all four strata. More importantly, the isolation coefficient for generalized floral accessibility remained positive in all four strata after conditioning on `selfing_core`, with the clearest support in northern high-latitude and tropical islands. Selfing-adjusted raw flower-colour composition changed in three of four contexts, and raw architecture conditional on raw colour also changed in several contexts: northern-high-latitude blue/purple flowers lost specialized/deep architecture, whereas tropical Direct evidence showed increased yellow/orange coupling to a deep-tube component. The selfing-independent component is therefore supported through observed display–architecture reorganization rather than a weighted pollinator score or identified visitor guild.

Independent experimental evidence then supported H3. Across 2,969 GloPL experiments from 1,248 sites and 919 publications, pollen limitation increased with distance from major continental landmasses (standardized slope `0.0794 ± 0.0377`; two-sided `p=0.0354`, one-sided positive `p=0.0177`). In H4, post-hoc exact-species triangulation showed the strongest association for autonomous selfing (`β=-0.4467`, `p=2.60×10^-8`), with additional but weaker global concordance for actinomorphy and generalized floral form.

Together, these results support a recurrent global island-syndrome direction containing reproductive, colour and structural components, but they do not support one compulsory serial pathway from isolation to selfing and only then to floral simplification. The strongest selfing-independent signal lies in floral accessibility rather than a named pollination syndrome. Pollen limitation provides an independent ecological-pressure pattern and current trait–pollen-limitation associations provide functional triangulation, while historical causal mediation remains unproven.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; pollination syndrome; GloPL

---

## Introduction

Islands provide a natural setting in which dispersal, reproduction and species interactions are simultaneously challenged. Geographic separation reduces connectivity with continental source regions, limits repeated immigration and can make mates and mutualists unreliable. These pressures motivated the expectation that island plants should increasingly possess traits that permit establishment and reproduction without dependable partners. Baker's law formalized one part of this expectation: uniparental reproduction can be advantageous when colonists experience mate limitation or uncertain pollination (Baker 1955; Pannell & Barrett 1998). Global comparisons subsequently showed that self-compatible plants are over-represented on islands and that breeding system interacts with other traits during island colonization (Grossenbacher et al. 2017; Zell et al. 2025).

The same logic extends beyond breeding system. Floral symmetry, tube depth, form, colour and display can increase compatibility with particular functional groups, but specialized architectures can also create dependence on the availability and effectiveness of suitable visitors. Pollination-syndrome theory therefore predicts correlated combinations of floral traits rather than one diagnostic character (Fenster et al. 2004; Rosas-Guerrero et al. 2014). If successful pollination becomes unreliable on increasingly isolated islands, at least two responses are plausible. Plants may gain from **reproductive assurance**, through self-compatibility, autonomous selfing or predominantly selfing mating systems. Independently, the benefit of maintaining restricted floral access may decline, favouring more open, generalized or radially accessible flowers.

These routes are related but need not form one obligatory sequence. The selfing syndrome includes recurrent reductions in floral investment associated with transitions toward selfing (Sicard & Lenhard 2011), but reproductive assurance and floral generalization are not synonymous. Self-compatible species need not self autonomously, autonomous selfing need not dominate realized mating, and floral architecture may respond to interaction environments without a change in breeding system. A macroecological island syndrome can therefore emerge from two partly independent pathways rather than a single chain in which isolation first causes selfing and selfing subsequently causes floral simplification.

A global claim also requires geographic replication. We therefore retain four predeclared island-flora strata—northern mid-latitude, northern high-latitude, tropical and southern extratropical—as independent checks on whether the same classic-positive orientation recurs in geographically distinct data. The submission does not require identical coefficients in every stratum and does not use a between-stratum difference as a primary result. The core question is recurrence of direction.

Direct evidence for the proposed ecological pressure is essential because floral traits alone cannot show that pollination is limiting. We therefore use GloPL, a global database of pollen-supplementation experiments (Bennett et al. 2018), as an independent functional evidence layer. Pollen-supplementation experiments quantify the reproductive deficit attributable to insufficient natural pollen receipt. They do not identify whether that deficit arises from visitor abundance, visitation frequency, pollen quality, mate availability or other components of pollination service, but they directly measure pollen limitation independently of the island floral-trait response.

Here we organize the evidence around four questions in biological order. **H1** asks whether isolation is associated with a recurrent island-syndrome direction spanning reproductive assurance, flower colour and floral accessibility. **H2** then asks how the floral component is generated: is it reducible to measured reproductive assurance, or do colour/accessibility shifts remain after conditioning on a strict `selfing_core` that excludes floral attraction traits? **H3** moves to an independent experimental dataset and asks whether pollen limitation itself increases with geographic isolation. **H4** finally asks whether frozen syndrome trait states are associated with lower current pollen limitation in exact-species GloPL matches.

This order distinguishes pattern, plant-side decomposition, independent ecological pressure and functional bridge. Named floral templates are retained only as secondary architecture concordances. A residual floral shift after selfing adjustment does not by itself identify direct pollinator selection or a realized pollinator guild.

---

## Materials and Methods

### Geographic universe and plant database

The plant analysis used a fixed geographic universe of 8,265 islands and a fixed denominator of 106,295 accepted angiosperm species. The final trait snapshot contained 222,688 resolved species-by-raw-axis cells among 318,885 possible cells (69.83%) across flower colour, floral structural complexity and reproductive assurance. Missing trait evidence was never converted to absence. Species-direct high/medium evidence was distinguished from validated lower-confidence evidence; the broader `all_analysis_eligible` scope was primary and the high/medium `direct_only` scope was retained as an evidence sensitivity.

The primary response describes contemporary observed island-flora composition. It does not by itself identify historical native colonization, species sorting or within-lineage evolution.

### Geographic exposure and covariates

The primary plant exposure was log distance to the nearest continental source boundary. Models included log island area and four climate principal components and used spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of geographic separation, connectivity and source supply rather than as a pure causal treatment. No post-hoc distance threshold was introduced.

For GloPL, valid georeferenced experimental sites were placed on a continuous distance axis to six seeded major continental landmasses used by the frozen global-distance extension. Distances were log1p-transformed and standardized. Mainland sites were retained at zero distance in the primary model; offshore-only shape analysis was post-hoc.

### H1: recurrent global floral/reproductive island syndrome

The v14 H1 response contains seven atomic probability outcomes coded so that a positive isolation coefficient represents the classic island-syndrome direction:

1. `self_compatibility`;
2. `selfing_mating_system`;
3. `autonomous_selfing`;
4. `plain_colour`;
5. `generalized_form`;
6. `actinomorphic_symmetry`;
7. `shallow_open_tube`.

The seven outcomes span three biological domains. Reproductive assurance comprises the first three outcomes. Colour dulling is represented by `plain_colour`, defined as white or green/brown-inconspicuous relative to yellow/orange, red/pink or blue/purple. Accessibility/generalization comprises generalized form, actinomorphy and shallow/open tube. `plain_colour` is a directional colour-composition contrast, not a direct measurement of pigment investment, animal-visible contrast or attraction intensity.

Island-level trait counts were fitted using the same beta-binomial logit framework as the frozen all-data analysis, with standardized log distance, log island area and climate PC1–PC4, and spatial-block cluster-robust covariance. Within each of four predeclared geographic replication strata, a seven-dimensional Wald test evaluates whether the isolation-response vector differs from zero.

Domain means are descriptive summaries of the fitted standardized isolation coefficients. Because the domains contain 3, 1 and 3 atomic outcomes, a single descriptive orientation is calculated as the arithmetic mean of the **three domain means**, not as the mean of all seven atomic coefficients. No new p-value is attached to this descriptive orientation.

### H2: selfing-syndrome versus pollinator-facing floral decomposition

H2 asks whether the floral component of H1 is fully reducible to reproductive assurance. `selfing_core` contains self-incompatibility/compatibility, mating-system and autonomous-selfing evidence and deliberately excludes colour, form, symmetry, tube-depth and flower-size traits.

H2a tests the reproductive-assurance route with:

`selfing_core ~ isolation + area + climate`.

H2b tests whether pollinator-facing floral responses remain after measured reproductive assurance is conditioned on. The two primary conditional responses are:

`plain_colour ~ isolation + selfing_core + area + climate`

and

`generalized_accessible ~ isolation + selfing_core + area + climate`.

`plain_colour` is fitted with the beta-binomial count model used in H1. `generalized_accessible` is an island-level continuous concordance score combining open/generalized form, actinomorphy and shallow/open tube states and is fitted with equal island weight and spatial-block cluster-robust covariance.

Pollination-syndrome concordance is then evaluated from the **reported raw trait states**, not from a weighted guild score. Five reported colour states (`white`, `red_pink`, `yellow_orange`, `blue_purple`, `green_brown_inconspicuous`) are fitted jointly and separately after `selfing_core` adjustment. Raw colour × raw floral-form and tube-depth combinations are then analysed in two ways: joint phenotype prevalence and, more stringently, the conditional frequency `P(raw architecture | raw colour, architecture resolved)`. The latter removes changes in the marginal frequency of the focal colour and tests whether colour–architecture coupling itself changes with isolation.

The raw architecture sets are compared with predeclared large-bee-, butterfly- and bird-associated trait combinations only as **concordance labels on the observed form/tube states**. The weighted `large_bee_like / butterfly_like / bird_like` scores are historical secondary summaries and do not define v14 H2 evidence.

Persistence of an isolation association after `selfing_core` adjustment means that the floral response is not reducible to the measured selfing core. Changes in raw colour–architecture coupling can be compared with known pollination-syndrome geometry, but these models remain conditional decompositions rather than causal mediation tests and do not identify realized pollinator identity.

### H3: global experimental pollen limitation

GloPL data were taken from the version-pinned public repository associated with Bennett et al. (2018), commit `abbe193770981b2c61713e75a0f01e004bdded8e`; canonical `Output/GloPL.csv` SHA-256 was `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`. The effect size is the log response ratio of pollen-supplemented to natural reproduction, so positive values indicate stronger pollen limitation.

The frozen full-global analysis contained 2,969 experimental rows, 1,248 sites and 919 publications. Duplicate measurements were summarized at publication-by-coordinate-by-measurement cells. Each publication received total analysis weight one. Models included broad-context intercepts and fixed effects for effect-size type, supplementation definition, constant addition and supplementation level, with publication-cluster robust covariance. The primary hypothesis was a positive global coefficient of standardized log distance. Supplemental-only and no-zero-constant subsets were frozen sensitivities.

### H4: post-hoc functional triangulation

The functional bridge was designed after the frozen Route A and Route B distance-by-trait moderation results had been inspected. It is therefore explicitly **post-hoc functional triangulation**, not a confirmatory mediation analysis. The design, code, parent artifacts and claim ceiling were frozen before the reproduction workflow was run.

The analysis reused exact-species matched effect rows generated by frozen reproductive-assurance and floral-architecture GloPL workflows. No synonym rescue, genus fallback, trait recoding or support-threshold relaxation was permitted. The evaluable pre-existing contrasts were self-compatibility, autonomous selfing, generalized floral form and actinomorphic symmetry; selfing mating system and shallow/open tube remained non-evaluable under parent support gates.

For each trait separately, duplicate effect rows were aggregated at the publication-by-site-by-species-by-measurement cell and each publication received total weight one. The primary model was:

`PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects`

with publication-cluster robust covariance. A negative trait-state coefficient means that the reproductive-assurance or generally accessible state is associated with lower current experimental pollen limitation after adjustment for distance and measured study structure.

Measurement sensitivities used supplemental effects only and excluded effects requiring a zero constant. Within-publication and publication-by-site comparisons retained groups containing both trait states and weighted-demeaned outcomes and covariates before fitting the remaining contrast with publication clustering. These checks reduce between-study confounding but remain observational trait-state comparisons.

The earlier Route A/B moderation tests remain separate. They asked whether trait state changed the *slope* of pollen limitation with isolation and were not supported under their frozen family rules. v14 does not reclassify those failures.

The v14 evidence hierarchy adds two validation layers without changing the inferential label of this discovery. First, PR #236 prospectively froze a post-2015 wild-plant temporal replication with two co-primary predictions: autonomous selfing should predict lower current pollen limitation, and an equal-weight accessibility/generalization score (generalized form + actinomorphy + shallow/open tube) should predict lower current pollen limitation. Outcome-blind support screening reached only 9 matched species for each co-primary family, below the frozen minimum of 30 species and 10 publications, so the wild temporal test stopped before outcome unblinding. This is a support failure, not a biological null.

Second, the same two directions were frozen for an independent PolLimCrop agricultural-domain transportability test. That analysis is support-gated and separate from the wild temporal attempt; it cannot convert the v13 discovery into a confirmatory wild-plant result or establish historical mediation. Because the H1/H2 colour response is context dependent, v14 does not add a global primary colour-to-pollen-limitation bridge.

### Supplementary pollination evidence

Named `large_bee_like`, `butterfly_like` and `bird_like` scores are multivariate floral-architecture concordances, not realized pollinator identities. Their substantial shared variance supports interpretation at the architecture level. GloBI functional-channel breadth is retained only as supplementary evidence because documented interactions depend strongly on study effort and source definition. Neither evidence family identifies global pollinator loss.

### Evidence roles and claim ceiling

The original plant and GloPL global-distance models retain their frozen inferential roles. The v13 exact-species H4 remains explicitly **post-hoc discovery**. The prospectively frozen post-2015 wild-plant validation remains separate and stopped before outcome unblinding because its support gate was not met; it therefore neither confirms nor biologically rejects H4. The independently frozen PolLimCrop analysis is a secondary cross-domain transportability test and cannot reclassify the wild discovery. No layer is interpreted as proving historical mediation from pollen limitation to trait evolution.

---

## Results

### H1: the seven-response island syndrome recurs across all four geographic strata

Adding plain colour to the six v13 functional responses did not weaken the multivariate result. In the primary all-analysis scope, the seven-response Wald test was supported in northern mid-latitude (`p=1.77×10^-6`, BH `q=2.36×10^-6`), northern high-latitude (`p=q=2.61×10^-6`), tropical (`p=1.89×10^-7`, `q=3.78×10^-7`) and southern extratropical (`p=4.01×10^-13`, `q=1.60×10^-12`) island floras. Direct-only evidence retained all four joint responses (`q=0.00461`, `9.98×10^-9`, `5.90×10^-6` and `1.09×10^-23`).

The equal-weight three-domain descriptive orientation was positive in every stratum in the primary scope: `0.0130`, `0.0900`, `0.1005` and `0.0824`. Direct-only orientations were also positive (`0.0159`, `0.0925`, `0.0854`, `0.0855`).

The component domains were not identical. Reproductive-assurance mean slopes were positive in all four primary strata (`0.0190`, `0.0730`, `0.1838`, `0.1490`), as were accessibility/generalization means (`0.0200`, `0.1783`, `0.0810`, `0.0175`). Colour dulling was approximately zero in northern mid-latitudes (`-0.00007`) but positive in northern high-latitude (`0.0188`), tropical (`0.0367`) and southern extratropical (`0.0808`) strata. Thus H1 supports recurrence of the multivariate syndrome direction, not a claim that every atomic trait changes identically in every region.

### H2: floral accessibility is not reducible to the measured selfing core

The H2a reproductive-assurance coefficient (`selfing_core`) was positive in all four geographic strata in both evidence scopes. In all-analysis data, estimates were `0.0105`, `0.0483`, `0.0621` and `0.1023`; the tropical and southern extratropical coefficients were individually supported, while the two northern estimates were positive but less precise. Direct-only estimates were `0.0252`, `0.0617`, `0.1152` and `0.1561`.

The stronger H2 result came from the selfing-adjusted structural response. After conditioning on `selfing_core`, the isolation coefficient for `generalized_accessible` remained positive in all four all-analysis strata: `0.0214`, `0.1248`, `0.0733` and `0.0536`. After BH correction across the eight primary H2b context-by-response tests, northern high-latitude (`q=0.00144`) and tropical (`q=0.00494`) accessibility effects remained supported. Direct-only estimates retained the same positive direction in all four strata, with northern high-latitude supported (`q=0.00806`) and tropical near the correction boundary (`q=0.0501`).

Selfing-adjusted colour change was more heterogeneous. `plain_colour` was approximately zero in northern mid-latitude and northern high-latitude strata, positive but borderline after multiplicity control in the tropical all-analysis scope (`β=0.0611`, `q=0.0501`), and strongly positive in the southern extratropical stratum in both all-analysis (`β=0.0803`, `q=9.24×10^-4`) and Direct-only (`β=0.0937`, `q=3.00×10^-6`) evidence.

The raw-state analyses provide the more direct pollination-syndrome concordance test. After conditioning on `selfing_core`, the five-colour vector changed with isolation in northern mid-latitude (all `p=3.14×10^-5`; Direct `p=3.38×10^-7`), tropical (all `p=2.70×10^-6`; Direct `p=3.49×10^-4`) and southern extratropical (all `p=0.00121`; Direct `p=1.23×10^-9`) floras, but not in northern high latitudes.

Colour-conditioned raw architecture then showed context-specific restructuring. In northern high latitudes, among species retaining `blue_purple`, butterfly-associated form declined with isolation in both evidence scopes (all `β=-0.20186, q=0.00144`; Direct `β=-0.19967, q=0.000681`), as did intermediate/deep large-bee-associated tube architecture (all `β=-0.25505, q=0.00144`; Direct `β=-0.24724, q=0.000266`). Tropical Direct evidence instead showed increased `yellow_orange × deep-tube` coupling (`β=+0.11703, q=0.0355`). Southern extratropical yellow/orange flowers showed mixed form and tube-depth restructuring, while northern mid-latitudes showed no FDR-supported colour-conditioned architecture despite a replicated raw `red_pink` decline.

H2 therefore supports a floral response beyond measured selfing and shows that raw display–architecture coupling can reorganize with isolation. The directions are context dependent and cannot be reduced to one named pollinator syndrome or realized visitor identity.

### H3: experimental pollen limitation increases with geographic isolation

Across 2,969 experimental rows, 1,248 sites and 919 publications, the standardized global distance coefficient was `0.07937` (SE `0.03773`; two-sided `p=0.03543`, one-sided positive `p=0.01772`). The no-zero-constant sensitivity retained a positive supported coefficient (`0.07789`, one-sided `p=0.02095`); the supplemental-only estimate was smaller but remained positive (`0.03349`).

A post-hoc shape audit showed that the signal was not explained solely by a mainland-versus-offshore level shift. The continuous gradient among sampled offshore sites was positive (`0.19983`, `p=0.01110`), as was the offshore-only fit (`0.23545`, `p=0.005865`; 260 sites, 158 publications). These shape analyses are descriptive and do not alter the inferential role of the parent global-distance result.

### H4: the functional bridge has a post-hoc discovery layer and support-gated validation layers

The reproduced post-hoc discovery yielded four evaluable frozen trait contrasts. Autonomous selfing provided the strongest functional bridge. Species coded as autonomous or delayed selfers had substantially lower current experimental pollen limitation after adjustment for distance, broad context and measurement definition (`β=-0.44672`, SE `0.08025`, two-sided `p=2.60×10^-8`). The association remained negative in supplemental-only (`β=-0.30372`) and no-zero-constant (`β=-0.45480`) sensitivities, within publications (`β=-0.30603`, two-sided `p=0.0295`) and within publication-by-site groups (`β=-0.18636`, one-sided negative `p=0.0417`).

Self-compatibility pointed in the same direction but was imprecise (`β=-0.11773`, two-sided `p=0.154`). Floral architecture provided weaker but concordant evidence. Actinomorphic species had lower pollen limitation in the primary model (`β=-0.38119`, `p=1.18×10^-5`) and both global measurement sensitivities retained negative direction. Generalized floral form was also negative in the primary model (`β=-0.18410`, `p=0.0449`) but was less stable across sensitivities and sparsely represented in within-study contrasts.

These results do not overturn the negative moderation results. The predeclared Route A/B analyses asked whether reproductive-assurance or accessibility states weakened the *increase of pollen limitation with isolation*; neither family met its promotion rule. H4 instead establishes functional association with the level of current pollen limitation, strongest for autonomous selfing.

The prospective validation layer does not yet upgrade that inference. The post-2015 wild-plant temporal replication stopped at its outcome-blind support gate: both co-primary families had only 9 matched species, below the frozen requirement of 30 species and 10 publications, so pollen-limitation outcomes were not opened for the primary validation. This is not evidence against either prediction. A separate PolLimCrop cross-domain test remains support-gated and cannot rescue the wild temporal support failure.

### Supplementary interaction evidence does not change the mechanism claim

GloBI documented interaction-channel breadth remains supplementary because its result is sensitive to source definition and recording effort. Named large-bee-like, butterfly-like and bird-like plant templates are also treated as floral-architecture concordance rather than pollinator observations. Neither layer is required for the v13 synthesis.

---

## Discussion

### The island syndrome recurs globally, but its components need not be uniform

The central H1 result is recurrence of a seven-response multivariate direction across geographically independent strata. Reproductive assurance, a directional plain-colour contrast and floral accessibility are now represented explicitly in the same formal response vector. All four strata retain supported joint responses in both evidence scopes, and the equal-weight three-domain orientation is positive everywhere.

The component pattern is nevertheless informative. Reproductive assurance and accessibility/generalization are positive as domain means in all four contexts, whereas plain colour is nearly flat in northern mid-latitudes and positive in the other three. The global result should therefore be described as a recurrent **multivariate syndrome direction**, not as a universal requirement that every island flora becomes visibly duller in the same way.

### The floral component is not one compulsory selfing syndrome

H2 addresses the key ambiguity in H1. One route is reproductive assurance: isolation is associated with increasing `selfing_core`, consistent with Baker-type establishment or reproductive-assurance processes. A strict serial selfing-syndrome model would then predict that the floral response is largely absorbed when reproductive assurance is conditioned on.

That prediction is too strong. Generalized/accessibility architecture remains positively associated with isolation after conditioning on `selfing_core` in every geographic stratum and both evidence scopes. The strongest corrected effects occur in northern high-latitude and tropical islands. This supports a second plant-side response dimension in which pollinator-facing floral architecture changes beyond the measured reproductive core.

Colour provides a useful qualification. Selfing-independent plain-colour change is strong in southern extratropical islands, weaker in tropical islands and essentially absent in the northern strata. Thus colour belongs in the H1 syndrome description, but it is not the strongest general marker of the H2 selfing-independent route.

The raw colour–architecture results set a more informative boundary. They show that some reported colour states change in frequency after `selfing_core` adjustment and that, in several regions, the association between colour and form/tube architecture itself changes with isolation. Northern high-latitude blue/purple flowers lose specialized/deep architecture, whereas tropical Direct evidence shows increased yellow/orange coupling to a deep-tube component. Southern results are internally mixed and northern mid-latitudes lack FDR-supported conditional architecture despite a robust raw-colour shift. The data therefore support context-dependent pollinator-facing reorganization without identifying one universal pollinator guild.

### Pollen limitation supplies an independent global ecological pressure

H3 is deliberately downstream of the plant-side decomposition. The GloPL result gives the syndrome an independent ecological counterpart: geographic separation is associated with stronger experimental pollen limitation across the full global sampling frame, and the positive relationship continues among sampled offshore sites.

This ordering matters. H2 does not infer pollen limitation from plant traits, and H3 does not assume that a trait shift proves pollinator decline. Pollen limitation can arise because visitation is infrequent, visitor identity is mismatched, pollen transfer is inefficient, compatible mates are scarce, or several processes covary. The current analysis therefore identifies a global **pollination-service constraint**, not a global decline in a named pollinator group.

### Functional triangulation is strongest for autonomous selfing

H4 asks a different question again: given the independent H3 pressure, are H1/H2 trait states associated with lower current pollen limitation? The strongest association is autonomous selfing. Species possessing autonomous or delayed selfing show substantially lower experimental pollen limitation after geography and measurement adjustment, and the association persists through global sensitivities and within-study comparisons.

This is functional compatibility rather than historical mediation. Current trait state and current pollen limitation can be related through evolutionary response, pre-existing differences, ecological sorting or unmeasured lineage structure. Self-compatibility alone is weaker, consistent with the fact that compatibility permits self-fertilization but does not guarantee autonomous pollen transfer.

Actinomorphy and generalized form provide weaker but directionally concordant H4 evidence. Their global associations with lower current pollen limitation are compatible with a broad accessibility interpretation, but within-study evidence is sparse or imprecise.

### Why failed slope moderation does not contradict the functional bridge

### Why failed slope moderation does not contradict the functional bridge

The frozen Route A/B moderation tests asked whether reproductive assurance or accessible architecture specifically flattened the positive relationship between isolation and pollen limitation. That interaction was not supported. H4 instead asks whether, at comparable geography and measurement structure, the protected trait state is associated with a lower level of pollen limitation.

Both statements can be true. A trait can reduce pollen limitation across the sampled gradient without changing how residual limitation scales with geographic distance. For transparency, the moderation failures remain negative results and are not “rescued” by the post-hoc analysis.

### Pollination syndromes are useful as trait geometry, not visitor labels

The study does not require flowers to be assigned to a realized “bee,” “butterfly” or “bird” pollinator. Colour, form, symmetry and tube depth overlap among functional groups. For v14, pollination-syndrome concordance is therefore evaluated from **raw reported colour and raw form/tube states**, especially the conditional frequency of architecture among plants retaining a focal colour, rather than by treating a weighted guild score as the biological outcome.

What remains useful is the underlying trait geometry. Restricted, deep or bilaterally specialized architectures can represent greater dependence on particular visitor access, whereas open, radial and generalized architectures represent broader accessibility. The northern-high-latitude loss of blue/purple coupling to specialized/deep architecture and the tropical increase in yellow/orange deep-tube coupling illustrate why the raw-state approach is preferable: different regions can reorganize along different pollinator-associated trait combinations without implying a direct census of the visitor community (Fenster et al. 2004; Rosas-Guerrero et al. 2014).

GloBI is similarly supplementary. Documented interaction breadth depends on study effort and source definition and does not carry the paper's mechanism claim.

### Limits of causal inference

The main limitation is temporal. The data triangulate three contemporaneous relationships: isolation with island-flora composition, isolation with experimental pollen limitation and trait state with current pollen limitation. They do not directly observe the historical sequence from pollination constraint through selection to trait evolution. The word “aligns” in the title is deliberate.

The primary flora response describes contemporary observed composition and should not be interpreted automatically as historical native colonization. Trait coverage is incomplete, and some key GloPL trait contrasts remain support-limited. Selfing mating system and shallow/open tube depth could not be evaluated in the functional bridge, while architecture within-study contrasts were sparse.

Finally, global distance is a composite exposure. It covaries with source separation, connectivity, colonization opportunity and potentially other environmental and historical factors. Models adjust for island area and climate in the plant analysis and broad context and measurement structure in GloPL, but residual confounding remains possible. A causal test of evolutionary sequence requires repeated lineage-level or local-system measurements of pollinator service, reproductive success and trait variation.

### A unified interpretation of the floral island syndrome

The evidence supports a four-step synthesis. First, island isolation is associated with a recurrent multivariate syndrome spanning reproductive assurance, a directional colour contrast and floral accessibility. Second, the floral component is not fully reducible to measured selfing: accessibility/generalization remains associated with isolation after conditioning on `selfing_core`, whereas colour contributes more heterogeneously. Third, independent pollen-supplementation experiments show that pollen limitation increases with isolation. Fourth, exact-species functional triangulation links autonomous selfing—and more weakly accessible floral architecture—to lower current pollen limitation.

The resulting model is not a single serial chain. Reproductive assurance and pollinator-facing floral architecture are partially separable response dimensions that can coexist and vary in relative importance across regions. The data support a common ecological pressure and recurrent plant response space, while leaving the historical causal route and pollinator identity unresolved.

### Conclusion

Across thousands of contemporary island floras, increasing geographic separation is associated with a recurrent seven-response floral/reproductive island-syndrome direction in all four geographic strata and both evidence scopes. The syndrome contains reproductive-assurance, colour and accessibility components, but the detailed colour response is not uniform among regions.

Plant-side decomposition shows that reproductive assurance increases with isolation and that generalized floral accessibility remains positively associated with isolation after conditioning on measured reproductive assurance. Raw colour composition and raw colour–architecture coupling add context-specific pollinator-syndrome concordance without reducing the result to a weighted guild score. This rejects a compulsory model in which floral change is only a downstream consequence of selfing, while the heterogeneous raw-state directions prevent attribution to one universal pollinator guild.

Independent pollen-supplementation experiments show that pollen limitation also increases with geographic isolation. The functional bridge is strongest for autonomous selfing in the post-hoc exact-species discovery, while a genuinely prospective wild-plant validation attempt remains non-evaluable because of frozen support limits and an independent crop-domain test remains separate. Together, H1–H4 establish pattern, decomposition, ecological pressure and functional compatibility without claiming historical mediation. The next decisive step is lineage-resolved or longitudinal evidence linking pollination service, reproductive success and trait change through time.

---

## Main figure legends

**Figure 1 | Reordered H1–H4 hypothesis architecture.** H1 tests a seven-response floral/reproductive island syndrome spanning reproductive assurance, plain colour and floral accessibility. H2 decomposes the floral response into a selfing-syndrome/reproductive-assurance route and a selfing-adjusted pollinator-facing route. H3 independently tests the global pollen-limitation gradient. H4 tests exact-species functional compatibility with current pollen limitation. Solid arrows denote supported associations; historical `pollen limitation -> selection -> trait evolution` remains dashed and unclaimed.

**Figure 2 | A recurrent seven-response island-syndrome direction across four geographic replication strata.** Standardized beta-binomial isolation slopes are shown for reproductive assurance, plain colour and accessibility/generalization outcomes in the primary and Direct-only evidence scopes. All four strata have supported seven-dimensional response vectors and positive equal-domain descriptive orientations. Plain colour is approximately flat in northern mid-latitudes and positive in the other three contexts, illustrating recurrence without identical component vectors.

**Figure 3 | H2 plant-side decomposition and raw pollination-syndrome concordance.** Reproductive assurance (`selfing_core`) is shown alongside selfing-adjusted isolation coefficients for `generalized_accessible` and `plain_colour`. Accessibility remains positive in all four strata after selfing adjustment. A second panel shows the five-colour joint response and FDR-supported `P(raw architecture | raw colour)` effects: northern-high-latitude blue/purple loses specialized/deep architecture, tropical Direct yellow/orange gains a deep-tube component, and southern yellow/orange shows mixed restructuring. These are raw trait concordances, not realized pollinator identities.

**Figure 4 | Independent pollen limitation and the tiered functional bridge.** The full GloPL global-distance model shows increasing pollen limitation with geographic separation across 2,969 experiments, 1,248 sites and 919 publications. The v13 exact-species discovery shows the strongest negative association for autonomous selfing, with weaker global concordance for actinomorphy and generalized form. A separate prospective post-2015 wild-plant replication is shown as support-limited and stopped before outcome unblinding (9 matched species for each co-primary family versus a 30-species gate). The independently frozen PolLimCrop transportability test is displayed separately and cannot reclassify the wild discovery.

---

## Claim ceiling for v14

The manuscript may state that geographic isolation is associated with a recurrent seven-response global island-syndrome direction replicated across four broad geographic strata; that the syndrome spans reproductive assurance, a directional plain-colour contrast and floral accessibility/generalization; that reproductive assurance increases with isolation; that accessibility/generalization remains positively associated with isolation after conditioning on measured `selfing_core`; that raw colour composition and raw colour–architecture coupling show context-specific selfing-independent reorganization; that experimental pollen limitation increases with geographic isolation globally; and that autonomous selfing is robustly associated with lower current experimental pollen limitation in the explicitly post-hoc exact-species discovery. It may additionally state that a prospectively frozen post-2015 wild validation was attempted but stopped before outcome unblinding because its frozen support gate was not met.

The manuscript must not state that H2 is causal mediation, that the selfing-adjusted floral residual proves direct pollinator selection, that a named bee/butterfly/bird syndrome or realized pollinator identity has been identified, that historical pollen limitation is proven to have selected observed traits, that pollen limitation statistically mediates the global syndrome, that pollinator abundance or visitation globally declines with isolation, that response vectors are identical across geographic strata, or that current cross-sectional data distinguish species sorting from within-lineage evolution.

---

## Literature cited

Baker, H. G. 1955. Self-compatibility and establishment after long-distance dispersal. *Evolution* 9:347–349.

Bennett, J. M., Steets, J. A., Durka, W., Burns, J. H., Vamosi, J. C., Arceo-Gómez, G., Burd, M., Burkle, L. A., Ellis, A. G., Freitas, L., Li, J. & Rodger, J. G. 2018. GloPL, a global data base on pollen limitation of plant reproduction. *Scientific Data* 5:180249.

Fenster, C. B. et al. 2004. Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics* 35:375–403.

Grossenbacher, D. L. et al. 2017. Self-compatibility is over-represented on islands. *New Phytologist* 215:469–478.

Hetherington-Rauth, M. C. & Johnson, M. T. J. 2020. Floral trait evolution of angiosperms on Pacific islands. *The American Naturalist* 196.

Pannell, J. R. & Barrett, S. C. H. 1998. Baker’s law revisited: reproductive assurance in a metapopulation. *Evolution* 52:657–668.

Rosas-Guerrero, V. et al. 2014. A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters* 17:388–400.

Sicard, A. & Lenhard, M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany* 107:1433–1443.

Zell, A. N. et al. 2025. Island colonization in flowering plants is determined by the interplay of breeding system, lifespan, floral symmetry, and arrival opportunity. *New Phytologist* 245:420–432.

---

## Reproducibility anchors

- v14 hypothesis architecture: `config/chapter1_v14_hypothesis_architecture.yml`;
- v14 H1 probability config: `config/chapter1_v14_all_data_probability.yml`;
- v14 H2 decomposition config: `config/chapter1_v14_h2_decomposition.yml`;
- v14 preflight result lock: `config/chapter1_v14_preflight_result_lock.json` (non-canonical until CI reproduction);
- v14 implementation PR: #237;
- v13 parent paper lock: `config/chapter1_v13_unified_island_syndrome_result_lock.json`;
- v13 functional bridge lock: `config/chapter1_v13_functional_bridge_result_lock.json`;
- v13 functional bridge run: `35141624253`;
- v13 functional bridge artifact: `10465048981`;
- v13 functional bridge digest: `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`;
- all-data primary probability run: `34961775336`;
- all-data primary artifact: `10394245237`;
- all-data primary digest: `sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b`;
- full-global GloPL run: `35090599662`;
- full-global GloPL artifact: `10444156159`;
- full-global GloPL digest: `sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623`.

The v13 result locks and v11/v12 evidence remain preserved unchanged as parent provenance. The v14 preflight lock is promoted to a canonical result lock only after the dedicated v14 workflow independently reproduces the frozen-input results.
