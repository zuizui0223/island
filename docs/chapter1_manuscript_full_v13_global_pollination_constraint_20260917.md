# A recurrent global floral island syndrome aligns with pollen limitation through two trait pathways

## Full working manuscript v13 — global-only submission draft — 2026-09-17

### Abstract

Island plants are often expected to become more self-reliant and less dependent on specialized pollination as geographic isolation increases, but whether this expectation forms a recurrent global syndrome remains unresolved. We combined a fixed universe of 8,265 islands and 106,295 accepted angiosperm species with an independent global database of pollen-supplementation experiments. We tested a unified hypothesis: increasing isolation is associated with a recurrent floral/reproductive island-syndrome direction, while successful pollination becomes increasingly limiting, and plants express the syndrome through partially separable reproductive-assurance and floral-accessibility pathways.

Across contemporary observed island floras, six atomic floral and reproductive responses showed a supported multivariate association with isolation in all four predeclared geographic replication strata. The descriptive classic-syndrome direction was positive in all four strata in both the primary and Direct-only evidence scopes. Reproductive assurance and floral accessibility/generalization were also positive as descriptive family-level summaries in all four strata, supporting recurrence of the same broad response across geographically distinct island floras.

Independent experimental evidence supported a global pollination constraint. Across 2,969 GloPL experiments from 1,248 sites and 919 publications, pollen limitation increased with distance from major continental landmasses (standardized slope `0.0794 ± 0.0377`; two-sided `p=0.0354`, one-sided positive `p=0.0177`). A post-hoc functional triangulation using frozen exact-species trait states then asked whether traits expected to reduce dependence on successful pollen delivery were associated with lower current pollen limitation. Autonomous selfing showed the strongest association (`β=-0.4467`, `p=2.60×10^-8`) and remained negative within publications (`β=-0.3060`, `p=0.0295`) and within publication-by-site groups (one-sided `p=0.0417`). Actinomorphy and generalized floral form also showed negative global associations, with weaker within-study support.

Together, these results support a recurrent global floral island syndrome characterized by greater reproductive assurance and greater floral accessibility, aligned with an independent increase in experimental pollen limitation. The evidence converges on a pollination-constraint hypothesis through two partially separable plant strategies, but it does not establish historical causal mediation from pollen limitation to trait evolution.

**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; pollination syndrome; GloPL

---

## Introduction

Islands provide a natural setting in which dispersal, reproduction and species interactions are simultaneously challenged. Geographic separation reduces connectivity with continental source regions, limits repeated immigration and can make mates and mutualists unreliable. These pressures motivated the expectation that island plants should increasingly possess traits that permit establishment and reproduction without dependable partners. Baker's law formalized one part of this expectation: uniparental reproduction can be advantageous when colonists experience mate limitation or uncertain pollination (Baker 1955; Pannell & Barrett 1998). Global comparisons subsequently showed that self-compatible plants are over-represented on islands and that breeding system interacts with other traits during island colonization (Grossenbacher et al. 2017; Zell et al. 2025).

The same logic extends beyond breeding system. Floral symmetry, tube depth, form, colour and display can increase compatibility with particular functional groups, but specialized architectures can also create dependence on the availability and effectiveness of suitable visitors. Pollination-syndrome theory therefore predicts correlated combinations of floral traits rather than one diagnostic character (Fenster et al. 2004; Rosas-Guerrero et al. 2014). If successful pollination becomes unreliable on increasingly isolated islands, at least two responses are plausible. Plants may gain from **reproductive assurance**, through self-compatibility, autonomous selfing or predominantly selfing mating systems. Independently, the benefit of maintaining restricted floral access may decline, favouring more open, generalized or radially accessible flowers.

These routes are related but need not form one obligatory sequence. The selfing syndrome includes recurrent reductions in floral investment associated with transitions toward selfing (Sicard & Lenhard 2011), but reproductive assurance and floral generalization are not synonymous. Self-compatible species need not self autonomously, autonomous selfing need not dominate realized mating, and floral architecture may respond to interaction environments without a change in breeding system. A macroecological island syndrome can therefore emerge from two partly independent pathways rather than a single chain in which isolation first causes selfing and selfing subsequently causes floral simplification.

A global claim also requires geographic replication. We therefore retain four predeclared island-flora strata—northern mid-latitude, northern high-latitude, tropical and southern extratropical—as independent checks on whether the same classic-positive orientation recurs in geographically distinct data. The submission does not require identical coefficients in every stratum and does not use a between-stratum difference as a primary result. The core question is recurrence of direction.

Direct evidence for the proposed ecological pressure is essential because floral traits alone cannot show that pollination is limiting. We therefore use GloPL, a global database of pollen-supplementation experiments (Bennett et al. 2018), as an independent functional evidence layer. Pollen-supplementation experiments quantify the reproductive deficit attributable to insufficient natural pollen receipt. They do not identify whether that deficit arises from visitor abundance, visitation frequency, pollen quality, mate availability or other components of pollination service, but they directly measure pollen limitation independently of the island floral-trait response.

Here we organize the frozen evidence around four hypotheses. **H1** asks whether a recurrent classic island-syndrome direction occurs across four geographic replication strata. **H2** asks whether experimental pollen limitation increases with geographic isolation globally. **H3** separates the plant response into reproductive-assurance and accessibility/generalization pathways and asks whether both recur globally. **H4** is an explicitly post-hoc functional triangulation: using trait states frozen before the new analysis, are reproductive-assurance or generally accessible states associated with lower current experimental pollen limitation after accounting for geography and measurement structure?

The central test is whether three independent evidence layers converge: isolation with a recurrent plant syndrome, isolation with experimental pollen limitation, and functional trait states with current pollen limitation. Named floral templates and GloBI interaction records remain supplementary compatibility evidence and are not required to identify a named pollinator mechanism.

---

## Materials and Methods

### Geographic universe and plant database

The plant analysis used a fixed geographic universe of 8,265 islands and a fixed denominator of 106,295 accepted angiosperm species. The final trait snapshot contained 222,688 resolved species-by-raw-axis cells among 318,885 possible cells (69.83%) across flower colour, floral structural complexity and reproductive assurance. Missing trait evidence was never converted to absence. Species-direct high/medium evidence was distinguished from validated lower-confidence evidence; the broader `all_analysis_eligible` scope was primary and the high/medium `direct_only` scope was retained as an evidence sensitivity.

The primary response describes contemporary observed island-flora composition. It does not by itself identify historical native colonization, species sorting or within-lineage evolution.

### Geographic exposure and covariates

The primary plant exposure was log distance to the nearest continental source boundary. Models included log island area and four climate principal components and used spatial-block cluster-robust covariance. Distance is interpreted as a composite gradient of geographic separation, connectivity and source supply rather than as a pure causal treatment. No post-hoc distance threshold was introduced.

For GloPL, valid georeferenced experimental sites were placed on a continuous distance axis to six seeded major continental landmasses used by the frozen global-distance extension. Distances were log1p-transformed and standardized. Mainland sites were retained at zero distance in the primary model; offshore-only shape analysis was post-hoc.

### H1: recurrent global island-syndrome direction

The primary all-observed response consisted of six atomic probability outcomes:

1. `generalized_form`;
2. `actinomorphic_symmetry`;
3. `shallow_open_tube`;
4. `self_compatibility`;
5. `selfing_mating_system`;
6. `autonomous_selfing`.

Each outcome was coded so that a positive isolation coefficient represents the classic island-syndrome direction: greater floral accessibility/generalization or greater reproductive assurance. Island-level trait counts were fitted with the frozen beta-binomial logit framework using distance, island area and climate PC1–PC4, with spatial-block cluster-robust covariance. A matched grouped-binomial analysis was retained as model-form sensitivity.

Within each of four predeclared geographic replication strata, the already locked six-dimensional omnibus test establishes whether the response vector differs from zero. Orientation is summarized descriptively as the arithmetic mean of the six identically coded distance slopes. No new p-value is attached to that arithmetic mean. A positive mean in all four strata, together with supported joint vectors, is used as evidence that the classic-positive direction recurs globally. Between-stratum contrasts are not used in the current submission.

### H2: global experimental pollen limitation

GloPL data were taken from the version-pinned public repository associated with Bennett et al. (2018), commit `abbe193770981b2c61713e75a0f01e004bdded8e`; canonical `Output/GloPL.csv` SHA-256 was `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`. The effect size is the log response ratio of pollen-supplemented to natural reproduction, so positive values indicate stronger pollen limitation.

The frozen full-global analysis contained 2,969 experimental rows, 1,248 sites and 919 publications. Duplicate measurements were summarized at publication-by-coordinate-by-measurement cells. Each publication received total analysis weight one. Models included broad-context intercepts and fixed effects for effect-size type, supplementation definition, constant addition and supplementation level, with publication-cluster robust covariance. The primary hypothesis was a positive global coefficient of standardized log distance. Supplemental-only and no-zero-constant subsets were frozen sensitivities.

### H3: reproductive-assurance and floral-accessibility pathways

The six atomic responses were grouped into two biological families without changing individual fits. The reproductive-assurance family comprised self-compatibility, selfing mating system and autonomous selfing. The floral-accessibility family comprised generalized form, actinomorphic symmetry and shallow/open tube. Family means are descriptive summaries of already fitted identically oriented atomic coefficients, not newly fitted composite traits.

A separate frozen plant-side analysis evaluated whether the floral-accessibility response was reducible to measured reproductive assurance. `selfing_core` combined strict reproductive evidence while excluding floral-attraction traits. An attraction/accessibility shift contrasted generally accessible architecture against large-bee-like restricted architecture and was fitted with distance, selfing core and covariates. These models are interpreted as conditional decomposition rather than causal mediation.

### H4: post-hoc functional triangulation

The functional bridge was designed after the frozen Route A and Route B distance-by-trait moderation results had been inspected. It is therefore explicitly **post-hoc functional triangulation**, not a confirmatory mediation analysis. The design, code, parent artifacts and claim ceiling were frozen before the reproduction workflow was run.

The analysis reused exact-species matched effect rows generated by frozen reproductive-assurance and floral-architecture GloPL workflows. No synonym rescue, genus fallback, trait recoding or support-threshold relaxation was permitted. The evaluable pre-existing contrasts were self-compatibility, autonomous selfing, generalized floral form and actinomorphic symmetry; selfing mating system and shallow/open tube remained non-evaluable under parent support gates.

For each trait separately, duplicate effect rows were aggregated at the publication-by-site-by-species-by-measurement cell and each publication received total weight one. The primary model was:

`PL ~ trait_state + z_distance + context_intercepts + measurement_fixed_effects`

with publication-cluster robust covariance. A negative trait-state coefficient means that the reproductive-assurance or generally accessible state is associated with lower current experimental pollen limitation after adjustment for distance and measured study structure.

Measurement sensitivities used supplemental effects only and excluded effects requiring a zero constant. Within-publication and publication-by-site comparisons retained groups containing both trait states and weighted-demeaned outcomes and covariates before fitting the remaining contrast with publication clustering. These checks reduce between-study confounding but remain observational trait-state comparisons.

The earlier Route A/B moderation tests remain separate. They asked whether trait state changed the *slope* of pollen limitation with isolation and were not supported under their frozen family rules. v13 does not reclassify those failures.

### Prospective post-2015 temporal validation audit

After the post-hoc H4 discovery, but before reading any new post-2015 validation outcomes, we froze a separate temporal-replication protocol. The validation sampling frame was restricted to publications dated 2016-01-01 through 2026-09-18 and was assembled from bibliographic metadata without using abstracts, Results text, reproductive-output values or pollen-limitation effect sizes. Frozen species-direct trait states were matched exactly; synonym rescue, genus fallback and post-support trait recoding were prohibited.

Two co-primary predictions were fixed in advance. H4a predicted lower current pollen limitation in species with autonomous selfing. H4b predicted lower current pollen limitation with a three-component accessibility/generalization score combining generalized form, actinomorphic symmetry and shallow/open tube depth; shallow/open tube was retained despite its weaker v13 recurrence and could not be dropped after support inspection. Each co-primary test had a one-sided alpha of 0.025.

Outcome-blind support gates were also fixed before validation outcomes. H4a required at least 30 matched species, 10 species per autonomous-selfing state, 10 publications and 5 publications per state. H4b required at least 30 matched species, 10 publications, accessibility-score SD >=0.10, and at least 8 species in each low- and high-score tail. If a gate failed, the corresponding validation outcomes were not authorized to be opened or analysed. This prospective audit is separate from the v13 H1-H4 result lock and cannot relabel the original H4 analysis as confirmatory.

### Supplementary pollination evidence

Named `large_bee_like`, `butterfly_like` and `bird_like` scores are multivariate floral-architecture concordances, not realized pollinator identities. Their substantial shared variance supports interpretation at the architecture level. GloBI functional-channel breadth is retained only as supplementary evidence because documented interactions depend strongly on study effort and source definition. Neither evidence family identifies global pollinator loss.

### Evidence roles and claim ceiling

The original plant and GloPL global-distance models retain their frozen inferential roles. The v13 recurrence summaries re-express locked coefficients, and H4 is explicitly post-hoc. No analysis is interpreted as proving historical mediation from pollen limitation to trait evolution. The strongest allowed conclusion is triangulation among a recurrent geographic plant pattern, an independent experimental pollen-limitation gradient and functional associations of frozen trait states with current pollen limitation.

---

## Results

### H1: a recurrent island-syndrome direction occurs in all four geographic replication strata

The six-atomic isolation-response vector was supported within every geographic replication stratum in both evidence scopes. In the primary all-analysis scope, omnibus p-values were `1.06×10^-6`, `1.10×10^-6`, `1.56×10^-5` and `2.09×10^-13` across the four strata. Direct-only evidence retained all four responses (`p=0.00231`, `1.79×10^-9`, `2.83×10^-5` and `1.14×10^-20`).

Under common classic-positive coding, the descriptive mean of the six atomic slopes was positive in all four strata: `0.0195`, `0.1256`, `0.1324` and `0.0832` in the primary scope. Direct-only means were also positive (`0.0228`, `0.1299`, `0.1113`, `0.0779`). In the primary scope, five, six, six and five of the six atomic slopes were positive; Direct-only evidence showed six, six, five and five positive slopes.

The submission therefore treats the four geographically distinct strata as independent replication of a recurrent classic-positive orientation. No between-stratum contrast is used to define the global syndrome.

### H2: experimental pollen limitation increases with geographic isolation

Across 2,969 experimental rows, 1,248 sites and 919 publications, the standardized global distance coefficient was `0.07937` (SE `0.03773`; two-sided `p=0.03543`, one-sided positive `p=0.01772`). The no-zero-constant sensitivity retained a positive supported coefficient (`0.07789`, one-sided `p=0.02095`); the supplemental-only estimate was smaller but remained positive (`0.03349`).

A post-hoc shape audit showed that the signal was not explained solely by a mainland-versus-offshore level shift. The continuous gradient among sampled offshore sites was positive (`0.19983`, `p=0.01110`), as was the offshore-only fit (`0.23545`, `p=0.005865`; 260 sites, 158 publications). These shape analyses are descriptive and do not alter the inferential role of the parent global-distance result.

### H3: reproductive assurance and floral accessibility form two recurrent response pathways

The reproductive-assurance family mean was positive in all four geographic replication strata in both evidence scopes. Primary means were `0.0190`, `0.0730`, `0.1838` and `0.1490`; Direct-only values were `0.0298`, `0.1003`, `0.1309` and `0.0967`.

The accessibility/generalization family mean was also positive in all four strata: `0.0200`, `0.1783`, `0.0810` and `0.0175` in the primary scope and `0.0158`, `0.1596`, `0.0916` and `0.0590` in Direct-only evidence. Shallow/open tube depth was not universally positive, so the global statement concerns multivariate accessibility/generalization rather than universal simplification of every floral trait.

The two routes were not reducible to one obligatory sequence. Frozen conditional analyses retained a distance-associated attraction/accessibility component after accounting for `selfing_core`, and the distance-by-selfing interaction did not require accessibility shift to occur only where reproductive assurance was high. The plant data therefore support partially separable reproductive-assurance and floral-accessibility responses.

### H4: current pollen limitation is lower in key island-syndrome trait states

The reproduced post-hoc functional analysis yielded four evaluable frozen trait contrasts. Autonomous selfing provided the strongest functional bridge. Species coded as autonomous or delayed selfers had substantially lower current experimental pollen limitation after adjustment for distance, broad context and measurement definition (`β=-0.44672`, SE `0.08025`, two-sided `p=2.60×10^-8`). The association remained negative in supplemental-only (`β=-0.30372`) and no-zero-constant (`β=-0.45480`) sensitivities, within publications (`β=-0.30603`, two-sided `p=0.0295`) and within publication-by-site groups (`β=-0.18636`, one-sided negative `p=0.0417`).

Self-compatibility pointed in the same direction but was imprecise (`β=-0.11773`, two-sided `p=0.154`). Floral architecture provided weaker but concordant evidence. Actinomorphic species had lower pollen limitation in the primary model (`β=-0.38119`, `p=1.18×10^-5`) and both global measurement sensitivities retained negative direction. Generalized floral form was also negative in the primary model (`β=-0.18410`, `p=0.0449`) but was less stable across sensitivities and sparsely represented in within-study contrasts.

These results do not overturn the negative moderation results. The predeclared Route A/B analyses asked whether reproductive-assurance or accessibility states weakened the *increase of pollen limitation with isolation*; neither family met its promotion rule. H4 instead establishes functional association with the level of current pollen limitation, strongest for autonomous selfing.

### Prospective temporal validation stopped at the outcome-blind support gate

The prospectively specified post-2015 wild-plant replication did not reach its minimum support requirements. Even under the most permissive predeclared core-design diagnostic ceiling, H4a contained 9 matched species across 6 publications (autonomous-selfing state 0: 7 species; state 1: 2), and H4b contained 9 matched species across 4 publications. Both were below the frozen minima of 30 matched species and 10 publications, with additional state- or score-balance requirements also unmet.

Accordingly, neither co-primary validation test was evaluable and the post-2015 primary pollen-limitation outcomes remained unopened. Thresholds were not relaxed and no substitute trait definition was introduced. This is a prospective **support insufficiency**, not evidence for absence of either biological association, and it does not change the post-hoc inferential role of H4.

### Supplementary interaction evidence does not change the mechanism claim

GloBI documented interaction-channel breadth remains supplementary because its result is sensitive to source definition and recording effort. Named large-bee-like, butterfly-like and bird-like plant templates are also treated as floral-architecture concordance rather than pollinator observations. Neither layer is required for the v13 synthesis.

---

## Discussion

### The floral island syndrome recurs globally

The central result is recurrence across geographically independent strata. Every stratum has a supported six-atomic isolation response, and the average direction of identically coded traits is positive in every stratum and both evidence scopes. This is stronger biologically than requiring one identical global coefficient vector: a global syndrome can recur when geographically distinct floras repeatedly move toward greater reproductive assurance and greater floral accessibility even though individual component magnitudes differ.

The result also avoids treating a syndrome as a rigid checklist. Shallow floral access is not uniformly positive, and individual reproductive components differ in strength and precision. The global inference rests on recurrent multivariate orientation and on both higher-level plant-response families pointing in the same classic direction across all four replication strata.

### Pollen limitation supplies an independent global ecological pressure

The GloPL result gives the plant pattern an independent ecological counterpart. Geographic separation is associated with stronger experimental pollen limitation across the full global sampling frame, and the positive relationship continues among sampled offshore sites. Pollen limitation is therefore measured independently by reproductive response to pollen supplementation rather than inferred from floral morphology.

The mechanism should be described carefully. Pollen limitation can arise because visitation is infrequent, visitor identity is mismatched, pollen transfer is inefficient, compatible mates are scarce, or several processes covary. The current analysis does not identify a global decline in pollinator abundance. “Pollination-service constraint” is therefore more defensible than a named global pollinator-loss mechanism.

### Reproductive assurance is the strongest functional bridge

The post-hoc functional triangulation is most compelling for autonomous selfing. Autonomous selfing is precisely the form of reproductive assurance expected to reduce the reproductive consequences of failed pollen delivery, and species possessing autonomous or delayed selfing show substantially lower experimental pollen limitation. The association persists after geography and measurement adjustment, through both measurement sensitivities and within publications and publication-by-site groups containing both trait states.

This does not prove that historical pollen limitation selected autonomous selfing on islands. Current trait state and current pollen limitation can be related through evolutionary response, pre-existing trait differences, ecological sorting or unmeasured lineage covariates. The appropriate inference is functional: autonomous selfing is associated with reduced current pollen limitation in an independent experimental database, making it a biologically credible response to unreliable pollen delivery.

Self-compatibility alone is a weaker bridge, which is biologically reasonable. Compatibility permits self-fertilization but does not guarantee autonomous pollen transfer, so a self-compatible species can remain strongly dependent on visitors.

### Floral accessibility is a second, more cautiously supported route

Open or generalized forms and actinomorphy can make flowers accessible to a wider set of visitors or reduce dependence on tightly matched interaction morphology. Globally, actinomorphy and generalized form are associated with lower current pollen limitation after adjustment for distance and measurement structure. The architecture evidence is weaker than the autonomous-selfing bridge because within-study contrasts are sparse or imprecise.

The plant-side conditional analysis nevertheless matters: accessibility/generalization persists as a distance-associated component after accounting for strict reproductive assurance. Together, plant and GloPL evidence support a **dual-pathway model** more strongly than a compulsory one-dimensional selfing syndrome.

### Why failed slope moderation does not contradict the functional bridge

The frozen Route A/B moderation tests asked whether reproductive assurance or accessible architecture specifically flattened the positive relationship between isolation and pollen limitation. That interaction was not supported. H4 instead asks whether, at comparable geography and measurement structure, the protected trait state is associated with a lower level of pollen limitation.

Both statements can be true. A trait can reduce pollen limitation across the sampled gradient without changing how residual limitation scales with geographic distance. For transparency, the moderation failures remain negative results and are not “rescued” by the post-hoc analysis.

### Pollination syndromes are useful as trait geometry, not visitor labels

The study does not require flowers to be assigned to a realized “bee,” “butterfly” or “bird” pollinator. Colour, form, symmetry and tube depth overlap among functional groups, and named template scores share substantial floral-architecture variance.

What remains useful is the underlying trait geometry. Restricted, deep or bilaterally specialized architectures can represent greater dependence on particular visitor access, whereas open, radial and generalized architectures represent broader accessibility. This is the level at which pollination-syndrome theory supports v13: floral architecture is biologically meaningful, but phenotype is not a direct census of the visitor community (Fenster et al. 2004; Rosas-Guerrero et al. 2014).

GloBI is similarly supplementary. Documented interaction breadth depends on study effort and source definition and does not carry the paper's mechanism claim.

### Limits of causal inference

The main limitation is temporal. The data triangulate three contemporaneous relationships: isolation with island-flora composition, isolation with experimental pollen limitation and trait state with current pollen limitation. They do not directly observe the historical sequence from pollination constraint through selection to trait evolution. The word “aligns” in the title is deliberate.

We therefore attempted a genuinely prospective temporal replication after H4 was discovered: the predictions, trait definitions, multiplicity rule and support thresholds were fixed before post-2015 validation outcomes. That audit stopped before outcome unblinding because the independent wild-plant cohort did not meet the frozen support gates. This strengthens the transparency of the evidence boundary but does not provide a confirmatory replication, and the failure to reach support is not a biological null.

The primary flora response describes contemporary observed composition and should not be interpreted automatically as historical native colonization. Trait coverage is incomplete, and some key GloPL trait contrasts remain support-limited. Selfing mating system and shallow/open tube depth could not be evaluated in the functional bridge, while architecture within-study contrasts were sparse.

Finally, global distance is a composite exposure. It covaries with source separation, connectivity, colonization opportunity and potentially other environmental and historical factors. Models adjust for island area and climate in the plant analysis and broad context and measurement structure in GloPL, but residual confounding remains possible. A causal test of evolutionary sequence requires repeated lineage-level or local-system measurements of pollinator service, reproductive success and trait variation.

### A unified interpretation of the floral island syndrome

The evidence supports a simple global synthesis. Island isolation is associated with a recurrent direction toward greater reproductive assurance and greater floral accessibility. Independent experiments show that pollen limitation increases with isolation. The functional trait analysis shows that autonomous selfing is associated with substantially lower current pollen limitation, with weaker but concordant evidence for generally accessible floral architecture.

The syndrome is best viewed as an **emergent solution space** to unreliable pollination. Reproductive assurance reduces dependence on pollen delivery. Generalized floral accessibility may reduce dependence on narrowly matched visitor access. The two strategies can contribute independently to the same recurrent global direction.

### Conclusion

Across thousands of contemporary island floras, increasing geographic separation is associated with a recurrent floral/reproductive island-syndrome direction, replicated across all four geographic strata and both evidence scopes. Independent pollen-supplementation experiments show that pollen limitation also increases with geographic isolation. The plant response contains at least two partially separable components—reproductive assurance and floral accessibility/generalization—and post-hoc functional triangulation provides the strongest bridge for autonomous selfing, which is associated with markedly lower current experimental pollen limitation even within studies.

These results support a unified pollination-constraint hypothesis while preserving a strict causal boundary. They show convergence among geography, experimental reproductive limitation and functional trait states, but they do not identify historical pollen limitation as the cause of observed trait evolution. The next decisive step is longitudinal or lineage-resolved evidence linking pollination service, reproductive success and trait change through time.

---

## Main figure legends

**Figure 1 | Unified global hypothesis and evidence hierarchy for a recurrent floral island syndrome.** Geographic isolation is associated independently with a global experimental pollen-limitation gradient and with a recurrent plant response. The plant response separates into reproductive-assurance and accessibility/generalization pathways. Solid arrows denote supported associations; the trait-state link to current experimental pollen limitation is post-hoc functional triangulation; the historical `pollen limitation -> selection -> trait evolution` arrow remains dashed and unclaimed.

**Figure 2 | A recurrent global island-syndrome direction across four geographic replication strata.** Six-atomic standardized isolation slopes are shown for four predeclared geographic strata in the primary and Direct-only evidence scopes. All four strata have supported multivariate response vectors and positive descriptive mean classic directions. Reproductive-assurance and accessibility/generalization family summaries show the two recurrent plant pathways. The strata are used as replication of recurrence rather than a primary between-stratum comparison.

**Figure 3 | Experimental pollen limitation and the post-hoc functional bridge.** The full GloPL global-distance model shows increasing pollen limitation with geographic separation across 2,969 experiments, 1,248 sites and 919 publications. A visually separate post-hoc offshore diagnostic shows that the association continues across sampled offshore sites. Frozen exact-species trait states are then compared with current pollen limitation: autonomous selfing shows the strongest negative association and remains negative within publications and publication-by-site groups; actinomorphy and generalized form provide additional global concordance with weaker within-study support. Frozen distance-by-trait moderation failures remain displayed separately as negative results.

---

## Claim ceiling for v13

The manuscript may state that geographic isolation is associated with a recurrent global classic island-syndrome direction replicated across four broad geographic strata; that experimental pollen limitation increases with geographic isolation globally; that reproductive assurance and floral accessibility/generalization are partially separable plant-response pathways; and that, in explicitly post-hoc exact-species triangulation, autonomous selfing is robustly associated with lower current experimental pollen limitation while actinomorphy and generalized floral form provide additional global functional concordance of differing robustness.

The manuscript must not state that historical pollen limitation is proven to have selected observed traits, that pollen limitation statistically mediates the global syndrome, that pollinator abundance or visitation globally declines with isolation, that a named pollinator group caused the pattern, that GloBI identifies the causal mechanism, that response vectors are identical across geographic strata, that between-stratum differences are a primary submission claim, or that current cross-sectional data distinguish species sorting from within-lineage evolution.

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

- v13 paper lock: `config/chapter1_v13_unified_island_syndrome_result_lock.json`;
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

Historical v11/v12 evidence remains preserved unchanged outside the current submission surface.
