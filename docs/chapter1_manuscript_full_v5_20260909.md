# The island floral syndrome is not one syndrome: source-distance filtering and biogeographic decoupling of reproductive assurance and floral architecture

## Full working manuscript v5 — 2026-09-09

### Abstract

Island floras are often expected to converge toward a coherent floral “island syndrome”, combining increased reproductive assurance with reduced floral specialization or attraction. Yet geographic isolation does not act directly on flowers. An oceanic barrier first changes which plant lineages from a regional source pool can disperse, arrive, establish and persist, and it can also change which pollination channels remain available after colonization. These two filters need not act at the same rate or in the same direction. We therefore asked whether increasing source separation produces one universal floral/reproductive response, or instead reorganizes partially separable reproductive and floral components in a biogeographically contingent manner.

We analysed a fixed universe of 8,265 islands and 106,295 accepted angiosperm species under a progressive, preregulated analysis contract. The final trait snapshot resolved 222,688 of 318,885 possible species-by-axis cells (69.83%) across flower colour, floral structural complexity and reproductive assurance. The geographic exposure was log-transformed mainland distance, interpreted as a composite gradient of geographic separation, connectivity and source accessibility rather than as a mechanistically pure causal treatment. The analysis proceeded sequentially: H1 tested a universal syndrome; H2 tested biogeographic branching of a pollinator-name-free two-axis plant response; H3 decomposed supported responses by source and taxonomic depth; H4 tested continuous area moderation; and H5 was withheld from Chapter 1 unless independent pollination-channel evidence was available.

A universal floral/reproductive response was not recovered. In the Palearctic, increasing separation was associated with greater accessibility/generalization and greater reproductive assurance. In all-analysis all-native assemblages, the corresponding slopes were 0.0795 and 0.0534 per unit of log-distance (both FDR q=0.00463); direct-only estimates were 0.0616 (q=0.0389) and 0.0966 (q=8.48×10^-14). The pattern persisted in native non-endemics. Tropical assemblages followed a different trajectory: in direct-only all-native data, accessibility/generalization declined with distance (-0.1005, q=0.00275) while reproductive assurance increased (0.1360, q=0.00428). Thus reproductive assurance and floral architecture did not form one obligatory serial syndrome.

Source-matched taxonomic decomposition showed that the broad Palearctic response survived family adjustment but disappeared after genus adjustment in all four evidence-scope-by-floristic-stratum combinations (4/4 → 4/4 → 0/4). Moreover, the Palearctic attraction/access shift remained positive after conditioning on selfing core across four source definitions (distance estimates 0.091–0.101 in all natives; q=0.0079–0.034), whereas tropical attraction/access shifts remained negative after the same decomposition. Secondary large-bee-like, butterfly-like and bird-like floral templates also diverged regionally, but a single source-trained factor explained 86.94% and 86.44% of template variance in all-analysis and direct-only evidence, respectively, showing that these labels mainly capture overlapping plant architecture rather than realized pollinator identity.

The results recast the floral island syndrome as an emergent assemblage-level outcome with at least two partially separable components: reproductive assurance and pollination-associated floral architecture. The strongest plant-side distance response is compatible with genus-level source filtering, while measured-climate-independent regional causation, an area/capacity mechanism and pollinator-channel causation remain beyond the Chapter 1 identification ceiling. Large additions of trait evidence no longer changed this claim structure, providing an inference-based stopping rule for the computational campaign. We propose that the next mechanistic layer should test the second half of a double geographic filter: whether pollination channels differ in their ability to cross and persist across the same oceanic barriers that filter plant source pools.

**Keywords:** island biogeography; floral traits; reproductive assurance; selfing syndrome; pollination syndrome; source pool; dispersal filtering; biogeographic contingency; mutualism; island syndrome

---

## Introduction

Geographic isolation is one of the defining processes of island biogeography. Classical island theory treats separation from a regional source pool as a constraint on colonization opportunity, while empirical island biology has repeatedly identified characteristic shifts in life history, growth form, dispersal and reproduction. In plants, these observations have motivated the idea of an “island syndrome”: a recurrent suite of traits favoured by long-distance colonization, small population size, novel environments and altered biotic interactions.

Floral versions of the island syndrome usually combine at least two biological arguments. The first is reproductive assurance. Long-distance colonists may have difficulty finding mates or reliable pollen vectors, creating an advantage for self-compatibility, autonomous selfing or other forms of uniparental reproduction. This argument is closely related to Baker’s law and its later formulation as a colonization filter: reproductive traits may be over-represented on islands because they increase establishment probability, not necessarily because the same reproductive transition repeatedly evolves after arrival. The second argument concerns pollination-associated floral architecture. If island pollinator assemblages are species-poor, functionally altered or dominated by different visitor groups, the returns to conspicuous displays, deep floral tubes, restricted access or precise mechanical matching may change. These processes can generate superficially similar phenotypes, but they are not equivalent mechanisms.

The distinction matters because reproductive assurance and floral architecture need not change together. The classical selfing syndrome often includes reduced floral investment, but self-compatibility is not autonomous selfing, autonomous selfing is not necessarily the dominant realized mating mode, and a plant can retain specialized floral architecture while gaining reproductive assurance. Conversely, floral access or attraction traits can change without a detectable shift in the measured selfing core. Treating all of these outcomes as one syndrome risks turning a broad pattern into a presumed causal sequence.

A second problem is that island floral patterns are assemblage patterns. A distant island can differ from a near island because different source-available lineages arrived, established or persisted, even if no lineage changed phenotype after colonization. Baker’s law itself can be understood this way: pre-existing reproductive variation is sorted during colonization. The same logic applies to floral form, symmetry and access. If lineages differ in both dispersal/establishment ability and floral strategy, geographic isolation can change community-level trait composition through differential lineage representation.

This source-pool perspective changes the interpretation of distance. We do not treat kilometres from a continent as a direct causal force on floral phenotype. Instead, distance is a composite axis of separation, connectivity and source accessibility. Conceptually, as distance approaches zero, an island approaches a high-accessibility boundary where geographic filtering relative to the regional source system should be weaker. As distance increases, the opportunity for source-to-island filtering increases. This limiting argument does not require a literal zero-distance island to equal mainland vegetation, nor does it imply that the nearest continent is the exact historical source for every island. Rather, it provides a continuous alternative to a binary island-versus-mainland comparison.

The same geographic barrier may also act on the other side of plant–pollinator interactions. A regional pollination channel must itself cross oceanic barriers, establish and persist before it can provide effective service on an island. Flight ability, long-distance dispersal, stepping-stone use, habitat requirements and establishment probability can therefore make a given ocean distance a different biological barrier for different pollinator functional groups. Structural absence of a pollinator from a regional source system is also fundamentally different from the loss of a channel that was available in the source region. This motivates a double geographic-filter framework: geographic separation can simultaneously filter plant source-pool assembly and pollination-channel continuity. In the present paper, however, only the plant side is directly tested; the pollinator-side filter is retained as a mechanistic hypothesis rather than inferred from floral phenotype.

A third challenge is analytical. Global island comparisons combine heterogeneous floristic histories, climates, island sizes, observation effort and taxonomic composition. A claim that one region is “significant” and another is not does not demonstrate that regional response vectors differ. Likewise, adjusting for lineage as a nuisance variable can remove the very assembly process that island biogeography seeks to explain. Missing trait evidence is also geographically and taxonomically structured, so increasing database coverage can change which biological questions are estimable.

We address these problems using a progressive analysis contract established in PR #142. Instead of allowing the model to change as trait coverage improved, we fixed the island universe, trait ontology, hypothesis order, support thresholds, multiple-testing families, source/lineage safeguards and causal claim ceiling. New trait snapshots were then rerun from the beginning. This creates a sequential hypothesis ladder in which each result determines the next biological question while preventing the question itself from drifting with the data.

We test four plant-data hypotheses and define a fifth mechanistic extension. **H1** treats a universal floral island syndrome as a rival hypothesis: does increasing source separation generate one coherent floral/reproductive response vector across contexts? **H2** asks whether isolation-associated response vectors instead branch among biogeographic contexts. **H3** asks at what source/taxonomic depth a supported response remains, distinguishing broad assemblage patterns from family- and genus-level representation. **H4** asks whether continuous island area modifies isolation-associated filtering. **H5**, which is deliberately not claimed from the present plant database, asks whether an independently measured pollination-channel deficit explains residual filtering after source, lineage, area, climate, observation and spatial safeguards.

The primary H2 response is deliberately pollinator-name-free. We separate (i) **accessibility/generalization**, representing increasingly open, shallow or radially accessible floral architecture, from (ii) **reproductive assurance**, representing strict reproductive-assurance evidence. Only after the primary response is established do we evaluate fixed large-bee-like, butterfly-like and bird-like floral templates. These are sampled pollination-associated architecture concordances, not pollinator classifiers.

Finally, we use the progressive trait campaign itself to decide when additional database acquisition ceases to be the highest-value scientific activity. The target is not complete global trait coverage. The relevant criterion is inference saturation: whether large additions of usable evidence continue to change the major hypothesis decisions, and whether remaining gaps can be recovered without disproportionate source-specific review or weaker evidential standards.

Our central prediction is therefore not that distant islands should be uniformly selfing, simple or inconspicuous. It is that increasing source separation should reorganize floral and reproductive assemblages, but the signs, coupling and taxonomic depth of those responses should depend on biogeographic context. If the floral island syndrome is an emergent outcome rather than one deterministic rule, reproductive assurance and floral architecture should sometimes decouple, and broad distance gradients should often resolve into lineage assembly rather than repeated within-lineage change.

---

## Materials and Methods

### Geographic and taxonomic universe

The fixed geographic universe contained 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. The denominator was not allowed to shrink between trait waves. Island occurrences, floristic-status information, island area, climate covariates, geographic exposure and source-pool infrastructure were held constant while trait evidence was allowed to improve.

The primary geographic exposure was `log1p_distance_to_continent_km`. We interpret this as a composite of geographic separation, connectivity and source accessibility rather than as a mechanistically pure isolation treatment. The conceptual near-source limit (`d → 0`) represents high source accessibility and relatively weak oceanic filtering, not a literal mainland observation or a claim that a zero-distance island should have mainland composition.

Analyses were conducted at two context layers. Broad **analysis regimes** included northern mid-latitude, northern high-latitude, tropical and southern extratropical strata. Formal **biogeographic realms** included Palearctic, Neotropical and other realm assignments where support permitted. We analysed `all_native`, `native_nonendemic` and `endemic` floristic strata where support thresholds were met. Persistence in native non-endemics was used to show that a pattern was not confined to island endemics; it was not interpreted as a causal endemicity contrast.

### Trait database and evidence hierarchy

Each accepted species could contribute to three primary axes:

1. `flower_colour`;
2. `floral_structural_complexity`;
3. `reproductive_assurance`.

The fixed denominator was therefore 318,885 species-by-axis cells. The final snapshot resolved 82,556 colour cells (77.67%), 91,635 structural cells (86.21%) and 48,497 reproductive-assurance cells (45.62%), for a total of 222,688 resolved cells (69.83%).

Evidence precedence was fixed. Species-direct High/Medium evidence preceded trait-specific validated Low evidence. Family-level trait imputation, global fallback and post hoc threshold relaxation were prohibited. Missing evidence was never converted to trait absence. Reproductive evidence retained distinctions among self-incompatibility, mating system and autonomous selfing where source methodology supported them. Self-compatibility alone was not mapped to autonomous selfing or realized selfing.

### Progressive analysis contract

All trait waves were materialized into the same `chapter1_trait_snapshot_v1` schema and reanalysed under the frozen `chapter1_progressive_analysis_v1` contract. The contract fixed scientific hypotheses, estimands, model order, support gates, multiplicity rules, source/lineage/area safeguards and the causal claim ceiling. Biological estimates and significance were allowed to change as evidence improved; a changed conclusion was treated as a scientific result rather than a pipeline failure.

#### H1: universal-syndrome rival

H1 asked whether increasing isolation produced one coherent floral/reproductive response across contexts. A universal syndrome was not inferred from a significant slope in one region and a non-significant slope in another. Evidence required multivariate within-context support and direct between-context response-vector comparisons.

#### H2: biogeographic branching

H2 tested whether response vectors differed among prespecified contexts. The primary plant response contained two axes:

- `accessibility_generalization`;
- `reproductive_assurance`.

Within-context tests used joint cluster-robust Wald tests of the two-axis distance-response vector. Direct between-context tests compared the vectors themselves. The same response had to satisfy the same support tier in both contexts before a pairwise comparison was promoted.

#### H3: source and lineage assembly

Supported H2 responses were decomposed against outcome-blind source expectations. GIFT source flora and frozen source definitions were used to estimate source-matched trait positions. The taxonomic-depth analysis evaluated the response at three stages:

1. observed assemblage response;
2. residual after family composition;
3. residual after source-matched genus composition.

Where supported, we additionally distinguished genus entry from within-represented-genus species loading. Family and genus were grouping structures, not trait-imputation devices.

Conceptually, if `T_source` is the regional source-pool expectation and `T_island(d)` is the island assemblage at separation `d`, the biological target is the change in `ΔT(d)=T_island(d)-T_source`. The primary H2 model does not fit this counterfactual literally. Instead, H2 estimates the geographic pattern and H3 asks how much of that pattern is compatible with source/lineage assembly.

#### H4: continuous area moderation

H4 tested whether isolation-associated filtering was stronger on smaller islands using continuous distance-by-continuous-area interactions. We did not introduce a post hoc small-versus-large island threshold. Promotion to a founder, capacity or pollinator-persistence mechanism required frozen equal-island, capped-information, direct-only, common-support and heteroskedastic-null safeguards.

#### H5: pollination-channel mechanism boundary

H5 was defined prospectively but was not used to generate the Chapter 1 headline. It requires an independently measured source pollination channel, retained/disrupted/structurally-absent state, realized visitation, single-visit effectiveness and effective service, followed by evidence that channel information explains residual plant filtering beyond H2–H4 safeguards. Floral phenotype, pollination-syndrome concordance and climatic suitability alone cannot satisfy H5. H5 is therefore treated here as an explicit identification boundary and handoff, not as a missing final model that is implicitly assumed to be true.

### Support gates and multiple testing

Global trait fill fraction was descriptive and was never used as a model admission threshold. Response-specific support was classified as:

- <30 supported islands: not promoted;
- 30–49: pilot;
- ≥50: confirmatory count support.

Pairwise context comparisons required the same response to meet the same tier in both contexts. BH-FDR correction was applied within prespecified test families. Multivariate gates preceded atomic trait interpretation.

### Separating reproductive assurance from floral architecture

To distinguish reproductive from floral components, `selfing_core` included only self-incompatibility, mating-system and autonomous-selfing information. Flower colour, size and access traits were excluded. A separate attraction/access pathway contrasted generalized-accessible and large-bee-like plant architecture.

We fitted an unconditional attraction/access distance response and a conditional model including `selfing_core`. Persistence of the distance coefficient after selfing-core adjustment was interpreted as statistical decomposition: the floral response could not be reduced to the measured selfing core. It was not interpreted as causal mediation.

### Pollination-associated floral architecture

Secondary concordance analyses used frozen `large_bee_like`, `butterfly_like` and `bird_like` templates. The templates combined prespecified aspects of colour, floral form, symmetry, tube depth and flower size. They were evaluated only after the primary plant response and were repeated in all-analysis and direct-only evidence as well as source-adjusted analyses.

Because the templates shared many traits, V4 estimated a source-trained first common factor and then examined template-specific residuals. This prevents a simultaneous increase in several guild-labelled scores from being interpreted as simultaneous increases in several realized pollinator guilds.

### Missingness and robustness

V5 evaluated trait-resolution missing-not-at-random sensitivity over a finite prespecified odds-ratio grid. V1 evaluated climate overlap using outcome-blind common-support weighting and distance-by-climate interactions. V2 evaluated family-to-genus taxonomic depth. V3 tested whether area patterns survived measurement and heteroskedastic-null safeguards. V4 decomposed shared versus template-specific floral architecture. These procedures were part of the claim boundary rather than optional post hoc robustness checks.

### Progressive trait-wave stopping rule

Wave comparisons were retained as prospective evidence history. The final computational freeze was not defined by a target percentage. We considered the campaign saturated when (i) a large increase in analysis-usable evidence no longer changed the main H1–H5 claim ceiling, (ii) realized yield from new strict sources declined sharply, and (iii) the major remaining mechanism question required a different data type rather than more records of the same traits.

---

## Results

### Trait evidence increased substantially while the latest hypothesis structure stabilized

The final snapshot contained 222,688 resolved species-by-axis cells (69.83%). Relative to Wave52, which contained 184,917 analysis-usable cells, 37,771 additional cells became usable. Reproductive-assurance coverage increased by 11,315 cells. The large gain was not equivalent to discovery of a second large source database: much of the improvement came from making previously reported but non-materializable evidence analytically usable under the fixed ontology and provenance rules.

New strict reproductive sources simultaneously entered strong diminishing returns. One high-yield Orchidaceae source contributed 192 new strict reproductive cells, whereas subsequent reviewed sources generally contributed only one to three cells each. A residual audit of 1,187 staging files and more than 2.25 million rows found no immediately strict-ready reproductive records under the frozen contract. Thus the late-stage bottleneck shifted from computation toward provenance and source-specific interpretation.

Importantly, the large Wave52-to-final increase did not reverse the current H1–H5 claim ceiling. The main biological conclusions below were therefore retained as the working Chapter 1 freeze.

### H1: one universal floral island syndrome was not recovered

The plant response did not converge on a common direction across context layers. Regional responses differed in both strength and sign, and direct between-context comparisons confirmed heterogeneity in supported contrasts.

For example, in direct-only evidence the northern-midlatitude versus tropical primary response vectors differed in both all-native assemblages (χ²=8.91, df=2, q=0.0116; 367 islands) and native non-endemics (χ²=11.62, df=2, q=0.00300; 367 islands). In all-analysis evidence, the same difference was supported in native non-endemics (χ²=14.63, df=2, q=0.000664) but not in all natives (q=0.197). These evidence-scope differences were retained rather than tuned away.

The absence of one universal vector means that “island syndrome” cannot be treated as a single global trajectory from isolation to selfing and floral simplification.

### H2: Palearctic and tropical assemblages followed different source-distance trajectories

The Palearctic showed the clearest broad primary response. In all-analysis all-native assemblages, accessibility/generalization increased with source separation (slope=0.0795, SE=0.0261, q=0.00463; 144 islands) and reproductive assurance also increased (slope=0.0534, SE=0.0173, q=0.00463; 145 islands). The joint two-axis vector was strongly supported (χ²=41.38, df=2, q=2.07×10^-9). Native non-endemics showed the same directions: accessibility/generalization 0.0671 (q=0.00279) and reproductive assurance 0.0513 (q=0.0122), with a supported joint vector (q=1.10×10^-8).

The result strengthened rather than disappeared in direct-only evidence. In all natives, accessibility/generalization was 0.0616 (SE=0.0264, q=0.0389; 139 islands) and reproductive assurance was 0.0966 (SE=0.0126, q=8.48×10^-14; 144 islands); the joint vector had q=2.44×10^-14. Native non-endemics again showed both positive axes (0.0649, q=0.0297; 0.0958, q=2.14×10^-12).

Tropical assemblages followed a different branch. In direct-only all natives, accessibility/generalization decreased with distance (slope=-0.1005, SE=0.0296, q=0.00275; 127 islands), while reproductive assurance increased (0.1360, SE=0.0443, q=0.00428; 124 islands). The joint vector was strongly supported (χ²=20.03, df=2, q=8.96×10^-5). Native non-endemics showed almost the same pattern: accessibility/generalization -0.0942 (q=0.00325) and reproductive assurance 0.1361 (q=0.00474), with joint q=5.99×10^-5.

All-analysis tropical estimates were less uniformly resolved at the all-native level: the two individual axes had q=0.0947 and the joint vector had q=0.108. However, native non-endemics retained a strong opposite-signed branch (accessibility/generalization -0.0856, q=0.00357; reproductive assurance 0.1466, q=6.71×10^-7; joint q=2.27×10^-6). Thus the tropical counter-pattern was not a simple artifact of one evidence scope, but its strongest all-native support came from direct evidence.

These results establish region-associated branching but not a causal effect of region identity. V1 climate-overlap tests retained a poor-overlap contrast, and the Palearctic–Neotropical primary vector difference was not supported in the unweighted realm comparison (all-native q=0.284; native-nonendemic q=0.396). We therefore treat biogeographic labels as context descriptors whose underlying climatic and historical mechanisms remain only partly identified.

### Reproductive assurance and floral architecture were partially separable

The Palearctic floral response persisted after conditioning on the measured selfing core. Across four source definitions, the all-native conditional attraction/access distance coefficient ranged from 0.0910 to 0.1007, with FDR q values from 0.00791 to 0.0340 (142 islands). The corresponding selfing-core coefficients were small and non-significant (P=0.77–0.95). Native non-endemics showed the same result: conditional distance coefficients 0.0825–0.0917 (q=0.00409–0.0270), while selfing-core coefficients remained non-significant.

Thus the Palearctic floral-architecture shift was not statistically reducible to the measured reproductive-assurance component. This does not identify relaxed pollinator-mediated selection, but it rejects a purely serial interpretation in which selfing alone accounts for the floral response.

Tropical source-adjusted attraction/access responses moved in the opposite direction. Across four source definitions, all-native conditional distance coefficients ranged from -0.0970 to -0.1022 (q=0.00102–0.00303; 126 islands). Native non-endemics ranged from -0.1097 to -0.1169 (q=0.000274–0.00169; 121 islands). In these conditional models, selfing-core coefficients were again not significant. Combined with the positive primary reproductive-assurance slopes, this shows that reproductive assurance can increase while floral architecture shifts toward the specialized/attractive side rather than toward generalized access.

The empirical pattern is therefore better represented by partially parallel branches:

```text
source separation
      |\
      | \----> reproductive assurance
      |
      \------> floral architecture

sign and coupling depend on context
```

### H3: the broad Palearctic response was compatible with genus-level lineage assembly

Taxonomic-depth decomposition produced the strongest change across progressive trait waves. At Wave36, the broad Palearctic two-axis response appeared to survive both family and genus adjustment. By Wave52, increased source and reproductive evidence changed that conclusion. The response remained after family composition but was no longer supported after source-matched genus composition.

The final snapshot reproduced the Wave52 classification exactly. In all-analysis/direct-only evidence and all-native/native-nonendemic strata, the Palearctic response was supported under all four source modes at the observed stage (4/4), retained after family adjustment (4/4), and unsupported after genus adjustment (0/4). The classification was therefore `compatible_with_genus_level_assembly_beyond_family` in all four combinations.

This result is central to interpretation. It means that the broad distance gradient can be biologically real at the assemblage level while being generated primarily by which source-available genera are represented on islands. The result does not support a claim that distance repeatedly transformed the same lineages into the same floral phenotype.

### Pollination-associated floral templates diverged strongly, but mostly along a shared architecture axis

Secondary template analyses connected the plant response to pollination theory without assigning realized pollinators. In the Palearctic all-analysis all-native dataset, butterfly-like concordance declined with distance (slope=-0.0603, q=0.00915) and large-bee-like concordance declined more strongly (-0.1037, q=0.00915), whereas bird-like concordance was unsupported. Native non-endemics showed similar declines. Direct-only data retained both butterfly-like (-0.0580, q=0.0147) and large-bee-like (-0.0685, q=0.0147) declines.

Tropical assemblages moved in the opposite direction. In all-analysis all-native data, bird-like (0.0588, q=7.2×10^-5), butterfly-like (0.0572, q=0.0113) and large-bee-like (0.0463, q=7.2×10^-5) scores all increased with distance. Native non-endemics showed the same pattern. Direct-only results again supported all three axes.

A hard pollinator-classification interpretation would be biologically incoherent: all three guild-labelled scores can rise together. V4 showed why. A source-trained first factor explained 86.94% of variation among the three templates in all-analysis evidence and 86.44% in direct-only evidence. The guild-labelled scores therefore primarily share a common plant architecture dimension. Tropical triple increases are best interpreted as maintenance or strengthening of a shared specialized/attractive architecture, not as evidence that birds, butterflies and large bees all increased in realized importance.

Template-specific residuals were retained only as secondary plant-side contrasts. Some tropical residual structure survived source/genus adjustment, including positive large-bee-like residuals in several source-mode combinations, but these cannot identify large-bee occurrence, visitation, effectiveness or replacement.

### H4: area remained a measurement-sensitive modifier

All 16 frozen primary V3 classifications remained `retain_area_as_measurement_sensitive_modifier_only`. Direct-only analyses did not contradict the main directional patterns, but zero of 16 heteroskedastic-null gates passed. We therefore did not promote distance-by-area patterns to a founder-filtering, habitat-capacity or pollinator-persistence mechanism.

This negative promotion decision is part of the result. Island area may alter precision or the apparent magnitude of filtering, but the current data do not identify the biological process responsible.

### Missingness and observation safeguards bounded but did not erase the broad conclusions

The broad Palearctic plant-response vector survived the finite predeclared V5 MNAR grid in both evidence scopes and both focal floristic strata, although individual reproductive details remained assumption-sensitive. The missingness analysis therefore supports a bounded claim: the broad response is not dependent on one simple missing-at-random assumption, but arbitrary unmeasured state-dependent resolution cannot be ruled out.

V1 climate-overlap analysis also prevented overinterpretation. Regional response heterogeneity is real at the plant-pattern level, but common measured-climate support is incomplete for some contrasts, and a measured-climate-independent categorical “realm effect” is not established.

### H5 remained outside the Chapter 1 claim

No Chapter 1 result identifies historical pollinator loss, pollinator replacement or effective service. The floral templates are plant traits, not visitor data. Climatic compatibility is not realized channel loss. Occurrence is not visitation, visitation is not per-visit effectiveness, and per-visit effectiveness is not rate-weighted effective service.

Accordingly, H5 is not used to rescue or complete the present paper. It remains the explicit mechanistic extension generated by the plant-side results.

---

## Discussion

### The floral island syndrome is an assemblage outcome, not one deterministic syndrome

The central result is not simply that a classical island syndrome failed. The data show a more informative structure: source separation reorganizes reproductive assurance and floral architecture along different trajectories, and those trajectories depend on biogeographic context.

The Palearctic approximates one familiar island-syndrome direction. More remote assemblages contain a stronger accessibility/generalization component and greater reproductive assurance. Yet even there, the floral shift is not reducible to measured selfing alone. In the tropical analysis regime, the two components decouple even more clearly: reproductive assurance increases while specialized/attractive architecture is retained or strengthened. Thus the simple sequence `isolation → selfing → reduced floral specialization` cannot serve as a general model.

This distinction helps reconcile apparently conflicting island studies. A local archipelago may show a coherent syndrome because multiple filters align in one direction. Another archipelago may show only a reproductive response, only a floral response, or opposite changes in the two components. The appropriate global object is therefore not a single syndrome score but a multivariate response architecture.

### Distance should be interpreted as a source-accessibility gradient

A useful consequence of the analysis is a shift away from binary island-versus-mainland thinking. The ecological reference is the regional source system, and distance indexes increasing separation from that reference. In a conceptual limiting sense, `d → 0` represents weak oceanic source filtering; greater distance creates more opportunity for differential arrival, establishment and persistence.

This does not imply that the formal `distance_to_continent` variable is the exact historical source distance for every species. The analysis deliberately separates the raw geographic exposure from H3 source matching. That separation is important: a distance slope can reflect both actual geographic filtering and pre-existing regional differences in lineage composition. H3 asks whether the gradient survives when source-matched taxonomic expectations are made explicit.

The `4/4 → 4/4 → 0/4` result shows that much of the strongest Palearctic gradient is compatible with genus-level source filtering. This is not a nuisance correction that weakens the ecological story; it identifies the level at which the story occurs. Island floral filtering can be real because isolation changes which floral lineages are represented.

### Source and lineage assembly provide a plant-side mechanism for the distance response

The strong genus-depth result places the study close to the logic of colonization filtering. Pre-existing combinations of floral access and reproductive strategy are sorted as source-available lineages cross geographic and ecological barriers. This is conceptually related to Baker’s law, but extends the assembly perspective from reproductive traits to a multivariate floral/reproductive phenotype.

The result also imposes a clear evolutionary boundary. Cross-sectional assemblage data cannot show that individual species evolved their current floral state after colonization. A source-matched genus signal can arise from dispersal, establishment, persistence or unmeasured traits correlated with genus membership. Repeated within-lineage evolution remains possible in particular clades, but it is not required to generate the broad macroecological pattern observed here.

This framing turns lineage from a confounder into a biological outcome. The next plant-side question is not simply whether genus should be “controlled”, but why particular source-available genera are disproportionately represented along the distance gradient and which dispersal, establishment or demographic traits generate that sorting.

### Selfing and pollination-associated floral responses can run in parallel

The decomposition of `selfing_core` and attraction/access architecture provides the paper’s clearest biological result beyond source assembly. In the Palearctic, attraction/access change remains after selfing-core adjustment. In the tropics, attraction/access moves toward the specialized side even while the primary reproductive-assurance axis increases.

This argues against using selfing syndrome as a universal explanation for floral simplification. Reproductive assurance and floral architecture are linked by ecology and evolution, but they are not the same axis. A colonist may gain or already possess reproductive assurance while retaining a flower that continues to depend on or benefit from effective animal visitation. Conversely, altered interaction structure may change floral access and attraction before—or without—a detectable change in reproductive mode.

The result is especially useful for interpreting island floras because assemblage turnover can mix these pathways. One lineage may be retained because it is reproductively assured, another because its pollination channel persists, and a third may disappear altogether. A community-level “syndrome” can emerge from the frequencies of these response modes rather than from a shared within-lineage trajectory.

### Pollination syndromes are informative as architectures, not as visitor identities

The guild-labelled analyses support a cautious return to pollination theory. Palearctic assemblages move away from large-bee-like and butterfly-like architecture with distance, whereas tropical assemblages move toward a shared architecture that scores positively on bird-, butterfly- and large-bee-like templates. These contrasts are ecologically meaningful, but the V4 decomposition shows why pollinator identity cannot be read directly from them.

The three templates share most of their variance. Their labels should therefore be interpreted as sampled multivariate floral architectures associated in the literature with functional groups, not as deterministic pollinator assignments. This is consistent with broader pollination-syndrome work showing that trait combinations can predict functional tendencies without perfectly identifying realized or effective visitors.

The tropical result is particularly informative. If bird-, butterfly- and large-bee-like scores all rise together, the useful conclusion is not that all three pollinator guilds increased. Rather, remote tropical assemblages retain or strengthen a shared specialized/attractive architecture. The residual template structure can generate hypotheses for future interaction work, but it cannot establish functional replacement.

### A double geographic filter provides a mechanistic hypothesis for regional branching

The plant-side result motivates, but does not prove, a broader double-filter model. The same oceanic barrier that filters plant propagules also acts on pollinator dispersal and establishment. For a plant lineage `z`, source continuity can be represented conceptually as

`E_plant(z,d) = f(source availability, distance, propagule dispersal, establishment, persistence)`.

For a pollination channel `g`, functional continuity can be represented as

`E_poll(g,d) = f(source availability, distance, flight/dispersal, establishment, habitat, realized community)`.

There is no reason for these functions to have identical shapes, either across plants or across pollinator groups. The same 100 km of ocean can be a severe barrier for one lineage or functional group and a weak barrier for another. Flight ability is only one component: successful channel continuity also requires arrival, establishment, suitable habitat and enough abundance to provide effective service.

This asymmetry offers a biological explanation for why reproductive and floral components branch among regions. Where source-pool filtering and pollination-channel disruption align, a classical generalization/assurance pattern may emerge. Where alternative effective channels persist, reproductive assurance can increase without floral simplification. Where a channel is structurally absent from the source region, its absence on islands is not a loss event at all.

We emphasize that this pollinator-side filter remains a hypothesis in Chapter 1. The current global plant database does not estimate guild-specific ocean-crossing ability, retained-versus-disrupted channel states, visitation, single-visit effectiveness or effective service. The value of the double-filter framework is that it specifies the missing links rather than silently treating floral phenotype as evidence that those links occurred.

### Climate and area are boundaries on interpretation, not failed covariates

The climate-overlap and area-falsification results are also biologically informative. Regional differences should not be narrated as intrinsic effects of biogeographic labels when measured climate support is insufficient to transport the contrast. Similarly, distance-by-area directionality should not be upgraded to founder filtering or pollinator persistence when the same pattern fails the prespecified heteroskedastic-null criterion.

These safeguards narrow the claim to what the data identify: source-distance-associated plant assemblage branching and its taxonomic depth. They leave open several biological explanations—climate, colonization history, demographic persistence and interaction structure—that can be addressed with additional data types rather than by changing the current statistical story.

### Progressive trait waves define a stopping rule rather than a novelty claim

The database campaign is methodologically useful because it shows that the biological story is no longer obviously an artifact of one incomplete snapshot. Wave52-to-final integration added 37,771 usable cells, including more than eleven thousand reproductive-assurance cells. This was large enough to have changed the result if missing coverage were the dominant explanation.

It did not. The current H1–H5 claim ceiling remained stable, including the genus-level interpretation that had already changed between earlier waves. At the same time, strict new source acquisition became increasingly expensive and low-yield. The appropriate conclusion is not that 69.83% represents a universal sufficiency threshold. It is that indiscriminate acquisition is no longer the most efficient way to change the inference for the currently supported contexts.

Poorly supported contexts remain legitimate targets for focused data collection. However, the main unresolved question of the present paper is no longer “what is the flower colour of another thousand species?” It is mechanistic: how plant-source filtering and interaction-channel filtering combine to generate the regional response architectures already observed.

### From WHEN/WHERE to HOW/WHY

Chapter 1 therefore closes at a useful inferential boundary. It identifies **when and where** source separation is associated with floral/reproductive assemblage filtering, demonstrates that the response components can decouple, and shows that the strongest broad Palearctic response is compatible with genus-level assembly. It does not identify the historical interaction pathway that generated the residual pattern.

The next mechanistic layer should test the pollinator side of the double filter without using floral phenotype as a proxy. The required chain is explicit:

```text
source pollination channel
        -> retained / disrupted / structurally absent
        -> realized visitation
        -> single-visit effectiveness
        -> rate-weighted effective service
        -> reproductive dependency / assurance
        -> plant response
```

That chain is deliberately outside the present paper’s causal claims. It provides the handoff to the companion mechanistic programme and prospective island field tests, where visitation, pollen deposition and reproductive treatments can be measured in the same linked populations.

---

## Conclusions

A universal floral island syndrome was not recovered. Instead, increasing separation from regional source systems was associated with biogeographically contingent assemblage trajectories whose reproductive-assurance and floral-architecture components could couple or decouple.

The Palearctic showed the clearest broad classical direction: more distant assemblages were more accessible/generalized and more reproductively assured. Yet the floral component persisted after conditioning on selfing core, and the broad response disappeared after source-matched genus adjustment. Tropical assemblages provided the complementary counter-pattern: reproductive assurance could increase while specialized/attractive floral architecture was maintained or strengthened.

These results support three conclusions. First, floral island syndromes are better treated as emergent assemblage outcomes than deterministic lineage-level rules. Second, source-pool and lineage assembly are part of the biological mechanism, not merely confounders. Third, pollination-associated floral architecture can generate mechanistic hypotheses but cannot identify realized pollinator loss or replacement without independent interaction data.

We therefore interpret oceanic distance as a source-accessibility gradient that creates a **plant-side geographic filter** demonstrable in the current data and potentially a **pollinator-side geographic filter** that remains to be tested. Because plant lineages and pollination channels differ in dispersal and establishment ability, the same geographic barrier need not generate the same island syndrome.

The progressive analysis also provides a practical stopping rule. Large late increases in trait coverage no longer changed the main inferential structure, while the remaining mechanistic uncertainty requires a different data type. Chapter 1 can therefore be closed as a global WHEN/WHERE and plant-assembly paper, with the pollinator-side causal extension reserved for direct mechanistic testing.

---

## Proposed main figures

### Figure 1. Double geographic-filter hypothesis and fixed PR142 inference ladder

Panel A: source-pool continuum from conceptual `d→0` to remote islands.  
Panel B: plant source filter versus candidate pollinator-channel filter.  
Panel C: H1 universal rival → H2 regional branching → H3 source/lineage assembly → H4 area moderation → H5 independent mechanism boundary.

**Key visual rule:** plant-side arrows supported by Chapter 1 should be visually distinct from the untested pollinator-side extension.

### Figure 2. Biogeographic branching of the primary two-axis response

Plot distance slopes ± 95% CI for `accessibility_generalization` and `reproductive_assurance` across analysis regimes and formal realms, with all-analysis and direct-only estimates side by side.

Emphasize:

- Palearctic: both axes positive;
- tropical: accessibility negative, reproductive assurance positive in direct-only and native-nonendemic all-analysis;
- unresolved contexts shown as unresolved, not zero.

### Figure 3. Decoupling of reproductive assurance and floral architecture

Panel A: Palearctic source-adjusted attraction/access slope before and after conditioning on selfing core across four source definitions.  
Panel B: tropical equivalent, showing the opposite floral direction.  
Panel C: conceptual parallel-pathway diagram replacing a compulsory `selfing → simplification` sequence.

### Figure 4. Taxonomic-depth decomposition of the Palearctic response

Four evidence-scope × stratum combinations, each showing:

`observed 4/4 → after family 4/4 → after genus 0/4`.

Include the Wave36 → Wave52 → final transition to show that the genus-level interpretation stabilized after evidence expansion.

### Figure 5. Pollination-associated architecture without pollinator assignment

Panel A: Palearctic large-bee-like / butterfly-like / bird-like distance slopes.  
Panel B: tropical slopes.  
Panel C: first-factor variance explained (86.94% all-analysis; 86.44% direct-only) and template residuals.

### Figure 6. Progressive evidence saturation and claim ceiling

Panel A: Wave52 versus final usable cells.  
Panel B: reproductive coverage gain.  
Panel C: gain per strict source batch.  
Panel D: H1–H5 decision matrix across waves, showing which conclusions changed and which stabilized.

---

## Proposed main tables

### Table 1. Data universe and final trait coverage

| Quantity | Final value |
|---|---:|
| Islands in fixed universe | 8,265 |
| Accepted angiosperm species | 106,295 |
| Species-axis denominator | 318,885 |
| Resolved species-axis cells | 222,688 (69.83%) |
| Flower colour | 82,556 (77.67%) |
| Floral structural complexity | 91,635 (86.21%) |
| Reproductive assurance | 48,497 (45.62%) |

### Table 2. Key primary distance responses

| Context / scope / stratum | Accessibility slope | q | Reproductive slope | q | Joint-vector q |
|---|---:|---:|---:|---:|---:|
| Palearctic, all-analysis, all-native | 0.0795 | 0.00463 | 0.0534 | 0.00463 | 2.07×10^-9 |
| Palearctic, direct-only, all-native | 0.0616 | 0.0389 | 0.0966 | 8.48×10^-14 | 2.44×10^-14 |
| Palearctic, direct-only, native-nonendemic | 0.0649 | 0.0297 | 0.0958 | 2.14×10^-12 | 2.00×10^-12 |
| Tropical, direct-only, all-native | -0.1005 | 0.00275 | 0.1360 | 0.00428 | 8.96×10^-5 |
| Tropical, direct-only, native-nonendemic | -0.0942 | 0.00325 | 0.1361 | 0.00474 | 5.99×10^-5 |
| Tropical, all-analysis, native-nonendemic | -0.0856 | 0.00357 | 0.1466 | 6.71×10^-7 | 2.27×10^-6 |

### Table 3. PR142 hypothesis decisions

| Hypothesis | Final Chapter 1 decision |
|---|---|
| H1 universal syndrome | not supported |
| H2 biogeographic branching | supported at plant-pattern level; climate-independent categorical causation not established |
| H3 source/lineage assembly | broad Palearctic response compatible with genus-level assembly beyond family |
| H4 area/capacity | measurement-sensitive modifier only |
| H5 independent pollination channel | not claimed in Chapter 1; requires independent channel chain |

---

## Working literature anchors

The reference list should be finalized against source-locked bibliographic records before submission. The current narrative is anchored by:

- Baker HG (1955). Self-compatibility and establishment after long-distance dispersal. *Evolution*.
- Pannell JR & Barrett SCH (1998). Baker’s law revisited: reproductive assurance in a metapopulation. *Evolution*.
- Pannell JR (2015). Evolution of the mating system in colonizing plants. *Molecular Ecology*.
- Grossenbacher DL et al. (2017). Self-compatibility is over-represented on islands. *New Phytologist*.
- Sicard A & Lenhard M (2011). The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany*.
- Fenster CB et al. (2004). Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics*.
- Rosas-Guerrero V et al. (2014). A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters*.
- Hetherington-Rauth MC & Johnson MTJ (2020). Floral trait evolution of angiosperms on Pacific islands. *The American Naturalist*.

Additional island source-pool, lineage-filtering and pollination-network references should be inserted from the existing source-audited project bibliography rather than reconstructed from memory.

---

## Reproducibility anchors

- Scientific contract: `config/chapter1_progressive_analysis.yml` (`chapter1_progressive_analysis_v1`).
- Latest trait integration: Run `34191508045`, 222,688 / 318,885 resolved cells.
- Latest PR142 progressive reanalysis: Run `34232450884`.
- Analysis artifact: `chapter1-progressive-analysis-34232450884`, artifact ID `10058653212`.
- Artifact digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- H5 causal extension is intentionally outside the present manuscript claim and is tracked separately under the independent channel-chain contract.
