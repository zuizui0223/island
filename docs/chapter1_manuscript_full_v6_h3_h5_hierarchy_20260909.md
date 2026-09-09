# The island floral syndrome is not one syndrome: source-distance filtering, lineage assembly, and a double geographic filter

## Full working manuscript v6 — 2026-09-09

### Abstract

Island floras are often expected to converge toward a coherent floral “island syndrome”, combining increased reproductive assurance with reduced floral specialization or attraction. Yet geographic isolation does not act directly on floral phenotype. An oceanic barrier first changes which plant lineages from a regional source pool can disperse, arrive, establish and persist; the same barrier can also change which regionally available pollination channels reach, establish and remain functionally available. Because plant lineages and pollinator functional groups differ in dispersal and establishment ability, the same geographic distance need not produce one floral or reproductive outcome.

We tested whether increasing source separation generates one universal floral/reproductive response or instead reorganizes partially separable reproductive and floral components in a biogeographically contingent manner. We analysed a fixed universe of 8,265 islands and 106,295 accepted angiosperm species using a progressive analysis contract in which hypotheses, support gates, source-pool safeguards, model order and causal claim ceilings were fixed while trait evidence improved. The final snapshot resolved 222,688 of 318,885 possible species-by-axis cells (69.83%) across flower colour, floral structural complexity and reproductive assurance. The geographic exposure was log-transformed mainland distance, interpreted as a composite gradient of geographic separation, connectivity and source accessibility rather than as a mechanistically pure causal treatment.

A universal floral/reproductive island syndrome was not recovered. In the Palearctic, increasing separation was associated with greater accessibility/generalization and greater reproductive assurance. In all-analysis all-native assemblages, slopes were 0.0795 and 0.0534 per unit of log-distance, respectively (both FDR q=0.00463); direct-only estimates were 0.0616 (q=0.0389) and 0.0966 (q=8.48×10^-14). The pattern persisted in native non-endemics. Tropical assemblages followed a different trajectory: in direct-only all-native data, accessibility/generalization declined with distance (-0.1005, q=0.00275) while reproductive assurance increased (0.1360, q=0.00428). Thus reproductive assurance and floral architecture did not behave as one obligatory serial syndrome.

Source-matched taxonomic decomposition showed that the broad Palearctic response survived family adjustment but disappeared after genus adjustment in all four evidence-scope-by-floristic-stratum combinations (4/4 → 4/4 → 0/4). The Palearctic attraction/access shift nevertheless remained positive after conditioning on selfing core across four source definitions (distance estimates 0.091–0.101 in all natives; q=0.0079–0.034). Secondary large-bee-like, butterfly-like and bird-like floral templates diverged among regions, but a source-trained common factor explained 86.94% and 86.44% of template variance in all-analysis and direct-only evidence, respectively, indicating that these labels mainly describe overlapping plant architecture rather than realized pollinator identity.

We interpret H3 and H5 as different levels of the same double-filter problem. H3 directly identifies the plant-side level at which the broad response is represented: differential source-matched genus assembly. H5 remains a candidate pollinator-side explanation for why such lineages are differentially filtered, as well as for any response that persists beyond genus. Consequently, disappearance after genus adjustment does not falsify a pollination mechanism; a pollination-channel disruption could itself alter lineage entry or persistence. Pollination-syndrome concordance is therefore treated as an intermediate compatibility layer, not as direct evidence of channel loss.

The results recast the floral island syndrome as an emergent assemblage-level outcome with at least two partially separable components: reproductive assurance and pollination-associated floral architecture. The strongest plant-side distance response is compatible with genus-level source filtering, while measured-climate-independent regional causation, an area/capacity mechanism and pollinator-channel causation remain beyond the present identification ceiling. Large additions of trait evidence no longer changed this claim structure, providing an inference-based stopping rule for the computational campaign.

**Keywords:** island biogeography; floral traits; reproductive assurance; selfing syndrome; pollination syndrome; source pool; dispersal filtering; lineage assembly; mutualism; biogeographic contingency

---

## Introduction

Geographic isolation is one of the defining processes of island biogeography. Classical island theory treats separation from a regional source pool as a constraint on colonization opportunity, while empirical island biology has repeatedly identified characteristic shifts in dispersal, growth form, reproductive systems and ecological interactions. In plants, these observations have motivated the idea of an “island syndrome”: a recurrent suite of traits associated with long-distance colonization, small population size, altered environments and changed biotic interactions.

Floral versions of this idea typically combine at least two biological arguments. The first is reproductive assurance. Long-distance colonists may experience mate limitation or unreliable pollen transfer, creating an advantage for self-compatibility, autonomous selfing or other forms of uniparental reproduction. This argument is closely related to Baker’s law and its modern interpretation as a colonization filter: reproductive traits may be over-represented on islands because they increase establishment probability, not because the same reproductive transition necessarily evolves repeatedly after arrival. The second argument concerns pollination-associated floral architecture. If island pollinator assemblages are functionally altered or dominated by different visitor groups, the returns to conspicuous displays, deep floral tubes, restricted access or precise mechanical matching may change.

These processes can generate superficially similar phenotypes, but they are not equivalent mechanisms. Self-compatibility is not autonomous selfing; autonomous selfing is not necessarily the dominant realized mating mode; and a plant can retain specialized floral architecture while gaining reproductive assurance. Conversely, floral attraction or access traits can shift without a corresponding change in the measured selfing core. Treating all of these outcomes as one syndrome therefore risks turning a broad pattern into a presumed causal sequence.

A second complication is that island floral patterns are assemblage patterns. A distant island can differ from a near island because different source-available lineages arrived, established or persisted, even if no lineage changed phenotype after colonization. Baker’s law itself can be understood this way: pre-existing reproductive variation is sorted during colonization. The same logic extends to floral form, symmetry, tube depth and access. If lineages differ in both dispersal/establishment ability and floral strategy, geographic isolation can change community-level trait composition through differential lineage representation.

This source-pool perspective changes the meaning of distance. We do not treat kilometres from a continent as a direct causal force on floral phenotype. Instead, distance is a composite axis of separation, connectivity and source accessibility. Conceptually, as distance approaches zero, an island approaches a high-accessibility boundary where geographic filtering relative to a regional source system should be weaker. As distance increases, the opportunity for source-to-island filtering increases. This limiting argument does not imply that a literal zero-distance island equals mainland vegetation, nor that the nearest continent is the exact historical source for every island. It provides a continuous alternative to a binary island-versus-mainland comparison.

The same geographic barrier can also act on the other side of a plant–pollinator interaction. A pollination channel that is regionally available must itself cross oceanic barriers, establish and persist before it can provide effective service on an island. Flight ability, long-distance dispersal, stepping-stone use, habitat requirements and establishment probability can therefore make a given ocean distance a different biological barrier for different pollinator functional groups. Structural absence of a pollinator from a regional source system is also fundamentally different from the loss of a channel that was available in that source system.

This motivates a **double geographic-filter framework**. Geographic separation can simultaneously filter plant source-pool assembly and pollination-channel continuity. Importantly, these two filters can interact. Pollinator-channel loss need not only change floral selection within surviving lineages; it can also alter which plant lineages are able to establish or persist. Thus a pollination mechanism may be expressed statistically as lineage assembly. This point becomes critical when interpreting source/genus adjustment: if a floral/reproductive gradient disappears after genus composition is controlled, that does not imply that pollinators were irrelevant. It means that any ecological mechanism producing the gradient may have acted at the level of lineage entry or persistence.

Global island comparisons also combine heterogeneous floristic histories, climates, island sizes, observation effort and taxonomic composition. A significant slope in one region and a non-significant slope in another does not demonstrate different regional responses. Similarly, lineage should not always be treated as a nuisance variable to be erased, because differential lineage representation is itself a possible biological outcome of island filtering. Missing trait evidence is also geographically and taxonomically structured, so increasing database coverage can change which questions are estimable.

We address these problems using the progressive analysis contract established in PR #142. Rather than allowing the model to change as trait coverage improved, we fixed the island universe, trait ontology, hypothesis order, support thresholds, multiple-testing families, source/lineage safeguards and causal claim ceiling. New trait snapshots were rerun from the beginning. This produced a sequential hypothesis ladder in which each result determines the next biological question while preventing the question itself from drifting with the data.

We test four plant-data hypotheses and define a fifth mechanistic extension. **H1** treats a universal floral island syndrome as a rival hypothesis: does increasing source separation generate one coherent floral/reproductive response vector across contexts? **H2** asks whether isolation-associated response vectors instead branch among biogeographic contexts. **H3** asks at what source/taxonomic depth a supported response remains, distinguishing broad assemblage patterns from family- and genus-level representation. **H4** asks whether continuous island area modifies isolation-associated filtering. **H5** asks whether independently measured pollination-channel change can explain the plant-side pattern after the source/lineage structure has been characterized.

We distinguish two future H5 routes. **H5a, channel-dependent lineage filtering**, asks whether source-channel retention or disruption changes which source-available lineages enter or persist, thereby generating H3. **H5b, beyond-genus channel response**, retains the original strict residual logic and asks whether independently measured channel change explains floral/reproductive response that persists after source/genus composition is held fixed. H5a and H5b are prospective causal extensions; neither is inferred from floral phenotype in the present paper.

The primary H2 response is deliberately pollinator-name-free. We separate (i) **accessibility/generalization**, representing increasingly open, shallow or radially accessible floral architecture, from (ii) **reproductive assurance**, representing strict reproductive-assurance evidence. Only after the primary response is established do we evaluate fixed large-bee-like, butterfly-like and bird-like floral templates. These are sampled pollination-associated architecture concordances, not pollinator classifiers.

Finally, we use the progressive trait campaign itself to decide when additional database acquisition ceases to be the highest-value scientific activity. The target is not complete global trait coverage. The relevant criterion is inference saturation: whether large additions of usable evidence continue to change the major hypothesis decisions, and whether remaining gaps can be recovered without disproportionate source-specific review or weaker evidential standards.

Our central prediction is therefore not that distant islands should be uniformly selfing, simple or inconspicuous. It is that increasing source separation should reorganize floral and reproductive assemblages, but the signs, coupling and taxonomic depth of those responses should depend on biogeographic context. If the floral island syndrome is an emergent assemblage outcome rather than one deterministic rule, reproductive assurance and floral architecture should sometimes decouple, and broad distance gradients should often resolve into lineage assembly.

---

## Materials and Methods

### Geographic and taxonomic universe

The fixed geographic universe contained 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. The denominator was not allowed to shrink between trait waves. Island occurrences, floristic-status information, island area, climate covariates, geographic exposure and source-pool infrastructure were held constant while trait evidence was allowed to improve.

The primary geographic exposure was `log1p_distance_to_continent_km`. We interpret this as a composite of geographic separation, connectivity and source accessibility rather than as a mechanistically pure isolation treatment. The conceptual near-source limit (`d → 0`) represents high source accessibility and relatively weak oceanic filtering, not a literal mainland observation or a claim that a zero-distance island should have mainland composition.

Analyses were conducted at two context layers. Broad analysis regimes included northern mid-latitude, northern high-latitude, tropical and southern extratropical strata. Formal biogeographic realms included Palearctic, Neotropical and other realm assignments where support permitted. We analysed `all_native`, `native_nonendemic` and `endemic` floristic strata where support thresholds were met. Persistence in native non-endemics was used to establish that a pattern was not confined to island endemics; it was not interpreted as a causal endemicity contrast.

### Trait evidence

Each accepted species could contribute to three primary axes: flower colour, floral structural complexity and reproductive assurance. The fixed denominator was therefore 318,885 species-by-axis cells. The final snapshot resolved 82,556 colour cells (77.67%), 91,635 structural cells (86.21%) and 48,497 reproductive-assurance cells (45.62%), for a total of 222,688 resolved cells (69.83%).

Evidence precedence was fixed. Species-direct High/Medium evidence preceded trait-specific validated Low evidence. Family-level trait imputation, global fallback and post hoc threshold relaxation were prohibited. Missing evidence was never converted to trait absence. Reproductive evidence retained distinctions among self-incompatibility, mating system and autonomous selfing where source methodology supported them. Self-compatibility alone was not mapped to autonomous selfing or realized selfing.

### Progressive hypothesis contract

All trait waves were materialized into the same `chapter1_trait_snapshot_v1` schema and reanalysed under `chapter1_progressive_analysis_v1`. The contract fixed hypotheses, estimands, model order, support gates, multiplicity rules, source/lineage/area safeguards and the causal claim ceiling.

H1 evaluated whether one coherent response vector recurred across contexts. H2 evaluated regional branching using multivariate within-context tests and direct between-context response-vector comparisons. H3 decomposed supported responses across observed, after-family and source-matched after-genus stages and evaluated lineage entry/loading where support permitted. H4 tested continuous distance-by-area moderation. H5 was withheld from primary Chapter 1 inference because independent channel exposure, retention/disruption, visitation and effective-service data were not present in the global plant database.

### Support gates and robustness

Global fill fraction was descriptive and never an analysis gate. Response-specific support was classified as <30 supported islands = not promoted, 30–49 = pilot and ≥50 = confirmatory count support. Pairwise context comparisons required the same response to meet the same support tier in both contexts. BH-FDR correction was applied within predeclared families. Multivariate response tests preceded interpretation of atomic trait states.

Mandatory robustness included all-analysis and direct-only evidence, equal-island and capped-information weighting, trait-resolution selection adjustment, outcome-blind source-pool sensitivity, alternative distance functional forms, leave-one-spatial-block-out analysis, genus-composition safeguards and continuous area moderation.

### Source and taxonomic-depth decomposition

Supported responses were decomposed against outcome-blind mainland source expectations using GIFT source flora and fixed source assignments. Family and genus were grouping structures, not trait-imputation devices.

Conceptually, if `T_source` is the expected regional source-pool trait composition and `T_island(d)` is the observed island composition at source separation `d`, the ecological deviation is `DeltaT(d) = T_island(d) - T_source`. The H2 model does not fit this expression literally. H2 estimates the geographic response, whereas H3 asks how much of that response is represented by source and lineage assembly.

### Decomposing reproductive and floral components

`selfing_core` included only self-incompatibility, mating-system and autonomous-selfing information and excluded floral attraction traits. A separate attraction/access pathway was constructed from floral architecture. Conditional models asked whether the isolation-associated attraction shift remained after conditioning on `selfing_core`. This analysis is a decomposition, not causal mediation.

Secondary pollination-associated concordance used fixed `large_bee_like`, `butterfly_like` and `bird_like` templates. These combined colour, floral form, symmetry, tube depth and flower size with prespecified weights. Because the templates share traits, V4 estimated a source-trained common factor before template-specific residuals were interpreted.

### H3–H5 causal hierarchy

H3 is an empirical decomposition of the plant assemblage. It identifies **where the distance response is represented**, not **why** that representation occurs. Genus-level assembly may reflect plant dispersal, habitat filtering, demographic persistence, pollination dependence or combinations of these.

We therefore distinguish two prospective H5 mechanisms:

- **H5a — channel-dependent lineage filtering:** `source channel -> retention/disruption -> effective service -> differential lineage establishment/persistence -> genus entry/loading -> assemblage trait pattern`.
- **H5b — beyond-genus channel response:** `source/genus composition fixed -> channel retention/disruption -> visitation -> per-visit effectiveness -> effective service -> residual floral/reproductive response`.

The present paper does not fit either H5a or H5b. This distinction is used only to prevent genus absorption from being misread as evidence against pollination-mediated filtering.

---

## Results

### H1 — no universal floral island syndrome

The primary two-axis response did not recur as one coherent direction across context layers, evidence scopes and floristic strata. Supported responses depended on biogeographic context. Thus increasing geographic separation did not universally push island floras toward the same combination of generalized floral access and increased reproductive assurance.

This conclusion was based on multivariate within-context response tests and direct between-context comparisons rather than on comparing significance labels across regions.

### H2 — source separation produced different regional trajectories

The Palearctic showed the clearest broad primary response. In all-analysis all-native assemblages, accessibility/generalization increased with log-distance (estimate 0.0795, q=0.00463) and reproductive assurance also increased (0.0534, q=0.00463). Direct-only evidence retained both directions: accessibility/generalization 0.0616 (q=0.0389) and reproductive assurance 0.0966 (q=8.48×10^-14). The response persisted in native non-endemic assemblages.

Tropical assemblages followed a different trajectory. In direct-only all-native data, accessibility/generalization declined with distance (-0.1005, q=0.00275), while reproductive assurance increased (0.1360, q=0.00428). Thus the same direction of source separation was associated with a different combination of floral and reproductive responses.

Regional heterogeneity did not establish a climate-independent causal effect of realm identity. North–Tropical comparisons failed the outcome-blind common-support positivity gate, and Palearctic–Neotropical climate-adjusted support did not replicate consistently across both evidence scopes and focal floristic strata. H2 is therefore interpreted as biogeographically contingent, region-associated filtering with a measured-climate identification ceiling.

### Reproductive assurance and floral architecture decoupled

The Palearctic attraction/access shift remained positive after conditioning on `selfing_core` across all four source definitions. In all-native assemblages, conditional distance estimates ranged from 0.091 to 0.101 with q-values from approximately 0.0079 to 0.034. Thus the observed floral-architecture shift could not be reduced to the measured selfing syndrome alone.

Tropical responses provided the complementary pattern. Reproductive assurance could increase while specialized or attractive floral architecture was maintained or strengthened. Together these results reject a compulsory serial model of `isolation -> selfing -> floral simplification` and support partially parallel response branches whose directions depend on context.

### H3 — the broad Palearctic response was represented at genus level

Taxonomic-depth decomposition produced the strongest mechanistic constraint available from the plant data. In Wave36, the Palearctic primary two-axis response appeared to persist beyond both family and genus composition. By Wave52, additional reproductive evidence changed that inference. The response persisted beyond family composition but disappeared after source-matched genus adjustment.

The final snapshot reproduced the Wave52 classification exactly. Across all-analysis/direct-only evidence and all-native/native-nonendemic strata, the response was supported at the observed stage (4/4 source modes), retained after family adjustment (4/4) and absent after genus adjustment (0/4). The number of scored GIFT source species increased to 7,750 in all-analysis and 4,552 in direct-only evidence without changing this classification.

The broad Palearctic gradient is therefore compatible with differential representation of source-available genera rather than a robust beyond-genus primary residual.

### H3 does not identify the cause of lineage filtering

The genus-depth result closes the plant-side question of **where** the broad response is expressed but not the ecological question of **why** those genera are differentially represented. Differential genus representation can arise from propagule dispersal, establishment, habitat compatibility, demographic persistence, dependence on particular interaction channels or combinations of these processes.

Consequently, disappearance after genus adjustment should not be interpreted as evidence against a pollination mechanism. A pollination-channel deficit could act upstream by changing which source-available plant lineages can establish or persist. This route corresponds to prospective H5a. Conversely, any floral/reproductive response that remains after source/genus composition is fixed would be eligible for the stricter within-lineage H5b route.

### Pollination-associated floral architecture provided an intermediate compatibility layer

Palearctic assemblages moved away from butterfly-like and large-bee-like floral architecture with increasing isolation, whereas tropical assemblages increased along bird-like, butterfly-like and large-bee-like templates in both evidence scopes and focal floristic strata. Neotropical guild-labelled responses were not robustly supported.

One source-trained factor explained 86.94% of variation among the three templates in all-analysis evidence and 86.44% in direct-only evidence. Thus simultaneous tropical increases cannot be read as evidence that birds, butterflies and large bees all became more important. They primarily indicate movement along a shared specialized/attractive plant-architecture dimension.

These concordance scores therefore bridge H2/H3 and future H5 rather than identify H5. They show which functional architectures are compatible with the observed plant-side response, but they do not measure visitor mobility, channel retention, visitation, single-visit effectiveness or effective service.

### H4 — area remained a modifier, not an identified mechanism

All 16 frozen V3 primary classifications remained `retain_area_as_measurement_sensitive_modifier_only`, and none passed the heteroskedastic-null promotion gate. Directional distance-by-area patterns therefore did not justify a founder-filter, habitat-capacity or pollinator-persistence causal interpretation.

### Robustness and inference saturation

The broad Palearctic response survived the predeclared finite MNAR grid in both evidence scopes and focal floristic strata. Major conclusions also remained stable under direct-only evidence and the fixed source/taxonomic decomposition.

Wave52 contained 184,917 analysis-usable species-by-axis cells. The final snapshot contained 222,688, a gain of 37,771 usable cells, including 11,315 reproductive-assurance cells. Despite that increase, the principal H1–H4 interpretation and current H5 claim ceiling remained stable. New strict-source acquisition also entered strong diminishing returns: one high-yield Orchidaceae source contributed 192 reproductive cells, while subsequent reviewed source packets generally added only one to three cells; a residual audit of 1,187 staging files and more than 2.25 million rows found no immediately strict-ready reproductive records under the fixed evidence contract.

We therefore treat the final snapshot as a defensible computational freeze based on inference stability and recoverability rather than an arbitrary global coverage percentage.

---

## Discussion

### The floral island syndrome is an emergent assemblage pattern, not one rule

The central result is not simply the absence of a universal syndrome. Increasing source separation produced reproducible floral and reproductive responses, but those responses branched among biogeographic contexts. In the Palearctic, increasing distance was associated with both greater accessibility/generalization and greater reproductive assurance. In tropical assemblages, reproductive assurance could increase while specialized or attractive floral architecture was maintained or strengthened.

This pattern is difficult to reconcile with a compulsory sequence in which island isolation causes reproductive assurance, which then necessarily causes floral simplification. Instead, reproductive and floral components behaved as partially separable response axes. Their coupling is itself context dependent.

This reframing connects two literatures that are often discussed separately. Baker-type colonization filtering predicts that reproductive strategies can be sorted during establishment. Pollination-syndrome theory predicts that floral architecture reflects multivariate interaction environments. Our results suggest that island assemblages can combine these components differently rather than move along one shared axis.

### Source pools convert distance from a direct effect into an assembly process

The genus-depth result changes the biological interpretation of the Palearctic distance gradient. The broad response was real at the assemblage level but was absorbed by source-matched genus composition. This is consistent with a source-filter view of island biogeography: geographic separation changes which source-available lineages are represented, and those lineages carry different floral and reproductive strategies.

The result should not be translated into repeated within-lineage evolution without additional evidence. Instead, it demonstrates that an island-syndrome pattern can emerge from differential lineage representation. This is not a statistical nuisance. It is a biologically meaningful component of the island filter.

### H3 and H5 are complementary, not rival explanations

A key consequence of the double-filter framework is that H3 and H5 operate at different levels. H3 answers **where in the plant hierarchy the distance response is expressed**. H5 asks **why that plant hierarchy is filtered in that way**.

If pollination-channel disruption changes which source-available genera can establish or persist, then the pollination mechanism will be visible statistically as H3. Genus adjustment would remove the downstream trait pattern because the mechanism acted through genus representation. Therefore a `4/4 -> 4/4 -> 0/4` result cannot by itself distinguish plant-propagule filtering from pollinator-mediated lineage filtering.

This motivates H5a: channel-dependent lineage filtering. A future test should combine independently measured source-channel availability and island retention/disruption with lineage-level functional dependence, then ask whether channel state predicts source-matched genus entry or persistence. Formally, the relevant family is closer to `genus entry ~ channel state × functional dependence` than to a trait residual regression.

H5b remains the stricter residual route. If a floral/reproductive signal persists after source/genus composition is fixed, independently measured visitation, single-visit effectiveness and effective service can test a within-lineage or below-genus ecological response. A null H5b result would not falsify H5a.

This distinction matters because it prevents the analysis from demanding a beyond-genus residual before considering interaction-mediated assembly. It also preserves the frozen PR142 contract: the existing H5 residual gate is not retroactively altered, while H5a is explicitly prospective.

### Pollination-syndrome concordance is a bridge, not an endpoint

The regional template results are biologically informative but cannot identify visitors. Palearctic declines in large-bee-like and butterfly-like architectures are compatible with reduced representation of those plant-side architectures. Tropical increases across all three sampled templates indicate movement toward a shared specialized/attractive architecture, not simultaneous increases of three pollinator guilds.

The ~87% shared-factor result makes this distinction unavoidable. Floral syndromes are useful here as multivariate compatibility hypotheses: they identify which plant architectures covary with source separation and therefore which functional-channel explanations are plausible enough to test. They cannot establish channel mobility, retention, loss, replacement or effectiveness.

### Pollinator movement belongs upstream of H5

The double-filter model predicts that different pollination channels can have different isolation-susceptibility curves. For a functional channel `g`, continuity can be expressed conceptually as `E_poll(g,d) = f(source availability, distance, dispersal/flight ability, ocean crossing, establishment, habitat suitability, realized community)`. The same distance can therefore represent a strong barrier for one channel and a weak barrier for another.

This idea is compatible with the observed biogeographic branching but remains untested in Chapter 1. A region with a source-available, isolation-sensitive channel may experience strong pollinator-mediated lineage filtering as distance increases. A region where alternative functional channels remain available may retain specialized floral architecture even while reproductive assurance also increases.

Structural absence must remain distinct from loss. If a channel was absent from the source region, its absence on an island cannot be interpreted as a distance-driven disruption. This is why visitor identity and mobility must be measured independently rather than inferred from floral traits.

### Why the tropical branch is particularly informative

The tropical result is not merely a failure to reproduce the Palearctic pattern. It demonstrates that increased reproductive assurance can coexist with stronger specialized/attractive architecture. This is a direct counterexample to treating selfing and floral simplification as one obligatory syndrome.

The shared architecture factor also shows why simple replacement narratives are premature. The tropical increase of bird-like, butterfly-like and large-bee-like scores cannot identify which pollination channel was retained or replaced. What it does establish is that the plant assemblage does not converge toward a single generalized architecture as isolation increases.

### Area and climate define boundaries on interpretation

The area analyses did not satisfy the frozen promotion gates needed to label distance-by-area patterns as founder filtering, habitat capacity or pollinator persistence. Area should therefore remain a modifier rather than a named mechanism.

Likewise, the regional branching cannot be promoted to a climate-independent effect of biogeographic realm. The analysis identifies region-associated response trajectories under measured climate and source constraints, not an intrinsic causal property of a realm label.

These negative promotions are part of the result. The PR142 architecture is designed to preserve a difference between a reproducible pattern and an identified mechanism.

### The progressive waves provide a stopping rule, not the paper's novelty

The trait campaign materially changed what could be estimated and, at Wave52, changed the taxonomic-depth interpretation from apparently beyond-genus to genus-level assembly. That history demonstrates why the progressive design was necessary.

After Wave52, however, a further 37,771 analysis-usable cells did not change the broad claim structure. At the same time, source acquisition entered strong diminishing returns. This combination provides a defensible stopping rule: additional global trait collection is unlikely to change the main Chapter 1 inference enough to justify remaining the dominant research activity.

The point is not that 69.83% coverage is universally sufficient. It is that the inferential bottleneck has changed. The next unresolved questions require a different data type: independent pollination-channel retention, visitation, per-visit pollen function, effective service and reproductive dependence.

### From Chapter 1 to the mechanistic programme

Chapter 1 therefore closes three questions. First, it establishes that increasing source separation does not generate one universal floral/reproductive vector. Second, it identifies how reproductive assurance and floral architecture branch by context. Third, it shows that the strongest broad Palearctic response is represented at genus level.

It does not identify whether plant dispersal, habitat, demography or pollination-channel continuity caused that lineage sorting. The next mechanistic layer should therefore work both **through** and **below** H3: test whether channel disruption predicts lineage entry/persistence (H5a), and separately whether effective service predicts any beyond-genus floral/reproductive response (H5b).

This creates a direct bridge to the `izu-core` programme and prospective Izu field work, where visitor identity, observation effort, single-visit pollen deposition, rate-weighted effective service, reproductive dependency and autonomous assurance can be measured rather than inferred from floral phenotype.

---

## Conclusion

Increasing geographic separation from a regional source system does not impose one universal floral island syndrome. Instead, it reorganizes plant assemblages along context-dependent trajectories in which reproductive assurance and pollination-associated floral architecture can couple or decouple. The strongest broad Palearctic distance response is compatible with source-matched genus-level assembly, demonstrating that an island-syndrome pattern can emerge through differential lineage representation.

This genus-level result does not exclude a pollination mechanism. It identifies the plant-side level at which such a mechanism may operate. Pollination-channel disruption could act upstream by changing lineage establishment and persistence, or downstream by affecting floral/reproductive responses that remain after genus composition is fixed. Pollination-syndrome concordance provides a plant-side bridge to those hypotheses but cannot identify the channel itself.

The resulting framework is a double geographic filter: oceanic separation filters plant source-pool assembly and may simultaneously filter pollination-channel continuity. Chapter 1 directly resolves the plant-side WHEN, WHERE and taxonomic depth. The pollinator-side causal extension remains the next empirical problem.

---

## Proposed main figures

### Figure 1 — Double geographic filter and H1–H5 hierarchy

Show distance acting on both plant source-pool assembly and pollinator-channel continuity. Place H3 at lineage assembly. Split future H5 into H5a (channel-dependent lineage filtering) and H5b (beyond-genus effective-service response). Place pollination-syndrome concordance between observed plant response and H5 as an interpretation bridge.

### Figure 2 — Regional distance-response vectors

Accessibility/generalization and reproductive-assurance slopes for Palearctic and tropical contexts, with all-analysis/direct-only and native/non-endemic sensitivities.

### Figure 3 — Decomposition of the floral island syndrome

Plot selfing-core and attraction/access responses separately. Highlight Palearctic attraction shift persisting after selfing-core conditioning and tropical reproductive assurance coexisting with specialized architecture.

### Figure 4 — Taxonomic-depth result

Observed -> after family -> after genus for the Palearctic primary response, emphasizing `4/4 -> 4/4 -> 0/4` across the four scope × stratum combinations.

### Figure 5 — Pollination-associated architecture

Large-bee-like, butterfly-like and bird-like slopes by region plus the ~87% shared-factor decomposition. Explicitly label these as plant architecture, not pollinator identity.

### Figure 6 — Progressive-wave stability and claim ceiling

Wave36 -> Wave52 -> final trait coverage, taxonomic-depth transition, and stable final H1–H5 claim boundary. Keep diminishing returns secondary to the ecological figures.

---

## Reproducibility anchors

- Scientific contract: `config/chapter1_progressive_analysis.yml`
- Final trait input: Run `34191508045`, 222,688 / 318,885 resolved cells
- Final PR142 reanalysis: Run `34232450884`
- Artifact: `chapter1-progressive-analysis-34232450884`, ID `10058653212`
- Digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- H3/H5 hierarchy note: `docs/chapter1_h3_h5_causal_hierarchy_20260909.md`
