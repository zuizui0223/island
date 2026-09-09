# The island floral syndrome is not one syndrome: source-distance filtering, lineage assembly, and biogeographic decoupling of reproductive assurance and floral architecture

## Full working manuscript v7 — submission-order draft — 2026-09-09

### Abstract

Island floras are often expected to converge toward a coherent floral “island syndrome”, combining increased reproductive assurance with reduced floral specialization or attraction. Yet geographic isolation does not act directly on floral phenotype. An oceanic barrier changes which plant lineages from a regional source pool can disperse, arrive, establish and persist, and it can also change which regionally available pollination channels reach, establish and remain functionally available. Because plants and pollinators differ in dispersal and establishment ability, the same geographic separation need not produce one floral or reproductive outcome.

We tested whether increasing source separation generates one universal floral/reproductive response or instead reorganizes partially separable reproductive and floral components in a biogeographically contingent manner. We analysed a fixed universe of 8,265 islands and 106,295 accepted angiosperm species under a progressive analysis contract in which hypotheses, support gates, source-pool safeguards, model order and claim ceilings were fixed while trait evidence improved. The final snapshot resolved 222,688 of 318,885 possible species-by-axis cells (69.83%) across flower colour, floral structural complexity and reproductive assurance. Mainland distance was treated as a composite gradient of separation, connectivity and source accessibility rather than as a mechanistically pure causal treatment.

A universal floral/reproductive island syndrome was not recovered. In the Palearctic, increasing separation was associated with greater accessibility/generalization and greater reproductive assurance. In all-analysis all-native assemblages, slopes were 0.0795 and 0.0534 per unit of log-distance, respectively (both FDR q=0.00463); direct-only estimates were 0.0616 (q=0.0389) and 0.0966 (q=8.48×10^-14). The pattern persisted in native non-endemics. Tropical assemblages followed a different trajectory: in direct-only all-native data, accessibility/generalization declined with distance (-0.1005, q=0.00275) while reproductive assurance increased (0.1360, q=0.00428). Thus reproductive assurance and floral architecture did not behave as one obligatory serial syndrome.

Source-matched taxonomic decomposition showed that the broad Palearctic response survived family adjustment but disappeared after genus adjustment in all four evidence-scope-by-floristic-stratum combinations (4/4 -> 4/4 -> 0/4). The Palearctic attraction/access shift nevertheless remained positive after conditioning on selfing core across four source definitions (distance estimates 0.091–0.101 in all natives; q=0.0079–0.034). Secondary large-bee-like, butterfly-like and bird-like floral templates diverged among regions, but a source-trained common factor explained 86.94% and 86.44% of template variance in all-analysis and direct-only evidence, respectively, indicating that these labels mainly describe overlapping plant architecture rather than realized pollinator identity.

We therefore interpret the island floral syndrome as an emergent assemblage-level outcome with at least two partially separable components: reproductive assurance and pollination-associated floral architecture. The strongest Palearctic distance response is expressed primarily through source-matched genus assembly. H3 identifies the plant-assemblage level at which this filter is represented; it does not identify why those genera are differentially represented. Pollination-channel retention or loss remains a candidate upstream mechanism that could itself generate lineage filtering, while independent visitation and effective-service data are required to test any mechanism that persists beyond genus. Large additions of trait evidence no longer changed this claim structure, providing an inference-based stopping rule for the computational campaign.

**Keywords:** island biogeography; floral traits; reproductive assurance; selfing syndrome; pollination syndrome; source pool; dispersal filtering; lineage assembly; mutualism; biogeographic contingency

---

## Introduction

Geographic isolation is one of the defining processes of island biogeography. Classical island theory treats separation from a regional source pool as a constraint on colonization opportunity, while empirical island biology has repeatedly identified characteristic changes in dispersal, life history, growth form and reproduction. In plants, these observations have motivated the idea of an “island syndrome”: a recurrent suite of traits associated with long-distance colonization, small population size, novel abiotic environments and altered biotic interactions.

Floral versions of the island syndrome commonly combine at least two biological arguments. The first is reproductive assurance. Long-distance colonists may have difficulty finding compatible mates or reliable pollen vectors, favouring self-compatibility, autonomous selfing or other forms of uniparental reproduction. This logic is closely related to Baker’s law and to later treatments of breeding system as a colonization filter rather than necessarily a repeated post-colonization evolutionary transition (Baker 1955; Pannell & Barrett 1998; Grossenbacher et al. 2017). The second argument concerns pollination-associated floral architecture. If island pollinator assemblages are species-poor, functionally altered or dominated by different visitor groups, the returns to conspicuous displays, restricted floral access, deep tubes or precise mechanical matching may change (Fenster et al. 2004; Rosas-Guerrero et al. 2014; Hetherington-Rauth & Johnson 2020).

These arguments are related but not equivalent. The selfing syndrome often includes reduced floral investment, yet self-compatibility is not autonomous selfing, autonomous selfing is not necessarily the dominant realized mating system, and a plant can gain reproductive assurance without losing specialized floral architecture (Sicard & Lenhard 2011). Conversely, floral attraction or accessibility traits can change without a detectable shift in reproductive assurance. Treating all of these outcomes as one syndrome risks converting a broad assemblage pattern into a presumed causal sequence.

A second difficulty is that island floral patterns are assemblage patterns. A distant island may differ from a near island because different source-available lineages arrive, establish or persist, even if no lineage changes phenotype after colonization. This is the ecological logic behind interpreting Baker’s law partly as filtering of pre-existing reproductive variation. Recent work has similarly shown that breeding system, floral symmetry and arrival opportunity jointly influence island colonization, and explicit source-pool frameworks have emphasized the need to distinguish trait filtering from generic island–mainland differences (Zell et al. 2025; Schrader et al. 2024). Thus, lineage composition should not automatically be treated as a nuisance to be removed: it can be a biological outcome of island assembly.

This source-pool perspective also changes the interpretation of geographic distance. We do not treat kilometres from a continent as a direct causal force on floral phenotype. Instead, distance is a composite axis of separation, connectivity and source accessibility. Conceptually, as distance approaches zero, an island approaches a high-accessibility boundary where geographic filtering relative to the regional source system should be weaker. As separation increases, the opportunity for differential source-to-island filtering increases. This limiting argument does not require a literal zero-distance island to equal mainland vegetation and does not imply that the nearest continent is the exact historical source of every island.

The same geographic barrier may also act on the other side of plant–pollinator interactions. A pollination channel that is available in a regional source system must itself cross oceanic barriers, establish and persist before it can provide effective service on an island. Flight capacity, long-distance dispersal, stepping-stone use, habitat requirements and establishment probability may therefore make a given ocean distance a very different biological barrier for different pollinator functional groups. Structural absence of a channel from the source region is also fundamentally different from loss of a channel that was initially available. These considerations motivate a double geographic-filter framework in which isolation jointly filters plant source-pool assembly and pollination-channel continuity. The present Chapter 1 directly evaluates the plant side; the pollinator side remains a mechanistic extension rather than an inference from floral phenotype.

A third challenge is analytical. Global island comparisons combine heterogeneous floristic histories, climates, island sizes, observation effort and taxonomic composition. A claim that one region is significant and another is not does not demonstrate that regional response vectors differ. Missing trait evidence is also geographically and taxonomically structured, and new trait evidence can change which biological questions are estimable. Consequently, both the biological hypotheses and the rules that determine when a claim is promoted must be fixed before outcome inspection.

We addressed these problems with a progressive analysis contract established in PR #142. The island universe, trait ontology, hypothesis order, support thresholds, multiple-testing families, source/lineage safeguards and causal claim ceiling were fixed, while trait evidence was allowed to improve. Each new trait snapshot was reanalysed from the beginning under the same contract. This created a sequential hypothesis ladder in which each result generated the next biological question without allowing the question itself to drift with the data.

**H1** treated a universal floral island syndrome as a rival hypothesis: does increasing source separation generate one coherent floral/reproductive response vector across contexts? **H2** asked whether response vectors instead branch among biogeographic contexts. **H3** asked at what source/taxonomic depth a supported response remains, distinguishing broad assemblage patterns from family- and genus-level representation. **H4** asked whether continuous island area modifies isolation-associated filtering. **H5**, which is deliberately not claimed from the present plant database, asks why the identified plant filter arises and whether independently measured pollination-channel change explains either lineage assembly itself or any response that persists beyond genus.

The primary H2 response is deliberately pollinator-name-free. We separate (i) **accessibility/generalization**, representing increasingly open, shallow or radially accessible floral architecture, from (ii) **reproductive assurance**, representing strict reproductive-assurance evidence. Only after the primary response is established do we evaluate fixed large-bee-like, butterfly-like and bird-like floral templates. These are sampled pollination-associated architecture concordances, not visitor observations or pollinator classifiers.

Our central prediction is therefore not that distant islands should be uniformly selfing, generalized or inconspicuous. Instead, increasing source separation should reorganize floral and reproductive assemblages, but the signs, coupling and taxonomic depth of those responses should depend on biogeographic context. If the floral island syndrome is an emergent assemblage outcome rather than one deterministic trait rule, reproductive assurance and floral architecture should sometimes decouple, and broad geographic gradients should often resolve into lineage assembly rather than repeated within-lineage change.

---

## Materials and Methods

### Geographic and taxonomic universe

The fixed geographic universe contained 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. The denominator was not allowed to shrink between trait waves. Island occurrences, floristic-status information, island area, climate covariates, geographic exposure and source-pool infrastructure were held constant while trait evidence was allowed to improve.

The primary geographic exposure was `log1p_distance_to_continent_km`. We interpret this as a composite of geographic separation, connectivity and source accessibility rather than as a mechanistically pure isolation treatment. The conceptual near-source limit (`d -> 0`) represents high source accessibility and relatively weak oceanic filtering, not a literal mainland observation or a claim that a zero-distance island should have mainland composition.

Analyses were conducted at two context layers. Broad analysis regimes included northern mid-latitude, northern high-latitude, tropical and southern extratropical strata. Formal biogeographic realms included Palearctic, Neotropical and other realm assignments where support permitted. We analysed `all_native`, `native_nonendemic` and `endemic` floristic strata where support thresholds were met. Persistence in native non-endemics was used to show that a pattern was not confined to island endemics; it was not interpreted as a causal endemicity contrast.

### Trait evidence and progressive wave integration

Each accepted species could contribute to three primary axes: flower colour, floral structural complexity and reproductive assurance. The fixed denominator was therefore 318,885 species-by-axis cells. The final trait snapshot resolved 82,556 colour cells, 91,635 structural cells and 48,497 reproductive-assurance cells, for 222,688 resolved cells overall (69.83%).

Evidence precedence was fixed. Species-direct High/Medium evidence preceded trait-specific validated Low evidence. Family-level imputation, global fallback and post hoc threshold relaxation were prohibited. Missing evidence was never converted to trait absence. Reproductive evidence retained distinctions among self-incompatibility, mating system and autonomous selfing when source methodology permitted; self-compatibility alone was not mapped to autonomous or realized selfing.

All trait waves were materialized into the same `chapter1_trait_snapshot_v1` schema and reanalysed under `chapter1_progressive_analysis_v1`. Biological estimates were allowed to change as evidence improved, but model definitions, support gates and claim ceilings were not.

### Support gates and multivariate priority

Global fill fraction was descriptive and was never an inferential gate. Response-specific support was classified as <30 supported islands = not promoted, 30–49 = pilot, and >=50 = count component of confirmatory support. Pairwise context comparisons required the same response to meet the same support tier in both contexts. BH-FDR correction was applied within predeclared test families.

Multivariate response-vector tests preceded interpretation of individual axes. A significant slope in one context and a non-significant slope in another was not treated as evidence of heterogeneity without a direct between-context comparison.

### H1 — universal-syndrome rival

H1 asked whether increasing source separation generated one coherent floral/reproductive response vector across contexts. The primary test combined within-context multivariate support with direct between-context vector comparison. A universal syndrome was therefore rejected only by evidence of response-vector heterogeneity, not by comparing significance labels across regions.

### H2 — biogeographic branching

H2 asked whether the primary plant response branched among prespecified biogeographic contexts. The pollinator-name-free response contained `accessibility_generalization` and `reproductive_assurance`. Analyses were repeated across evidence scopes and floristic strata. V1 climate-overlap analyses separately evaluated whether regional contrasts could be transported over common measured-climate support.

### Decomposing reproductive and floral components

To distinguish reproductive assurance from pollination-associated floral change, `selfing_core` included only self-incompatibility, mating-system and autonomous-selfing information and excluded floral attraction traits. A separate attraction/access pathway was constructed from floral architecture. Conditional models asked whether the isolation-associated attraction shift remained after conditioning on `selfing_core`; these analyses were treated as decomposition rather than causal mediation.

### H3 — source and lineage assembly

H3 asked at what taxonomic depth a supported H2 response remained. Outcome-blind source expectations were constructed using GIFT source flora and fixed source assignments. Taxonomic decomposition evaluated the response at three stages: observed assemblage response, residual after family composition, and residual after source-matched genus composition. Genus entry and within-represented-genus species loading were evaluated where support permitted.

Family and genus were grouping structures rather than trait-imputation devices. Disappearance after genus adjustment was interpreted as evidence that the broad response is represented by genus-level assembly; it was not interpreted as proof that plant dispersal alone caused the filter.

### H4 — continuous area moderation

H4 used continuous distance-by-continuous-area models. No post hoc small/large island threshold was introduced. Equal-island, capped-information, direct-only, common-support and heteroskedastic-null safeguards were used to decide whether a directional interaction could be promoted beyond a measurement-sensitive modifier.

### Pollination-associated floral-architecture concordance

Only after the H2 plant response was established were fixed `large_bee_like`, `butterfly_like` and `bird_like` templates evaluated. These combined colour, form, symmetry, tube depth and size with prespecified weights. They were treated as floral-architecture concordances rather than pollinator observations. V4 source-trained factor decomposition quantified how much of the three-template covariance was shared plant architecture, and template-specific residuals were interpreted only after this common factor and source/taxonomic expectations were considered.

### H5 — mechanistic extension and claim boundary

The frozen PR142 H5 is a strict post-H2–H4 residual gate: independently measured pollination-channel deficit must add explanatory information beyond source, lineage, area, climate, observation and spatial safeguards. In the double-filter interpretation, however, an additional prospective route is biologically possible. Pollination-channel retention or loss may affect which plant genera establish or persist, in which case the channel mechanism generates H3 itself and is absorbed by genus composition.

We therefore distinguish two causal routes in interpretation while leaving the fitted PR142 analysis unchanged. **H5a** is channel-dependent lineage filtering: source channel -> retention/disruption -> effective service -> differential plant establishment/persistence -> genus entry/representation. **H5b** is the original residual route: source/genus composition fixed -> retention/disruption -> visitation -> single-visit effectiveness -> effective service -> beyond-genus floral/reproductive response. Neither route is claimed from the Chapter 1 plant database.

### Robustness and missingness

V5 used a finite predeclared MNAR grid for trait-resolution sensitivity rather than assuming missingness at random. V1 evaluated common measured-climate support. V2 evaluated taxonomic depth. V3 applied area-support falsification. V4 separated shared versus template-specific floral architecture. The workflow also retained direct-only evidence, information-weight sensitivities and other predeclared support checks.

### Computational stopping rule

The final computational freeze was not based on reaching a target coverage percentage. We required (i) a large increase in usable evidence without a corresponding change in the principal H1–H5 claim ceiling, (ii) declining realized yield from newly reviewed sources, and (iii) exhaustion of immediately promotable committed reproductive evidence under the strict evidence rules.

---

## Results

### H1 — one universal floral island syndrome was not recovered

The primary two-axis response did not recur as one coherent direction across context layers, evidence scopes and floristic strata. Response vectors differed among contexts, rejecting a simple model in which increasing source separation universally pushes island floras toward the same combination of greater reproductive assurance and greater floral generalization.

This conclusion was based on the predeclared multivariate architecture rather than on comparing significance labels among regions. The principal result was therefore not that one trait “worked” in one region and failed in another, but that the multivariate response to source separation was context dependent.

### H2 — source separation produced biogeographically contingent response branches

The Palearctic showed the clearest broad primary response. In all-analysis all-native assemblages, the accessibility/generalization slope was 0.0795 and the reproductive-assurance slope was 0.0534 per unit of log-distance; both were FDR-supported at q=0.00463. In direct-only evidence, the corresponding estimates were 0.0616 (q=0.0389) and 0.0966 (q=8.48×10^-14). The broad response persisted in native non-endemic assemblages, indicating that it was not confined to island endemics.

Tropical assemblages followed a contrasting trajectory. In direct-only all-native data, accessibility/generalization decreased with distance (-0.1005, q=0.00275) while reproductive assurance increased (0.1360, q=0.00428). Thus increasing isolation was compatible with stronger reproductive assurance without the floral-generalization response observed in the Palearctic.

These patterns establish region-associated branching, not a measured-climate-independent causal effect of realm identity. North–Tropical contrasts did not satisfy the outcome-blind common-support positivity requirement, and Palearctic–Neotropical climate-adjusted support did not reproduce consistently across both evidence scopes and both focal floristic strata.

### Reproductive assurance and floral architecture were partially separable

The Palearctic attraction/access shift remained positive after conditioning on `selfing_core` across four source definitions. In all-native assemblages, conditional distance estimates ranged from approximately 0.091 to 0.101, with q values from 0.0079 to 0.034. Thus the floral-architecture response was not reducible to the measured reproductive-assurance component alone.

Tropical assemblages showed the complementary pattern: reproductive assurance could increase while specialized or attractive floral architecture was maintained or strengthened. Together these results reject a compulsory serial model of `isolation -> selfing -> floral simplification` and instead support partially parallel response components whose coupling depends on context.

### H3 — the broad Palearctic response was represented by genus-level assembly

Taxonomic-depth decomposition produced the strongest plant-side mechanistic result. Across all-analysis/direct-only evidence and all-native/native-nonendemic strata, the broad Palearctic response was supported at the observed stage in all four source modes, remained supported after family composition, and disappeared after source-matched genus composition: `4/4 -> 4/4 -> 0/4`.

The same classification persisted as the number of scored GIFT source species increased to 7,750 in all-analysis evidence and 4,552 in direct-only evidence. Thus the strongest broad source-distance response is compatible with differential representation of source-available genera rather than a robust beyond-genus primary residual.

This result identifies the plant-assemblage depth of the filter but does not identify its ultimate cause. Plant propagule dispersal, establishment, habitat filtering, demographic persistence and interaction dependence remain possible contributors to differential genus representation.

### Pollination-associated floral concordance linked H2/H3 patterns to candidate functional channels

Palearctic assemblages moved away from butterfly-like and large-bee-like floral architecture with increasing separation, whereas tropical assemblages increased along bird-like, butterfly-like and large-bee-like templates in both evidence scopes and both focal floristic strata. However, one source-trained factor explained 86.94% of variation among the three templates in all-analysis evidence and 86.44% in direct-only evidence.

The simultaneous tropical increase therefore cannot be interpreted as evidence that three realized pollinator guilds all became more important. Instead, the templates largely track a shared specialized/attractive plant-architecture dimension. Some template-specific tropical residual structure remained after source/genus adjustment, including positive large-bee-like residuals in robust comparisons, but these remain plant-side architecture signals rather than measurements of visitor identity, mobility, effectiveness or replacement.

### H4 — area remained a modifier rather than an identified mechanism

All 16 frozen primary area classifications remained `retain_area_as_measurement_sensitive_modifier_only`, and none passed the heteroskedastic-null promotion gate. Directional distance-by-area patterns therefore did not justify claims of founder filtering, habitat capacity, demographic persistence or pollinator persistence.

### H5 remained beyond the Chapter 1 identification ceiling

The Chapter 1 plant database contained no complete independent source-channel -> retention/disruption -> visitation -> single-visit effectiveness -> effective-service chain. Consequently, no floral architecture or genus-assembly result was promoted to historical pollinator loss, functional replacement or effective-service causation.

Importantly, the H3 genus result does not falsify a future H5 mechanism. A pollination-channel filter could act upstream by changing which source-available genera establish or persist, producing an H3 signal that is then absorbed by genus composition. Conversely, any response that remains beyond genus can be tested through the original strict H5 residual route. These possibilities are mechanistic interpretations generated by the Chapter 1 result, not additional fitted Chapter 1 effects.

### Progressive trait waves reached practical inference saturation

Wave52 contained 184,917 analysis-usable species-by-axis cells (57.99%). The final integrated snapshot contained 222,688 cells (69.83%), a gain of 37,771 usable cells, including 11,315 additional reproductive-assurance cells. Despite this large increase, the current H1–H5 claim structure remained stable.

Late source acquisition simultaneously entered strong diminishing returns. One high-yield Orchidaceae source contributed 192 strict reproductive cells, whereas later reviewed source packets generally added one to three cells. A committed residual audit scanned 1,187 staging files and more than 2.25 million rows and found no immediately strict-ready reproductive records under the fixed evidence rules. We therefore treat the final snapshot as a defensible computational freeze based on inference stability and declining recoverability rather than on declaring 69.83% universally sufficient.

---

## Discussion

### The floral island syndrome is not one syndrome

Our results do not support a single universal floral/reproductive trajectory with increasing island isolation. Instead, source separation is associated with context-dependent assemblage responses in which reproductive assurance and pollination-associated floral architecture can couple or decouple. The strongest contrast is not simply a difference in effect size. In the Palearctic, accessibility/generalization and reproductive assurance increase together with separation. In tropical assemblages, reproductive assurance can increase while specialized or attractive architecture is maintained or strengthened.

This distinction matters because reproductive assurance and floral simplification are often bundled into one island or selfing syndrome. The Palearctic conditional analyses show that the floral component cannot be reduced to the measured selfing core alone, while the tropical result demonstrates that increased reproductive assurance need not imply floral generalization. The data therefore favour a parallel-component interpretation over a compulsory serial pathway.

### Source and lineage assembly are ecological results, not nuisance structure

The strongest plant-side result is the transition from an apparently broad response to a genus-level assembly signal. The Palearctic response remains after family composition but disappears after source-matched genus composition across both evidence scopes and both focal floristic strata. This pattern is consistent with a source-pool view of island assembly in which geographic separation changes which lineages are represented rather than repeatedly transforming the same lineages after colonization.

This interpretation complements previous evidence that self-compatibility and other traits can influence island colonization (Grossenbacher et al. 2017; Zell et al. 2025) and explicit calls to compare island trait distributions against source-pool expectations (Schrader et al. 2024). The contribution here is to show that a multivariate floral/reproductive distance response can itself resolve to a particular taxonomic depth, and that this depth is stable to substantial additional trait evidence.

The genus result also places an important limit on evolutionary interpretation. We do not infer repeated within-lineage evolution from the cross-sectional assemblage. At the same time, genus composition should not be read as proving that plant propagule dispersal alone generated the pattern. Lineage representation integrates arrival, establishment, persistence, habitat and interaction dependence.

### H3 and H5 are different levels of the same double-filter problem

A purely sequential reading of H3 followed by H5 could imply that pollination mechanisms matter only if a floral response remains after genus composition is removed. Biologically, this is too narrow. Pollination-channel retention or loss can act upstream of lineage assembly. If a source-available lineage depends strongly on a pollination channel that fails to reach or persist on remote islands, that lineage may fail to establish or may disappear. The resulting pollination effect would then be expressed as H3 genus composition rather than as a beyond-genus residual.

We therefore distinguish two future mechanistic routes. H5a asks whether independently measured pollination-channel geography predicts plant lineage entry or persistence conditional on lineage dependency. H5b retains the original strict PR142 residual logic and asks whether visitation and effective service explain floral/reproductive response that remains after source/genus composition is fixed. A null H5b result would not falsify H5a.

This hierarchy also clarifies the role of pollination-syndrome concordance. The large-bee-like, butterfly-like and bird-like scores do not identify realized pollinators. They provide an intermediate compatibility layer between the plant-side response and candidate functional channels. Their strong shared factor is itself a warning against assigning visitor identity from floral phenotype.

### The same oceanic barrier can filter both sides of a mutualism

The double geographic-filter interpretation follows naturally from the combination of H2 and H3. Distance constrains access to plant source pools, and the strongest broad Palearctic response is expressed through source-matched genus assembly. The same distance can also constrain pollination channels, but the relevant susceptibility curve is likely to differ among functional groups because flight, long-distance dispersal, ocean crossing, habitat requirements and establishment differ.

Conceptually, for a plant lineage `z`, source-to-island continuity can be written as `E_plant(z,d) = f(source availability, distance, propagule dispersal, establishment, persistence)`. For a pollination channel `g`, continuity can be written as `E_poll(g,d) = f(source availability, distance, flight/dispersal ability, ocean crossing, establishment, habitat, realized community)`. These are conceptual functions rather than fitted Chapter 1 models, but they explain why equal geographic distances need not be biologically equivalent across lineages, pollination channels or biogeographic contexts.

This framework also prevents a common categorical error. Absence of a pollinator channel from a tropical source system is not the same event as loss of that channel from an island whose source system contains it. Structural absence must therefore be distinguished from island-driven disruption before any mechanistic interpretation is promoted.

### Regional contrast is biologically informative even with a climate ceiling

The failure to establish measured-climate-independent realm causation does not erase H2. It narrows its interpretation. The plant response is clearly context dependent, but the current global data do not identify whether categorical biogeographic history, measured climate, source composition or unmeasured environmental structure is ultimately responsible for each branch.

This distinction is important because the tropical result is not merely a weaker version of the Palearctic syndrome. It can point in a different floral direction while retaining increased reproductive assurance. This qualitative difference motivates mechanistic comparisons of source pools and interaction channels rather than a search for one universal coefficient.

### Area defines a boundary condition, not a closed causal mechanism

The area analyses show how the progressive contract prevents directional patterns from becoming over-interpreted mechanisms. Although some distance-by-area estimates were compatible with stronger filtering on small islands, the frozen V3 safeguards did not support promotion to a founder, habitat-capacity or pollinator-persistence mechanism. Area therefore remains a modifier whose biological basis must be resolved with data that directly distinguish these alternatives.

### Progressive waves provide an inference-based stopping rule

The trait-acquisition campaign is methodologically useful because it shows when additional database work stops being the main scientific bottleneck. The final snapshot added 37,771 usable species-by-axis cells relative to Wave52, including more than 11,000 reproductive-assurance cells, without changing the principal H1–H5 claim structure. At the same time, new strict-source acquisition fell from one relatively large reproductive source to isolated one-to-three-cell gains, while a large residual audit found no immediately promotable reproductive records.

The stopping rule is therefore not that 69.83% coverage is intrinsically sufficient. It is that substantial evidence growth no longer changes the major inference, while additional recovery becomes increasingly source-specific and costly. Remaining uncertainty is now dominated by identification of mechanism rather than indiscriminate trait filling.

### From global WHEN/WHERE to mechanistic HOW/WHY

Chapter 1 establishes when and where source separation is associated with floral/reproductive assemblage response, demonstrates that response components can decouple, and identifies genus-level assembly as the dominant plant-side depth of the strongest Palearctic pattern. It does not identify the pollinator-side causal chain.

The next mechanistic layer must therefore measure pollination channels independently of floral phenotype: source availability, retention or disruption, realized visitation, per-visit pollen transfer and rate-weighted effective service. Such data can test whether pollinator geography generates lineage filtering itself (H5a) or explains any response remaining beyond genus (H5b). This is the point at which the global pattern becomes a mechanistic question rather than another database-completion problem.

---

## Conclusion

Increasing geographic separation does not impose one universal floral island syndrome. Instead, island floras follow biogeographically contingent assemblage trajectories in which reproductive assurance and pollination-associated floral architecture can be partly decoupled. The strongest Palearctic response combines greater accessibility/generalization with greater reproductive assurance, but its broad multivariate signal is expressed primarily through source-matched genus assembly rather than a robust beyond-genus residual. Tropical assemblages demonstrate that increased reproductive assurance can coexist with maintained or strengthened specialized floral architecture.

These results support a source-pool view of floral island assembly in which distance filters lineage representation rather than acting directly on floral phenotype. They also motivate a broader double geographic-filter hypothesis: the same oceanic barrier can filter both plant colonists and their pollination channels, with different lineages and functional groups experiencing different isolation susceptibilities. Chapter 1 directly closes the plant-side WHEN/WHERE problem and identifies the taxonomic depth of the strongest filter. The remaining WHY question requires independent pollination-channel and effective-service evidence rather than further indiscriminate trait acquisition.

---

## Proposed main figures

**Figure 1. Hypothesis ladder and double geographic filter.** H1 -> H2 -> H3/H4 -> pollination-associated concordance -> H5a/H5b, with plant-side evaluated paths solid and pollinator-side mechanistic extensions dashed.

**Figure 2. Progressive wave convergence and stopping rule.** Wave36, Wave52 and final usable coverage by trait axis; source-reported versus materialized evidence; late source yield; H1–H5 decision stability.

**Figure 3. H2 regional source-distance response vectors.** Accessibility/generalization and reproductive-assurance slopes by context, evidence scope and floristic stratum.

**Figure 4. H3 taxonomic-depth decomposition.** Palearctic observed -> after family -> after source-matched genus, emphasizing `4/4 -> 4/4 -> 0/4`.

**Figure 5. Reproductive versus pollination-associated floral components.** Selfing-core conditional attraction shift, regional guild-template concordance and the approximately 87% shared architecture factor.

**Figure 6. Identification ceiling and mechanistic handoff.** Plant-side H3 result -> H5a lineage-filter mechanism / H5b beyond-genus effective-service mechanism -> prospective island field tests.

---

## Working references

Baker, H. G. 1955. Self-compatibility and establishment after long-distance dispersal. *Evolution* 9:347–349.

Fenster, C. B. et al. 2004. Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics* 35:375–403.

Grossenbacher, D. L. et al. 2017. Self-compatibility is over-represented on islands. *New Phytologist* 215:469–478.

Hetherington-Rauth, M. C. & Johnson, M. T. J. 2020. Floral trait evolution of angiosperms on Pacific islands. *The American Naturalist* 196.

Pannell, J. R. & Barrett, S. C. H. 1998. Baker’s law revisited: reproductive assurance in a metapopulation. *Evolution* 52:657–668.

Rosas-Guerrero, V. et al. 2014. A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters* 17:388–400.

Schrader, J. et al. 2024. Trait filtering in island floras: A conceptual framework. *Journal of Biogeography* 51:1596–1606.

Sicard, A. & Lenhard, M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany* 107:1433–1443.

Zell, A. N. et al. 2025. Island colonization in flowering plants is determined by the interplay of breeding system, lifespan, floral symmetry, and arrival opportunity. *New Phytologist* 245:420–432.

---

## Reproducibility anchors

- frozen scientific contract: `chapter1_progressive_analysis_v1`;
- final trait integration: Run `34191508045`;
- final PR142 progressive reanalysis: Run `34232450884`;
- artifact: `chapter1-progressive-analysis-34232450884`;
- artifact ID: `10058653212`;
- artifact digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
