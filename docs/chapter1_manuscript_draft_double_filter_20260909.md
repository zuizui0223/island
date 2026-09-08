# The island floral syndrome is not one syndrome: a double geographic filter on plant assembly and pollination channels

## Working manuscript draft — 2026-09-09

### Abstract

Geographic isolation is often treated as if it directly imposes a single floral “island syndrome”, including reproductive assurance, reduced floral specialization, and less conspicuous flowers. Yet an oceanic barrier acts simultaneously on two biological systems. It filters which plant lineages from a regional source pool can disperse, arrive, establish and persist, and it also filters which pollination channels can cross the same barrier, establish, persist and continue to provide effective service. Because plant lineages and pollinator functional groups differ in dispersal and establishment ability, the same geographic distance need not generate the same floral or reproductive outcome.

We tested this problem with a progressive global analysis in which the island universe, trait ontology, model order, support thresholds, multiplicity rules, source/lineage safeguards, and causal claim ceiling were fixed while species-level trait evidence improved across acquisition waves. The final snapshot contained 106,295 accepted angiosperm species scored on three primary axes—flower colour, floral structural complexity, and reproductive assurance—yielding 222,688 resolved species-by-axis cells out of 318,885 (69.83%). The primary exposure was a log-transformed mainland-distance gradient interpreted as a composite of geographic separation, connectivity, and source accessibility rather than as a mechanistically pure causal variable.

A universal floral/reproductive island syndrome was not recovered. Instead, isolation-associated response vectors differed among biogeographic contexts. The clearest broad response occurred in the Palearctic, where increasing source separation was associated with greater floral accessibility/generalization and greater reproductive assurance across all-analysis and direct-only evidence and in native non-endemic assemblages. Tropical assemblages could instead show increased reproductive assurance while specialized or attractive floral architecture was maintained or strengthened. Thus reproductive assurance and pollination-associated floral architecture behaved as partially separable response components rather than one obligatory serial syndrome.

Source-matched taxonomic decomposition showed that the primary Palearctic response survived family adjustment but not genus adjustment in all four evidence-scope-by-stratum combinations (`4/4 -> 4/4 -> 0/4`), making genus-level lineage assembly the strongest supported taxonomic depth. Secondary large-bee-like, butterfly-like, and bird-like floral templates showed strong regional contrasts, but one source-trained factor explained 86.94% and 86.44% of template variation in all-analysis and direct-only data, respectively. These templates therefore captured overlapping floral architecture rather than realized pollinator identity. Measured-climate-independent categorical biogeography and an area/capacity mechanism were not established, and no direct pollination-channel exposure, disruption, visitation or effective-service data entered Chapter 1.

The progressive wave series also provided an inferential stopping rule. Wave52 contained 184,917 analysis-usable cells; the final snapshot contained 222,688, including 11,315 additional reproductive-assurance cells, yet the major H1–H5 claim ceiling remained stable. New source acquisition simultaneously entered strong diminishing returns. We therefore treat the final snapshot as a defensible computational freeze based on inference stability and recoverability, not on an arbitrary coverage percentage.

Together, the results support a double geographic-filter view of island floral assembly. Increasing separation from a regional source system filters both plant source-pool representation and pollination-channel persistence. The plant-side filter is directly supported by the genus-depth result; the pollinator-side filter remains a mechanistic hypothesis requiring independent data. Chapter 1 therefore establishes WHEN and WHERE assemblage-level responses branch, while motivating a mechanistic handoff to interaction models and prospective island field tests.

---

## Introduction

Geographic isolation is one of the defining processes of island biogeography. Classical theory emphasizes the decline in colonization opportunity with increasing separation from a regional species pool, but plant island syndromes are often discussed as if isolation itself directly selects a characteristic phenotype. Floral versions of the syndrome commonly combine two expectations. First, mate and pollinator limitation may favour reproductive assurance, including self-compatibility or autonomous selfing. Second, altered pollination environments may change the returns to floral attraction and specialization, potentially favouring more accessible, less specialized, or less conspicuous flowers.

These expectations are biologically plausible, but they are not the same mechanism. A plant can possess reproductive assurance while retaining specialized floral architecture. Floral access or attraction traits can shift without a corresponding change in the measured selfing core. More fundamentally, island-level trait composition can change because different lineages reach and establish, not because the same lineages repeatedly evolve the same phenotype after colonization.

The geographic barrier itself also has two sides. The same oceanic separation filters plants and their interaction partners. On the plant side, source-pool membership is filtered by propagule dispersal, arrival probability, establishment, and persistence. On the pollinator side, a regionally available functional channel is filtered by flight or dispersal capacity, ocean crossing, establishment, habitat suitability, and realized community context. A kilometre of ocean is therefore not necessarily an equivalent biological barrier for all plant lineages or pollination channels.

This motivates a double geographic-filter framework:

```text
regional plant source pool                     regional pollinator pool
          |                                             |
          | increasing oceanic separation              |
          v                                             v
plant source-to-island filter                  pollination-channel filter
(dispersal, arrival, establishment)             (flight/dispersal, establishment)
          |                                             |
          +--------------------+------------------------+
                               v
                    realized island interaction context
                               |
                  +------------+------------+
                  v                         v
        reproductive assurance      floral architecture
                               |
                               v
                    observed island assemblage
```

Under this framework, no universal island floral syndrome is expected unless plant and pollinator filters repeatedly align in the same direction across regions. Instead, the syndrome can emerge as an assemblage-level property whose components couple or decouple depending on source composition and interaction context.

We test this logic using a progressive global analysis. Rather than freezing conclusions at one incomplete trait snapshot, we fixed the scientific questions and repeatedly re-estimated the same analysis as evidence improved. The resulting PR142 contract is a sequential hypothesis ladder. H1 treats a universal floral island syndrome as a rival. H2 asks where isolation-associated multivariate response vectors branch. H3 decomposes supported responses into source and lineage assembly. H4 asks whether continuous island area modifies the filtering. H5 is deliberately withheld from the plant-only data and requires independent pollination-channel evidence after the H2–H4 safeguards.

The primary plant response is deliberately pollinator-name-free and separates two components: accessibility/generalization of floral architecture and reproductive assurance. Only after this pattern is established do we examine fixed large-bee-like, butterfly-like, and bird-like floral-architecture templates. These secondary templates are concordance scores, not visitor observations or pollinator classifiers.

Finally, we ask when additional database acquisition stops being the main scientific bottleneck. The relevant stopping criterion is not complete trait coverage. It is whether large additions of usable evidence continue to change the inferential structure, and whether remaining gaps can be recovered at reasonable evidential cost without weakening the prespecified evidence rules.

---

## Materials and Methods

### Island and species universe

The fixed geographic universe comprised 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. Island occurrences, floristic-status information, geographic isolation, island area, climate covariates, and source-pool infrastructure were held constant across trait waves.

The formal exposure was `log1p_distance_to_continent_km`. We did not interpret this as a pure causal measure of one process. Instead, it represents a composite gradient of geographic separation, connectivity, and source accessibility. Conceptually, as distance approaches zero, an island approaches a high-accessibility boundary in which geographic filtering relative to a regional source system should be weaker. This is a limiting interpretation, not a claim that a zero-distance island is literally equivalent to mainland vegetation.

Analyses were stratified by floristic status (`all_native`, `native_nonendemic`, and `endemic` where support permitted) and evaluated at two context layers: broad analysis regimes and formal biogeographic realms. Native non-endemic persistence was used to establish that a pattern was not confined to endemic taxa; it was not treated as a causal endemicity contrast.

### Trait axes and evidence hierarchy

Each species could contribute to three primary axes:

1. `flower_colour`;
2. `floral_structural_complexity`;
3. `reproductive_assurance`.

The final snapshot therefore contained 318,885 possible species-by-axis cells. Evidence precedence was fixed: species-direct High/Medium evidence preceded trait-specific validated Low evidence. Family inference, global fallback, and post hoc threshold relaxation were prohibited. Missing trait evidence was never converted to trait absence.

Reproductive assurance retained distinctions among self-incompatibility, mating system, and autonomous selfing where source methodology allowed. Self-compatibility alone was not equated with realized selfing, and pollination/autogamy descriptions were not used as shortcuts for autonomous selfing.

### Progressive analysis contract

All trait waves were materialized into the same `chapter1_trait_snapshot_v1` schema and analysed under the frozen `chapter1_progressive_analysis_v1` contract. The contract fixed hypotheses, estimands, island universe, model order, support gates, multiplicity rules, source/lineage/area safeguards, and the causal claim ceiling. Biological estimates were explicitly allowed to change when trait evidence changed.

The analysis followed five linked questions.

**H1 — universal-syndrome rival.** Does increasing isolation generate one coherent floral/reproductive response vector across contexts?

**H2 — biogeographic branching.** If H1 fails, where do isolation-associated response vectors differ? The primary response contains `accessibility_generalization` and `reproductive_assurance` and is evaluated with multivariate within-context tests and direct between-context vector tests.

**H3 — source/lineage assembly.** At what taxonomic depth does a supported response remain? We compare the observed response, the response after family composition, the response after source-matched genus composition, and lineage entry/within-genus species loading where evaluable.

**H4 — area/capacity moderation.** Does continuous island area modify the isolation response? We do not introduce a post hoc small/large threshold.

**H5 — channel-gated residual mechanism.** After H2–H4, does an independently measured deficit in a regionally available pollination channel explain residual floral/reproductive filtering? H5 requires independent channel evidence and is not inferred from floral phenotype.

### Support gates and multivariate priority

Global fill fraction was descriptive and never an analysis gate. Response-specific support was classified as <30 supported islands = not promoted, 30–49 = pilot, and >=50 = count component of confirmatory support. Pairwise context comparisons required the same response to meet the same support tier in both contexts. BH-FDR correction was applied within prespecified test families.

Multivariate context-vector tests preceded interpretation of individual axes. A significant result in one context and a non-significant result in another was never treated as evidence of heterogeneity without a direct between-context test.

### Source-pool and taxonomic-depth decomposition

Supported responses were decomposed against outcome-blind mainland source expectations using GIFT source flora and fixed source assignments. Family and genus were used as grouping structures, not as trait-imputation devices.

Conceptually, if `T_source` is the expected regional source-pool trait composition and `T_island(d)` is the observed island composition at source separation `d`, the ecological deviation is

`DeltaT(d) = T_island(d) - T_source`.

The primary H2 model does not fit this expression literally. Instead, H2 estimates the geographic response and H3 determines how much of that response is generated by source and lineage assembly. This separation avoids pretending that `distance_to_continent` is the exact historical source distance for every island.

### Area/capacity moderation

H4 used continuous distance-by-continuous-area models. Equal-island, capped-information, direct-only, outcome-blind common-support, and heteroskedastic-null safeguards were applied. A directional interaction was not promoted to founder filtering, habitat capacity, or pollinator persistence unless the prespecified null and measurement-sensitivity gates passed.

### Reproductive versus pollination-associated floral components

`selfing_core` included only self-incompatibility, mating-system, and autonomous-selfing information and excluded floral attraction traits. A separate attraction/access pathway was constructed from floral architecture. Conditional models asked whether the isolation-associated attraction shift remained after conditioning on `selfing_core`. This was interpreted as decomposition, not causal mediation.

Secondary pollination-associated concordance used fixed `large_bee_like`, `butterfly_like`, and `bird_like` templates. These combined colour, floral form, symmetry, tube depth, and flower size with prespecified weights and were analysed only after the primary plant response. The same templates were rerun in direct-only and source-adjusted sensitivities.

Because the templates share many floral traits, V4 estimated a source-trained common factor before examining template-specific residuals. Simultaneous movement of several guild-labelled scores was therefore not interpreted as simultaneous change in several realized pollinator guilds.

### Pollination-channel mobility as a mechanistic extension

Chapter 1 does not estimate pollinator dispersal directly, but its ecological interpretation distinguishes three channel states: retained, disrupted/deficient, and structurally absent. Structural absence in the source region is not treated as island-driven loss.

For a pollination channel `g`, the unobserved continuity function can be written conceptually as

`E_poll(g,d) = f(source availability, distance, flight/dispersal ability, establishment ability, habitat suitability, realized community context)`.

For a plant lineage `z`, the corresponding source-to-island assembly function is

`E_plant(z,d) = f(source availability, distance, propagule dispersal ability, establishment ability, persistence conditions)`.

These functions are conceptual, not fitted Chapter 1 models. They formalize the double-filter mechanism passed to H5 and Chapter 2.

### Missingness and robustness

Trait-resolution MNAR sensitivity used a finite prespecified grid rather than assuming missingness at random. Climate-overlap analysis separated reference climate adjustment, outcome-blind overlap weighting, and distance-by-climate interaction tests. The progressive workflow also retained direct-only evidence, information-weight sensitivity, source adjustment, taxonomic-depth decomposition, and area-support falsification.

---

## Results

### H1 — a universal floral island syndrome was not recovered

The primary two-axis response did not recur as one coherent direction across context layers, evidence scopes, and floristic strata. Supported responses depended on biogeographic context. The result rejects a simple model in which increasing geographic separation universally pushes island floras toward the same combination of generalized floral access and increased reproductive assurance.

This is not based on comparing significance labels across regions. The analysis prioritized multivariate within-context response vectors and direct between-context comparisons.

### H2 — source separation was associated with different regional response trajectories

The Palearctic showed the clearest broad primary response. Increasing separation was associated with greater accessibility/generalization and greater reproductive assurance across all-analysis and direct-only evidence and in both all-native and native-nonendemic strata.

Tropical analyses provided a contrasting branch. Reproductive assurance could increase while accessibility/generalization decreased, corresponding to maintenance or strengthening of specialized-access architecture. In other words, the same direction of geographic separation was compatible with a different combination of floral and reproductive change.

These regional differences did not establish a purely categorical biogeographic mechanism independent of measured climate. North–Tropical contrasts failed the outcome-blind common-support positivity gate, while Palearctic–Neotropical climate-adjusted support did not replicate consistently across both evidence scopes and both primary floristic strata. The defensible conclusion is therefore biogeographically contingent, region-associated filtering with a measured-climate identification ceiling.

### The classic syndrome decomposed into partially independent reproductive and floral components

The Palearctic attraction/access shift persisted after conditioning on `selfing_core`. Thus the floral response could not be reduced to the measured selfing syndrome alone.

Tropical responses showed the complementary pattern: reproductive assurance could increase while specialized or attractive floral architecture was maintained or strengthened. These results reject a compulsory serial model of `isolation -> selfing -> floral simplification` and instead support partially parallel branches whose signs and coupling depend on context.

### H3 — the strongest broad distance response resolved to genus-level lineage assembly

Taxonomic-depth decomposition produced the most consequential conclusion across trait waves. In Wave36, the Palearctic primary response appeared to persist beyond both family and genus composition. By Wave52, additional evidence changed that inference. The response persisted beyond family composition but disappeared after source-matched genus adjustment.

The final snapshot reproduced the Wave52 classification exactly. Across all-analysis/direct-only evidence and all-native/native-nonendemic strata, the response was supported at the observed stage (`4/4` source modes), retained after family adjustment (`4/4`), and absent after genus adjustment (`0/4`). The number of scored GIFT source species increased to 7,750 in all-analysis and 4,552 in direct-only evidence without changing this classification.

This result supports the plant side of the double geographic filter. Much of the broad Palearctic source-distance response is compatible with differential representation of source-available genera rather than a robust beyond-genus primary residual.

### Pollination-associated floral architecture showed opposite regional trajectories

Secondary floral templates showed strong geographic contrasts. Palearctic assemblages moved away from butterfly-like and large-bee-like floral architecture with increasing isolation, whereas tropical assemblages increased along bird-like, butterfly-like, and large-bee-like templates in both evidence scopes and both primary floristic strata.

These guild-labelled results cannot be interpreted as realized visitor identities. A source-trained first factor explained 86.94% of variation among the three templates in all-analysis evidence and 86.44% in direct-only evidence. The simultaneous tropical increase therefore primarily represents a shared specialized/attractive floral-architecture dimension.

Some tropical template-specific residual structure remained after source/genus adjustment, including a positive large-bee-like residual in several robust comparisons. This is a plant-side residual and does not establish large-bee occurrence, effectiveness, or functional replacement.

### H4 — area remained a modifier rather than an identified mechanism

All 16 frozen primary V3 classifications remained `retain_area_as_measurement_sensitive_modifier_only`, with zero heteroskedastic-null passes. Directional distance-by-area patterns therefore did not justify promotion to founder filtering, habitat-capacity, or pollinator-persistence mechanisms.

### H5 — the pollinator side of the double filter remained unidentified

No independent pollination-channel exposure, disruption, visitation, single-visit effectiveness, or effective-service data entered the Chapter 1 primary analysis. The global plant database therefore cannot establish historical Bombus loss, bird/butterfly replacement, or another realized interaction mechanism.

The plant-side pattern nevertheless sharpens the mechanism that must be tested. A regionally available channel would need to show retention or disruption along the isolation gradient, and that measured disruption would need to explain residual plant filtering beyond source, lineage, area, climate, observation, and spatial safeguards.

### Progressive waves reached an inferential and acquisition saturation point

Wave52 contained 184,917 analysis-usable species-by-axis cells (57.99%). The final integrated snapshot contained 222,688 cells (69.83%), a gain of 37,771 usable cells, including 11,315 reproductive-assurance cells. Despite this large increase, the current H1–H5 claim ceiling remained stable.

The acquisition process simultaneously entered strong diminishing returns. One high-yield Orchidaceae source contributed 192 strict reproductive cells. Later reviewed source packets usually added only one to three cells. A residual audit scanned 1,187 staging files and more than 2.25 million rows, including approximately 99,000 reproductive rows, and found no immediately strict-ready reproductive records under the frozen evidence contract.

Thus the stopping point is inferential rather than numerical. Further targeted acquisition may still matter for poorly supported contexts, but indiscriminate global trait searching is no longer the highest-value route to changing the main Chapter 1 conclusion.

---

## Discussion

### Oceanic isolation is a double geographic filter

The central interpretation of the results is that distance should not be read as a direct pressure on floral phenotype. The same geographic barrier filters two coupled systems. It changes which plant lineages remain represented from the regional source pool and may also change which pollination channels remain functionally available after ocean crossing and establishment.

The strong H3 genus-depth result directly supports the plant side of this framework. Increasing source separation can generate an assemblage-level floral/reproductive gradient by changing which source-available genera arrive, establish, and persist. A real island syndrome can therefore emerge without repeated parallel within-lineage evolution.

The pollinator side remains prospective. Different functional groups may experience the same ocean distance differently because flight ability, long-distance dispersal, stepping-stone use, propagule number, habitat requirements, and establishment probability differ. The relevant concept is therefore **channel isolation susceptibility**, not a simple mobile-versus-immobile dichotomy.

This distinction also clarifies why structural absence is not loss. If a pollination channel is absent from the regional source system, its absence on an island cannot be interpreted as an isolation-driven disruption. H5 requires source availability to be defined independently of floral outcomes before retention or deficit is assessed.

### Distance is best interpreted as separation from a source system

The design is not fundamentally an island-versus-mainland comparison. It is a continuous source-accessibility analysis. In a conceptual limiting sense, `d -> 0` represents a high-accessibility boundary in which source-to-island filtering should be relatively weak. As separation increases, both plant assembly and interaction continuity can become increasingly filtered.

This limiting interpretation should not be confused with the formal exposure. PR142 uses distance to the continent in the primary model, while source-pool assignments and source-matched expectations enter explicitly in H3. We therefore interpret distance as a source-accessibility gradient rather than claiming that it equals exact historical source distance for every island.

The combination of H2 and H3 nevertheless provides the key biological decomposition. H2 identifies how island assemblages move along the geographic gradient; H3 asks how much of that movement is generated by source and lineage assembly.

### The island floral syndrome contains at least two response components

The traditional syndrome often bundles reproductive assurance with floral simplification. Our results show that these components can be partially independent. In the Palearctic, attraction/access architecture shifted with isolation even after conditioning on `selfing_core`. In tropical assemblages, reproductive assurance could increase while specialized or attractive architecture was maintained or strengthened.

The ecological implication is that reproductive assurance and floral architecture are better treated as parallel response branches whose coupling depends on biogeographic and assembly context. This explains why local studies may recover familiar syndrome-like combinations while global comparisons contain reversals and exceptions.

### Pollination-syndrome concordance links phenotype to candidate channels without identifying them

The secondary guild-labelled templates remain useful because they provide a multivariate connection to pollination theory. But their shared-factor structure makes hard visitor assignment inappropriate. The Palearctic decline in butterfly-like and large-bee-like architecture and the tropical increase of all three templates are best understood as movements along overlapping plant-side architecture dimensions.

The double-filter framework adds a testable mechanism to this pattern. If a source-available pollination channel is especially susceptible to oceanic separation, increasing distance could weaken its functional continuity. If a functionally comparable alternative has lower barrier susceptibility, specialized floral architecture may be maintained despite distance. Chapter 1 does not determine which guild has which susceptibility; that is the mechanistic hypothesis generated by the plant-side pattern.

### The progressive wave series strengthens the biological result rather than defining it

The database campaign is important because it shows that the main conclusion is no longer easily dismissed as one incomplete trait snapshot. Wave52-to-final integration added 37,771 analysis-usable cells, yet the main H1–H5 ceiling remained stable. At the same time, new strict-source yield collapsed to small increments.

The practical freeze therefore follows from a combination of inferential stability and acquisition saturation. The ecological novelty is the conditional, decomposable island response—not the final coverage percentage.

### From WHEN/WHERE to HOW/WHY

Chapter 1 ends at an explicit identification boundary. It establishes where source separation is associated with assemblage-level floral/reproductive filtering, where response components decouple, and at what taxonomic depth the strongest response is generated. It does not identify the realized pollination process linking geographic separation to effective service.

The companion `izu-core` programme addresses this missing layer. Its conditional-response framework separates partner loss/arrival balance, plant starting functional state, realized pollinator community, functional matching, effective service, local filtering, and autonomous assurance. This is the natural mechanistic extension of the double geographic filter: the same barrier changes partner arrival and plant assembly, but the realized downstream response depends on the particular interaction system that forms.

The Chapter 2 model is not used as retroactive validation of the Chapter 1 regional vectors. It instead supplies a mechanistic class capable of generating different response branches under a common broad perturbation.

### Returning to islands and zooming to Izu

The prospective Izu field programme measures the channels that remain absent from the global database in the same linked populations: observation effort and visit rate, visitor identity/contact, background-adjusted single-visit pollen deposition, rate-weighted effective pollen service, open-pollinated reproduction, bagged autonomous reproduction, supplemental-outcross reproduction, and mature fruit/seed outcomes.

This allows direct discrimination among service limitation, autonomous-assurance buffering, network/service allocation, and non-assurance buffering without using floral phenotype as a proxy for dependency. Izu is therefore a depth system selected for measurement continuity, not a confirmation case chosen because it resembles the global pattern.

---

## Conclusion

Increasing oceanic separation does not impose one universal floral island syndrome. It acts as a double geographic filter on plant source-pool assembly and, potentially, on pollination-channel persistence. The strongest broad Palearctic response combines greater accessibility/generalization and reproductive assurance, but its primary taxonomic depth is compatible with source-matched genus-level assembly. Tropical assemblages demonstrate that reproductive assurance can coexist with maintained or strengthened specialized floral architecture. Thus the classical floral island syndrome decomposes into partially independent response components whose coupling is biogeographically contingent.

The plant side of the double filter is supported directly by the source/genus decomposition. The pollinator side remains an identification problem: guild-labelled floral architectures show regional contrasts, but they do not establish realized visitor identity, dispersal, establishment, effective service, loss, or replacement.

Large late gains in usable trait coverage did not alter this claim ceiling, while new evidence acquisition entered strong diminishing returns. The global computational campaign can therefore be frozen defensibly and the main scientific effort moved to the missing data type.

The resulting research sequence is explicit: Chapter 1 uses global databases to establish WHEN and WHERE the double geographic filter is expressed; Chapter 2 asks HOW partner loss/arrival and realized interaction structure generate different response branches; comparative island systems confront those mechanisms with independent evidence; and Izu provides the linked depth test from visitor structure through effective pollen service and reproductive dependency to realized phenotype.

---

## Proposed figures

### Figure 1 — Double geographic filter and PR142 hypothesis ladder
Left: source distance acting on plant source-pool assembly and pollination-channel continuity. Right: H1 universal rival -> H2 regional branching -> H3 source/genus decomposition -> H4 area moderation -> H5 independent channel mechanism.

### Figure 2 — Distance as a source-accessibility continuum
Conceptual `d -> 0` high-accessibility boundary, near-source islands, intermediate islands, and remote islands. Overlay the distinction between the formal distance exposure and H3 source-matched expectations.

### Figure 3 — WHERE: primary response branches
Regional/realm slopes for accessibility/generalization and reproductive assurance, separated by evidence scope and floristic stratum.

### Figure 4 — Plant-side filter: source and taxonomic depth
Palearctic observed -> after family -> after genus (`4/4 -> 4/4 -> 0/4`) with scored source-species support.

### Figure 5 — Decomposing the floral island syndrome
`selfing_core` versus attraction/access architecture, highlighting Palearctic persistence of floral shift after selfing-core conditioning and tropical coexistence of reproductive assurance with specialized architecture.

### Figure 6 — Pollination-associated architecture and channel susceptibility hypothesis
Observed large-bee/butterfly/bird concordances, shared-factor decomposition, and a conceptual panel showing different channel continuity curves across source separation without assigning empirical guild-specific dispersal rates.

### Figure 7 — Wave convergence and computational stopping rule
Wave36/Wave52/final usable coverage, source-scale reproductive gain, residual-audit exhaustion, and H1–H5 decision stability.

### Figure 8 — Mechanistic handoff and Izu zoom
Chapter 1 global pattern -> Chapter 2 partner arrival/loss and response geometry -> comparative island confrontation -> Izu visitation -> SVD -> effective service -> dependency/assurance -> reproduction -> phenotype.

---

## Reproducibility anchors

- Progressive scientific contract: `chapter1_progressive_analysis_v1`.
- Wave52 reanalysis: Run `33587356935`, artifact `chapter1-progressive-wave52-33587356935`.
- Final trait integration: Run `34191508045`, `222,688 / 318,885` resolved cells.
- Final PR142 progressive reanalysis: Run `34232450884`, artifact `chapter1-progressive-analysis-34232450884`, artifact ID `10058653212`, digest `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
