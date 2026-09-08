# Decomposing the island floral syndrome: biogeographically contingent assemblage filtering, reproductive assurance, and pollination-associated floral architecture

## Working manuscript draft — 2026-09-08

### Abstract

Island floras are often expected to converge toward a coherent floral “island syndrome”, including simplified or less conspicuous flowers and increased reproductive assurance. Yet the same observed pattern can arise through colonization filtering, lineage composition, reproductive system differences, or changes in pollination interactions, making a universal syndrome difficult to distinguish from context-dependent assembly. We developed a progressive global analysis in which the island universe, trait ontology, model order, support thresholds, source-pool safeguards, and causal claim ceiling were fixed while species-level trait evidence was allowed to improve across acquisition waves. The final computational snapshot comprised 106,295 accepted angiosperm species scored on three primary axes—flower colour, floral structural complexity, and reproductive assurance—yielding 222,688 resolved species-by-axis cells out of 318,885 (69.83%). Relative to Wave52, 37,771 previously unresolved analysis cells became usable, including 11,315 reproductive-assurance cells. Despite this large increase, the principal hypothesis decisions remained stable.

We found no evidence for one universal floral/reproductive response to island isolation. The most reproducible broad response occurred in the Palearctic, where increasing isolation was associated with both greater floral accessibility/generalization and greater reproductive assurance across all-analysis and direct-only evidence and in native non-endemic assemblages. However, source-matched taxonomic decomposition showed that this primary response survived family adjustment but not genus adjustment in all four evidence-scope-by-stratum combinations, making genus-level lineage assembly the strongest supported depth. Measured-climate-independent categorical biogeography was not established, and area remained a measurement-sensitive modifier rather than an identified capacity or founder-filter mechanism. Secondary pollination-associated templates showed regionally contrasting floral architectures, but a single shared factor explained 86.9% and 86.4% of variation among large-bee-like, butterfly-like, and bird-like scores in all-analysis and direct-only data, respectively, demonstrating that these labels primarily represent overlapping floral architecture rather than realized pollinator identity.

The acquisition campaign itself approached an empirical saturation point. Wave52 contained 184,917 analysis-usable cells, whereas the final integrated snapshot contained 222,688; yet source-reported coverage increased by only 362 cells over the same comparison, showing that most of the late gain came from making existing evidence analytically materializable rather than discovering large new evidence pools. After one high-yield reproductive source (+192 cells), subsequent reviewed source batches yielded only small gains, while a committed-evidence audit of more than 2.25 million rows found no immediately strict-ready reproductive records. We therefore treat the final snapshot as a defensible computational stopping point based on recoverability and inference stability, not on an arbitrary global coverage percentage.

Together, the results recast the island floral syndrome as an assemblage-level outcome composed of at least two partially separable response components: reproductive assurance and pollination-associated floral architecture. Chapter 1 identifies when and where these components are expressed and where source/lineage assembly constrains interpretation. It does not identify historical pollinator loss, effective pollination service, or functional replacement. Those mechanisms require an explicit handoff to independent interaction data and prospective field tests.

---

## Introduction

Geographic isolation is one of the oldest natural experiments in ecology and evolution. Island floras repeatedly differ from mainland floras in dispersal, growth form, reproductive systems, and floral traits, motivating the idea of a plant “island syndrome”: a suite of predictable changes associated with insularity. Floral versions of this idea commonly combine two expectations. First, mate and pollinator limitation may favour reproductive assurance, consistent with Baker’s law and subsequent work on self-compatibility and colonization. Second, altered pollinator communities may change the returns to floral attraction or specialization, potentially favouring smaller, less conspicuous, or more accessible flowers. These expectations are biologically plausible but are not equivalent mechanisms.

This distinction matters because selfing syndrome and pollination-associated floral change can produce superficially similar phenotypes. Transitions toward selfing often involve reduced floral investment and changes in mating and compatibility traits, whereas pollinator-mediated selection acts on multivariate combinations of colour, form, symmetry, tube depth, and size. Floral morphology can therefore become simpler because reproductive assurance has increased, because the effective pollination environment has changed, because the source flora already differs in lineage composition, or because several processes act in parallel. Conversely, reproductive assurance can increase without loss of specialized floral architecture. Treating all of these outcomes as one island syndrome risks turning a broad pattern into a presumed mechanism.

A second difficulty is scale. Most evidence for plant island syndromes comes from particular archipelagos, taxonomic groups, or island–mainland contrasts. Global comparisons must also confront heterogeneous species pools, climates, island sizes, historical connectivity, observation effort, and lineage structure. Recent studies have questioned whether canonical floral island rules are globally consistent; for example, island flowers are not universally smaller than mainland relatives, and reviews of regional island syndromes reveal a mixture of supported, tentative, and unsupported components. A global test should therefore ask not only whether an average syndrome exists, but where a response is detectable, whether response vectors differ among biogeographic contexts, and how much of any signal is attributable to source-pool and lineage assembly.

Here we implement that strategy using a progressive analysis contract. Rather than freezing biological conclusions at an early, incomplete trait snapshot, we froze the hypotheses, island universe, trait ontology, support thresholds, model order, multiple-testing families, and claim ceiling. New trait waves were then reanalysed from the beginning under the same contract. This design allowed better data to change the biological conclusion without changing the question.

We test five linked hypotheses. H1 treats a universal island syndrome as a rival: if isolation imposes one coherent floral/reproductive response, similar multivariate vectors should recur across contexts. H2 tests biogeographic branching: isolation-associated response vectors may differ among prespecified regions or realms. H3 decomposes supported patterns into source availability, family and genus composition, lineage entry, and within-genus loading. H4 asks whether continuous island area modifies isolation effects strongly enough to support an area/capacity interpretation. H5 is deliberately outside the direct Chapter 1 evidence: an independently measured pollination-channel deficit should explain residual filtering only after source, lineage, area, climate, observation, and spatial safeguards are satisfied.

We further separate the traditional floral island syndrome into two interpretive components. The first is a reproductive-assurance component based on compatibility, mating-system, and autonomous-selfing evidence. The second is a pollination-associated floral component based on accessibility/generalization and fixed, non-exhaustive floral-architecture templates associated in the literature with large-bee, butterfly, and bird functional groups. These guild-labelled scores are not pollinator observations. Their role is to test concordance between plant assemblage responses and sampled pollination-associated architectures after the primary pollinator-name-free response has been estimated.

Finally, we use the progressive trait campaign itself to determine when additional computation and database acquisition cease to be the main scientific bottleneck. The goal is not complete trait coverage. It is to reach a point where large increases in usable evidence no longer change the major inference, while remaining unresolved cells require disproportionately expensive or methodologically weaker recovery. This makes the stopping decision part of the reproducibility record rather than an arbitrary endpoint.

---

## Materials and Methods

### Island and species universe

The fixed geographic universe comprised 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. Island occurrences, floristic-status information, geographic isolation, island area, climate covariates, and source-pool infrastructure were held constant across trait waves. Isolation was represented primarily by a log-transformed mainland-distance/source-accessibility gradient and interpreted as a composite of isolation, connectivity, and source supply rather than as a single causal mechanism.

Analyses were stratified by floristic status (`all_native`, `native_nonendemic`, and `endemic` where support permitted) and evaluated at two context layers: broad analysis regimes (including northern mid-latitude and tropical contexts) and formal biogeographic realms. Native non-endemic persistence was used to distinguish broad regional filtering from patterns confined to island endemics; it was not treated as a causal endemic-versus-nonendemic contrast.

### Trait axes and evidence hierarchy

Each species could contribute to three primary axes:

1. `flower_colour`;
2. `floral_structural_complexity`;
3. `reproductive_assurance`.

The final snapshot contained 318,885 possible species-by-axis cells. Evidence precedence was fixed: species-direct High/Medium evidence preceded trait-specific validated Low evidence. Family inference, global fallback, and post hoc threshold relaxation were prohibited. Missing trait evidence was never converted to trait absence.

Reproductive assurance retained distinctions among self-incompatibility, mating system, and autonomous selfing where source methodology allowed. Self-compatibility alone was not equated with realized selfing, and pollination/autogamy descriptions were not used as shortcuts for autonomous selfing.

### Progressive wave contract and computational stopping point

All trait waves were materialized into the same `chapter1_trait_snapshot_v1` schema and analysed under `chapter1_progressive_analysis_v1`. The contract fixed scientific hypotheses, estimands, support gates, multiplicity rules, source/lineage/area safeguards, and the causal claim ceiling. Biological estimates and significance were explicitly allowed to change as evidence improved.

Wave36 provided an early reference point with 57.14% analysis-usable coverage. Wave52 increased usable coverage to 184,917/318,885 cells (57.99%), including 37,182 reproductive-assurance species. However, Wave52 also contained 37,409 source-reported Low cells that lacked a materializable trait composition. Recovery and wave integration subsequently restored these evidence states without changing their evidence tier. In the final snapshot, 222,688 cells were analysis-usable (69.83%), including 82,556 colour, 91,635 structural, and 48,497 reproductive-assurance cells.

The final computational stopping point was not defined by a target percentage. Instead, we required evidence of practical saturation: (i) a large increase in analysis-usable coverage without a corresponding change in the major H1–H5 decisions; (ii) declining realized yield from new reviewed sources; and (iii) exhaustion of immediately promotable committed reproductive evidence under the strict evidence contract. The committed residual audit scanned 1,187 staging files and more than 2.25 million rows, including approximately 99,000 reproductive rows, and found no immediately strict-ready reproductive records. Remaining records required source-specific methodological review, resolution of conflicts, or new external evidence.

### Support gates and multivariate priority

Global fill fraction was descriptive and never an analysis gate. Response-specific support was classified as:

- <30 supported islands: not promoted;
- 30–49: pilot;
- >=50: count component of confirmatory support.

Pairwise context comparisons required the same response to meet the same support tier in both contexts. BH-FDR correction was applied within prespecified test families. Multivariate context-vector tests preceded interpretation of individual axes; significance in one context and non-significance in another was not treated as evidence of heterogeneity.

### Primary plant response and biogeographic branching

The primary plant-side response was pollinator-name-free and contained two axes:

- `accessibility_generalization`: increasing open, shallow, radially accessible floral architecture;
- `reproductive_assurance`: increasing strict reproductive-assurance traits.

Within-context multivariate tests asked whether these responses changed with isolation after area and climate adjustment. Direct between-context tests evaluated whether response vectors differed among prespecified contexts. Observation-selection sensitivity and direct-only evidence were retained as mandatory robustness layers.

### Source-pool and taxonomic-depth decomposition

Supported responses were next decomposed against outcome-blind mainland source expectations. GIFT source flora and fixed source assignments were used to estimate source-matched trait positions. Taxonomic decomposition was performed on common observed species and evaluated the response at three stages:

1. observed score;
2. residual after family composition;
3. residual after genus composition.

Family and genus were grouping structures, not trait-imputation devices. The decomposition does not identify historical ancestry or in-situ evolution; it only determines the taxonomic depth at which the contemporary assemblage pattern remains.

### Area/capacity moderation

H4 used continuous distance-by-continuous-area models; no post hoc small/large island threshold was introduced. Equal-island, capped-information, direct-only, outcome-blind common-support, and heteroskedastic-null safeguards were applied. A directional distance-by-area pattern was not promoted to founder filtering, habitat-capacity, or pollinator-persistence mechanism unless the prespecified null and measurement-sensitivity gates were passed.

### Decomposing reproductive and pollination-associated components

Reproductive and floral components were separated before interpretation. `selfing_core` included only self-incompatibility, mating-system, and autonomous-selfing information and excluded floral attraction traits. A separate attraction/access pathway was constructed from floral architecture. Conditional models tested whether an isolation-associated attraction shift remained after conditioning on `selfing_core`; this was interpreted as a decomposition, not causal mediation.

Secondary pollination-associated concordance used fixed `large_bee_like`, `butterfly_like`, and `bird_like` templates. Template definitions were frozen before the latest trait snapshot and combined colour, floral form, symmetry, tube depth, and flower size with prespecified weights. These scores were analysed only after the primary plant response and were repeated in direct-only and source-adjusted sensitivities.

Because the templates share many floral features, we additionally estimated a source-trained common factor and examined template-specific residuals after source and taxonomic adjustment. This prevented simultaneous movement of several guild-labelled scores from being interpreted as simultaneous changes in several realized pollinator guilds.

### Missingness and robustness

Trait-resolution MNAR sensitivity used a finite prespecified grid rather than assuming missingness at random. Climate-overlap analysis separated reference climate adjustment, outcome-blind overlap weighting, and distance-by-climate interaction tests. The progressive contract prohibited lowering support gates or selecting only the evidence scope producing the preferred result.

---

## Results

### Wave integration greatly increased usable evidence, but late source acquisition showed diminishing returns

Wave52 contained 184,917 analysis-usable species-by-axis cells (57.99%). The final integrated snapshot contained 222,688 cells (69.83%), a gain of 37,771 usable cells. The increase was concentrated in flower colour (+26,434 cells) and reproductive assurance (+11,315), with a smaller increase in structural complexity (+22).

Importantly, this was not equivalent to discovering 37,771 new source-reported cells. Source-reported coverage increased from 222,326 cells at Wave52 to 222,688 at the final snapshot, only +362 cells. The large analytical gain therefore came primarily from resolving the composition/provenance needed to materialize previously labelled evidence, especially validated-Low evidence, rather than from discovering a new database of comparable scale.

The source-scale reproductive campaign showed a second form of saturation. The first high-yield source contributed 192 strict new SI/SC cells. Subsequent reviewed sources produced much smaller increments, typically one to three cells per source packet. By the seven-batch checkpoint, a residual audit of 1,187 files and 2,253,340 rows found zero immediately strict-ready reproductive rows; later source packets continued to yield only isolated additions. Thus, the limiting step shifted from computation to provenance and primary-source interpretation.

Despite the 11.84 percentage-point increase in global analysis-usable coverage from Wave52 to the final snapshot, the main H1–H5 claim ceiling remained stable. We therefore treat the final snapshot as the working computational saturation point for Chapter 1.

### H1: a universal floral island syndrome was not recovered

The primary two-axis response did not recur as one coherent direction across all context layers, evidence scopes, and floristic strata. Instead, the supported response depended on biogeographic context and evidence scope. This rejects a simple model in which isolation universally pushes island floras toward the same combination of generalized floral access and increased reproductive assurance.

The result also argues against treating a classic syndrome projection as the primary statistic. Individual traits and secondary floral templates sometimes moved in the expected direction, but those movements were not globally coherent enough to define one universal island rule.

### H2: the clearest broad response was Palearctic, while regional branching remained climate-bounded

At the formal biogeographic-realm level, the Palearctic showed the most reproducible primary response. In all-analysis data, both all-native and native-nonendemic assemblages showed supported multivariate response vectors, with positive accessibility/generalization and reproductive-assurance slopes. The pattern was retained in direct-only evidence. Thus, the Palearctic signal was not dependent on validated-Low evidence and was not confined to endemic species.

Other contexts showed different combinations. Tropical direct-only analyses supported a branch combining decreased accessibility/generalization (maintenance or increase of specialized-access architecture) with increased reproductive assurance. In all-analysis data, the same tropical branch was strongest in native non-endemics. Neotropical responses were weaker and more evidence-scope-sensitive: some all-analysis vectors remained multivariate-supported without a stable single-axis interpretation, whereas direct-only support was reduced.

These differences did not establish a purely categorical biogeographic mechanism independent of measured climate. North–Tropical contrasts failed the outcome-blind common-support positivity gate. Palearctic–Neotropical comparisons had adequate overlap, but climate-adjusted support did not replicate consistently across both evidence scopes and both primary floristic strata. The defensible H2 conclusion is therefore region-associated, biogeographically contingent assemblage filtering with a measured-climate identification ceiling, not a causal effect of realm identity itself.

### H3: the primary Palearctic response resolved to genus-level lineage assembly

Taxonomic-depth decomposition provided the strongest change across earlier waves and the most stable conclusion in the final snapshot. In Wave36, the Palearctic primary two-axis response remained after both family and genus adjustment. By Wave52, additional reproductive evidence changed this result: the response persisted beyond family composition but disappeared after genus adjustment.

The final high-coverage snapshot reproduced the Wave52 classification exactly. Across all-analysis/direct-only evidence and all-native/native-nonendemic strata, the primary Palearctic response was supported at the observed stage (4/4 source modes), retained after family adjustment (4/4), and absent after genus adjustment (0/4). The number of scored source species increased to 7,750 in all-analysis and 4,552 in direct-only evidence, yet the taxonomic-depth conclusion did not revert.

Thus, the strongest current interpretation is not a robust beyond-genus trait shift. Instead, the broad Palearctic gradient is compatible with non-random representation of source-available genera beyond family-level composition. This makes lineage assembly part of the biological result rather than a nuisance correction.

### H4: island area modified precision but did not identify a capacity mechanism

All 16 frozen primary area classifications remained `retain_area_as_measurement_sensitive_modifier_only`. Several fitted cells showed the expected direction of stronger isolation effects on smaller islands, and direct-only evidence did not contradict the direction. However, none of the 16 classifications passed the heteroskedastic-null gate required for mechanistic promotion.

Accordingly, Chapter 1 supports area as a potentially important modifier of the observed gradient but does not identify target size, founder filtering, habitat capacity, pollinator persistence, or another specific area-mediated mechanism.

### Reproductive assurance and pollination-associated floral architecture partially decoupled

The two-component decomposition clarified why one island syndrome is insufficient. In the Palearctic source-adjusted analysis, attraction/access shifts remained positive after conditioning on `selfing_core` across all four source definitions. For example, conditional distance effects remained supported in both all-native and native-nonendemic strata under multiple source definitions, while the contemporaneous `selfing_core` coefficient in the attraction model was not supported. Thus, the Palearctic floral-architecture shift was not reducible to the measured reproductive selfing core.

Tropical assemblages showed a contrasting structure. Source-adjusted attraction shifts were negative, indicating maintenance or increase of specialized/attractive architecture with isolation. In native non-endemics, reproductive assurance could increase at the same time. Therefore increased reproductive assurance and maintenance of specialized floral architecture are not mutually exclusive responses.

This decoupling is central to the interpretation: an island-associated reproductive response need not imply a selfing-syndrome floral reduction, and a floral architecture shift need not be a by-product of reproductive assurance.

### Secondary pollination-associated concordance was regional but largely shared architecture

The secondary floral templates showed clear regional contrasts. In all-analysis data, northern mid-latitude all-native assemblages showed decreasing butterfly-like concordance with isolation, while large-bee-like and bird-like axes were not FDR-supported. In direct-only data, both butterfly-like and large-bee-like concordance declined. In the Palearctic, butterfly-like and large-bee-like concordance declined consistently across evidence scopes and primary floristic strata, whereas bird-like decline was unsupported.

Tropical assemblages showed the opposite broad pattern: bird-like, butterfly-like, and large-bee-like concordances all increased with isolation in both evidence scopes and in all-native and native-nonendemic strata. Neotropical guild-labelled concordances were not supported.

However, the three sampled templates were strongly correlated. A source-trained first factor explained 86.94% of their variation in all-analysis evidence and 86.44% in direct-only evidence. Simultaneous tropical increases therefore cannot be interpreted as evidence that birds, butterflies, and large bees all became more important. Instead, the result indicates a shared floral-architecture dimension, with template-specific residuals used only as secondary contrasts.

Some tropical template-specific residuals remained beyond genus after source adjustment, including a positive large-bee-like residual across several source-mode/evidence-scope combinations. This does not identify large-bee presence, pollinator effectiveness, or replacement; it shows only that a residual component of floral architecture is not absorbed by the tested shared factor and taxonomic expectations.

### H5 remained unidentified by the global plant database

No independent pollination-channel exposure, disruption, visitation, single-visit effectiveness, or effective-service data entered the Chapter 1 primary analysis. Floral phenotype therefore cannot establish historical Bombus loss, bird/butterfly replacement, or another realized pollination mechanism. H5 remained not evaluable.

The data nevertheless sharpened the mechanistic question. Chapter 1 identifies regional combinations of reproductive assurance and floral architecture, determines which patterns are source/genus-associated, and identifies where area/climate/observation explanations remain unresolved. The remaining gap is not primarily another global trait query. It is independent measurement of the interaction channel linking pollinator environment to effective service and reproduction.

---

## Discussion

### The island floral syndrome is better treated as a decomposable assemblage outcome

Our results do not support one universal floral island syndrome. They instead support a more conditional view in which isolation changes the representation of floral and reproductive strategies, but the resulting assemblage response depends on biogeographic context and source-lineage composition. This helps reconcile why local island studies can recover syndrome-like patterns while global studies often find exceptions or reversals.

The key conceptual move is to separate two response components that are often bundled together. Reproductive assurance concerns the capacity to reproduce when mates or effective pollen service are limited. Pollination-associated floral architecture concerns how plant traits related to attraction, access, and functional matching are represented. Selfing syndrome predicts that these components may sometimes move together, but they need not. Our data show both possibilities: Palearctic attraction/access shifts persisted after conditioning on selfing core, whereas tropical reproductive assurance could increase while specialized floral architecture was maintained or strengthened.

Thus, “island syndrome” need not describe a deterministic lineage-level rule. It can instead describe an emergent assemblage-level tendency generated by different mixtures of source filtering, lineage representation, reproductive filters, and interaction change.

### Source and lineage assembly are not confounders to be erased

The most consequential result of the progressive analysis was the Wave36-to-Wave52 reversal in taxonomic depth. With better trait evidence and broader source coverage, the primary Palearctic response changed from apparently beyond-genus to compatible with genus-level assembly. The final snapshot retained that result despite a large additional increase in usable trait coverage.

This illustrates why source and lineage composition should be treated as part of the island assembly process. If particular source-available genera differ in floral access or reproductive strategies, isolation can change community-level trait distributions by altering which lineages reach, establish, and persist. A global trait gradient can therefore be biologically real without representing repeated within-lineage evolution of the same phenotype.

This is also the reason Chapter 1 cannot infer historical evolutionary transitions from the cross-sectional assemblage alone. The present analysis identifies the contemporary taxonomic depth of the pattern, not whether each lineage changed after colonization.

### Pollination-syndrome concordance is useful precisely because it is not pollinator identification

Pollination syndromes remain useful as multivariate hypotheses about floral architecture, but their strongest use here is comparative rather than classificatory. The three sampled guild templates moved differently among regions, yet most of their variance lay on one shared plant-architecture factor. This makes hard visitor assignment from floral phenotype especially inappropriate.

The Palearctic decline in butterfly-like and large-bee-like concordance and the tropical increase of all three sampled templates should therefore be interpreted as opposite movements in overlapping floral architecture, not as direct evidence for large-bee loss in the north or bird/butterfly replacement in the tropics. The fact that all three tropical scores can rise together is itself evidence that the templates share plant-side trait structure.

This interpretation preserves a useful biological connection to pollination theory while respecting the causal boundary: floral traits can be concordant with architectures associated with functional pollinator groups, but realized visitor identity, abundance, pollen transfer, and effective service require independent data.

### Trait acquisition reached a practical information ceiling before complete coverage

The progressive campaign also provides a methodological lesson. A trait database does not become scientifically complete when it reaches a particular global percentage. Its value depends on whether added evidence changes the estimand, support tier, uncertainty, or biological conclusion in the strata that matter.

The final campaign achieved a large increase in analysis-usable evidence: 37,771 additional cells relative to Wave52, including more than eleven thousand reproductive-assurance cells. Yet the final H1–H5 claim ceiling remained stable. Meanwhile, the remaining acquisition problem became qualitatively harder. The large stock of previously labelled but nonmaterializable evidence was exhausted; high-yield reproductive sources became rare; later reviewed source packets produced one-to-three-cell gains; and more than 2.25 million committed rows contained no immediately strict-ready reproductive records under the fixed evidence rules.

This is genuine diminishing return. Further searching can still recover individual traits, and targeted additions may matter for poorly supported regions. But indiscriminate global acquisition is no longer the highest-value way to change the Chapter 1 inference. The working freeze therefore rests on recoverability and inference stability rather than on claiming that 69.83% is “enough” in the abstract.

### From WHEN/WHERE to HOW/WHY: the mechanistic handoff

Chapter 1 ends with an identification problem rather than a complete mechanism. It establishes where isolation-associated floral/reproductive assemblage structure is detectable, where response architectures differ, which components are compatible with source/genus assembly, and which interpretations fail the current climate/area/source safeguards. It cannot determine why the same broad island-associated perturbation produces different downstream branches.

That question motivates a separate mechanistic layer. The companion `izu-core` programme formalizes a conditional post-establishment response geometry in which partner loss/arrival balance, plant starting functional state, and the realized pollinator community jointly determine functional matching and effective service. Local interaction filtering then changes the realized response branch, while autonomous assurance alters downstream reproductive consequences. In that framework, syndrome-like shifts can emerge at the assemblage level even when individual lineages occupy opposing response branches.

This is a conceptual handoff, not retroactive validation of Chapter 1. The global regional vectors are not assigned to particular synthetic parameter regimes, and the model does not prove that any observed Chapter 1 gradient was caused by pollinator change. Instead, it supplies a mechanistic class that explains why a common broad interaction perturbation need not map to one floral or reproductive response.

### Returning from mechanism to islands: why Izu is the next empirical scale

The final step is to project the mechanistic uncertainty back onto a linked island system where the missing channels can be measured in the same populations. The Izu programme is suited to this because historical focal-lineage responses, repeated contemporary visitor records, pollinator functional traits, and prospective reproductive experiments can be connected along one island series.

The prospective field design in `izu-core` deliberately separates quantities that Chapter 1 cannot distinguish from phenotype:

1. usable observation effort and visit rate;
2. visitor identity/contact;
3. background-adjusted single-visit pollen deposition;
4. rate-weighted effective pollen service;
5. open-pollinated reproduction;
6. bagged autonomous reproduction;
7. supplemental-outcross reproduction;
8. mature fruit/seed outcomes.

The critical tests are therefore not “does the flower look bee-pollinated?” but: does lower effective service reduce open reproduction where dependency is high; does autonomous assurance buffer that loss; can similar visit rates produce different service because per-visit effectiveness differs; and can reproduction remain high under low service without autonomous assurance? These tests directly target the mechanism that remains unidentified in Chapter 1.

Izu should therefore be viewed as a zoom from global pattern to linked mechanism, not as a cherry-picked confirmation case. A null result would be informative: it would narrow the mechanisms capable of explaining the global assemblage pattern without invalidating the WHEN/WHERE result.

---

## Conclusion

A universal island floral syndrome is not recovered from the final global trait snapshot. Instead, island isolation is associated with biogeographically contingent assemblage structure whose reproductive and floral components can couple or decouple. The strongest broad Palearctic response combines accessibility/generalization and reproductive assurance, but its taxonomic depth is compatible with source-matched genus-level assembly rather than a robust beyond-genus residual. Tropical assemblages can show increased reproductive assurance while maintaining or strengthening specialized floral architecture. Pollination-associated guild templates capture real regional floral-architecture contrasts, but most of their variation is shared and they do not identify realized pollinators.

The progressive wave analysis also shows when global computation should stop being the default next step. Large late gains in usable trait coverage did not alter the major claim ceiling, while remaining recoverable evidence became increasingly source-specific and costly. The main unresolved question is now mechanistic rather than database-limited.

We therefore propose a staged island-syndrome research programme: use global databases to establish WHEN and WHERE assemblage-level components occur; use mechanistic models to determine HOW the same broad interaction perturbation can branch into different responses; return to comparative island systems to test whether those mechanisms survive contact with independent evidence; and finally zoom into Izu to link visitor structure, effective pollen service, reproductive dependency, autonomous assurance, and realized phenotype in the same populations.

---

## Proposed figures

### Figure 1 — Fixed hypothesis ladder and inferential ceiling
`H1 universal syndrome -> H2 biogeographic branching -> H3 source/lineage assembly -> H4 area moderation -> H5 independent channel mechanism`, with the boundary between Chapter 1 plant data and independent pollination evidence.

### Figure 2 — Wave convergence and diminishing returns
Panel A: Wave36, Wave52, and final analysis-usable coverage by axis.  
Panel B: source-reported versus analysis-usable coverage, highlighting the 37,409 Wave52 nonmaterializable Low cells and their later recovery.  
Panel C: realized reproductive gain per source-scale batch (large initial Orchidaceae gain followed by small increments).  
Panel D: H1–H5 decision matrix across Wave36, Wave52, and final snapshot.

### Figure 3 — WHERE: primary plant response branches
Regional/realm response vectors for accessibility/generalization and reproductive assurance, shown separately for all-analysis and direct-only evidence and native/non-endemic strata.

### Figure 4 — Source and taxonomic-depth decomposition
For the Palearctic primary response: observed -> after-family -> after-genus, emphasizing `4/4 -> 4/4 -> 0/4` across evidence scopes/strata and the increase in scored source species.

### Figure 5 — Decomposing the island floral syndrome
Left: reproductive `selfing_core`.  
Middle: attraction/access floral architecture.  
Right: sampled guild-labelled concordance and the shared-factor decomposition.  
Bottom: Palearctic attraction shift persists after selfing-core conditioning; tropical reproductive assurance can coexist with specialized floral architecture.

### Figure 6 — Mechanistic handoff and Izu zoom
Global database (WHEN/WHERE) -> conditional mechanism (`izu-core`, HOW/WHY) -> comparative island projection -> Izu linked field chain (visitation -> SVD -> effective service -> dependency/assurance -> reproduction) -> focal phenotype.

---

## Working literature anchors

- Baker, H. G. 1955. Self-compatibility and establishment after long-distance dispersal. *Evolution* 9:347–349.
- Pannell, J. R. & Barrett, S. C. H. 1998. Baker’s law revisited: reproductive assurance in a metapopulation. *Evolution* 52:657–668.
- Pannell, J. R. 2015. Evolution of the mating system in colonizing plants. *Molecular Ecology* 24:2018–2037.
- Grossenbacher, D. L. et al. 2017. Self-compatibility is over-represented on islands. *New Phytologist*. DOI: 10.1111/nph.14534.
- Sicard, A. & Lenhard, M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. *Annals of Botany* 107:1433–1443.
- Fenster, C. B. et al. 2004. Pollination syndromes and floral specialization. *Annual Review of Ecology, Evolution, and Systematics* 35:375–403.
- Rosas-Guerrero, V. et al. 2014. A quantitative review of pollination syndromes: do floral traits predict effective pollinators? *Ecology Letters* 17:388–400.
- Hetherington-Rauth, M. C. & Johnson, M. T. J. 2020. Floral trait evolution of angiosperms on Pacific islands. *The American Naturalist* 196.

## Reproducibility anchors

- Progressive scientific contract: `chapter1_progressive_analysis_v1`.
- Wave52 reanalysis: Run `33587356935`, artifact `chapter1-progressive-wave52-33587356935`.
- Final trait integration: Run `34191508045`, `222,688 / 318,885` resolved cells.
- Final PR142 progressive reanalysis: Run `34232450884`, artifact `chapter1-progressive-analysis-34232450884`, artifact ID `10058653212`, digest `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
