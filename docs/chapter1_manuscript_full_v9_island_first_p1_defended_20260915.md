# A floral island syndrome can emerge from hierarchical lineage assembly: biogeographic contingency across 8,265 islands

## Full working manuscript v9 — island-first, P1-defended draft — 2026-09-15

### Abstract

Why should island flowers become less conspicuous, more generalized, or more capable of reproductive assurance as islands become isolated? Floral island syndromes are often discussed as repeated adaptive responses of plants to altered mates and pollinators. Yet the same assemblage-level pattern could arise because geographic isolation changes which lineages from a regional source pool arrive, establish, and persist. Distinguishing those possibilities requires asking not only whether traits change with isolation, but also whether responses are shared across biogeographic contexts and where in the lineage hierarchy the response is represented.

We analysed a fixed universe of 8,265 islands and 106,295 accepted angiosperm species under a progressive, predeclared analysis contract. The final trait snapshot resolved 222,688 of 318,885 possible species-by-raw-axis cells (69.83%) across flower colour, floral structural complexity, and reproductive assurance. Primary fitted responses were pollinator-name-free accessibility/generalization and reproductive assurance; mainland distance was treated as a composite source-separation/connectivity gradient rather than as a mechanistically pure treatment.

A universal floral island syndrome was not recovered. Palearctic assemblages became more accessible/generalized and more reproductively assured with increasing separation, whereas tropical assemblages could show increasing reproductive assurance alongside decreasing accessibility/generalization. The broad Palearctic response persisted in native non-endemics and through family adjustment, but the predeclared multivariate gate failed after source-matched genus adjustment in all four evidence-scope-by-floristic-stratum combinations (`4/4 -> 4/4 -> 0/4`). Genus adjustment attenuated 78.8–85.9% of the observed response-vector magnitude. A post-baseline, outcome-blind matched-complexity test then showed that true genus boundaries produced stronger conditional attenuation than arbitrary within-family partitions preserving the exact genus-count and genus-size structure (observed median 0.721; null median 0.224; one-sided randomization `p=0.029`, 2,000 permutations). However, paired spatial-block bootstrap intervals for the *additional* family-to-genus attenuation included zero in all eight direct-only primary profiles, so the exact incremental depth is not precisely localized.

The Palearctic result survived a finite MNAR trait-missingness stress test and a separate species-list detection tipping analysis: 99/100 baseline-supported Palearctic accessibility surfaces remained supported, including all 80/80 scenarios in which remote islands were less complete. By contrast, tropical accessibility was more sensitive. Neither continuous island area, a common nonlinear breakpoint, coarse pollinator-channel heterogeneity, sampled source interaction breadth, nor an independently defined biotic-versus-wind pollination contrast supplied a promoted global mechanism. A prospective distributed-threshold audit further showed that the present macroecological design cannot reliably distinguish lineage-specific thresholds from heterogeneous smooth clines.

We therefore interpret the strongest floral island pattern not as evidence for one repeated adaptive trajectory, but as a **genus-structured assemblage syndrome**: a visible trait pattern whose direction and taxonomic representation depend on biogeographic context. The data support non-random genus-level structure beyond matched arbitrary grouping complexity, but not a precisely estimated taxonomic breakpoint or a causal assembly mechanism. More generally, trait syndromes observed across environmental gradients need not represent repeated organismal adaptation; identifying where they are represented in community lineage structure is a prerequisite for mechanistic interpretation.

**Keywords:** island biogeography; community assembly; floral traits; reproductive assurance; source pool; lineage sorting; taxonomic depth; trait syndrome; biogeographic contingency

---

## Introduction

Why are island flowers often expected to become duller, less specialized, or more self-reliant? The question is appealing because islands alter several ecological conditions at once. Colonists cross geographic barriers, mate availability may become unreliable, pollinator assemblages may differ from those of source regions, and small or isolated populations may experience strong demographic filtering. These observations have motivated the idea of a floral island syndrome in which increasing isolation favours reproductive assurance while reducing investment in specialized or conspicuous pollination-associated floral architecture.

That interpretation, however, contains two logically distinct claims. One is about **organisms**: lineages arriving on islands repeatedly change their floral or reproductive traits after colonization. The other is about **assemblages**: isolation changes which already-different lineages from the source pool can arrive, establish, and persist. Both processes can produce the same cross-sectional association between island isolation and mean trait composition. A distant island with more self-compatible or more generalized flowers need not contain lineages that evolved those traits on the island; it may instead contain a different subset of source-available lineages.

This distinction matters particularly for floral traits because reproductive assurance and pollination-associated architecture need not form one obligatory causal sequence. Baker's law and later work on breeding systems motivate a colonization filter favouring uniparental reproduction under mate or pollinator limitation. Separately, altered interaction environments may change the value of restricted floral access, deep tubes, bilateral symmetry, conspicuous colours, or large floral displays. Self-compatibility is not autonomous selfing, autonomous selfing is not necessarily the realized mating system, and reproductive assurance can increase without floral simplification. Conversely, floral architecture can change without a detectable shift in reproductive assurance. A community-level syndrome may therefore be a composite of partially independent biological channels.

A second complication is that a geographic barrier can filter both sides of a plant–pollinator interaction. Plants differ in propagule dispersal, establishment, habitat requirements, and persistence. Pollination channels differ in flight capacity, long-distance dispersal, stepping-stone use, nesting requirements, and persistence. In principle, isolation can therefore affect plant lineage assembly directly, alter interaction continuity, or do both. But floral phenotype alone cannot reveal which pollinator channel was historically lost or retained. A credible global analysis must first identify the plant-side response and its taxonomic representation, then ask whether independent interaction data add information beyond that structure.

A third complication is biogeographic context. If the same isolation gradient acts on source pools that differ in lineage composition, climate, interaction networks, and historical connectivity, a universal floral syndrome need not exist even when strong regional responses do. Comparing whether individual coefficients are significant in one region and not another is insufficient; the response vectors themselves must be compared directly. Likewise, a regional pattern cannot be interpreted as repeated adaptation until its dependence on lineage composition is known.

We therefore recast the island-flower question as a sequential inferential problem. **First**, does increasing source separation generate one coherent floral/reproductive syndrome, or do response vectors branch among biogeographic contexts? **Second**, when a branch is supported, how strongly is it represented by family and genus composition rather than by a residual beyond genus? **Third**, is the apparent genus localization biologically structured, or could it arise mechanically from sample loss or from inserting any sufficiently fine grouping factor? **Fourth**, can simpler alternatives such as area, a common nonlinear breakpoint, or observation bias explain the result? **Fifth**, after those layers are fixed, does independent pollination information identify a mechanism? This inference sequence is summarized in Fig. 1.

We implemented the primary sequence under a progressive analysis contract. The geographic universe, trait ontology, hypothesis order, support gates, multiple-testing families, source-pool safeguards, and claim ceilings were fixed while trait evidence improved. This allowed evidence coverage to increase without allowing the biological question to drift with each new result. The final database contains three raw measurement axes—flower colour, floral structural complexity, and reproductive assurance—but the primary hypothesis test does not collapse these into three arbitrary scalar responses. Instead, it uses a pollinator-name-free primary response comprising **accessibility/generalization** and **reproductive assurance**, with broader colour, structure, and self-compatibility contrasts and named floral templates treated as secondary diagnostics.

The hypothesis ladder was as follows. **H1** treated a universal floral island syndrome as a rival hypothesis. **H2** asked whether the primary response branches among biogeographic contexts. **H3** asked how a supported branch changes after source-matched family and genus composition are accounted for. **H4** asked whether island area modifies the isolation response strongly enough to support a capacity or founder interpretation. **H5** was reserved for mechanism: independent pollination information had to explain either lineage filtering itself or a residual response that survived the earlier composition controls.

The central prediction was not that distant islands should always be duller, more generalized, or more selfing. It was that, if floral island syndromes are emergent assemblage properties, their **direction, component coupling, and taxonomic representation should vary among biogeographic contexts**. A stronger assembly prediction is that biologically real taxonomic partitions should absorb the response more strongly than arbitrary partitions of matched complexity. Conversely, a repeated beyond-genus interpretation would be more plausible if a robust response persisted after genus composition was accounted for.

This framing creates a route from an island-specific question to a general macroecological problem. Repeated trait syndromes are often interpreted as repeated adaptation across urban, alpine, arid, fragmented, or otherwise filtered environments. But cross-sectional community data mix organismal change with changes in community membership. Islands allow us to ask explicitly when a visible syndrome represents a response of lineages and when it represents **which lineages are present**. The general inference problem is reached from the island question, rather than imposed on it in advance.

---

## Materials and Methods

### Geographic universe and exposure

The fixed geographic universe contained 8,265 islands. Species-level analyses used a fixed denominator of 106,295 accepted angiosperm species. The denominator, island universe, floristic-status information, source-pool infrastructure, island area, climate covariates, and primary geographic exposure were held constant while trait evidence was allowed to improve.

The primary exposure was log-transformed distance to the nearest continental source boundary. We interpret this variable as a composite gradient of source separation, connectivity, and accessibility, not as a physically pure ocean-crossing treatment and not as a historical-source assignment for every island. Consequently, any breakpoint on this axis would be an assemblage breakpoint in composite geographic separation, not a literal physiological ocean-crossing threshold.

Analyses were repeated across all-analysis-eligible and direct-only evidence scopes, and across all-native and native-nonendemic floristic strata where support permitted. Endemic analyses were retained as secondary diagnostics but were not interpreted as a time axis or as direct evidence of in-situ evolution.

### Trait evidence and response construction

Each accepted species could contribute evidence to three raw measurement axes: flower colour, floral structural complexity, and reproductive assurance. The final snapshot resolved 222,688 of 318,885 possible species-by-axis cells (69.83%). Missing evidence was never converted to trait absence. Species-direct high/medium evidence took precedence over validated lower-confidence evidence, and family-level imputation or unrestricted fallback was prohibited in the primary analysis.

The primary H2 response was deliberately pollinator-name-free. `accessibility_generalization` represented increasingly open, shallow, or radially accessible floral architecture. `reproductive_assurance` was derived from strict reproductive evidence. Broader atomic contrasts (`plain_colour`, `generalized_form`, and `self_compatibility`) were used in secondary response-geometry analyses. Named large-bee-like, butterfly-like, and bird-like templates were treated only as floral-architecture concordances, never as observations of realized pollinator identity.

### Progressive analysis and support gates

The analysis was prospectively ordered. Multivariate response-vector tests preceded individual-axis interpretation. Direct between-context comparisons were required before claiming biogeographic heterogeneity; a significant coefficient in one region and a non-significant coefficient in another was not sufficient. Response-specific island support thresholds and BH-FDR families were fixed before each outcome opening. Evidence-scope replication and native-nonendemic persistence were treated as robustness dimensions rather than as opportunities to select the most favourable result.

Post-baseline P1 robustness analyses were explicitly labelled as such. Their purpose was to attack the already-known H3 interpretation, not to replace historical gates or to search for a more favourable endpoint. The P1a, P1c and P1d designs were frozen before their respective new outputs were inspected, and failure was required to narrow the manuscript claim rather than trigger retuning.

### H1 and H2: universal versus branching response

H1 asked whether source separation generated one coherent floral/reproductive response across contexts. H2 asked whether the response vector branched among predeclared biogeographic contexts. The principal H2 comparison separated accessibility/generalization from reproductive assurance, allowing the two components to move together, decouple, or move in opposing directions.

To distinguish reproductive assurance from a broader floral-architecture response, `selfing_core` included self-incompatibility, mating-system, and autonomous-selfing evidence but excluded attraction traits. Conditional attraction/access models then asked whether the floral-architecture association with distance remained after conditioning on measured reproductive assurance. These models were interpreted as decomposition, not as causal mediation.

### H3: source-matched taxonomic structure and P1 defense

H3 treated taxonomic composition as a biological level of response rather than a nuisance correction. Outcome-blind source expectations were constructed from fixed source-flora assignments. The broad response was evaluated at three stages: observed assemblage response, residual after family composition, and residual after source-matched genus composition. The historical H3 estimands were the attenuation in response-vector magnitude between these stages and the predeclared support status at each depth.

Because genus adjustment could appear strong for purely analytical reasons, P1 attacked three rival explanations. **P1a** verified exact paired support: observed, family-adjusted and genus-adjusted stages had to use the same island–species observations and the same information weights. **P1c** tested matched grouping complexity. Within every family, species were randomly assigned to pseudo-genera while preserving the exact number of genera and the exact multiset of genus group sizes. The primary statistic was the median conditional attenuation from family to genus across the eight direct-only Palearctic source-mode-by-stratum profiles. Two thousand deterministic permutations were generated from a frozen seed schedule; the one-sided randomization p-value was `(1 + number(null >= observed)) / (1 + valid permutations)`. **P1d** propagated spatial uncertainty using 2,000 paired spatial-block bootstrap draws on the direct-only profiles, estimating uncertainty for total genus attenuation and for the additional family-to-genus attenuation.

P1c and P1d answer different questions. P1c asks whether *true genus membership* contains information beyond arbitrary fine grouping of identical complexity. P1d asks how precisely the *incremental amount* of family-to-genus attenuation is localized across spatial blocks. A positive P1c result therefore cannot override an imprecise P1d increment, and neither test identifies a causal assembly process.

### H4: continuous area moderation

H4 tested continuous distance-by-area interactions. No post hoc small/large island cutoff was introduced. Equal-island, capped-information, direct-only, common-support, and heteroskedastic-null safeguards determined whether an area interaction could be promoted beyond a measurement-sensitive modifier.

### Response geometry

Because a smooth assemblage trend can conceal local nonlinear responses, we prospectively evaluated five response geometries: flat, cline, step, hinge, and reversal. Initial naive information-criterion selection produced excessive false nonlinear calls under spatially clustered cline simulations. We therefore calibrated a nonlinear evidence statistic separately within each realized design cell and required held-out false-nonlinear rates <=0.10 before opening observed response geometry. Step detection was identifiable in nearly all design cells at the predeclared target effect, whereas hinge responses were not.

The observed geometry test used the frozen calibrated thresholds and required agreement across all-analysis and direct-only scopes before promoting an identified shape. Because distance is a composite exposure, even an identified step would not have been interpreted as a physical ocean-crossing threshold.

### Missingness and observation-process stress tests

Two distinct missingness problems were treated separately. **V5** varied trait-state-dependent resolution using a finite MNAR tipping grid. This asked how strongly trait resolution would have to depend on the unobserved state before the primary result was erased. **V6** targeted species-list incompleteness. It jointly varied distance-dependent flora-list completeness and trait-state-dependent recording odds while preserving the original information weights; hypothetical missing species did not increase regression precision.

V6 was designed as a tipping analysis rather than an occupancy estimate. It therefore does not estimate the true completeness of GBIF-derived or compiled island floras. It asks which combinations of distance-dependent under-survey and state-dependent recording would be required to erase the frozen biological results.

### H5: independent mechanistic tests and claim ceiling

The primary plant database cannot identify historical pollinator loss from floral phenotype. We therefore treated H5 as a set of independently gated tests rather than as an interpretation attached to H3.

The first global test (N1) evaluated whether isolation-associated occurrence deficits differed among independently measured pollinator channels after effort correction. A second source-side analysis used effort-matched sampled interaction breadth. Neither was allowed to rescue the other after a failed gate.

We then added two prospective tests after the primary manuscript architecture had been frozen. **H5c** used external GIFT `pollen_vector_mode` to compare the isolation response of biotically pollinated and wind-pollinated plants. Pollination mode was never inferred from the floral traits used as the response. Support and power were simulated on the realized island, species-support, covariate, and spatial-block design before the observed interaction was opened; only cells passing the predeclared false-positive and recovery gates could be inspected. **H5d** asked whether the realized assemblage design could distinguish a generator with lineage-specific thresholds from a generator with heterogeneous smooth clines. Observed genus-level threshold distributions were prohibited unless identifiability gates passed.

### Computational stopping rule

Trait acquisition stopped because the inference stabilized, not because a target percentage was reached. Large additions of evidence no longer changed the principal H1–H5 claim structure, later reviewed sources yielded sharply diminishing numbers of strict-ready cells, and a large residual audit found no immediately promotable reproductive evidence under the fixed rules.

---

## Results

### H1: increasing isolation did not produce one universal floral island syndrome

The primary multivariate response did not recur in one direction across biogeographic contexts. The principal result was therefore not a collection of region-specific significance labels but a direct failure of the universal-response model. Source separation was associated with different combinations of reproductive assurance and floral accessibility in different contexts (Fig. 2A,B).

### H2: reproductive assurance and floral architecture branched by biogeographic context

The Palearctic showed the clearest broad primary response. In all-analysis all-native assemblages, accessibility/generalization increased with source separation (slope 0.0795, FDR `q=0.00463`), as did reproductive assurance (0.0534, `q=0.00463`). Direct-only estimates pointed in the same broad direction: 0.0616 for accessibility/generalization (`q=0.0389`) and 0.0966 for reproductive assurance (`q=8.48 x 10^-14`). The response persisted in native non-endemics, showing that it was not confined to island endemics.

Tropical assemblages did not follow a weaker copy of the Palearctic branch. In direct-only all-native data, accessibility/generalization declined with distance (-0.1005, `q=0.00275`) while reproductive assurance increased (0.1360, `q=0.00428`). Thus increasing reproductive assurance did not imply increasing floral generalization. The two response components can be decoupled at macroecological scale (Fig. 2A,B).

The Palearctic attraction/access component also remained positive after conditioning on `selfing_core` across four fixed source definitions. In all-native assemblages, conditional distance estimates were approximately 0.091–0.101 with `q=0.0079–0.034`. This further argues against a compulsory serial model in which isolation first increases selfing and floral architecture then changes only as a secondary selfing syndrome (Fig. 2C,D).

Secondary large-bee-like, butterfly-like, and bird-like templates differed among regions, but a source-trained common factor explained 86.85% of their variance in all-analysis evidence and 86.44% in direct-only evidence. These named templates therefore mainly captured overlapping plant architecture. They were not treated as evidence for realized visitor identity.

### H3: the Palearctic floral-island response is genus-structured beyond matched grouping complexity

The historical taxonomic-depth decomposition showed that the broad Palearctic response passed the predeclared gate at the observed assemblage stage, remained supported after family adjustment, and failed after source-matched genus adjustment in all four evidence-scope-by-floristic-stratum combinations: `4/4 -> 4/4 -> 0/4`. Across the broader frozen profiles, family adjustment removed approximately 19.6–33.4% of the observed vector magnitude, whereas genus adjustment removed 78.8–85.9%.

P1a showed that this attenuation was not produced by stage-specific sample loss. Observed, family-adjusted and genus-adjusted stages used the same frozen island–species support and the same `n_species` information weights in all focal audit cells; family and genus entered as grouping variables rather than trait-imputation devices.

P1c then tested whether the genus result was merely a consequence of inserting a finer categorical structure. True genera were compared with 2,000 pseudo-genus partitions generated within family while preserving the exact genus count and the exact genus group-size multiset of every family. The true-genus median conditional attenuation was 0.7207, compared with a null median of 0.2243. Fifty-seven of 2,000 pseudo-genus permutations equalled or exceeded the observed statistic, giving a one-sided randomization `p=0.02899` (Fig. 3B). Thus the biological membership of real genera contains information relevant to attenuation beyond generic fine-grouping complexity.

P1d supplied the counterweight. In eight direct-only primary profiles, observed total genus attenuation remained large (0.788–0.802), and the lower bounds of the paired spatial-block bootstrap intervals for total genus attenuation remained positive (0.277–0.402). But the additional family-to-genus attenuation was not precisely localized: all eight 95% intervals for that incremental quantity included zero (Fig. 3C). We therefore do not claim a precisely estimated taxonomic breakpoint at the family-to-genus transition.

The combined inference is narrower and stronger than the original H3 wording. The Palearctic floral-island response is **strongly genus-structured**, and real genus boundaries outperform arbitrary within-family partitions of matched complexity. At the same time, the exact amount of attenuation attributable specifically to the family-to-genus increment is spatially imprecise. Genus structure localizes the response in taxonomy; it does not by itself identify dispersal, establishment, habitat filtering, persistence, mutualist dependence, or evolution as the causal process (Fig. 3A–D).

A separate hierarchical-depth audit showed that taxonomic attenuation differs among contexts within widespread native non-endemic assemblages. This supports the broader claim that biogeographic contingency is not limited to the sign of a trait coefficient; the degree to which broad responses are represented in lineage composition can also differ among contexts.

### H4: island area did not yield a promotable capacity mechanism

Across 16 predeclared area-sensitivity cells, none passed the full promotion rule for a directional area mechanism. Heteroskedastic-null safeguards also failed to establish a robust mechanistic interaction. Apparent small-island amplification therefore remained a measurement-sensitive modifier rather than evidence for a founder, habitat-capacity, or pollinator-persistence mechanism (Fig. 1C).

### A common nonlinear isolation threshold was not recovered

The response-geometry audit first showed why an unrestricted search for breakpoints would be misleading: spatially clustered cline simulations frequently selected nonlinear models. After calibration, held-out false nonlinear rates were controlled below 0.10. A true step of the target magnitude was recoverable in 11/12 cross-scope design cells.

When the observed outcomes were finally opened under those frozen gates, all 12 broad atomic response cells were classified as `monotonic_or_unresolved`; no step, hinge, or reversal was promoted. The global assemblage evidence therefore does not support one common isolation breakpoint. This result does not rule out lineage-specific or local-system thresholds (Fig. 4B).

### V5 and V6 constrained two different observation-bias explanations

The V5 MNAR analysis showed that the broad Palearctic response was not dependent on assuming trait missingness at random. The primary vector survived the finite predeclared state-dependent-resolution grid, although reproductive details and some regional contrasts remained less robust at the most extreme boundaries. Arbitrary missingness mechanisms were not ruled out.

V6 addressed the separate problem of species-list incompleteness. Baseline `OR_D=1` scenarios reproduced the frozen H2 results to numerical precision. Among baseline-supported Palearctic accessibility surfaces, 99/100 survived the complete frozen bias grid. More importantly, all 80/80 scenarios in the biologically concerning direction—remote islands less complete, combined with under-recording of the focal generalized/accessibility state—retained support. The specified detection bias therefore acts as a poor explanation for the positive Palearctic pattern.

The tropical negative accessibility component was more fragile: 40/75 baseline-supported surfaces crossed a tipping boundary under combinations of distance-dependent incompleteness and state-dependent recording. The North–Tropical multivariate contrast was more robust, breaking in only 5/75 baseline-supported surfaces, all under the strongest completeness decline. We therefore retain the broad biogeographic branching claim while treating the tropical single-axis accessibility result as more observation-sensitive than the Palearctic branch (Fig. 4A).

### H5: independent global pollination tests did not identify the upstream mechanism

The first independent channel analysis found no promoted isolation-by-channel heterogeneity (`W=1.6187`, df=3, `p=0.65516`). Simulation showed limited power for modest true differences, so this is non-identification rather than evidence that channels are equivalent. An effort-matched GloBI source-breadth analysis likewise promoted no context-by-stratum cell (Fig. 1C).

H5c provided a stronger negative-control test because pollination mode came from an external source and was not inferred from floral phenotype. GIFT v3.2 yielded 5,771 unambiguous mode assignments after exclusions (4,538 biotic; 1,233 wind). Before outcome inspection, only one of eight planned design cells met both support and power gates: direct-only, native-nonendemic, Palearctic. In that cell, the observed distance-by-biotic interaction was +0.06495 (95% CI -0.09030 to 0.22020; `p=0.41221`). The wind slope was -0.05393 (95% CI -0.20672 to 0.09886), and the biotic slope was +0.01102 (95% CI -0.01042 to 0.03246). We therefore found no independent evidence that the Palearctic accessibility response is stronger among biotically pollinated plants than among wind-pollinated plants (Fig. 4C).

H5d tested the attractive alternative that smooth global gradients could be generated by many lineage-specific thresholds. The idea is biologically plausible, but the current design could not identify it. Across eight source-mode-by-stratum design cells, 0/8 passed the prospective gate. Classification accuracy ranged from approximately 0.733 to 0.778, while smooth heterogeneous clines were falsely classified as distributed-threshold generators 19.0–25.5% of the time. Observed genus-level threshold distributions therefore remained unopened (Fig. 4D).

Taken together, H5c and H5d sharpen rather than fill the causal gap. The strongest plant-side pattern is real and strongly genus-structured, but the present global data do not identify a pollinator-specific mechanism or a threshold generator that can be distinguished from heterogeneous smooth alternatives.

### Trait acquisition reached an inference-based stopping point

The final trait snapshot contained 222,688 resolved cells, but the stopping decision was not based on 69.83% as a universal adequacy threshold. The principal H1–H5 claim structure had stabilized while later source review produced sharply diminishing strict-ready yield. The evidence campaign was therefore frozen when additional acquisition ceased to change the main inference more than it changed nominal coverage.

---

## Discussion

### From “why are island flowers dull?” to “where is the syndrome represented?”

The original question behind this analysis was simple: why should island flowers become duller, more generalized, or more self-reliant with isolation? The data do not support answering that question with one adaptive story. Instead, they change the object of explanation. The first question is no longer *which mechanism makes island flowers dull?* but *at what biological level is the apparent syndrome represented?*

Three results drive that shift. First, the response is not globally uniform: Palearctic and tropical assemblages can combine reproductive and floral changes differently. Second, the strongest Palearctic response is highly sensitive to real genus composition, and that genus structure is stronger than arbitrary within-family grouping of matched complexity. Third, several simple mechanistic stories fail explicit tests: continuous area does not promote, a common breakpoint is not recovered despite adequate step-detection ability, coarse pollinator-channel heterogeneity is not supported, sampled interaction breadth does not promote, and the independent biotic-versus-wind contrast is null in the only cell qualified for inspection (Figs. 2–4).

The parsimonious description is therefore a **genus-structured assemblage syndrome**. The visible floral/reproductive gradient is an assemblage property whose components can be generated by changing community membership, changing traits within represented lineages, or both. In the strongest Palearctic branch, real genus membership contains substantial information about the broad response. What P1 does *not* establish is a sharply measured family-to-genus breakpoint: spatial-block uncertainty around the incremental family-to-genus attenuation is wide.

### Taxonomic localization is an estimand, not merely a correction

Community-level trait studies often control for taxonomy to remove non-independence. Here the attenuation trajectory itself is biologically informative, but its interpretation must separate localization from precision. The historical decomposition shows that the broad response remains after family adjustment and is greatly reduced after genus adjustment. P1a shows that this is not a sample-loss artifact. P1c shows that true genus membership absorbs more response than arbitrary within-family partitions with the same grouping complexity. P1d shows that the exact *incremental amount* assigned to the family-to-genus step is not precisely estimated across spatial blocks.

Together these results make taxonomic localization an estimand without pretending taxonomy is a causal mechanism. The useful question is **where in lineage structure does the syndrome cease to behave like a robust beyond-group response?** Real genus boundaries matter for that answer. But the data do not justify a literal claim that a fixed percentage of the syndrome is generated at one sharply resolved taxonomic transition.

This framing avoids two opposite errors. One is to interpret every assemblage-level trait gradient as repeated adaptation. The other is to treat taxonomy as an inconvenient confounder that should simply be regressed away. Lineage composition is itself a potential outcome of geographic filtering. A genus-structured signal can reflect propagule movement, habitat compatibility, establishment, persistence, mutualist dependence, historical source structure, within-lineage change, or combinations of these processes. The correct inference is not that “taxonomy explains the biology,” but that the biology is substantially represented in non-random lineage composition at genus scale.

### A syndrome can be assembled without being a single syndrome

The decoupling of reproductive assurance and floral accessibility is equally important. In the Palearctic, both components increase with separation. In tropical assemblages, reproductive assurance can increase while floral accessibility decreases. This means the classical verbal package—less reliable pollination, more selfing, simpler flowers—cannot be assumed to operate as one serial pathway globally.

The source-trained pollination-template factor reinforces this point. Large-bee-like, butterfly-like, and bird-like labels share most of their variance because they reuse overlapping floral architecture. Without independent visitor or service data, naming a template after a pollinator does not identify a pollinator mechanism. The primary plant response is therefore intentionally pollinator-name-free.

### Robustness is asymmetric across contexts

The observation-bias analyses also change how strongly different parts of the story should be stated. The positive Palearctic accessibility branch is difficult to erase using the specified trait-resolution or species-detection mechanisms, and the direction in which remote floras are less complete is especially conservative for that result. By contrast, the tropical negative accessibility component is materially more sensitive to species-list incompleteness. The tropical result remains useful because the multivariate North–Tropical contrast is more robust than the single axis, but the two branches should not be presented as equally defended.

This asymmetry is a strength rather than a defect of the analysis. A robustness framework should reveal where claims are weak, not make every result look equally certain. The court-style logic therefore distinguishes a strongly defended Palearctic, genus-structured assemblage response from a more sensitivity-bounded tropical component.

### Why the global pollinator mechanism remains unidentified

It is biologically plausible that altered pollination contributes to lineage sorting on islands. Indeed, an interaction filter could act upstream of H3: if a lineage depends on a service that fails to persist on remote islands, differential lineage persistence would appear as genus structure. But plausibility is not identification.

The independent H5 tests did not promote a global pollinator mechanism. N1 did not establish channel heterogeneity. Sampled GloBI breadth did not predict the expected source-side filtering pattern. Most importantly, the independent biotic-versus-wind negative control did not show a stronger Palearctic accessibility response among biotically pollinated plants in the only prospectively qualified cell. These results do not demonstrate that pollinators are irrelevant. They show that the present global proxies are insufficient to name pollination as the causal generator of H3.

This distinction protects the main result from overreach. Chapter 1 identifies **where the syndrome is represented** and which simple explanations fail. It does not infer an unmeasured historical interaction mechanism merely because the traits are floral.

### Smooth global gradients do not resolve local response geometry

The distributed-threshold audit provides a second important boundary. Local ecological systems can respond nonlinearly even when a global assemblage average looks smooth. If different lineages or interactions have different thresholds, averaging them can produce an apparently gradual macroecological gradient. But the converse problem is equally serious: heterogeneous smooth clines can mimic a distribution of thresholds.

Our prospective simulation showed that the current assemblage design cannot reliably distinguish those generators. We therefore do not fit and interpret observed genus-specific threshold distributions in Chapter 1. This negative result creates a clear division of labour with local mechanistic work: threshold-versus-cline questions require systems in which interaction state, effective service, reproduction, and phenotype can be measured within the same biological chain.

### From the floral island syndrome to a broader problem of assembled trait syndromes

The broader implication is reached only after following the island problem to its limit. The classical island question asks why flowers should become more self-reliant or generalized as isolation increases. Our results first show that this expectation is not globally uniform, then that its components can decouple, and finally that the strongest regional branch is strongly represented in real genus composition. The general issue emerges from that sequence: a cross-sectional trait syndrome can look like repeated organismal adaptation even when much of its visible structure is carried by which lineages make up the assemblage.

The same inference problem can arise in urban, alpine, arid, fragmented, disturbed, or otherwise environmentally filtered communities. We therefore suggest a general diagnostic sequence for cross-sectional trait syndromes: **Which components move? In what direction? Does the response branch among contexts? How much of it is represented by lineage composition? Do biologically real taxonomic partitions outperform matched arbitrary partitions? Is the apparent response geometry identifiable? Which observation and process alternatives survive explicit falsification?** Only after those questions are answered should a named mechanism be promoted.

In this view, the contribution of the island case is not simply that some flowers become more generalized or more reproductively assured. It is that a long-standing island syndrome, pursued as a biological question, exposes a general ambiguity between repeated change of organisms and reassembly of lineages. The data support a genus-structured response, while also showing why `assembly depth` should be treated as a localization concept with uncertainty rather than as a perfectly sharp taxonomic breakpoint.

### Chapter 1 and Chapter 2 resolve different scales of the same problem

The global analysis mixes many lineages and many island histories. Its strength is breadth: it can identify component direction, biogeographic branching, taxonomic structuring, and failures of universal explanations. Its weakness is mechanism: it cannot reconstruct realized interaction service for each lineage.

A deeply resolved island system can do the opposite. Within a single biological system, interaction regime, effective pollination service, reproductive outcome, and phenotype can be measured along the same geographic sequence. Such a system can test whether one channel changes gradually while another crosses a threshold, whether multiple responses share a breakpoint, and whether a floral transition follows rather than merely covaries with a change in service.

The two scales are therefore complementary. **Chapter 1 asks where and at what lineage-assembly level the syndrome is represented. Chapter 2 asks how a local interaction–reproduction–phenotype chain generates a particular response geometry.** A local threshold and a smooth global assemblage gradient are logically compatible, but the present study demonstrates that the global data alone cannot identify a distributed-threshold generator against heterogeneous smooth alternatives (Fig. 1D).

### Conclusion

Geographic isolation does not impose one floral island syndrome. It reorganizes floral and reproductive composition in ways that depend on biogeographic context, and the strongest Palearctic response is strongly structured by real genus membership. That genus structure is not explained by stage-specific sample loss and exceeds matched arbitrary within-family grouping of the same complexity, although the precise incremental family-to-genus attenuation remains spatially imprecise. The result also survives substantial trait-missingness and species-detection challenges, while several simpler explanations fail prospective tests. At the same time, the global data do not identify a pollinator-specific causal mechanism and do not justify interpreting a smooth assemblage gradient as evidence for or against local thresholds.

The island-flower question therefore leads to a broader inference problem: **when does a macroecological trait syndrome represent repeated adaptation of organisms, and when does it emerge from which lineages can assemble?** Our results show that taxonomic localization should be established—and its uncertainty exposed—before mechanism is assigned to a community-level syndrome.

---

## Main figure legends

**Figure 1 | From an island floral syndrome to a hierarchy-of-assembly test.** An isolation-associated shift in floral or reproductive composition can arise through repeated within-lineage response, differential assembly of source-available lineages, or both (A). The global analysis rejects one universal floral/reproductive response, identifies biogeographic branching, and then asks how strongly the Palearctic branch is represented in taxonomic composition (B). The broad Palearctic response remains supported after family adjustment but fails the predeclared vector gate after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`); P1 subsequently tests whether this genus sensitivity is a sample-loss or generic grouping artifact. Independent robustness and falsification analyses define the claim boundary rather than supplying a post hoc mechanism (C): the Palearctic branch survives strong specified missingness challenges, while area, a common global nonlinear shape, coarse pollinator-channel heterogeneity, sampled interaction breadth and independent biotic-versus-wind specificity do not yield a promoted global mechanism; distributed thresholds are not identifiable against heterogeneous smooth clines in the present design. Chapter 1 therefore identifies response direction, biogeographic context and taxonomic structuring without identifying a pollinator-specific causal chain. A deeply resolved local system is required to distinguish interaction, effective-service, reproductive and phenotypic response geometries (D). Dashed connections indicate mechanistic links not identified by the global data.

**Figure 2 | Floral and reproductive responses branch among biogeographic contexts.** Primary response vectors are shown as the isolation-associated slopes of accessibility/generalization and reproductive assurance in all-native (A) and native-nonendemic (B) assemblages. Open and filled symbols distinguish all-analysis and direct-only evidence rather than separate hypotheses. The Palearctic branch is reproduced across evidence scopes and floristic strata (C), whereas tropical assemblages can combine increasing reproductive assurance with decreasing accessibility/generalization, so the tropical pattern is not a weaker scalar copy of the Palearctic response. The Palearctic floral-architecture association also remains positive after conditioning on the measured `selfing_core` across four frozen source definitions (D), showing component decoupling rather than causal mediation. The figure visualizes frozen H1/H2 estimates; direct between-context multivariate tests remain the inferential basis for biogeographic heterogeneity.

**Figure 3 | The Palearctic floral-island response is genus-structured beyond matched grouping complexity, but the exact taxonomic increment is imprecise.** Frozen taxonomic-depth trajectories show strong attenuation after source-matched genus adjustment across the Palearctic profiles (A). A matched-complexity randomization then compares true genera with pseudo-genera formed within family while preserving the exact genus-count and genus-size multiset: true genus membership produces a median conditional attenuation of `0.721`, compared with a pseudo-genus null median of `0.224`; 57/2,000 null permutations equal or exceed the true-genus statistic (`p=0.02899`; B). Paired spatial-block bootstrap intervals provide the counterweight: although total genus attenuation remains large, all eight direct-only 95% intervals for the *additional* family-to-genus attenuation include zero (C). The integrated decision is therefore genus-specific taxonomic structure without a precisely estimated family-to-genus breakpoint or a causal assembly mechanism (D).

**Figure 4 | Cross-examination narrows the interpretation without erasing the plant-side result.** Species-list detection sensitivity is asymmetric across contexts (A): the Palearctic accessibility branch survives `99/100` baseline-supported bias surfaces, whereas the tropical accessibility branch survives `35/75`; the direct North–Tropical multivariate contrast is more robust (`70/75`). Calibrated response-geometry tests require nonlinear evidence to exceed a design-specific threshold in both all-analysis and direct-only scopes; none of the 12 observed cells meets that promotion rule (B). An independent biotic-versus-wind specificity test in the sole prospectively qualified cell yields a positive but imprecise interaction (`+0.06495`, 95% CI `-0.09030` to `0.22020`, `p=0.41221`) and does not identify a pollinator-specific global mechanism (C). Finally, distributed lineage thresholds cannot be distinguished reliably from heterogeneous smooth clines in the realized global design (`0/8` qualified; D). Negative or non-identified results are not inverted into evidence that pollinators, local thresholds or other mechanisms are irrelevant.

---

## Claim ceiling for v9

The manuscript may state that the strongest Palearctic floral/reproductive syndrome is strongly structured by source-matched genus composition, that true genus boundaries attenuate the response more than matched arbitrary within-family partitions of identical grouping complexity, and that trait-syndrome direction and taxonomic representation depend on biogeographic context. It may state that several simple alternatives fail explicit falsification tests.

It must also state that the exact incremental family-to-genus attenuation is spatially imprecise. It must not state that the family-to-genus transition is a precisely estimated taxonomic breakpoint, that pollinator loss caused the Palearctic response, that genus attenuation proves dispersal alone, that the H5c null proves pollinators do not matter, that tropical accessibility is as observation-robust as the Palearctic branch, that endemicity establishes within-lineage evolution, or that smooth global responses imply the absence of local thresholds.

## Literature cited

Full reference formatting remains inherited from the v7 canonical bibliography and the frozen literature-positioning note; no new citation should be added to the submission version without source verification.
