# Chapter 1 v12 — manuscript-ready Discussion replacement for pollination and two pathways

Status: candidate prose for the exploratory v12 branch. This does **not** overwrite the frozen v11 manuscript.

## Why one floral island syndrome fragments into two plant response components

The broad island pattern is inconsistent with a compulsory serial model in which isolation first increases selfing and floral simplification follows automatically. Reproductive assurance and pollinator-facing floral architecture are statistically separable components of the response. In the defended Palearctic analyses, the attraction/accessibility shift remains associated with isolation after conditioning on the measured `selfing_core` across four frozen source definitions (conditional distance estimates approximately 0.091–0.101; q=0.0079–0.034). The tropical branch provides a complementary counterexample: reproductive assurance can increase while accessibility/generalization does not move in the same simplifying direction. Thus reproductive assurance is one plausible response to unreliable pollen delivery, but it is not sufficient to account for all isolation-associated floral change.

This distinction changes how the classical selfing-syndrome argument should be used. A reproductive-assurance route can reduce dependence on outcross pollen delivery and thereby lower the fitness value of some attraction investments, but that is not the only route available to an island flora. Floral architecture can also respond to the identity, continuity and effectiveness of animal pollination opportunities. The present data therefore support two partially separable plant-side components—reproductive assurance and pollination-associated floral architecture—without requiring one to be downstream of the other.

## Regional floral responses are consistent with pollination syndromes, but syndrome labels are not visitor observations

The predeclared syndrome templates provide a biological interpretation of the regional floral vectors. In northern mid-latitude native assemblages, large-bee-like architecture declines with isolation. In tropical assemblages, butterfly-like architecture increases with isolation, and warm colours are disproportionately coupled to tubular floral forms relative to northern mid-latitudes. These combinations are qualitatively consistent with different pollination-associated architectures being maintained or filtered in different biogeographic settings rather than with one universal movement toward floral simplification.

However, the named templates cannot be treated as realized visitor classifications. More than 86% of the covariance among large-bee-like, butterfly-like and bird-like scores is captured by a shared floral-architecture factor. A flower can therefore score highly on multiple named templates for the same structural reasons. The appropriate inference is that isolation is associated with different **pollination-associated floral architectures** among contexts, not that the analysis has identified which pollinator visited each species or island.

## Experimental pollen limitation provides a positive global ecological signal

The strongest new H5 result comes from GloPL, which measures experimental pollen limitation directly rather than inferring pollination pressure from floral phenotype or visitor occurrence. Across the full georeferenced sampling frame, the analysis included 2,969 experiments from 1,248 sites and 919 publications. Pollen limitation increased with standardized distance from the seeded major continental landmasses (slope=+0.07937; frozen one-sided positive p=0.01772). Thus geographic isolation is associated with a measurable decline in realized reproductive service at global scale.

A post-hoc shape audit suggests that this association is not explained only by a mainland-versus-island level shift. The mainland-to-offshore step was positive but unsupported (+0.13745, p=0.1084), whereas the within-offshore slope was positive (+0.19983, p=0.01110) and the offshore-only slope was also positive (+0.23545, p=0.005865; 260 sites, 158 publications). Because these shape tests were conducted after the parent global result, they remain descriptive rather than confirmatory. Even so, they make a simple binary mainland/island interpretation an incomplete description of the sampled pattern.

This result materially changes the pollination discussion. The question is no longer whether *any* pollination-related ecological pressure covaries with isolation. A general isolation-associated pollen-limitation gradient is supported. The unresolved question is how that common pressure is translated into the different floral, reproductive and taxonomic responses observed among biogeographic contexts.

## A common service pressure does not reproduce the North–Tropical plant branching

The global pollen-limitation result does not by itself explain H2. In the predeclared North–Tropical GloPL comparison, the northern slope was +0.06170 (one-sided p=0.1073), while the Tropical−North distance interaction was +0.14300. The predeclared prediction was a negative interaction, but its one-sided p-value was 0.9119; the implied tropical slope was +0.20470. Thus the point-estimated service gradient is stronger in the tropics rather than in the North.

This mismatch is informative. It rules out the simplest explanation in which the stronger northern floral-island response is merely the consequence of a stronger northern increase in pollen limitation. Instead, similar or even stronger reproductive service limitation can be associated with different trait outcomes in different biogeographic settings. The ecological pressure may therefore be broadly shared while its biological translation depends on lineage composition, pollinator identity and compensation, mating system, or other regional properties.

## Reproductive assurance does not yet provide the functional bridge from pollen limitation to H2

If reproductive assurance buffers the fitness consequences of increasing pollen limitation, species with stronger assurance should show a weaker pollen-limitation increase with isolation. We tested that prediction using exact species overlap between GloPL and the frozen trait database.

The frozen reproductive-assurance family did not pass. Self-compatible species showed an interaction of +0.0351 with distance (one-sided buffering p=0.6790), opposite the buffering prediction. The selfing-mating-system contrast was not evaluable under the frozen support gate because only 11 matched selfing species and one offshore protected-class site were available. Autonomous/delayed selfing showed the expected negative interaction (−0.0736), and both frozen sensitivities retained the negative direction, but the primary p-value was 0.1543. Under the predeclared family rule, 0/3 reproductive-assurance contrasts were supported.

This failure should not be read as evidence that reproductive assurance is biologically irrelevant. The plant-side results still show a distinct reproductive-assurance component, and the autonomous-selfing point estimate is directionally compatible with buffering. The narrower conclusion is that the current exact species-matched public data do not establish reproductive assurance as the functional mediator linking the global pollen-limitation gradient to the Chapter 1 H2 pattern.

## Floral accessibility also does not yet provide the functional bridge

We tested the second route with the three predeclared atomic architecture variables rather than named pollinator templates. If generally accessible flowers buffer increasing pollen limitation, the distance gradient should be weaker in generalized/open architectures after accounting for the frozen model structure.

For `generalized_form`, the restricted-architecture slope was +0.07409 and the generalized slope +0.00735, giving a distance×generalized interaction of −0.06674. This is directionally compatible with buffering, but the one-sided p-value was 0.2441 and the supplemental-only sensitivity reversed sign (+0.06455). For `actinomorphic_symmetry`, the interaction was +0.03496 (p=0.6690), opposite the buffering prediction. `shallow_open_tube` was non-evaluable because exact matched support was too sparse offshore. Under the frozen family rule, 0/3 architecture contrasts were supported.

Again, this is not a biological null for floral architecture. The plant-side floral response remains partly independent of measured selfing, and `generalized_form` remains one of the clearest components after source-free genus residualization in H3A. The GloPL result instead says that the currently measured atomic accessibility traits do not provide a supported species-level buffering bridge from the general pollen-limitation gradient to H2.

## Visitor occurrence explains even less than experimental service limitation

The GloPL result helps reinterpret the earlier occurrence and interaction analyses. Exact-island visitor states are much coarser than experimental pollen limitation. No single functional channel—Bombus, non-Bombus bees, Lepidoptera, flower-visiting birds or Diptera—provided adequate retained-versus-disrupted overlap in both northern-midlatitude and tropical contexts. Pooling the five channels solved that support problem, but total channel disruption still did not predict either higher `selfing_core` or greater accessibility after conditioning on selfing. Identity-aware diagnostics produced a nominal tropical reproductive-assurance signal and the expected sign for matched floral architecture, but neither survived the full inferential gate.

GloBI adds a different cross-examination: interaction breadth varies among contexts in a source-definition-sensitive way, but it does not yield a robust global filter and does not reproduce the plant-side area pattern. Taken together, these results suggest that coarse channel presence, channel count and source interaction breadth are poor substitutes for realized pollination service. GloPL shows that the service outcome itself carries a global isolation signal even when those coarse proxies do not.

## Relation to taxonomic assembly: common pressure, different representation depth

The pollination result becomes most informative when read alongside H3. The broad contemporary-flora North–Tropical response is not erased by source-free genus residualization, whereas the narrower defended Palearctic native response loses support after source-matched genus adjustment.

This creates a useful asymmetry. In the defended Palearctic layer, a general isolation-associated pollen-limitation pressure could in principle contribute through lineage sorting: genera differing in reproductive assurance, floral architecture, demographic sensitivity or dependence on particular pollinators may differ in colonisation, establishment or persistence. The observed genus structuring is compatible with such an assembly-mediated response, but the middle steps are not measured, so pollination cannot be named as the causal generator of H3B.

In the broad contemporary-flora layer, raw genus composition does not exhaust the context difference. If pollen limitation contributes there, it must be translated through finer taxonomic structure, genus-internal species sorting, within-lineage change, introduced/non-native assembly, or another ecological filter. The failed GloPL buffering families show that neither the measured reproductive-assurance traits nor the three atomic accessibility traits currently close that bridge.

The combined interpretation is therefore **common pressure, divergent translation, different taxonomic depth**. Geographic isolation is associated with increasing pollen limitation globally, but plant communities do not all resolve that pressure in the same way. Some responses are strongly represented by lineage assembly; others remain visible after raw genus composition is removed.

## What this adds to the island-syndrome problem

The classical floral island-syndrome narrative often assumes that remote islands experience reduced pollination, which then produces a predictable combination of selfing, reduced attraction and floral generalization. The present results split that narrative into empirically separable pieces.

First, a global ecological premise is supported: pollen limitation increases with geographic isolation. Second, the plant response is not globally uniform: reproductive and floral components branch among biogeographic contexts and can decouple. Third, the taxonomic depth of the response also differs: the broad contemporary pattern persists after raw genus residualization, whereas the defended native Palearctic response is strongly genus-structured. Finally, the two most direct species-level buffering bridges currently available do not explain the connection between the service gradient and H2.

The resulting picture is not the absence of an island syndrome. It is a **context-dependent assembled syndrome** in which a broadly shared ecological pressure is translated through different reproductive, floral and lineage-assembly routes.

## Revised H5 conclusion for the manuscript

> Geographic isolation is associated with increasing experimental pollen limitation at global scale, providing direct evidence that reproductive service becomes more limiting with isolation. Yet this common ecological pressure does not reproduce the North–Tropical plant branching, and exact species-matched tests do not establish functional buffering through either reproductive assurance or the three predeclared atomic generally accessible floral traits. Together with the H3 taxonomic-depth contrast, the results indicate that isolation-associated pollen limitation is translated differently among biogeographic contexts and assemblage layers rather than generating one universal floral island syndrome.

## Claim ceiling

Do not convert this section into any of the following claims:

- the global GloPL slope is a causal island effect;
- pollen limitation proves pollinator abundance or visitation decline;
- GloPL proves historical pollinator loss;
- the global GloPL gradient mediates H2;
- the post-hoc offshore gradient is confirmatory;
- unsupported reproductive-assurance buffering means selfing is irrelevant;
- unsupported architecture buffering means floral architecture is irrelevant;
- support-limited trait contrasts are biological nulls;
- genus structuring proves pollinator filtering;
- post-genus residuals prove within-lineage evolution.
