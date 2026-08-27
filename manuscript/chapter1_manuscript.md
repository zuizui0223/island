# Biogeographic contingency and island area structure floral and reproductive assemblage responses to mainland isolation

## Abstract

Island isolation can filter both the plants available to assemble a flora and the
ecological interactions that permit those plants to persist. A universal floral island
syndrome is therefore unlikely if source pools and interaction channels differ among
biogeographic regions. We tested when and where floral and reproductive assemblage
composition changes along a mainland-distance/source-accessibility gradient in a frozen
global universe of 8,265 islands. Raw occurrence multiplicity was collapsed to island x
species incidence, missing observations were not coded as biological absence, and primary
inference used multivariate response-vector tests with explicit support thresholds,
floristic-status strata, observation-process sensitivities, and spatial-block robust
uncertainty. Northern mid-latitude and tropical assemblages both changed along the
geographic gradient, but their response vectors differed in all-native and native
non-endemic floras. Outcome-blind source and observation-selection corrections retained
this direct regional difference in 24 of 24 frozen scenarios, while the independent broad
northern branch weakened. The clearest positive response was Palearctic. A source-matched
lineage decomposition showed that its broad floral-accessibility signal was primarily
associated with which source-available genera were represented rather than additional
species loading after genus entry. Continuous area moderation further showed that the
Palearctic genus-entry association with distance was stronger on smaller islands in all
eight source-mode x floristic-stratum scenarios. Neotropical reproductive assurance also
showed bounded smaller-island amplification, but with 59 supported islands and 14 spatial
blocks. These results reject a single universal floral island syndrome and instead support
biogeographically contingent, area-conditioned lineage assembly. Pollinator-mediated
filtering remains a compatible upstream hypothesis, but pollinator mobility, loss,
replacement, and effective service were not measured.

## Introduction

Distance from a continent and island area can influence island floras through more than one
route. Distance changes opportunities for plant dispersal and access to source floras,
whereas area combines target size, habitat and resource diversity, population capacity,
and local-extinction risk (MacArthur & Wilson 1963; Whittaker et al. 2008). The same
geographic variables can also filter pollinators and
other interaction partners through arrival, establishment, and persistence. Floral and
reproductive assemblages may therefore change because different plant lineages reach and
persist on islands, because reproductive assurance favors establishment, because the
benefits of floral attraction or access change, or through combinations of these routes.

This branching structure does not predict one obligatory direction. A decline in effective
pollination could favor self-compatibility or autonomous reproduction and secondarily
reduce floral investment (Baker 1955; Grossenbacher et al. 2017; Razanajatovo et al.
2019). Arrival opportunity and breeding system can contribute separately without a
detectable distance x breeding-system interaction, so distance should not be relabelled as
a Baker gradient (Zell et al. 2025). Pollinator decline could also reduce the return to attraction, increase
competition for scarce service, or leave floral traits intact when another functional
guild maintains service. Floral phenotypes can be compared with literature-defined
pollination syndromes, but phenotype alone neither identifies the visitor nor measures
effective pollen transfer (Ollerton et al. 2009; Rosas-Guerrero et al. 2014; Wang et al.
2020). Pollination syndromes are therefore used here only to interpret
an independently frozen plant-side pattern.

We tested a dispersal-filtered assembly hypothesis at the boundary observable with global
plant data. First, we asked whether a common floral/reproductive response vector emerges
worldwide along a mainland-distance/source-accessibility gradient. Second, we tested where
regional vectors occur and whether they differ directly. Third, we separated source-matched
genus entry from extra species loading within represented genera. Finally, we tested
whether island area continuously moderates these distance associations. We expected
regional branching rather than a universal syndrome, persistence beyond endemic-only
assemblages, and stronger filtering on smaller islands if limited capacity intensifies
either source-pool or interaction-channel constraints.

## Materials and methods

### Island universe and observation model

The frozen geographic universe contains 8,265 islands. Of these, 3,760 have no realized
flora record, 4,505 have at least one record, 425 have direct evidence in at least one
Chapter 1 trait domain, and 405 contain direct evidence from all three core domains
somewhere in the observed flora. All islands remain in the observation-process layer;
ecological trait models use only response-supported islands. Raw records were collapsed to
unique island x accepted-species incidence before calculating trait composition. An
unrecorded species, trait, or island was never converted to biological absence.

The primary estimand is the probability or assemblage score of a floral/reproductive state
conditional on the directly trait-resolved observed realized flora and geography. It is not
the probability that a trait evolved on an island or a census-level colonisation
probability.

### Traits, evidence, and support

Species-direct evidence was prioritized. Validated genus consensus could fill declared
species x trait gaps in broader sensitivity layers, without counting a species x trait
twice. The atomic when/where layer retained floral colour, floral form/access, and
reproductive states. The branch synthesis used two pollinator-name-free responses:
`accessibility_generalization` and strict `reproductive_assurance`; attraction traits were
excluded from the latter.

Within a vector test, atomic responses were selected by support before model fitting.
Responses observed on at least 50 islands met the count component of confirmatory support;
30-49 islands were pilot only. Direct comparisons retained only responses meeting the same
threshold in both contexts. False-discovery-rate correction was applied within frozen test
families, and an axis was interpreted only after its multivariate vector gate passed.

### Geographic and floristic comparisons

Primary operational contexts distinguished northern mid-latitude, tropical, northern
high-latitude, and southern-extratropical islands. Formal realm sensitivity separated the
Palearctic, Nearctic, Neotropical, and other supported realms. Models were stratified by
all native and native non-endemic floras; endemic-only strata remained under-supported for
confirmatory inference.

Distance to the nearest continent was log transformed in the canonical models and treated
as a composite of isolation and changing source accessibility, not a pure causal treatment.
Climate principal components and log island area were retained as prespecified controls.
Spatial-block cluster-robust covariance represented repeated information within geographic
blocks.

### Observation and spatial robustness

The when/where result was repeated under canonical species-count weighting, effective-trial
caps of 100, 50, and 20, and equal-island weighting. Trait-resolution fractions were
audited by island, stratum, and domain and entered as response-specific selection
adjustments. Raw, square-root, and log-distance forms retained the full island universe.
Each spatial block was deleted in turn. A later outcome-blind sensitivity combined four
frozen GIFT source definitions with unweighted and two stabilized inverse-probability
weighting specifications. Neither source assignment nor selection fitting used response
scores, effect sizes, or p-values.

### Source-matched lineage representation

For species with both frozen floral syndrome scores, functional position was defined as
`(-large_bee_like + generalized_accessible) / 2`. Genus source position was the mean among
scored species occurring in candidate GIFT mainland floras. For each island and frozen
source definition, genera were matched by source-region prevalence and source-species
richness class. The observed number of represented genera and island species was held fixed
within each availability class.

`entry_enrichment` gave one vote to each represented genus. `species_enrichment` weighted
each island species by its genus position. `loading_increment` was species enrichment minus
entry enrichment and therefore measured additional species weighting after entry. These are
assembly diagnostics; they do not reconstruct historical source ancestry.

### Continuous area moderation

Area-moderation models used continuous standardized distance and area without a small/large
threshold. Each response-specific stacked model contained distance, area,
distance x area, climate PCs 1-4, and spatial-block cluster-robust covariance. Direct
regional tests evaluated the corresponding response x distance x area x context
interaction vector. A negative distance x area term with a positive mean-area distance
slope was classified as a stronger distance association on smaller islands only after the
joint interaction-vector gate passed.

## Results

### Regional response vectors, not a universal syndrome

Northern mid-latitude and tropical assemblages each showed confirmatory multivariate
responses to the geographic gradient in all-native and native non-endemic floras. Their
vectors differed directly (all native q = 2.35 x 10^-8; native non-endemic q = 7.13 x
10^-7). Equal-island weighting, trait-resolution adjustment, three alternative distance
forms, and every single-block deletion retained the headline. Northern high-latitude and
southern-extratropical floras remained confirmatorily unresolved rather than demonstrated
null regions.

Self-compatibility alone did not explain the multivariate pattern. Its direct-evidence
distance association was weak in northern mid-latitude and tropical assemblages, and the
regional vectors combined floral-access and reproductive components differently.

### The clearest robust branch is Palearctic

After joint source-pool and observation-selection adjustment, the direct
northern-midlatitude versus tropical two-axis difference persisted in all 24 frozen
source-mode x selection-mode x stratum scenarios. The independent broad northern vector,
however, was unsupported in all 16 source-mode x IPW-mode x stratum scenarios. Palearctic
accessibility generalization plus reproductive assurance remained supported in 24 of 24
scenarios. Nearctic inference remained pilot and did not reproduce a universal northern
direction.

### The Palearctic floral branch is lineage-associated

The exact-self-incompatible Palearctic attraction/accessibility slope remained positive
after source adjustment and extensive block deletion, showing that measured selfing was
not necessary for that assemblage association. However, the observed slope did not exceed
a genus-composition-preserving source-pool null. Thus the data did not identify residual
within-genus floral filtering beyond measured genus composition.

In the source-matched broad flora, Palearctic genus-entry enrichment was positive and
FDR-supported in all eight source-mode x stratum scenarios. Species enrichment was also
positive in all eight, but loading increment was supported in zero of eight. The broad
result was therefore mainly associated with which source-available genera were represented.
Within exact-SI taxa, the surviving decomposition shifted toward loading among represented
genera, emphasizing that entry and loading need not respond identically.

### Smaller area amplifies broad Palearctic genus entry

The continuous distance x area extension retained positive mean-area Palearctic genus-entry
slopes and negative interaction terms in all eight broad source-mode x stratum scenarios.
Each scenario passed the hierarchical area-moderation gate, indicating a stronger distance
association on smaller islands. No equally stable area law was recovered for species
loading.

Neotropical reproductive assurance was also more strongly distance-associated on smaller
islands in both primary strata. This result was based on 59 response-supported islands in
14 spatial blocks, and the direct Palearctic-Neotropical reproductive-assurance interaction
difference was unsupported. It is therefore Baker-compatible assemblage evidence, not a
direct observation of establishment after a founding event.

### Area conditions a regional accessibility contrast

Within broad northern-midlatitude and tropical contexts, plant vectors did not show a
simple shared smaller-island amplification. Nevertheless, the direct tropical-minus-north
accessibility interaction was supported in both strata (q = 0.00780 and 0.0134). Area thus
conditioned a regional difference in accessibility geography without identifying the
pollinator process that generated it.

## Discussion

### Island floral responses are biogeographically contingent

The central result is a branching pattern rather than a universal floral island syndrome,
consistent with direct Pacific island-mainland comparisons that found lineage- and
archipelago-dependent rather than uniformly smaller island flowers
(Hetherington-Rauth & Johnson 2020).
Isolation/source accessibility is associated with plant traits in both northern
mid-latitude and tropical assemblages, but not through the same multivariate response.
The contrast survives observation and source-pool sensitivities, whereas a broad
latitude-defined northern branch does not. Formal biogeographic structure therefore matters
for deciding where a pattern occurs and prevents a Palearctic result from being generalized
to the Nearctic or all northern temperate islands.

### Source-pool filtering precedes a pollinator mechanism claim

Lineage analysis materially changes the interpretation. Broad Palearctic functional
reorganization is expressed primarily through genus entry, and the exact-SI floral slope is
fully reproducible under a genus-fixed null. The strongest supported sequence is therefore
`geography -> source-available lineage assembly -> floral/reproductive assemblage
composition`. Pollination-channel filtering could act on arrival or persistence of whole
lineages, but the present data cannot separate that process from non-pollination
biogeographic assembly.

### Reproductive assurance and attraction are distinct routes

Pollinator limitation could produce floral simplification through at least two routes. A
reproductive-assurance route can favor self-compatible or autonomous reproduction and lead
secondarily to reduced floral investment. An attraction/access route can reduce the benefit
of signaling, increase competition for scarce service, or preserve/redesign floral traits
when other visitors remain effective. The weak direct self-compatibility slope and the
positive exact-SI Palearctic floral association show that a simple selfing-syndrome sequence
is insufficient as the sole description. Conversely, the Neotropical area interaction is
compatible with reproductive assurance but cannot identify Baker's law, founding events,
or post-colonisation evolution.

### Pollination syndromes remain a bounded interpretation

Large-bee-like, bird-like, butterfly-like, and other syndrome scores summarize floral-trait
concordance, not realized visitor identity or effectiveness. The Palearctic shift away from
large-bee-like architecture and tropical maintenance of specialized-access traits can be
compared with hypotheses about differences in guild dispersal, establishment, and
preference. Yet neither Bombus loss nor bird/Lepidoptera replacement was observed. The
required causal test must measure visitors, single-visit pollen delivery, reproductive
outcome, and local persistence independently of floral phenotype.

### Island area is a capacity proxy, not a measured mechanism

The uniform smaller-island amplification of source-matched Palearctic genus entry is the
most stable area result. It is consistent with stronger assembly constraints when target
size, habitat/resource diversity, or population capacity are limited. The analysis does
not distinguish those components, measure founder number or drift, or resolve stepping-stone
connectivity. Area should therefore be interpreted as a shared capacity proxy rather than
a direct estimate of any one process.

## Conclusions

Global island floral/reproductive composition does not converge on one syndrome. Regional
vectors differ, the most robust positive branch is Palearctic, and its broad area-dependent
signal is primarily source-matched genus entry. These results resolve substantial parts of
the Chapter 1 when/where question while lowering the mechanism claim: pollination-channel
filtering is a biologically coherent explanation to test, not a process identified by the
current global data. The Izu *Campanula* studies provide a concrete local example in which
Bombus absence, smaller flowers, and breeding-system change covaried, but they do not test
the global Palearctic assemblage (Inoue & Amano 1986; Inoue 1988). Caribbean radiation data
likewise show that alternative guilds, generalization, and autonomous selfing can coexist
without defining one global replacement route (Martén-Rodríguez et al. 2010). The next step
is direct interaction evidence linking plant functional
position to visitor exposure, pollination effectiveness, and reproductive outcome.

## Data and code availability

The submission-candidate analysis is locked by
`config/chapter1_submission_freeze_lock.json`. The GitHub Actions artifact is checksum
verified but temporarily retained; replace this paragraph with the durable release DOI or
permanent repository URL before submission.

## References

Baker, H. G. (1955). Self-compatibility and establishment after “long-distance”
dispersal. *Evolution*, 9, 347-349. https://doi.org/10.1111/j.1558-5646.1955.tb01544.x

Grossenbacher, D. L., et al. (2017). Self-compatibility is over-represented on islands.
*New Phytologist*, 215, 469-478. https://doi.org/10.1111/nph.14534

Hetherington-Rauth, M. C., & Johnson, M. T. J. (2020). Floral trait evolution of
angiosperms on Pacific islands. *The American Naturalist*, 196, 87-100.
https://doi.org/10.1086/709018

Inoue, K. (1988). Pattern of breeding-system change in the Izu Islands in *Campanula
punctata*: bumblebee-absence hypothesis. *Plant Species Biology*, 3, 125-128.
https://doi.org/10.1111/j.1442-1984.1988.tb00178.x

Inoue, K., & Amano, M. (1986). Evolution of *Campanula punctata* Lam. in the Izu
Islands: changes of pollinators and evolution of breeding systems. *Plant Species
Biology*, 1, 89-97. https://doi.org/10.1111/j.1442-1984.1986.tb00018.x

MacArthur, R. H., & Wilson, E. O. (1963). An equilibrium theory of insular
zoogeography. *Evolution*, 17, 373-387.
https://doi.org/10.1111/j.1558-5646.1963.tb03295.x

Martén-Rodríguez, S., Fenster, C. B., Agnarsson, I., Skog, L. E., & Zimmer, E. A.
(2010). Evolutionary breakdown of pollination specialization in a Caribbean plant
radiation. *New Phytologist*, 188, 403-417.
https://doi.org/10.1111/j.1469-8137.2010.03330.x

Ollerton, J., et al. (2009). A global test of the pollination syndrome hypothesis.
*Annals of Botany*, 103, 1471-1480. https://doi.org/10.1093/aob/mcp031

Razanajatovo, M., et al. (2019). Autofertility and self-compatibility moderately benefit
island colonization of plants. *Global Ecology and Biogeography*, 28, 341-352.
https://doi.org/10.1111/geb.12854

Rosas-Guerrero, V., et al. (2014). A quantitative review of pollination syndromes: do
floral traits predict effective pollinators? *Ecology Letters*, 17, 388-400.
https://doi.org/10.1111/ele.12224

Wang, X., Wen, M., Qian, X., Pei, N., & Zhang, D. (2020). Plants are visited by more
pollinator species than pollination syndromes predicted in an oceanic island community.
*Scientific Reports*, 10, 13918. https://doi.org/10.1038/s41598-020-70954-7

Whittaker, R. J., Triantis, K. A., & Ladle, R. J. (2008). A general dynamic theory of
oceanic island biogeography. *Journal of Biogeography*, 35, 977-994.
https://doi.org/10.1111/j.1365-2699.2008.01892.x

Zell, A. N., Miranda, C. H., Grady, E. L., Grossenbacher, D. L., & Igić, B. (2025).
Island colonization in flowering plants is determined by the interplay of breeding system,
lifespan, floral symmetry, and arrival opportunity. *New Phytologist*, 245, 420-432.
https://doi.org/10.1111/nph.20234

Full author lists and claim-specific boundaries are audited in
`docs/chapter1_submission_literature_evidence_20260827.md`. Journal-style export remains a
pre-submission formatting gate.
