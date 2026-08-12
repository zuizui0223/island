# Rule-unlock wave 2 checkpoint (2026-08-12)

This checkpoint freezes individually reviewed species-direct evidence selected
for information gain in the current `genus x trait_name` queue. It contains
seventy-five reproductive records and thirty visible-morphology records.
It never emits genus,
family, or global inference itself.

Reproductive concepts remain separate: failure of unassisted bagged flowers is
stored as `autonomous_selfing_capacity=absent`, while successful hand-selfing is
stored independently as `self_incompatibility=SC`. Flower symmetry is accepted
only from explicit symmetry terminology or a directly described bilabiate
corolla; corolla-tube depth remains a separate trait.

The UNESP thesis and publisher/government/source pages were retrieved and
fingerprinted. Two otherwise inaccessible original full-text statements use a
hash of the verified exact excerpt and say so explicitly in
`content_sha256_basis`; search snippets are not represented as evidence.

The University of Calicut doctoral thesis contributes two distinct records per
studied species: controlled self-compatibility and the explicitly reported mixed
mating system. The downloaded original repository PDF is fingerprinted; neither
record is inferred from the other reproductive concept.

The Tristaniopsis record keeps a primary bagging result (zero fruit from 189
bagged flowers, versus 19.1% under open pollination) as autonomous-selfing
capacity only. It is not converted into self-incompatibility. A separate 2019
field study supplies explicit actinomorphic flower morphology. The PROTA Adenia
species account supplies a regular-flower statement, normalized to actinomorphic.

The latest morphology additions use a Botanical Survey of India district flora
for the corymbose cyme of *Polycarpaea corymbosa*, a non-cultivar Australian
National Botanic Gardens observation for the star-like flower of *Boronia
muelleri*, and an expert-attributed eFlora of India description of the branched
spadix and male catkins of *Benstonea foetida*. The latter retains the credited
India Biodiversity Portal treatment as its source lineage so redistribution is
not counted as independent evidence.

The next high-yield additions retain three different evidence tiers without
weakening the identity or trait gates. A peer-reviewed New Caledonian taxonomic
revision supplies the cream/white-edged corolla of *Pleioluma balansana*; cream
remains in the ontology's white state. A species-only New Zealand nursery page
supplies only the explicitly star-shaped flower form of *Kunzea ericoides* and
is retained as Medium Tier C morphology, not reproductive evidence. A published
species morphology account explicitly identifies *Corchorus olitorius* in
Malvaceae and describes its flower as actinomorphic; because it is a secondary
morphology account, it is conservatively retained as Medium rather than High.

The third morphology increment adds three source-fingerprinted, species-direct
statements. NatureSpot's specialist regional account identifies the rayed
flowerheads of *Jacobaea aquatica* and is retained as Medium. The National
Institute of Biological Resources' *Flora of Korea* treatment and downloaded page
image identify pistillate catkins in *Carpinus laxiflora*. The original article
abstract for *Citharexylum myrianthum*, retrieved from the FAO AGRIS record and
deduplicated under its DOI lineage, states that flowers occur in raceme-like
inflorescences. The two institutional/original-publication records are High.

The fourth increment targets two independent third-species rule unlocks. Gary
Starr's University of Adelaide doctoral thesis reports fruit after pre-anthesis
insect exclusion in two wild *Hakea carinata* populations and explicitly concludes
that the species can produce seed by autogamy; this is retained only as
`autonomous_selfing_capacity=autonomous`. The Chinese University of Hong Kong's
Shiu Ying Hu Herbarium Pro-Factsheet identifies *Celtis sinensis* and records its
flower symmetry as radial (actinomorphic). The thesis experiment is High. The
herbarium factsheet is conservatively Medium because it does not attribute that
field to an individual source or reviewer and categorizes the species as
horticultural, although no cultivar or hybrid is named. Both records are
source-fingerprinted and kept as distinct reproductive and morphology traits.

The fifth increment adds two non-overlapping direct records discovered by
credential-free institutional-site traversal. The Chinese University of Hong
Kong factsheet for *Phoenix roebelenii* explicitly records a panicle, but the
page is categorized as horticultural and has no field-level reviewer, so this
structure-only record remains Medium Tier C. Singapore NParks' government
species record for *Gonystylus confusus* gives the accepted species, matching
family, and flower colour as yellow/golden; it is retained as High Tier A. No
cultivar statement, search snippet, or cross-trait mapping is used.

The sixth increment adds eight source-fingerprinted records without forcing
genus consistency. CIRAD's species treatment for *Turraea cadetii* explicitly
records axillary cymes and white petals, while the original Willdenowia
description of *Mosiera yamaniguensis* records solitary axillary flowers and
white petals. These are the third species for the existing trait-specific
Turraea inflorescence and Mosiera inflorescence agreements; formal Low remains
conditional on the shared all-evidence masked and lineage validation.

The same increment also records counterevidence rather than manufacturing rule
unlocks. Vatke's original description gives yellow and brownish flowers for
*Dichaetanthera rutenbergiana*, contradicting the queued two-species red/pink
agreement. Munzinger & Gateble's peer-reviewed comparison preserves explicit
colour state sets for *Acropogon bullatus*, *A. mesophilus*, and *A. veillonii*.
These rows intentionally prevent an undifferentiated `multicolored_variable`
label from hiding biologically different colour combinations.

The seventh increment uses a peer-reviewed Galapagos species table for thirteen
previously unresolved, unflagged compatibility or autonomous-selfing cells. Its
numbered original references are preserved as lineage identifiers. Two
*Chiococca alba* cells are deliberately excluded because the table's dagger
footnote contradicts the displayed compatibility state. Repeated artificial
selfing of *Encyclia phoenicea* and *E. plicata* is retained as Medium evidence
of compatibility under a single living-collection lineage. A separate
pollinator-exclusion experiment supports only absence of autonomous selfing in
*E. tampensis* and is not converted to self-incompatibility.

The eighth increment adds thirty compatibility cells from Lord's Southern Ocean
Islands supplementary species table. The paper explicitly cautions that many
compatibility entries are statements rather than experiments, so all remain
Medium and share one conservative compilation lineage until their numbered
original sources are retrieved. `SCp` is retained as `mixed_or_variable`, not
collapsed to SC. Five High records come from Anderson et al.'s Juan Fernandez
bagging and controlled-cross experiments. They keep SI, autonomous or delayed
selfing, and mixed mating as distinct traits; the conflicting Escallonia
self-cross assays are not used to overwrite its compatibility state.

The ninth increment targets three visible traits whose current genus rules had
two agreeing species. An official Academie de La Reunion arboretum sheet names
*Pyrostria commersonii* and describes its corolla as yellowish white; the colour
is retained as the explicit state set `white|yellow_orange`. A non-cultivar
*Quintinia acutifolia* nursery species page explicitly reports white flowers and
is limited to Medium Tier C morphology. The Bangladesh Forest Information
System's government species record explicitly identifies *Suregada lanceolata*
and records actinomorphic floral symmetry. These records remain direct Medium
evidence; any resulting Low cells must still pass the shared all-evidence,
masked, and leave-one-lineage-out rule validation.

The combined files supersede the previous curated checkpoint as workflow input.
Every new row has an `accept` audit decision, exact source URL, supporting
excerpt, retrieval timestamp, content fingerprint, and source lineage.
