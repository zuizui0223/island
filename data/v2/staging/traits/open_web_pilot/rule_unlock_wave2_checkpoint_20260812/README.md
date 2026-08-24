# Rule-unlock wave 2 checkpoint (2026-08-12)

This checkpoint freezes individually reviewed species-direct evidence selected
for information gain in the current `genus x trait_name` queue. It contains
seventy-five reproductive records and 399 visible-morphology records.
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

The tenth increment freezes a credential-free bulk traversal of the Bangladesh
Forest Information System government index. The downloaded HTML snapshot
(SHA-256 `80f4d97075ebb66d41b5707b1ae4bc330dd41ea4df567b4ea19f63531aa88b4f`)
contains 281 species treatments. Exact accepted-name matching found 119 target
species, of which 108 had an explicit floral-symmetry field. Forty-nine records
without prior direct floral-symmetry evidence are added as Medium, each with a
species-specific provider-treatment lineage and quote. Existing direct values
were never overwritten; four BFIS statements that disagree with an existing
direct symmetry record were therefore not ingested. Historical family labels in
the source quote are retained, while the identity gate checks each species
against its current family in the fixed target master.

The eleventh increment traverses the Indian Institute of Science's India Flora
Online index rather than issuing one search request per species. Of 14,982
indexed angiosperm treatments, 5,607 names exactly matched the fixed Island
master; only 1,626 treatments with an unresolved colour or structure axis were
downloaded. The complete species page, stable URL, page SHA-256 and exact `Key
identification features` excerpt are preserved. A full audit reviewed all 262
extracted candidates: 250 were accepted (precision 95.42%) and 12 were rejected.
Two rejected horticultural-hybrid treatments yield cultivar contamination of
0.76%, below the 2% ceiling. The other rejections prevent fruit, pod, stem,
involucre or spathe states from being transferred to flowers. Seventeen
accepted colour records were corrected to preserve the exact floral state set
while dropping colours belonging only to non-floral organs. All 250 accepted
records remain species-direct Medium evidence despite the university Tier A
source because the individual morphology fields are compiled treatments rather
than independently attributed primary measurements.

The twelfth increment reuses the same downloaded IISc index but does not relax
the name gate. WFO June 2026 first connects historical treatment names to the
fixed accepted master, after which the GBIF species-match API must independently
return the same accepted species, family, species rank and exact match with at
least 95 confidence. This reduced 605 morphology-relevant synonym candidates to
373 treatments with two-backbone agreement. Only those pages were downloaded;
71 trait candidates were extracted and all were reviewed. Seventy were accepted
(precision 98.59%) across 64 accepted species, with no cultivar contamination.
One apparent flower-colour record was rejected because red referred only to a
ripe drupe. Three accepted state sets were corrected to remove leaf or fruit
colours while retaining the explicit floral colours. The searched historical
name, WFO and GBIF identifiers, full page hash, exact quote and one lineage per
IISc treatment are preserved; neither one-backbone nor fuzzy matches enter the
ledger.

The thirteenth increment replaces per-species search with the official Pl@ntUse
MediaWiki API for Plant Resources of Tropical Africa (PROTA). The public PROTA
category contains 3,331 species monographs; 1,702 exactly match the fixed
Island master, and 1,287 have at least one unresolved strict axis. Fifty-two
batched API requests captured the exact page revision, timestamp and wikitext
hash. Redirects were held back, leaving 699 species pages with an authored
Description section. The complete audit reviewed 469 trait candidates and
accepted 448 across 268 species (precision 95.52%, cultivar contamination 0%).
The 21 rejections prevent mature fig colours, fruits, cyathial involucres,
bracts, non-floral structures, reference-title matches and statements about a
different species from entering the ledger. The accepted Tier A records retain
149 colour, 126 symmetry, 159 inflorescence-display, 9 floral-form and 5
trait-specific reproductive statements as species-direct High evidence. Each
row stores its exact quote, stable monograph URL, authored citation, revision
ID, revision timestamp, source-treatment lineage and revision-content SHA-256.
Self-fertility is mapped only to compatibility; self-incompatibility,
autonomous selfing and mating system remain separate traits.

The combined files supersede the previous curated checkpoint as workflow input.
Every new row has an `accept` audit decision, exact source URL, supporting
excerpt, retrieval timestamp, content fingerprint, and source lineage.
