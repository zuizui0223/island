# Rule-unlock wave 2 checkpoint (2026-08-12)

This checkpoint freezes individually reviewed species-direct evidence selected
for information gain in the current `genus x trait_name` queue. It adds twenty-four
reproductive records and seventeen visible-morphology records. It never emits genus,
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

The combined files supersede the previous curated checkpoint as workflow input.
Every new row has an `accept` audit decision, exact source URL, supporting
excerpt, retrieval timestamp, content fingerprint, and source lineage.
