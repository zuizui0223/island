# Rule-unlock wave 2 checkpoint (2026-08-12)

This checkpoint freezes individually reviewed species-direct evidence selected
for information gain in the current `genus x trait_name` queue. It adds twenty-two
reproductive records and five visible-morphology records. It never emits genus,
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

The combined files supersede the previous curated checkpoint as workflow input.
Every new row has an `accept` audit decision, exact source URL, supporting
excerpt, retrieval timestamp, content fingerprint, and source lineage.
