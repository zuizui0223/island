# Rule-unlock wave 2 checkpoint (2026-08-12)

This checkpoint freezes individually reviewed species-direct evidence selected
for information gain in the current `genus x trait_name` queue. It adds twelve
reproductive records and two visible-morphology records. It never emits genus,
family, or global inference itself.

Reproductive concepts remain separate: failure of unassisted bagged flowers is
stored as `autonomous_selfing_capacity=absent`, while successful hand-selfing is
stored independently as `self_incompatibility=SC`. Flower symmetry is accepted
only where the fetched species page explicitly says `actinomorphic` or
`actinomorfas`.

The UNESP thesis and publisher/government/source pages were retrieved and
fingerprinted. Two otherwise inaccessible original full-text statements use a
hash of the verified exact excerpt and say so explicitly in
`content_sha256_basis`; search snippets are not represented as evidence.

The combined files supersede the previous curated checkpoint as workflow input.
Every new row has an `accept` audit decision, exact source URL, supporting
excerpt, retrieval timestamp, content fingerprint, and source lineage.
