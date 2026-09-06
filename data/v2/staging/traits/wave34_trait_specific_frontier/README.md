# Wave34 trait-specific frontier

This bundle contains additive Low candidates and the acquisition queue derived
from the frozen Wave33 secondary baseline. Every accepted rule retains
`genus x axis x trait_name`; application never joins on `genus x axis` alone.

Eligibility requires at least three directly supported species, species
leave-one-out accuracy >= 0.8, and source-lineage leave-one-out accuracy >= 0.8.
Predictions retain one to three states. Family inference and global fallback are
absent. The candidates fill only previously unresolved cells, so the materialized
frontier has zero loss from the frozen baseline.
