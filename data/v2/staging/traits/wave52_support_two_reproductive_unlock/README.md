# Wave52 support-two reproductive unlock

Wave52 adds one independently sourced species-level self-compatibility record to
each of three `genus × axis × trait_name` rules that had exactly two agreeing
species in the formal Wave51 artifact. The frozen targets are `Helichrysum`,
`Mammillaria`, and `Triumfetta`.

The packet keeps `self_incompatibility`, `autonomous_selfing_capacity`, and
`mating_system` separate. “Self-fertile”, “autogame”, and “self-pollinated” are
accepted only as Medium `self_incompatibility=SC`; none is silently promoted to
autonomous selfing. Family inference and global fallback remain disabled.

`Helichrysum italicum` already had direct evidence for another reproductive
trait. Its new `self_incompatibility` row therefore enriches the species-trait
ledger without adding a second species-axis cell. This is the regression case
that guards against the former axis-only acquisition gate.
