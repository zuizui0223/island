# Wave58 reproductive-assurance recovery

Wave58 keeps the canonical checkpoint unchanged until collision audit against private TRY plus Waves53-55. It adds two source-reviewed direct candidates (`Cinnamomum verum`, `Neolitsea sericea`), blocks the ambiguous `Ourisia poeppigii` third-species shortcut, validates the public-Wave52 `Schoenoplectiella` SC unlock separately, and searches exact unresolved species only in coherent support-two genus/trait pairs.

## Schoenoplectiella

Wave52 has two independent species-direct SC records: `S. mucronata` and `S. supina`. Wave56 adds independently sourced `S. juncoides` SC (Sada et al. 2013, DOI 10.1016/j.pestbp.2013.05.013). The Wave58 audit requires 3 species, one state, unique source lineages, reproductive dominance >=0.95, species LOO >=0.85 and source-lineage LOO >=0.85. On the immutable Wave52 baseline this produces one direct candidate cell plus 34 Validated-Low candidate cells = 35 public-baseline candidate cells. Every row remains promotion_allowed=false because the current private-plus-public collision audit is not present in public artifacts.

## Search strategy

The broad Wave57 Crossref query was too noisy. Wave58 therefore uses exact-species Europe PMC queries. It selects only genera where Wave52 already has exactly two direct species for the same reproductive trait and those two agree on one normalized state. Conflicted or already-handled genera are blocked. At most 8 unresolved species per genus are queried, preventing one large genus from monopolizing the wave. Search hits are review leads only: no value, quality grade, genus-rule vote, or coverage promotion is inferred from abstract co-occurrence.

Canonical remains reproductive assurance 48,527 / 106,295 and total Chapter 1 222,759 / 318,885 = 69.8556% until current-ledger net gains are verified.
