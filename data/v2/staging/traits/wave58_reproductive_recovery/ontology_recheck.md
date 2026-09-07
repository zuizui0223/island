# Wave58 trait-specific recheck

Reviewer: Codex (source-text review); date: 2026-09-07.

Original source: https://www.nparks.gov.sg/florafaunaweb/flora/3/0/3086

The page identifies Pouteria campechiana and lists self-pollination together with fauna-mediated pollination. Neither term demonstrates autonomous selfing in the absence of a pollen vector. Their coexistence also does not demonstrate a mixed or variable autonomous-selfing phenotype.

The retained record is an explicitly rejected mapping, with empty normalized value and quality, promotion_allowed=false and genus_rule_training_allowed=false. It must not resolve a species-trait or species-axis cell or count as a genus-rule vote. Its raw source quote and lineage remain available for audit. It is not a confirmed counterexample to absence of autonomous selfing. The precautionary Pouteria acquisition block remains pending trait-specific re-review.

This correction changes no accepted coverage: the earlier candidate already prohibited promotion. It does not claim an independently recomputed canonical coverage or loss count. Wave58 Run 34085090235 contains discovery and public-Wave52 unlock candidates, not a private-plus-public collision audit.

## Schoenoplectiella source independence remains unverified

The downloaded support CSV contains three distinct lineage strings, but two are species-keyed BiolFlor/FloraWeb identifiers. Their distinctness does not establish distinct original studies. The current executable's `source_lineage.nunique()==3` and lineage-LOO score of 1.0 therefore verify only the supplied identifier partition. They do not certify original-source independence. Recover the provider bibliographies and reconcile against the juncoides DOI before promotion. Do not call the 34 Low candidates independently validated on this evidence alone. No claim is made that the two provider records necessarily share an original study either.

The per-record review is saved in `source_independence_recheck.csv`. This is an additional gate, alongside current-ledger collisions and shared trait-specific rule reconstruction; all candidate promotion remains disabled.
