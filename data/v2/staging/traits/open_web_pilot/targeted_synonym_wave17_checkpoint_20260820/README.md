# Wave 17 strict synonym recovery

This checkpoint recovers already-published BSdb reproductive records that the
accepted-name join had missed.  It queried 1,383 unresolved BSdb names against
the World Flora Online 2026-06 matcher and the GBIF species matcher.  A mapping
was accepted only when both services returned an exact species-rank match to
the same fixed-master species and source, GBIF, and master families agreed.
The 1,383 frozen response records are stored in
`bsdb_name_resolution_responses_20260820.jsonl.gz`; no API credential or live
request is required to reproduce the 21 accepted mappings.

The pinned BSdb release is `dirtyplants/BSdb@9e87946d1e3121d39e657b702cf9f92ccc10936e`
(`Zell_df_12_29_23.csv`, SHA-256
`f935c8150b3d719aba5f62f14335c1a185304155403bf50db1e3ef1393fc55f4`).
Underlying study keys, rather than BSdb itself, remain the source lineages.

Ten new BSdb evidence rows cover eight previously missing species x trait
pairs.  Three independent rows for *Doona cordifolia* disagree (SC versus SI),
so the checkpoint preserves the rows as an unresolved direct conflict.  It
does not choose the majority value.  Two additional reviewed rows come from a
newly fetched exact-species Flora of Zambia treatment of *Monopsis
stellarioides*; the quote describes the corolla directly and supports both a
bilabiate form and a variable colour set.

`combined_curated_evidence_20260820.csv` and its audit append these 12 rows to
Wave 16.  The checkpoint contains no genus, family, or global inference.
Validated Low and realized strict coverage are rebuilt only by the shared
all-evidence workflow.
