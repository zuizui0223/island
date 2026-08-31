# Wave41 external-congener reproduction checkpoint

This checkpoint adds species-direct reproductive evidence for species outside
the fixed 106,295-species island universe. The evidence may train only the
shared trait-specific `genus x axis x trait_name` Validated Low rules. It is
never counted as confirmatory direct coverage, family inference, or global
fallback.

## Frozen acquisition

- Retrieval date: 2026-08-28 UTC
- Queries: 2,001 exact names against WFO and GBIF, 4,002 requests total
- Query cost: USD 0
- Frozen WFO/GBIF response files: all 2,001 names have both responses
- Strict identity gate: exact species rank in both backbones, identical
  accepted species, and family agreement
- Fail-closed exclusions: ambiguous/fuzzy matches, one-backbone matches,
  family conflicts, infraspecific records, and all species mapping into the
  fixed target universe

## Sources

- Zell et al. BSDB source, pinned repository commit
  `9e87946d1e3121d39e657b702cf9f92ccc10936e`, raw file
  `Zell_df_12_29_23.csv`, SHA-256
  `f935c8150b3d719aba5f62f14335c1a185304155403bf50db1e3ef1393fc55f4`
- Rodger et al. 2021 pollinator-contribution dataset, Figshare v1 snapshot
  `S3_pc_publish.csv`, committed gzip SHA-256
  `fa745c578f3537933fafedc1d36b4ea266348cd83d7f6cbb231c253b0f348d3f`

Original publication or dataset citations are retained as source lineages.
Provider redistribution is therefore not counted as independent evidence.

## Frozen result

- 1,591 evidence rows
- 1,391 external species
- 1,541 external species x trait cells
- 726 underlying source lineages
- 1,152 self-incompatibility rows and 439 autonomous-selfing rows
- 0 target-universe species entered confirmatory direct coverage

`wave41_external_congener_checkpoint_summary.json` records every input/output
hash, mapping and selection reason, request count, cost, and fail-closed check.
The formal GitHub workflow regenerates the checkpoint from the pinned raw
source plus frozen responses before rebuilding all-evidence genus rules.
