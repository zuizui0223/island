# Restart source screen

Date: 2026-09-07. Reviewer: Codex (automated screening, not independent human review). Baseline Run 34093932899. No promoted rows; no claimed coverage gain.

Two discovery queries were executed using the available Web search tool: `Spermacoce self compatible breeding system pollination` and `Dicliptera self compatibility breeding system`. Search API billing is not exposed; cost is unknown, not asserted zero.

The original publisher page was fetched at https://revistas.um.es/analesbio/article/view/analesbio.39.13 (Universidad de Murcia; Raju and Radha Krishna 2017, Anales de Biologia 39, 111-126), following AGRIS discovery. The page describes Spermacoce hispida, S. articularis and S. pusilla. Exact short excerpt: "Las flores son débilmente protándrica, nectaríferas y auto-polinizadoras."

Lineage: citation:Raju-Radha-Krishna-2017-Anales-Biologia-39-111-126. AGRIS and publisher are not independent studies. Status: original full-text experimental-methods review pending; promotion=false; genus_rule_training=false. The abstract does not establish a controlled self-compatibility test or justify translating self-pollination into autonomous selfing. Check the species-specific results and prior ledger before adopting any trait. This screen reopens a source lane, not a new accepted observation.

## Full-text follow-up and deduplication

The publisher PDF was subsequently retrieved from https://revistas.um.es/analesbio/article/download/analesbio.39.13/281161/1418601. SHA256: `18c86920287c9057dd98f8c9b990e0e5c1a556b83637972efff8aa6bc03d71ff`. Methods (printed pp.113-114) measure open-pollinated fruit/seed set; discussion (p.124) attributes reproduction partly to autonomous pollination, without a controlled compatibility experiment in these methods. The recovered direct ledger already contains autonomous_selfing_capacity for all three species under `doi:10.6018/analesbio.39.13`. Therefore this paper is a duplicate, not three new records. Do not revisit this source for that same trait. Citation alias above resolves to this DOI, not another independent lineage.

## Refreshed search queue

`scripts/export_restart_support_two_queue.py` was executed against the recovered GitHub artifact. It found 248 genus-trait pairs with exactly two agreeing direct species and at least one unresolved reproductive-axis congener. The CSV includes species lists and lineage strings for audit, not inferred rules. The larger prefilter count includes pairs without unresolved axes and is not the acquisition queue size. Missing individual traits in already resolved axes are outside this queue.

Follow-up queries: `Lithocarpus self incompatibility breeding system pollination`; `Callicarpa autonomous selfing bagging breeding system`. Two new source leads were checked: PubMed 37341801 returned no readable body through the Web tool, and SSRN 7038775 returned HTTP403. Neither is accepted evidence from snippets. Lithocarpus's two supporting species share one dataset DOI and require upstream-source reconciliation; Callicarpa SC/SI variation must not be substituted for its autonomous-selfing trait. Four discovery queries have been issued across the restart screen in total; accepted coverage gain remains zero. This is ongoing screening, not a completed acquisition wave.
