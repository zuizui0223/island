# Prioritize source permissions for the H1 reproductive core

26 November 2026 JST. The preceding overlap-aware cell audit was independently reproduced in [run 45718126892](https://github.com/zuizui0223/island/actions/runs/45718126892); see [anonymous receipt](RIGHTS_PRIORITY_VERIFIED.json).

## Exact-state crosswalk, not exact rights adjudication

Four original archives were SHA256-verified in the private working environment. The pinned H1 Direct state audit contains 12,401 resolved records with an informative single state for self-incompatibility, mating system or autonomous capacity, corresponding to the final H1 binary state definitions. Exact accepted-species × trait × canonical-state matching to the v1 Direct ledger retains 10,654 rows (7,852 species). The H1 state audit does not retain original source lineage; matching does not establish identical upstream evidence. 1,747 unmatched rows remain outside this crosswalk. Join the matched records to the **v1 reproductive axis-cell** rights audit (many trait records can share one axis cell). This is a conservative prioritization proxy, not a release audit of the H1 rows.

| Trait | H1 informative rows | Matched v1 state | Matched records in authorized v1 cells | Matched records in review-required cells |
|---|---:|---:|---:|---:|
| Self-compatibility | 6671 | 6556 | 1135 | 5421 |
| Mating system | 2054 | 2047 | 56 | 1991 |
| Autonomous capacity | 3676 | 2051 | 87 | 1964 |
| Total | 12401 | 10654 | 1278 | 9376 |

## Targeted direct evidence

| Source | v1 Direct reproductive trait rows (any state) | Exact-state overlap, three core traits | SC/SI | Mating | Autonomy |
|---|---:|---:|---:|---:|---:|
| FloraWeb/BiolFlor | 5420 | 3869 | 1807 | 1464 | 598 |
| Pladias | 1656 | 1102 | 134 | 968 | 0 |
| Europe PMC | 1364 | 1092 | 629 | 77 | 386 |
| Ken Fern | 235 | 235 | 235 | 0 | 0 |
| Go Botany | 1233 | 0 | 0 | 0 | 0 |

These source associations overlap; do not add counts as independent species. Go Botany's 1,233 v1 Direct reproductive records are cleistogamy, not any of the three H1 core traits. Its family can nevertheless block an axis cell supporting a matched core record. Ken Fern family normalization collapses PFAF, Useful Tropical Plants and Useful Temperate Plants in the acquisition script; the original domain and release are not uniformly resolved. Do not assign one website's licence to the family.

## Revised order

For a planning objective based on matched H1 core trait records, a whole-blocker-set greedy heuristic selects BiolFlor (923 additional matched records), Pladias (1,399), Europe PMC (1,016), Go Botany (836), unresolved other DOI (558), Ken Fern (274), Ecoflora (241), unresolved citation hash (278). The first three would add 3,338 records to the 1,278 baseline if every affected blocker were cleared in full. Europe PMC and unresolved source labels are provenance-resolution tasks, not single licensors. These counts are conditional ceilings under the conservative v1 axis gate, not actual permission, island-model sample sizes, or scientifically independent confirmations.

Decision: BiolFlor/FloraWeb and Pladias first for provider contact; Europe PMC article-level licence verification in parallel. Keep the broader five-family inquiry register for whole-axis reproducibility. Do not send an attachment before its specific source scope and sharing rights are reviewed by the authors.

## Artifacts

- [Aggregate result](CORE_RIGHTS_SCOPE_RESULT.json) and [reproducible local audit](core_rights_scope.py). Local calculation completed on original pinned archives; the private-scope option is deliberately not run in public CI.
- [Targeted request drafts](TARGETED_PERMISSION_REQUESTS_UNSENT.md).
- Five provider-specific private scope CSVs, containing accepted name, trait name, quality and only the named provider's lineage references (no trait values or other providers' lineage records), are stored in the already authorized owner-only private Drive. Metadata readback verified private access. Do not commit the CSVs or private Drive identifiers.
- No new legal permission, no edit to frozen results or scientific selector, no public record release or external correspondence.
