# Source-scale trait acquisition, 2026-08-07

## Scope

This acquisition round uses the fixed denominator of 106,295 accepted species
and three strict axes (318,885 species × axis cells). Family inference and
global fallback are excluded. Direct evidence is resolved before
trait-specific `genus × trait_name` Validated Low inference.

The round has two credential-free source-scale lanes:

1. official World Flora Online flora archives for floral descriptions; and
2. Europe PMC official REST full-text XML for reproduction traits.

Search results and snippets are discovery only. Every promoted record retains
an accepted species, exact source quote, stable URL, source lineage, retrieval
manifest, and review decision.

## WFO official flora archives

The frozen candidate artifact contains 6,696 rows from 24 official archive
downloads and 18 flora providers. The second fail-closed organ/context gate
accepted 6,445 rows (4,348 species; 5,894 species × trait) and rejected 251.
The reproducible stratified context audit reviewed 185 rows, of which 181 were
correct (97.84%). This is an agent context review, not an independent human
audit.

After direct conflict resolution, source-lineage deduplication, and
trait-specific Low rebuilding:

- direct species × trait: 44,248 → 48,570;
- Validated Low species × trait: 61,021 → 66,747;
- strict filled species × axis: 90,743 → 96,490;
- net coverage gain: 5,747 cells;
- strict coverage: 28.4563% → 30.2586%.

The final WFO parent artifact has 318,885 rows, 106,295 species, zero duplicate
species × axis keys, zero direct/Low species × trait overlaps, and zero
High > Medium > Low precedence mismatches.

## Europe PMC reproduction acquisition

Thirty high-information `genus × reproduction trait` queries returned 753
hits and 602 unique open-access PMCIDs. All 602 full-text XML documents were
fetched successfully at zero API cost. The acquisition is resumable by cached
query JSON and PMCID XML; the final source manifest hashes 636 files.

The exact-name/full-text gate emitted 43 candidates. Full context review found
17 correct species-level statements (39.53% candidate precision), demonstrating
that exact names plus nearby trait terms are insufficient without semantic
attribution review. Twenty-six rows were rejected for table-column bleed,
comparison-species attribution, or genus/parent/hybrid mismatch. Nine repeated
Panicum virgatum statements were retained as corroboration but not counted as
independent source origins. Eight species × trait rows were promoted.

The review also found one pre-existing direct error: a sentence saying that
switchgrass (*Panicum virgatum*), unlike *P. hallii*, is self-incompatible had
previously been assigned to *Panicum hallii*. It was explicitly invalidated,
and new source-backed evidence changed *P. hallii* from SI to SC.

After the reviewed update:

- direct species × trait: 48,570 → 48,574;
- new direct species × trait: 4;
- existing direct values corrected: 1;
- Low upgraded to direct: 3;
- Low invalidated: 245;
- newly filled species × axis: 1 (*Rumex bucephalophorus*, reproduction);
- removed noisy Low species × axis: 242, all Poa reproduction cells;
- strict filled species × axis: 96,490 → 96,249;
- final strict coverage: 30.1830%.

The negative net change is intentional evidence correction, not acquisition
failure: new Poa SC and SI counterexamples invalidate an old genus-wide SI
rule. Retaining those Low cells would increase coverage by preserving a rule
that now fails its direct-evidence contract.

## Reproduction commands

The GitHub workflow
`.github/workflows/acquire-source-scale-trait-evidence.yml` downloads formal
run `30440345597`, fetches Europe PMC XML, performs the complete 43-row review
contract, updates only affected `genus × trait_name` partitions, independently
checks denominator, precedence, uniqueness and hashes, and uploads the source
and integrated artifacts.

Local source and result directories:

- `data/v2/staging/traits/direct_llm_pilot/20260807_wfo_bulk_trait_candidates_round2`
- `data/v2/staging/traits/direct_llm_pilot/20260807_wfo_round2_context_gated_integrated`
- `data/v2/staging/traits/direct_llm_pilot/20260807_europe_pmc_reproduction_scale`
- `data/v2/staging/traits/direct_llm_pilot/20260807_europe_pmc_reproduction_integrated`

The raw Europe PMC XML is kept in the workflow artifact rather than committed
to Git. The committed article and query manifests, exact reviewed candidates,
source hashes, and deterministic review mapping are sufficient to audit and
rerun the integration.
