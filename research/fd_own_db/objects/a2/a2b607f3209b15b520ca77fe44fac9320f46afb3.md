# Functional-trait turnover under island taxonomic turnover

**Status (2026-11-10): source-pinned private-input empirical prevalence/coverage and sequential equal-carrier rarefaction diagnostics completed; these are exploratory descriptions, NOT independent functional-diversity significance, geological causality, or evolutionary findings.** See [restricted-source aggregate results](../results/20261110_functional_turnover_readiness/private_functional_rarefaction_diagnostic.md).

## Scientific question and independence from island Chapter 1

The existing fd category-profile analyses measure **recorded** colour/form configurations, not actual pollination performance. The upstream island Chapter 1 already tests reproductive-assurance, floral architecture, and pollen-limitation associations with isolation. This independent fd extension asks:

**When different islands replace constituent plant taxa, is reproductive-assurance or floral-access structure retained more than lineage membership, geographic source availability, and observation coverage predict?**

That question is distinct from reprinting the Chapter 1 isolation slopes. Functional similarity among different species is a hypothesis, not a demonstrated biological equivalence.

## Input dimensions

Use seven source-verified, binary plant traits. Analyse each separately; do not invent a universal mixture of colour, attraction and selfing:

| Domain | Input columns | Interpretation |
| --- | --- | --- |
| Reproductive assurance | self_compatibility, selfing_mating_system, autonomous_selfing | Capacity to reduce dependence on mates/pollen vectors; **not** directly measured seed set |
| Floral accessibility | generalized_form, actinomorphic_symmetry, shallow_open_tube | Visitor access proxies; **not** realized visitor identity |
| Signal proxy | plain_colour | Visual trait proxy; **not** measured pollinator response |

Independent GloPL pollen-addition experiments may be used for **external, species-matched functional interpretation** but are not contemporaneous fd island-pair reproductive outcomes.

The upstream corrected island database includes 106,295 accepted species and 222,688/318,885 resolved species-axis cells. The raw reproductive-assurance axis resolves only 48,497/106,295 species. These are upstream denominators, not proof of usable complete support in fd near-island matched comparisons. Avoid global complete-case requirements that leave tiny, unrepresentative support.

## Rights-protected staging contract

A local, approved private CSV must have **one row per island_id x taxon_name** in an identical, explicitly documented All *or* Direct evidence scope. Required fields are:

- island_id, taxon_name (accepted species/name identifier), genus;
- seven listed binary traits (0/1/empty). Empty means **unknown**, not absent or zero.

The paired-island CSV has island_a, island_b, comparison, with donor categories defined *without* looking at floral outcomes. Never combine All/Direct as independent replications. A shared taxon must not have inconsistent trait values or genus between islands within one staging scope.

A runnable diagnostic that preserves these rules is in src/fd/functional_turnover.py. It measures observed trait prevalence gaps for **all taxa**, **each island's species-exclusive taxa**, and **each island's genus-exclusive taxa**. The last two explicitly remove shared taxonomic identities. Every species-level datum stays in local memory; by default each trait requires >=20 observed name carriers per island. Unresolved pairs remain missing, never scored zero.

Example:

    python -m fd.functional_turnover \
      --private-traits /protected/accepted_all_traits.csv \
      --private-pairs /protected/outcome_blind_island_pairs.csv \
      --output work/private_fd_functional_audit \
      --min-carriers 20 --minimum-pairs 10

Output is only a grouped trait-scope x source-comparison table and status flag in an **ignored local work/ directory**. The code never exports island names, species names, genus names, or raw pair-level data. Do not add private input to any public workflow or release; conduct rights review before releasing any further aggregates.

## Three gates to a genuinely new macroecological finding

1. **Real coverage and attrition first.** Check how many physically distinct islands, named archipelagos, spatial blocks, unique species and genera actually survive for each trait *after* shared-taxon/genus removal. Confirm identical outcome-blind pairs and All/Direct scope sensitivity; do not assume the earlier 606 colour/form islands all have seven traits.
2. **Geography/source negative controls.** Compare local islands with geographically distance-matched other-archipelago islands *at equal taxonomic and trait-observation opportunity*, controlling area/climate and regional native source pools. Absolute trait-prevalence differences are descriptive; a negative gap is not by itself evidence of stabilizing selection or functional retention.
3. **Conditional inference, not volume chasing.** Hold observed island x taxon membership and shared-taxon structure fixed. Reassign *whole observed trait vectors* within defensible genus x source/region opportunity strata, retaining missingness, and use spatial components/source archipelagos as inferential units. A genuine functional-retention result must survive species identity removal, source/lineage control, and geographic independence. Document support loss and preregister tests before examining outcomes.

A color/form PCoA area or 17D KDE Gaussian volume is **not** a substitute for reproductive function. Pollinator guilds inferred from colour/shape and GloBI interaction occurrence are not measured pollination effectiveness. Chapter 1's already reported H3/H4 results cannot count as independent fd confirmations.

## Provenance and inference ceiling

- island Chapter 1 current: submission/chapter1_current/H1_H4_EVIDENCE_MAP.md, docs/DATABASE_BUILD.md; self-compatibility recurrent, corrected experimental pollen-limitation association, post-hoc function-trait consistency.
- fd's local floral profile result: results/20261109_archipelago_flower_transfer/README.md. Its 606-island local prediction does **not** exclude the same plant taxa.
- The continuing source/genus and joint-colour/form prediction work is in fd PR #11; its small coupling residual is a **phenotype prediction**, not measured functional replacement.
- The bundled synthetic tests validate correct missingness handling, species/genus exclusion, shared taxon consistency and only aggregated output. **They do not establish any real-data significance.**

**Stop rule:** No language about ecological niche replacement, pollination equivalence, adaptive syndrome evolution or universal function conservation without an independently verified biological functional endpoint and spatially/source-robust real-data result.

## 2026-11-10 real public taxonomic recording-frequency result

A new source-pinned **GBIF observation-frequency audit** was completed on the `island` public `data/v2/staging/gbif/collected/island_taxa.csv`, source commit `3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5`, blob `5ae92313517118f0dbb29d8863aa2e379ce88747`. It was verified against `collection_status.json` (blob `9b186799a43f9108cf62f6a54eef074f8d42d3c3`).

- **115,328 recorded candidate taxon names**, **1,039,757 island×taxon occurrences**, **4,549 islands with records** from 103 collected source blocks.
- **58,660 taxa (50.86%) were recorded on exactly one sampled island**; another 15,981 on exactly two, 17,184 on 3–5, 9,015 on 6–10, and 14,488 on 11 or more islands.
- The median number of sampled islands per taxon is 1, the 90th percentile 14, the 99th percentile 148, and the maximum 993.

These are **recording frequencies**, NOT evidence that half the taxa are true biological island endemics. Opportunistic GBIF sampling and unequal islands' taxonomic coverage could produce much of the concentration. The earlier 115,328 candidate-name GBIF campaign is **not** the corrected Chapter 1's 106,295 trait-analysis-applicable species. This public table contains **no per-island functional trait outcome**.

Aggregate-only output: [public_island_taxon_frequency.json](../results/20261110_functional_turnover_readiness/public_island_taxon_frequency.json). Reproduction script: [audit_public_island_taxon_frequency.py](../scripts/audit_public_island_taxon_frequency.py); the dedicated [public-only CI](../.github/workflows/functional-readiness.yml) checks the frozen values against both pinned files. The source status and the raw record lists stay in the upstream `island` repository.

## Matched local-versus-distance-control functional comparison (new executable module)

The [functional_matched.py](../src/fd/functional_matched.py) extension takes one outcome-blind **matched two-donor triplet per focal island**. A row of the private matched-pairs input must include:

- `match_id` and `island_a`: unique match and focal island identifier;
- `island_b`: donor island, with two different donor islands per focal;
- `source_group`: known focal archipelago;
- `donor_group`: known donor archipelago;
- `role`: exactly one `local` donor from the same named archipelago and one `distance_control` donor from another;
- `distance_km`: positive actual focal–donor distance; donor-distance ratio must be <=3, preassigned from geography, not floral outcomes.

No focal island may appear in more than one match, and all labels must be known. This is a **stricter and explicitly archipelago-identity-conditioned** test than a simple near-versus-far comparison; it does not by itself estimate the main geographic distance slope. All/Direct inputs are run **separately**, not treated as independent replications.

For each of the seven traits, score the same focal island against both donors under three scopes: all taxa, remove plant species shared in each donor comparison, remove genera shared in each donor comparison. Each scored trait needs **>=20 recorded carriers on EACH side of EACH donor comparison**. Define:

`local_advantage = |prevalence_focal - prevalence_control| - |prevalence_focal - prevalence_local|`.

A positive value means the local donor has a smaller observed prevalence gap than an otherwise distance-matched alternative. For the two primary domains, take the equal-weight mean over at least two of the three component traits **available on both donor comparisons**. The signal-only `plain_colour` is secondary. This measures traits of **different sampled taxa**, not equivalent pollination performance or selection. Deleting shared genera can create different selected taxon universes on the local and control comparisons; outcome-blind support diagnostics and balanced-source checks remain essential.

Example local rights-authorized execution:

    python -m fd.functional_matched \
      --private-traits /protected/accepted_all_traits.csv \
      --private-matches /protected/outcome_blind_triplets.csv \
      --output work/private_fd_functional_matched \
      --min-carriers 20 --minimum-matches 10 --minimum-groups 5

The CLI exports ONLY `functional_matched_grouped.csv` and `status.json`. The mean advantage weights known source groups equally; a group must have eligible matches, >=10 independent focal matches and >=5 named source groups are required to print a row. Raw ids, taxa, pair-level differences and donor names are never written by this CLI. The output is **descriptive and post-hoc**, with **no p-values**; source archipelagos and plant phylogenies may remain dependent.

The older colour/form geographic control retained only 51 focal islands in 18 named archipelagos after matching donor distances. This does **not** imply the same 51 have adequate reproductive or floral-access trait evidence; eligibility must be newly measured. Underlying Chapter 1 trait coverage includes just 48,497 resolved reproductive-assurance species-axis cells among 106,295 analysis-applicable species. A rights-protected seven-trait join and exploratory sampled-pair audit HAS now been run; only coarse summaries are released in [the 2026-11-10 functional readiness report](../results/20261110_functional_turnover_readiness/README.md). The distinct same-focal local-versus-distance-matched different-archipelago contrast remains **unrun** on the functional seven-trait input.

**Empirical next gate:** produce a privacy-safe table of supported focal matches and named source groups for each primary domain after genus exclusion, from an authorized private checkpoint; if too few remain, do not claim functional retention. If support remains, freeze a whole-trait-vector genus × source-availability conditional null and inference groups before fitting, with true ecological service interpretation kept separate from plant-side trait proxies.

## New 2026-11-10 empirical research gate — sampled genus-linked floral access

Using rights-checked archived original plant-trait/GBIF checkpoints and the **7,618** fd G0 sampled pairs, we successfully evaluated the Chapter 1 frozen seven binary traits in All/Direct scopes. Taxon/genus exclusions use **whole floras**; missing trait records are never false negatives. After requiring >=20 verified observed trait states on both sides, shared-genus-exclusive coverage is **5,054/4,414 reproductive-domain** and **5,121/4,753 floral accessibility** pairs (All/Direct), averaging at least two of three domain components. These pairs reuse islands and do not represent independently evolved lineages.

The interpretable **additional genus** diagnostic uses a sequential counterfactual. First remove identical species on both islands, then compare the observed gap after further removal of all shared genera to 399 draws that *randomly thin the already species-exclusive known trait carriers* to exactly the same final support sizes. The mean excess gap is **floral accessibility All +0.0182 / Direct +0.0166**, while reproductive assurance is **All +0.0103 / Direct +0.0024**. Positive excess means genus-linked taxonomic removal changes prevalence similarity more than random sampling loss alone. It is a **post-hoc sample-size control**, not an independent source/genus evolutionary null, geologic-origin test, or measured pollination service.

The full nonidentifying [source SHA/aggregate report](../results/20261110_functional_turnover_readiness/README.md), [machine-readable rounded regional results](../results/20261110_functional_turnover_readiness/functional_stagewise_real_aggregate.json), [sequential traitwise rarefaction helper](../src/fd/functional_rarefaction.py), and [domain summary generator](../src/fd/functional_support_null.py) are the reviewable output. No original source species/island identities, named private botanical data, or individual pair-level results are committed. Exact reproducibility requires user-authorized original source archives at their pinned workflow SHA.

**Revised priority:** Stop adding general KDE- or geology-only p-values. The next actual discriminating ecological test is same-target local-vs-distance-matched donors with shared species/genus exclusions, then source-availability and block dependence controls, *if* sufficient functional-trait carriers survive. The observed genus-linked accessibility tendency itself should not be called conserved realized pollination function.
