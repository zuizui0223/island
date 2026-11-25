# Effort-normalized positive visitor-guild portfolios and floral motifs — 2026-11-17

> **Stricter actual insect-species-degree source null (18 Nov 2026):** [Original 4505 GIFT × DoPI degree-fixed bipartite network result](../20261118_degree_fixed_dopi_visitor_network_null/README.md), original source [run #44213550467](https://github.com/zuizui0223/island/actions/runs/44213550467) SUCCESS. The new null preserves both every real insect species' observed plant degree and every plant's observed real visitor degree, beyond the looser guild-label marginal benchmark below. Open radial source guild-documentation residual remains negative **−0.188 [−0.305,−0.071]** but on 63 eligible source studies rather than 73, so the two samples cannot be treated as independent replications. The apparent insect guild breadth is not an effective pollen-guild or island-specialization estimate.

**Completed original source numerical analysis. Exploratory European positive flower-visitation data, not observed pollination on GIFT islands, effective pollen deposition, or island-specific evolution.**

## Original source and reproducibility

- Source-protected [original island Actions #44102655957](https://github.com/zuizui0223/island/actions/runs/44102655957) **SUCCESS**, source workflow commit be9ff9df4488ef098fa2f2c7fa9a2570dd9c6f10 and tested pinned fd model a68970e15c2598c68fa3a611de31c9de0c929a8b.
- Original anonymous artifact ID **14335684406**, extracted JSON SHA256 **1a45cd7bad541f9e8bba6feef2f911f34ca3215a4db2e144f2da6ed21a6c804f** (2,659,332 bytes). Original exact plant species, botanical island IDs, private site and publication IDs are NOT exported.
- Original botanical GIFT **4,505 island floras** retained; original accepted flora SHA256 6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a, species-Direct flower ledger SHA256 9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745, [DoPI v2 DOI](https://doi.org/10.5281/zenodo.22876365) frozen CSV SHA256 b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8, original **676,900** records.
- In this supplementary source analysis, **490,162** DoPI rows refer to a GIFT-island-listed accepted plant with at least one Direct phenotype; **150,368** lack the necessary original study/site/year/month; **32,703** source visitor labels are not genuine resolved species binomials and were rejected (for example Bombus sp.); **268,406** grade 1–3/site–month qualified rows remain before exact plant×visitor-species duplication collapse. **193** source paper IDs, **16,757** source-positive site/month units, and **1,611** source site/months with >=4 plants documenting >=2 distinct real visitor species per plant.

## Effort-conditioned estimand

Each source study × published coordinate/grid × year × month is analyzed separately. Use species-wide original GIFT Direct botanical categories for **4 one-axis, 11 two-axis, 9 three-axis, 3 four-axis motifs**. Only positive source flower visits (grade 1/2/3) from five named insect guild channels are used; grade 4 and unresolved binomials are excluded. A named insect species with multiple authority spellings contributes one positive species, not multiple guild records.

For each documented botanical plant i, n_i is its number of DISTINCT positively recorded visitor species, K_g the study/site/month number of positive plant×visitor-species edges in guild g and N the site's total edges. The exact finite hypergeometric **documentation-opportunity baseline** is:

E[G_i | n_i, K_g] = SUM_g [1 − choose(N−K_g,n_i) / choose(N,n_i)].

Observed G_i is the number of distinct positively recorded guild categories for the plant. Compare observed versus baseline guild counts between documented motif-positive and motif-negative plants from the SAME original source site/year/month. Require >=2 distinct visitor species on each analyzed plant, >=2 plant species in each phenotypic class, >=10 independent papers and >=10 real source papers containing positive multi-guild records. Average source sites within paper; bootstrap paper keys 999 times. **It is a combinatorial documentation/record-count standardization, NOT a biologically realistic rewiring null, true ecological preference or estimate of unrecorded guild absence.** Insect identity/degree, plant abundance and actual observation effort are not fixed by this null.

## Original actual source numerical outcomes (number of documented insect guilds)

All effects are motif-positive minus motif-negative botanical plant records, in **guild count units**, not percent or pollen transfer. Conditional 95% intervals use whole-source-paper resampling, no multiplicity control (27 overlapping morphology tests).

| Genuine Direct flower motif | Source papers / eligible source site-months | Observed difference | Expected difference from source positive record count | **Excess** [95% bootstrap] |
|---|---:|---:|---:|---|
| **Open-radial form** | **73 / 419** | −0.16225 | +0.01481 | **−0.17706 [−0.27080,−0.08511]** |
| Tubular / bell | 56 / 336 | +0.03956 | −0.03210 | +0.07166 [−0.02521,+0.16334] |
| Zygomorphic | 85 / 647 | +0.04485 | +0.00251 | +0.04234 [−0.03091,+0.12339] |
| Deep tube | 18 / 27 | −0.15772 | −0.16197 | +0.00424 [−0.12715,+0.13012] |
| Actinomorphic × tubular/bell | 42 / 154 | +0.02888 | −0.05131 | +0.08019 [−0.02875,+0.18822] |
| Zygomorphic × tubular/bell | 20 / 79 | +0.16108 | +0.04861 | +0.11247 [−0.02259,+0.29082] |

Only **6/27** motif panels reach the original source opportunity support criteria; **0/27** reach independent within-genus-variable source-study support. **No three- or four-axis true syndrome profile** reaches the primary source gate. Open radial source subset includes **536** distinct accepted plant names, **3,358** accepted plant×site-month positive source units (at least two visitor species), **1,969** with at least two positively recorded guilds, and **419** site/months from 73 original source papers. Within these positive observations, **only 1 source site/month** documents enough within-genus flower morphological variation to attempt a lineage-adjusted contrast: the effect is **NOT_IDENTIFIABLE** within genus.

### Critical real visitor taxonomic-quality falsification

A first successful but pre-final original source run (#44102023953) counted unresolved visitor labels as if they were independent named visitor species. This produced apparently positive tubular/bell source guild-breadth excess **+0.12610 [+0.03740,+0.21504]**, but once requiring REAL visitor binomials and normalizing author-suffixed duplicates, it fell to **+0.07166 [−0.02521,+0.16334]**, crossing zero. This is evidence that an unqualified positive tubular/bell conclusion was sensitive to taxon-resolution quality. Open radial remained negative (old −0.19386, corrected −0.17706), still only an exploratory source association and not proof of ecological specialization.

## Relation to island pollination syndromes

- **Not supported as a universal rule:** open radial or shallow flowers necessarily carry a larger realized pollinator-guild repertoire on islands. The external European positive-only records actually show lower observed/source-standardized GUILD *DOCUMENTATION* among open radial flowers, and lineage independence is not identifiable.
- **Supported limited observation:** some individual accepted plant species have actual positive records of visits from multiple guilds at matched source sites/periods (see [real co-visitation audit](POSITIVE_MULTI_GUILD_COACCESS.md)); those records do not establish effective pollen service or historical guild replacement.
- **Primary global botanical results retained:** original GIFT 4,505 floras, regional source-relative real same-species floral morphology pairs, taxonomic turnover, functional overlap and exact signed genus accounting. These are flora-composition results, not island-specific local evolution.
- The DoPI source is principally EUROPEAN, not a measured tropical bird/butterfly/hawkmoth guild community. Different visitor and botanical lineage/documentation opportunities remain confounded; positive guild labels cannot justify treating undocumented animal guilds as absent.

**Manuscript position:** auxiliary effort/documentation stress-test of a proposed broad-access functional mechanism, with transparent limits and negative/counterexample results. Do not turn this source visitation test into a causal or confirmatory island adaptation claim. Original PR #16 remains DRAFT.

[Reproducible finite-margins source model](../../src/fd/dopi_guild_portfolio_effort_null.py) · [Tests](../../tests/test_dopi_guild_portfolio_effort_null.py) · [Protected source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_dopi_gift_guild_motifs.py)
