# 14 November 2026 — Original global GIFT floral source coverage versus GBIF background-record effort

**Original rights-protected source workflow completed successfully:** [island Actions #43593501472](https://github.com/zuizui0223/island/actions/runs/43593501472), anonymous source artifact **#14113536684**. Every original GIFT accepted plant/island source ID, original independent GBIF animal occurrence, background observation geocoordinate and GBIF source dataset identity remains inside the protected original CI. Only anonymous region/geology/guild and trait recording counts are exported.

## Input and scope

All **4,505 original island botanical floras** are retained. Four separate effort metadata columns from the original independently NEE-qualified GBIF channel evidence are used: background animal occurrence record count, independent broad-background spatial units, temporal units, and distinct GBIF datasets. These are documentary **coverage proxies**, not equal-effort formal pollinator surveys, true guild abundance, effective pollen delivery or historical declines.

The source Direct botanical trait evidence is unambiguous high/medium species-specific, original source species-wide data; **five separately recorded** trait categories are assessed without requiring any triple-complete floral evidence, plus independently paired colour×form (only TWO Direct flower traits required). Original GIFT native-status documented plants and all originally recorded flora are separate source tiers. No missing GBIF effort is set to zero; no guild absent is inferred from nonobservations.

The original private join has **4,505 original islands × five source-qualified guild channels = 22,525** island×GBIF source cells; 75 anonymous global/region/geology×guild summary panels are public.

## Original qualified GBIF state–effort imbalance

| Original five-channel GBIF guild evidence | n GBIF-detected source islands | n adequately nondetected source islands | Median background RECORD count on detected islands | Median background RECORD count on adequate nondetections | Median background spatial units detected / nondetected | Median temporal units detected / nondetected |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Bombus | **817** | **23** | **14** | **107** | 2 / 5 | 4 / 14 |
| Non-Bombus bees | **1,083** | **11** | **5** | **98** | 1 / 5 | 2 / 10 |
| Lepidoptera | 2,051 | **0** | 36 | no qualified nondetections | 2 / NA | 5 / NA |
| Diptera | 1,256 | **7** | **20** | **83** | 2 / 3 | 4 / 6 |
| Flower-visiting bird catalog | 2,718 | **41** | **292** | **147** | 2 / 4 | 4 / 10 |

**This is a mechanical selection issue intrinsic to the qualified GBIF nondetection scheme:** observing ONE focal guild record can certify presence at low background effort, whereas certifying adequate nondetection requires independent broad-background records spread across multiple survey space/time units. It is therefore not valid to treat retained vs qualified nondetection as if it were a controlled treatment or random ecological absence. The bird background record COUNT exception (292 vs147) reinforces that source effort is multidimensional—spatial/temporal/dataset diversity also matters. The source bird record numbers may reflect archive query limits/censoring, not true bird community magnitude, so keep them as reported documentary counts.

Discrete background record count bins were frozen **[0,10), [10,50), [50,200), [200,1000), [1000,+∞)**. Only TWO of these strata contain BOTH qualified observation statuses for Bombus, non-Bombus bees and bird-catalog guilds (Diptera only ONE, Lepidoptera ZERO because no adequately nondetected sources). Their overlap is at most a source-opportunity preflight; it does NOT ensure overlap within spatial block, area, candidate GIFT donor source or climate.

## Floral trait recording vs GBIF background source documentation

Descriptive Spearman rank correlation, using *only islands with both nonmissing original source-available/evaluable GBIF background records and original GIFT botanical-flora Direct trait coverage fractions*. All records are observational, correlated within source species, geography, guild background and active biological communities. No p values are promoted.

| Original guild / number of islands with evaluable GBIF background | Flower primary colour source-Direct recording fraction | Floral form source-Direct recording fraction | Autonomously selfing trait source-Direct recording fraction | Joint source-Direct flower colour × form recording fraction |
| --- | ---: | ---: | ---: | ---: |
| Bombus / **840** | +0.0299 | +0.0434 | **−0.2106** | **+0.1376** |
| Non-Bombus bees / **1,094** | +0.0102 | +0.0867 | −0.0414 | +0.0691 |
| Lepidoptera / **2,051** | +0.1425 | **+0.3269** | +0.1390 | **+0.3518** |
| Diptera / **1,263** | +0.0925 | +0.1355 | −0.0876 | +0.2031 |
| Bird-catalogued channel / **2,759** | +0.1233 | **+0.3284** | +0.2627 | **+0.3310** |

These are original GIFT plant trait **recording fractions**, not true native plant trait diversity and not confirmed flower visitors. A positive correlation may reflect greater zoological and botanical documentation intensity, biogeographic patterns, true biodiversity, or catalog construction. No within-geographic-block causal control has been achieved.

## Interpretation after the independent same-island/mixed-block guild FD result

The separate original same-cohort study [#43592661994](https://github.com/zuizui0223/island/actions/runs/43592661994), [full report](ORIGINAL_GUILD_SAME_COHORT_WITHIN_BLOCK_MULTIPLICITY_AUDIT.md), shows that Bombus presence vs adequate nondetection has almost no documented within-spatial-block contrast on source-comparable original island floras. The global raw colour×form separately-block contrast **−0.02404** cannot be interpreted as a guild effect; the same 50/16 source-comparable islands have raw contrast **−0.01696** and GIFT source-relative **−0.00538**, with ZERO blocks containing both states, while the FULL 760/23 sample has only FOUR mixed blocks. At six mixed first-GIFT source groups of the raw full cohort, raw colour×form effect **−0.00257**, cluster-sign descriptive p **0.78125**.

This independent GBIF effort source audit gives a direct observable reason to **downgrade pollinator guild ecology to a carefully qualified correlation and observation-process layer**. It does not prove all trait associations are artifacts; it demonstrates that their observation process and source-group compositional explanations cannot be ignored.

**Global manuscript main population remains 4,505 island floras**, 10 high/medium direct single traits and three true two-trait functional diversity endpoints. The strict 928 triple-complete flora islands are only a historical sensitivity, not a new primary sampling bottleneck. Geologically typed independently oceanic/continental island floras (425 with botanical data among 434 typed geometric islands) are secondary analytical strata; no independent geological moderation of the guild association has yet been identified from source-qualified nondetection.

## Reproducibility

[fd global GBIF background-effort × floral Direct documentation auditor](../../src/fd/global_guild_effort_trait_coverage.py), [unit tests](../../tests/test_global_guild_effort_trait_coverage.py), [island original-source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_global_guild_effort_trait_coverage.py), [original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-global-guild-effort-coverage.yml). The source original GBIF NEE projected channel has categorical source_state and channel_state PLUS all four independent background effort variables; prior categorical-only guild analysis discarded the latter. This audit restores them without falsely filling unknown with zero. Keep manuscript PRs DRAFT pending independent phylogenetic and actual visitor-efficacy evidence.
