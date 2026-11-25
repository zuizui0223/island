# 13 November 2026: actual European flower-visitor network guild-removal sensitivity

**Executed original EuPPollNet v1.3.0 source data; simulated visitor network perturbation, NOT observed historical pollinator decline or measured effective pollination.** Source experiment: successful [island Actions #43356298751](https://github.com/zuizui0223/island/actions/runs/43356298751), anonymous original aggregate [artifact #14005374360](https://github.com/zuizui0223/island/actions/runs/43356298751/artifacts/14005374360). Parent source-geography and plant-flora gate: [EUPPOLLNET_V130_ORIGINAL_NETWORK_OVERLAP.md](EUPPOLLNET_V130_ORIGINAL_NETWORK_OVERLAP.md).

## Identification and contrast

From 623,476 source-corrected original EuPPollNet flower-visitor interaction rows, 38,612 date/site network interactions occur on 11 of the **928 directly flower-trait-complete original fd island floras** (all 11 northern midlatitudes, only 4 original geographic blocks). Reuse *only* source study×network×year×coordinate networks with (i) >=10 observed unique plant species, (ii) >=2 distinct observed visitor orders and (iii) study metadata explicitly reporting that both bees/Hymenoptera and an alternative Diptera/syrphid group were in scope. Exclude studies that looked only for bees: a missing alternative category would not prove ecological absence.

Two **separate in silico removal** scenarios:
1. Remove **all Bombus visitor records** for the network, retaining every other documented visitor link.
2. Remove **all Bombus plus other documented bee taxa**, retaining explicitly documented nonbee links.

Outcome = proportion of original network plant names retaining *any* documented other group interaction, not actual seed set, pollination efficiency, evolutionary buffering or a field-measured guild shift.

**Negative-control stress:** for each network independently, remove the **same number of unique originally documented visitor species** as in the targeted guild, selected uniformly from all its documented visitor species, 199 deterministic random subsets. This is a descriptive network-specific *equal-taxon* contrast, not a randomized ecological experiment; visitor abundances, detection and species-degree distributions may differ.

## Direct empirical source-backed results

| Scenario | Networks / islands / original studies / spatial blocks | Median exposed plants with alternative recorded visitor | Median **all** visited plants losing all their documented links | Island-equal mean loss | Study-equal mean loss |
| --- | --- | ---: | ---: | ---: | ---: |
| **Bombus loss simulation** | **57 / 9 / 11 / 4** | **61.54%** | **18.75%** | **9.63%** | **17.98%** |
| **All bees loss simulation** | **76 / 11 / 11 / 4** | **58.95%** | **26.14%** | **32.58%** | **31.79%** |

These weighting estimands differ deliberately: networks, plants, islands and studies have very unequal sampling representation. The two scenarios are nested/overlapping networks and therefore cannot be treated as independent replications. Four original spatial blocks are insufficient for a confirmatory geography-disjoint effect.

### Is target-group dependency stronger than generic documented visitor removal?

Compared with a same-network **199-draw random deletion of an equal number of visitor SPECIES**, the target guild-specific fraction of documented plant names losing every alternative link is higher:

| Network stress | Median observed-minus-random excess plant disconnection | Equal-study mean observed-minus-random |
| --- | ---: | ---: |
| Bombus removal, 57 networks | **+0.03769** (3.77 percentage points) | **+0.05491** |
| All-bee removal, 76 networks | **+0.09360** (9.36 points) | **+0.05965** |

**Positive descriptive network result:** bees (and specifically Bombus in the overlapping selected subset) occupy a nonrandom share of the networks' documented plant visitor connections. Their hypothetical loss disconnects more observed plant links than a random removal of the same number of other documented visitor names in this source-specific cohort. This is network *structural sensitivity* and **potential visitor-link dependence**, NOT the effect of historical insect declines or fitness loss. There is no blocked/spatially independent CI and no causal p-value.

### Which alternatives were actually documented?

Nonexclusive medians of the fraction of *removed-guild-exposed plant names* that also have a documented visitor of each remaining category within the same source study×site×year network:

| Alternative visitor group | After Bombus removal | After all-bee removal |
| --- | ---: | ---: |
| Other bees | **33.3%** | — (removed) |
| Diptera | **38.1%** | **40.0%** |
| Coleoptera | 0.0% (network median) | **24.6%** |
| Lepidoptera (moths+butterflies pooled) | 0.0% | **11.3%** |

These shares can overlap within the same plant and are not additive. Even if another group visited a flower, no source-specific pollen deposition, fruit set, temporal compensation or successful guild replacement has been verified.

## Critical falsification of the desired floral functional-diversity impact

On the exact same island/plant morphology data, **0/76** all-bee removal networks and **1/57** Bombus removal networks retain eight directly colour/form/symmetry-complete plant species both before AND after hypothetical visitor-guild loss. Therefore **joint-profile Rao Q, 17D hypervolume, flower-pairing MI, reproductive-assurance change and simulated pollination service cannot be meaningfully estimated from this matched subset**. Lowering the 8-plant requirement posthoc to obtain a favourable result would not establish functional diversity. Direct plant source records are structurally sparse at original site×network×year grain even though 74 accepted directly measured species join across all 11 islands.

**Scientific inference boundary:** the robust current island plant-side result remains **source-relative recorded morphological/reproductive trait composition**, but its causal connection to declining/shifted local pollinator communities is not identified. The new network stress test advances the *plausible functional partner-redundancy mechanism* without confusing documentation with ecology.

## Reproducibility

- [Original fd network module](../../src/fd/euppollnet_guild_loss_simulation.py), [pure synthetic tests](../../tests/test_euppollnet_guild_loss_simulation.py), [source-pinned island execution](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/scripts/run_fd_euppollnet_v130_overlap.py), [successful Actions #43356298751](https://github.com/zuizui0223/island/actions/runs/43356298751).
- Original source Lanuza et al. *Global Ecology and Biogeography* (2025), https://doi.org/10.1111/geb.70000; v1.3.0 Zenodo https://doi.org/10.5281/zenodo.15183272 (CC BY 4.0). No raw source reference, island IDs, coordinates, or species identities are in the fd aggregate.
- Original multi-guild metadata pre-screen excludes bee-only original studies, but **cannot** equalize all observation durations, flower detectability or original fieldwork methods; source researcher and source geography effects remain.

**Status: validated exploratory in silico bee-guild structural-dependency endpoint, NOT confirmed pollinator-induced change in island floral FD.**
