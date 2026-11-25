# 2026-11-13 EuPPollNet v1.3.0 actual island plant–pollinator network overlap

**Executed original-data geospatial/taxonomic audit, not demonstrated pollinator decline, true shift or pollen-transfer service.**

Source: source-corrected [EuPPollNet v1.3.0](https://zenodo.org/records/15183272), CC BY 4.0, [published study](https://doi.org/10.1111/geb.70000), exact original GitHub source `joseblanuza/EuPPollNet` tag v1.3.0. Original `Interaction_data.csv.gz` blob SHA `6e639d4af81feee76eefa6804622f8b901c39feb`, byte SHA256 `cde17692dc66946002e19026433d9db429357d28873a369d884941f71131db61`. Original `Metadata.csv` blob SHA `f9708e2f2b142ab2437a8be9d16e33aaefdef0ed`, SHA256 `adb396f31e7e7e3257425872b92d56c1051ca05849c6228a39675a71364d4833`.

The actual [original-source `island` Actions #43354857729](https://github.com/zuizui0223/island/actions/runs/43354857729) **SUCCEEDED**, releasing only [anonymous summary artifact #14005368043](https://github.com/zuizui0223/island/actions/runs/43354857729/artifacts/14005368043). Script: [original-source runner](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/scripts/run_fd_euppollnet_v130_overlap.py), reusable [fd join and original-guild checks](../../src/fd/euppollnet_island_overlap.py).

## Original source provenance and non-fabricated units

| Actual original input/evidence type | Count |
| --- | ---: |
| **EuPPollNet v1.3.0 source interaction rows** (duplicate-corrected; not old pre-correction 1.16m) | **623,476** |
| Legacy invalid-UTF8 source rows removed; taxonomic names NOT guessed | **8** |
| Dated/georeferenced source interaction rows passing original source quality | **623,468** |
| Unique **study × network × coordinates** source site keys, NOT all unique geographic localities | **1,638** |
| Source original research studies | **54** |
| Such network-site keys intersecting any original verified GSHHG island polygon | **373 keys on 13 distinct islands** |
| Same-site original network keys on **eligible 928 fd Direct floral islands** | **11 islands**, all northern midlatitude |
| Distinct spatial blocks and studies among these 11 | **4 blocks, 12 studies** |
| Observed visitor–plant interaction **rows** inside these eligible island polygons | **38,612** |
| Original source study × network × year cells inside eligible islands | **453** |
| Interaction rows with exact plant species name also recorded on same original fd island flora | **2,730** |
| Of these, rows with strict source-native original floral species | **2,170** |
| Distinct actual visited species matching original fd floral trait-complete flora names | **74** |
| Original eligible island sites with observations in >=2 calendar years | **1 island** |
| Source study×network×coordinate keys repeated across years on eligible islands | **90** |
| Geography/source independence support rule >=10 eligible islands, >=5 spatial blocks, >=3 studies | **0/4 regions pass** |

A single observational row is **not** one independently sampled pollinator or confirmed effective pollen deposition. The 1,638 site *keys* are not in contradiction with upstream publication's 1,214 geographic *locations*: a geographic location may have multiple study/network keys. The 90 repeated network-coordinate keys cluster in only **one eligible island**, not 90 independent island histories.

The 38,612 island source visit rows include **5,021 documented Bombus relations**, **11,668 other-bee relations**, **6,080 Diptera**, **1,118 Lepidoptera (butterflies/moths still unresolved)** and **12,833 Coleoptera**; other/uncertain groups remain separate. These are interaction records from nonuniform original observational protocols. **They cannot be used as five independent guild sampling-effort-standardized island abundance proportions.**

EuPPollNet is a superior *actual documented flower-visit* evidence class to global GloBI genus-level documentary guild references. But that does not confer measured pollen-transfer effectiveness. The source genus/functional-guild determination is constrained by recorded visit morphology/taxonomy, not inferred visitor assignments from flower colour.

## Ecological question and exact limits

Can guild loss cause decline in floral functional diversity? **Not identified by this original network cohort.** There are 11 eligible islands from 4 blocks, not 11 spatial/source-independent tests, and a single island with source-repeated years. Different studies also have distinct pollinator order targeting and survey protocols. The original fd floral colour, form and symmetry labels are one fixed species-level state, not island-population phenotypes measured in each visit study.

An exploratory, **explicitly simulated guild-removal** network resilience diagnostic can be built on these observed visitor links, but it must distinguish `would retain documented non-removed visitor links` from `would retain effective pollination or seed set`. No network guild-loss simulation should be presented as an observed historical visitor switch or evolutionary response.

### Next source opportunities

1. Check strictly comparable within-study sampled channels and network×site×year richness before any guild-removal simulations. Do **not** confuse a study deliberately sampling only bees with actual local lack of flies or butterflies.
2. Evaluate actual sampled visits to **identical island flora species** versus nearby non-island network sites only if they were surveyed with comparable original methods, and independently establish which are true mainland source habitats.
3. Obtain independent real same-plant visit and pollen-success surveys across >=5 geographically separate source blocks before claiming generality. In the currently source-verified full 928 island flora panel, no historical Bombus decline has been observed; the separately executed 80-island original GBIF year-window preflight yielded just **3** minimally comparable islands.

## Publication and privacy boundary

Public data attribution: Lanuza et al., *Global Ecology and Biogeography*, DOI 10.1111/geb.70000, source Zenodo 10.5281/zenodo.15183272 under CC BY 4.0. Original restricted `island` flora records, island IDs/polygons, raw third-party network study IDs, species names, original fine-scale sites and coordinates were kept in the original CI source job. Only anonymous region summary in fd GitHub; no republishing original protected species-identity rows.

**Main FD paper remains a source- and nativeness-conditional island floral morphology/reproductive repertoire biogeography analysis, not proven pollinator-loss-induced functional diversification.**
