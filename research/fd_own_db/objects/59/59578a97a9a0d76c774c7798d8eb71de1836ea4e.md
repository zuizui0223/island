# Source-pool botanical sorting versus an actual island pollinator shift

**12 November 2026 | fd draft PR #14 | source-verified PRIVATE original species/island analysis, public aggregates only.**

## One question, three fundamentally different measurements

**Question:** When a candidate island adequately lacks a Bombus GBIF record but has a flower-visiting bird occurrence, is there evidence that the same plant switched effective pollinators, or only that the *island plant flora* differs from potential mainland GIFT sources?

The data contain three non-interchangeable layers:
1. **Direct single-state species flora traits** (flower colour, form, symmetry) for 2,234 globally documented plant species. The state is attached to **one accepted plant name worldwide**, so it CANNOT measure a mainland-to-island phenotypic change within that species.
2. **GIFT five-nearest native mainland source checklists**, `geo_k5`; potential dispersal opportunities, not a reconstruction of historical colonization/ancestor or effective pollinator community. A match is a source-list **name match**, not source-population fieldwork.
3. **GBIF exact-island guild status**: `retained` = recorded animal occurrence, `disrupted` = source-qualified adequate GBIF nondetection, `unresolved` = not assessable. Even retained bird occurrence does not mean a bird visited any recorded plant. Global **GloBI genus** documentation shares come from completely different source lists and are NOT local guild abundances.

Every statement below is about descriptive opportunity or plant-community assembly. No historical pollinator shift is certified.

## Focal matched cohort and actual source-name overlap

The frozen `island` original progressive checkpoint (SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`) and original pollinator GBIF evidence qualification (artifact `10314984302`) yield **928 islands** with ≥30 strictly measured Direct species colour/form/symmetry. All five pollinator group states are available as retained / qualified nondetection / unresolved.

- **22/928** islands adequately lack Bombus records. **21/22** have flower-visiting bird occurrence; the same 21 also have other-bee, Lepidoptera (butterfly + moth pooled) and Diptera occurrences. **All four alternative guilds are recorded**, so a one-for-one bee-to-bird switch cannot be inferred.
- The 21 candidate islands represent **14 spatial blocks**: 15 tropical islands/8 blocks, 5 northern midlatitude/5 blocks, and 1 southern extratropical/1 block.
- Candidate island flora include **1,125 island×plant species occurrences**, **375 unique species / 277 genera** with all three Direct morphology states. **251** occurrences have strictly recorded island-native status; **793** have *unresolved* origin status (not all 1,125 can be asserted native).
- **586/1,125** direct species×island occurrences match an identical accepted plant name in the candidate five-nearest **native GIFT mainland** lists; among island-native records **197/251** also match the same species name in those source lists. Of the 21 focal islands, 13 have ≥1 strictly native species sharing the same mainland source-list name.
- Across the 21 candidate islands, the median per island is **53** Direct complete-trait flora species, of which **30** match a candidate mainland source-list species; the corresponding potential GIFT source has **163** complete Direct trait species at median. The source list is a geographic opportunity, not guaranteed actual maternal ancestry or visitation.

**Feasible future fieldwork:** The **197** source-matched, strictly native island×species occurrence records collapse to **54 unique candidate native plant species in 49 genera** across **13 islands**. Of those 13 islands, **11** have ≥5 and **8** have ≥10 such measured, source-matched native species. At the GENUS level only, 33/54 target species have a matching original GloBI genus with ≥3 independent references, and 15/54 have a genus with globally documented bird references (18 Bombus, 30 non-Bombus bees, 32 Lepidoptera, 24 Diptera; categories overlap). These are an informative prioritized *sampling universe*, **not** 54 known bird-pollinated native island species. Counts of 197 or 54 should never be read as verified original pollinators, independent source publications or geographic flower phenotype contrasts. The [anonymous fieldwork-feasibility receipt](native_same_species_GIFT_source_fieldwork_feasibility.json) documents these denominators.

Those same-name matches identify a manageable, evidence-backed *candidate pool* for same-species mainland/island visits and flower measurements. They do not represent 197 distinct taxa (counts are island×species occurrences), verified original pollinators, or geographic flower phenotype contrasts.

## Source-pool floral morphology contrasts: a floristic-composition result

For each island, use exactly the original Direct species-level categorical state labels, aggregate among observed island floral taxa, and subtract the average state proportion in the union of its five GIFT native mainland source checklists. This is a **species composition** comparison: the same species has precisely the same assigned flower morphology in both locations.

Among the 21 Bombus-adequate-nondetection/bird-recorded candidate islands:

| Morphological trait | Island median minus own five-source pool | Islands with positive source-relative difference |
| --- | ---: | ---: |
| Open-radial floral form | **+0.0431** (4.31 percentage points) | **16/21** |
| Actinomorphic symmetry | **+0.0908** (9.08 points) | **15/21** |
| Green/brown inconspicuous colour | +0.0041 (0.41 points) | 11/21 |

Source-verified native-only floral species are too sparse for independent generality: only 13 focal islands retain ≥5 strictly native direct-trait species. Of those, 9/13 show positive open-form source-relative difference and 11/13 positive actinomorphy difference, but their spatial block coverage is only six blocks.

The *14-block-equal* candidate open-form excess is **+0.0107** with conditional 4,999 block-resampling 95% range **[−0.0507,+0.0722]**; actinomorphy +0.0264 [−0.0149,+0.0671]. **Neither is robust against spatial-block weighting**. Similar positive source-relative open-form differences occur among many islands where both Bombus and birds are recorded, so the 21-island direction is **not a specific pollinator-loss effect**. The strict tropical detected/undetected contrast lacks positivity (only 2 Bombus-detected tropical islands versus 16 adequately nondetected); northern midlatitudes have 480 detected versus only 5 adequately nondetected. The data do not identify a causal local shift.

## Is marginal floral diversity maintained while recorded guild channels change?

To connect functional-trait diversity to the potential mainland source, the same **21 candidate islands** and five GIFT source pools were additionally rarefied to **exactly 30 Direct fully trait-complete species each**. Perform **120 independently seeded species subsamples** per island and per source pool; Rao Q averages equal-weight Gini–Simpson diversity across colour/form/symmetry, while colour–architecture coordination uses the chance-corrected mutual information `I(colour;[form,symmetry])` with **29 within-subsample colour-label shuffles**. This is an exploratory island-flora-versus-potential-source botanical assemblage comparison, not an independent source-population flower measurement.

| Source-relative quantity, n=21 islands / 14 blocks | Island-equal median (island−source) | Islands above source | Equal-spatial-block mean and conditional 95% range |
| --- | ---: | ---: | --- |
| **Morphological Rao Q** | **−0.018858** | **7/21** | −0.001569 **[−0.018059,+0.014363]** |
| **Excess colour–architecture mutual information**, nats | **−0.026626** | **9/21** | −0.015501 **[−0.049419,+0.018282]** |

Neither functional-related floral disparity nor multi-trait coordination is robustly higher in the candidate island flora than in its potential mainland source pools. With a different unbiased full-composition Rao Q estimator on all 928 eligible islands, tropical Bombus-undetected/bird-recorded candidate islands (15/8 blocks) have median source-relative Q **−0.0376**, while tropical observation-unresolved islands (88/40 blocks) have median **−0.0104**. This is a **context comparison only**: the unresolved group is NOT a Bombus-present biological control, local effort differs, and the 5-region source pool is generally richer than the recorded island flora.

Hence neither increased flower morphological hypervolume nor increased syndrome coordination provides a demonstrated compensatory signature for the hypothetical change in pollinator channel. Source-pool guild-documentation contrasts in the next section can coexist with lower morphological Q because **which species are present** and **how diverse their trait states are** are different assemblage properties. In particular, no island-versus-mainland *same-species* trait evolution can be estimated from globally fixed accepted-species trait categories.

[Anonymized 21-island matched rarefaction](source_matched_21_island_RaoQ_trait_coordination.csv); source-covering reference measurements are not the actual ancestral community or observed pollination effectiveness.

## Global GloBI genus reference composition: descriptive evidence of source-pool floristic sorting

Use the original GloBI 2026-06-04 genus breadth v2 source archive (SHA256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`), requiring ≥3 independent references per botanical genus. Each island and mainland source has its own **plant species weighted mean of the genus's worldwide GloBI documentary five-guild reference shares**. These are NOT the local visitor composition of those islands, nor an independent source-provider-controlled survey.

Among the 21 focal islands, comparing mean documented genus reference fractions with the corresponding potential mainland GIFT source flora:

| Documented reference channel | Median island minus source pool | Positive focal islands | 14-block-equal mean and exploratory 95% interval |
| --- | ---: | ---: | --- |
| Bombus | **−0.0294** | **3/21** | **−0.0579 [−0.0908,−0.0265]** |
| Other bees | **+0.0398** | **17/21** | **+0.0302 [+0.0065,+0.0527]** |
| Lepidoptera (combined moths and butterflies) | +0.0157 | 12/21 | +0.0292 [−0.0055,+0.0594] |
| Flower-visiting birds | +0.0017 | 11/21 | +0.0037 [−0.0146,+0.0208] |
| Diptera | −0.0005 | 10/21 | −0.0052 [−0.0221,+0.0116] |

Conditional block intervals leave **source-identical species, repeated genera, overlapping GIFT donor lists and globally shared GloBI literature providers** uncontrolled. The underlying v0.9 GloBI source audit also found **95.96%** of genus×guild×citation units concentrated in one archive-first source label, with substantial change on exclusion. Thus these literature profiles are **supporting descriptive composition proxies only**, not generalizable local pollinator matching or functional fitness.

### Why this is *not* a bee→bird switch

Break the candidate plant flora into **shared species** in both island/source lists, **island-only recorded** species, and **source-pool-only** documented species. With globally fixed source genus GloBI reference shares, the 21-island averages are:

| Plant-list portion | GloBI Bombus reference fraction | Other-bee fraction | Bird fraction |
| --- | ---: | ---: | ---: |
| Shared species | 0.078 | 0.323 | 0.057 |
| Island-only list species | 0.115 | 0.330 | 0.056 |
| Source-pool-only species | **0.153** | 0.287 | 0.062 |

The source-only species have, on average, a more Bombus-weighted *global reference profile*; both shared and island-only plant species have more other-bee-weighted global documentation. The documented bird fraction does **not** systematically increase. The shift in these composition means results from which plant taxa occur in each list and their weighting, **not from measuring a within-species change in flower shape or an animal visitor**.

## Two biological hypotheses with nonoverlapping evidence demands

**(A) Community sorting / channel gating:** a regional plant pool yields a different floral and genus-associated reference profile on islands, even though each species retains its globally recorded trait state. The source-pool descriptive comparisons above are compatible with this; their ecological mechanism and independence remain untested.

**(B) Real pollinator shift / compensatory generalization:** the SAME plant species, or a phylogenetically justified island sister lineage, has an actual flower-linked pollinator guild composition/effectiveness that differs on the island versus mainland. A bird may replace bees, or multiple flies, bees, moths and birds may compensate without strong morphology evolution. The existing GBIF guild record and GloBI global documentation **cannot test this directly**.

To support (B), a row must identify the same island/native plant and mainland source plant population, original visitor species and guild, actual visitation to THAT flower, standardized survey effort including adequate nondetection, and pollen deposition/seed set (or mark effectiveness unknown). The flower itself must be independently measured on both populations to claim floral trait evolution. `source_state=available` does not establish ancestral bee pollination.

## Decision: keep the positives without overstating pollination syndromes

The fd **structural floral results remain positive**: flower trait diversity and nonrandom colour×architecture pairing are distinct; across 928 Direct complete-flora islands, genus-linked floral form structure is widespread. The new candidate study finds **same-species mainland source opportunities** and descriptive flora-composition/visitor-documentation differences, but **none is a verified pollinator shift, local ecological functional diversity outcome or selection signal**. The earlier family/source-omitted GloBI prediction work also shows that genus-specific source documentation is highly source-dependent.

The next genuine macroecological test should compare **floristic community assembly/lineage sorting** versus **actual island population-level visitor replacement** with the same genus and repeated source regions, using independent island-local flower-specific observations and pollen effectiveness. No amount of re-shuffling existing island-wide GBIF animal occurrence statuses can substitute for that.

## Code and privacy

- `src/fd/pollinator_shift_source_pool.py` performs original five-nearest native GIFT source-list matching with *indexed region→species set*, **not** a giant island×mainland species join, accepts only original species-direct high/medium single-state categorical evidence, and returns ONLY aggregated region/evidence classes with ≥10 islands / ≥5 distinct spatial blocks, plus anonymous candidate cohort receipt.
- `tests/test_pollinator_shift_source_pool.py` synthetic five-condition tests: no species-level-shift claim, qualified GBIF nondetection, correct five GIFT sources, ignore imputed genus-consensus traits, aggregate privacy, duplicate GloBI genus rejection.
- Reproduced original matched candidate counts **928 islands / 21 candidates / 375 species / 277 genera / 1,125 island-plant rows / 251 native rows / 197 native source-matched island-species rows** in a second aggregate-only implementation. Original island IDs, exact plant/source genus labels and GloBI reference IDs are never pushed to public fd.

**Boundary:** do not merge into a paper with a claimed bird-pollination shift without a source/provider-disjoint plant×animal interaction and flower phenotype measurement.
