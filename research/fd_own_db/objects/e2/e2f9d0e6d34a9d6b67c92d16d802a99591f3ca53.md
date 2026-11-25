# Same-species pollinator-shift feasibility: original GloBI v0.9 visitor evidence

**12 November 2026 · fd PR #14 · original-source, evidence-gated extension; ecological guild switching NOT established.**

## Why this is stricter than the original 928-island comparison

A geographical overlap between island plant species and a candidate GIFT native mainland list establishes a *source-opportunity NAME MATCH*, not that a documented island plant was visited by an independently identified pollinator or that a mainland population was ancestral. Similarly, a globally documented GloBI visitor guild associated with the **plant GENUS** is not direct evidence that an identified plant species or island population received pollen from that visitor.

We therefore rebuilt the **exact source-matched species** candidate set using ONLY:
- original Direct high/medium evidence with complete colour, form and symmetry;
- original source-verified native island flora, retaining only species present in one of the five geographically closest candidate GIFT native mainland source checklists;
- island record cohort with ≥30 directly trait-complete species;
- **qualified GBIF Bombus nondetection** and independent recorded flower-visiting bird animal occurrence.

This deterministic original frozen source selection reproduces **928 original eligible islands**, **21** Bombus-adequately-nondetected/bird-occurred islands, **197** native source-matched island×species rows, **54 accepted species / 49 genera**, on **13 focal islands**. Individual species, original island/source IDs and original GIFT checklist identities remain in the protected study workspace.

## Baseline direct floral / reproductive evidence among all 54 candidates

| Trait | Original high/medium Direct known species | Missing | Comment |
| --- | ---: | ---: | --- |
| Flower colour | **54** | 0 | Species-wide categorical record only, not island population photo phenotype |
| Floral form | **54** | 0 | Open radial in 9/54, not universal across candidate pool |
| Floral symmetry | **54** | 0 | Actinomorphic in 38/54 |
| Flower size class | 30 | 24 | Direct categorical size, not flower-specific tube/fit measure |
| Tube depth class | 22 | 32 | Not all source-island pairs can compare potential access morphology |
| Self-compatibility / incompatibility | **17** | 37 | Only 17/54 direct SC/SI records; **not selfing rate or reproductive assurance** |
| Autonomous selfing capacity | 12 | 42 | One original state is mixed/variable, thus 11 unambiguous single categories |
| Mating system | **8** | 46 | Too sparse to estimate macro-scale selfing-shift coevolution |
| Inflorescence display | 42 | 12 | One direct state is mixed, not force-coded |
| Cleistogamy | 7 | 47 | Incomplete |

Using the previously frozen **GloBI global plant-genus** reference breadth data only, 33/54 candidates have a genus with ≥3 independent documented visitor-interaction references, and 15/54 have a genus with some bird documentary records. **Only 3 of those 15** have original open-radial floral form. These are taxonomic/documentation overlaps, not observed visitors of the 54 individual plant species.

### Cross-evidence attrition inside the native 54-species cohort

The old *globally documented GENUS-level* GloBI visit/breadth ledger yields ≥3 reference support for **33/54** species' genera and any global bird-channel reference for **15/54** species' genera. Among those **15** genus-bird-documentation candidates, only **6** have Direct species-level SC/SI evidence and **5** have Direct floral tube-depth evidence. Of the full 54, **9** have original Direct open-radial form (just **3** of the 15 genus-bird documentation candidates). Strict Direct high/medium reproductive record coverage is therefore inadequate for a bird-versus-bee pollination-service or mating-system adaptation model. These statements do **not** indicate actual bird interactions with any of the 54 named species.

The three older morphology-template candidates are poorly discriminating in exactly this 54-species cohort: only **5** have ≥3 Direct informative traits for `bird_like` and **4** have `bird_like` concordance ≥0.75 on those ≥3 traits. Only **3** candidates have ≥3 Direct informative traits for EACH of bee-like, butterfly-like and bird-like, and all three have high score in **all three templates**. These are overlapping morphology templates, not verified guild membership. Anonymous cohort aggregate: [`native_54_morphology_template_ambiguity.csv`](native_54_morphology_template_ambiguity.csv) and [`native_54_species_direct_functional_trait_coverage.csv`](native_54_species_direct_functional_trait_coverage.csv).

## Exact plant species × visitor reference retrieval implemented

To avoid projecting genus labels onto species, we implemented original GloBI v0.9 source-stream parsing in [island Draft PR #323](https://github.com/zuizui0223/island/pull/323): `src/island_v2/fd_globi_species_shift_v09.py`, `tests/test_fd_globi_species_shift_v09.py`, and the existing v0.9 SHA-pinned workflow.

The scanner:
1. Independently reconstructs and freezes the **same 54 native source-matched accepted plant names** from the original progressive island source and archived GBIF observation-state files; fails closed on sample drift.
2. Reads the exact **GloBI 2026-06-04 v0.9 `interactions.tsv.gz`** (SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`) in the SAME source pass used by the established provider audit.
3. Matches the **full plant species name**, not merely genus, on the correctly oriented source or target side of each of the four frozen flower-visitation/pollination RO predicates. Rejects records lacking a real plant kingdom, an identified visitor species, or one of the five preselected visitor guilds (Bombus, other bees, undivided Lepidoptera, birds, Diptera).
4. Separately counts raw exact-species hits, available original **provider AND bibliographic reference** evidence, duplicate references shared across provider labels, distinct cited interactions by species and guild, records explicitly labeled `pollination` vs `flower_visit`, plus **availability of record coordinates/event times**. A coordinate is NOT treated as exact island polygon membership.
5. Returns ONLY anonymous counts: number of the 54 species with any source-cited interaction, ≥3 distinct cited references, ≥2 recorded guilds, ≥2 source-archive labels, any georeference/date, and each of five guilds' count. **Zero documented matching records means GloBI evidence absent, not true pollinator absence**.

Source/role/missing-documentation/reference reuse and candidate selection synthetic tests all passed locally; the GitHub workflow performs source-pinned original-data tests before publishing a restricted-ID-free receipt. At the time of this report, **the original v0.9 stream was still running; no empirical exact-species GloBI support count is claimed yet**. A positive global species-linked record, when found, is STILL NOT proof it belongs to the focal island or the candidate GIFT mainland source population.

## Biological question and decision conditions

The paper should treat three nonidentical hypotheses separately:
- **Flora sorting / source-pool filtering:** existing floral trait states and taxon frequencies differ between the island flora and geographically plausible mainland plant pools. Already observable, but no ecological mechanism established.
- **Real ecological pollinator switching:** the *same native plant species* or a supported island sister lineage has standardized and *location-specific* pollen-effective visits from different guilds in the source mainland and on the focal island.
- **Functional compensation / generalized visitation:** source pollinator reduction does not induce a one-for-one switch; multiple visitors could share pollination service on the island, possibly with no large change in floral structure.

The exact 54-species GloBI scan is **evidence feasibility**, not the endpoint: an independently verified island and mainland coordinate assignment or direct population interaction dataset is required to claim switching. Additional original flower phenotype **by population**, effective pollen deposition/fruit set or controlled pollen-limitation tests are required to make a claim about syndrome evolution or service compensation.

Only the first hypothesis currently has source-pinned descriptive results; the other two are not rejected, simply **unidentified in existing input schemas**. We keep fd PR #14 Draft and do not count globally genus-linked bird references as 15 actual bird-pollinated island species.

## References and provenance
- Frozen `island` progressive source archive `10058653212`, SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Exact-island GBIF qualified channel evidence artifact `10314984302`.
- GloBI source genus v2 artifact `10328696801`, SHA256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`.
- Original GloBI v0.9 pinned DOI `10.5281/zenodo.20546682`, SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`.
- Original GloBI documented interaction records include species identity, source/provider and bibliographic citation keys and optional `decimalLatitude`, `decimalLongitude`, `eventDateUnixEpoch`. Cite original GloBI contributing datasets and do not redistribute protected raw input; source-specific anonymity thresholds apply to future releases.


## Verified original-data result: 27 source-cited exact species

This gate **has now been passed at the GLOBAL species-name level**. Original source [island run #43144106425](https://github.com/zuizui0223/island/actions/runs/43144106425) **completed successfully**, with [anonymous artifact #13911307901](https://github.com/zuizui0223/island/actions/runs/43144106425/artifacts/13911307901) and SHA-verified exact GloBI June 4 2026 v0.9 `interactions.tsv.gz` (`095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`).

The frozen **54 native plant species whose accepted names match candidate mainland GIFT native source checklists**, from 13 of 21 candidate islands, had the following actual source-reconstructed evidence:

| Real source-v0.9 identity and provenance check | Observed value |
| --- | ---: |
| Original GloBI rows streamed | 24,577,183 |
| Flower-interaction predicate rows with exact plant-species match | **1,123** |
| Retained records among five chosen guilds, correct plant kingdom and named visitor species | **792** |
| Distinct exact plant species with ≥1 real bibliography-backed guild interaction | **27/54** |
| Plant species with ≥3 distinct bibliography keys | **12/54** |
| Plant species with two or more globally documented visitor guilds | **11/54** |
| Species with at least one explicitly classified `pollination` interaction | **7/54** |
| Species with explicitly classified `flower_visit` interaction | **25/54** |
| Species represented under ≥2 archive-first original provider labels | **5/54** |
| Species with any source coordinates | **27/54** |
| Source-eligible georeferenced interaction rows | **642/792** |
| Plant species with coordinates in at least two distinct rounded 0.1-degree grid cells | **13/54** |
| Distinct original species×guild×provider×bibliography×evidence-class units | **622** |
| Original v0.9 `eventDateUnixEpoch` available | **NO** (missing field, dating unassessable) |

**Exact species-based global documented visitor guild evidence** (guilds overlap, so counts cannot be added):

| Global documented channel | Species with ≥1 bibliography-backed record | Species with ≥3 bibliography keys |
| --- | ---: | ---: |
| Bombus | 4 | 3 |
| Other bees | 14 | 7 |
| Lepidoptera (moth and butterfly NOT separated) | 20 | 11 |
| **Flower-visiting birds** | **5** | **2** |
| Diptera | 5 | 1 |

The previously reported **15 bird-associated *genera*** among these 54 plant species were a different and weaker proxy; only **5 plants have independently retrievable exact plant-species–bird-channel GloBI records**, and 2 have ≥3 distinct bird bibliography keys. This is a major difference between genus-level extrapolation and actual plant identity. The **11 species with 2+ globally documented guilds** are candidates for **pollinator generalization/interaction rewiring**, but **not** evidence that distinct guilds visited the same location, that both effectively transferred pollen, or that a change occurred after island colonization.

Neither 622 source-documentation units nor 642 georeferenced relation rows establishes the number of independent field studies, and the earlier source v0.9 concentration audit found one archive-first provider label disproportionately dominates global reference units. Source crossprovider duplicates are excluded by the original audit; only **5/54 plant species have evidence associated with two or more provider labels**. Provenance unit definitions remain stricter than raw interaction row counts.

### True island-versus-mainland switch remains geographically unassessed

The **642** nonempty-coordinate records have NOT been spatially joined to the original frozen GSHHG island polygons and relevant candidate GIFT mainland regions; a globally georeferenced species record cannot be called an observation on any of the 13 focal islands by species name alone. The original frozen geometry artifact #8066083419 was checked and is **no longer accessible (expired 2026-10-01)**. Do **not** replace it with an unverified Natural Earth geometry, approximate island centroid radius or arbitrary map rectangle and label that a local visitor record.

If the original GSHHG 2.3.7 high-resolution polygon layer is reconstructed from its independently pinned archive SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`, island polygon equality and exact-island plant species overlap must be revalidated BEFORE any genuinely local interaction evidence can be counted. Even then, the v0.9 original has no event-time field and the `pollination` label may reflect an observer interpretation, not pollen transfer/fruit set.

### Current scientific conclusion

- **Confirmed source-backed result:** 27 native island/GIFT-source candidate *plant species* have original GloBI five-guild *documented global* flower interaction evidence, 11 have multiple guilds documented, 7 have explicit pollination-labeled rows, and 5 have exact species-associated bird-channel references.
- **Not confirmed:** a single 13-island plant species with both independently sampled island and mainland flower visitation guild composition or effectiveness; a historical bird-for-bee, moth-for-bee or multiple-guild compensatory **pollinator shift**; island-local floral form evolution. Data absence is **not proof of no true shift**.
- **New study design implication:** test the original documented **5 bird species**, **11 multiguild species** and **12 bibliography-rich species** for actual spatially linked source/island visitor information; do not select them because of favorable sign in a guild model and do not turn genus template scores into known pollinator identities.

[Aggregate-only five-guild data](native_54_exact_species_globi_v09_guild_reference_coverage.csv) · [frozen source receipt](native_54_exact_species_globi_v09_coverage_receipt.json) · [original code and linked Source PR #323](https://github.com/zuizui0223/island/pull/323). No exact plant taxon, coordinate, locality, original provider name or bibliography identifier was committed.


## Dominant GloBI archive omission: the exact species-level bird signal largely disappears

A SOURCE-v0.9 original [Actions #43160907960](https://github.com/zuizui0223/island/actions/runs/43160907960) **SUCCEEDED**, yielding anonymous [artifact #13917963603](https://github.com/zuizui0223/island/actions/runs/43160907960/artifacts/13917963603) with an exact-species (NOT genus) **original provider-label omission audit**. The **same 54 candidate plants** were retained as the candidate universe; no observation or flower form was used to choose a favorable source subset. The largest archive-FIRST label accounts for **95.016%** of the **622** exact species×visitor guild×provider×reference×evidence-class units. Only **8** original archive-first source labels were represented for these 54 species. Removing the dominant source label leaves:

| Named-plant GloBI original evidence | Before source-label deletion | After |
| --- | ---: | ---: |
| Species with one or more globally cited flower-interaction records | **27/54** | **7/54** |
| Species with 3 or more bibliography keys | 12/54 | **3/54** |
| Species documented globally with 2+ visitor guilds | 11/54 | **4/54** |
| Species with ≥1 explicitly `pollination`-classified record | **7/54** | **7/54** |
| Documented Bombus guild species | 4 | **0** |
| Documented other-bee guild species | 14 | **7** |
| Documented Lepidoptera guild species | 20 | **4** |
| Documented bird guild species | 5 | **1** |
| Documented Diptera guild species | 5 | **2** |
| Species with ≥3 bird-channel bibliography keys | 2 | **0** |

This is a decisive **data-source concentration / evidence-attrition finding**, not evidence that a biological pollinator guild was extirpated. The persistently documented seven species with explicitly labeled pollination claims are potentially useful follow-up cases but original `pollination` record types are not independently verified seed set, and one archive may group several contributing studies. The **genus-level** positive bird-reference prediction after dominant-source omission has a much weaker **exact species-level** supporting evidence base; these two distinct analyses should never be presented as equivalent.

*Research-safe conclusion:* The original data can identify some species with references to more than one visitor group, but **source-disjoint local visitor turnover is not measured**. The narrow observational claim `11 distinct plant species with ≥2 globally documented guilds` drops to `4` after removing a coarsely defined dominant archive-label group. It is more accurate to interpret this as **global documented multi-guild interaction opportunities** than an evolutionary pollinator shift or field-level functional compensation.

[Original five-group before/after counts](native_54_exact_species_globi_v09_provider_omission.csv) · [54-species source citation coverage](native_54_exact_species_globi_v09_guild_reference_coverage.csv). No original plant/visitor/provider identities have been exported.

## 12 November 2026 confirmed spatial overlap: exact native island plant×GloBI cited visit references

The original frozen GSHHG global 8,265-island prepared artifact #8066083419 expired on 1 October 2026. Rebuilding from the EXACT same GSHHG 2.3.7 high-resolution ZIP (SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`) yielded **8,260** geometry IDs, with **9 historical IDs not reconstructed and 4 additional IDs not historically present**. This is a real full-universe geometry migration mismatch; no full 8,265-island equivalence is claimed.

However, the strict [original focal polygon recovery #43162949729](https://github.com/zuizui0223/island/actions/runs/43162949729) **SUCCEEDED**: all **21 prior, source-selected** Bombus-qualified-nondetection/bird-occurrence candidate island polygon IDs were reconstructed *identically*. Only geometries carrying exact original IDs could be used; source-area centroids and alternate Natural Earth coastlines were prohibited.

The completed [SHA-verified original GloBI v0.9 source-and-island overlay #43163702599](https://github.com/zuizui0223/island/actions/runs/43163702599), [anonymous artifact #13919492018](https://github.com/zuizui0223/island/actions/runs/43163702599/artifacts/13919492018), linked the **54 accepted native source-matched flower species**, original **197 native island×plant-species name matches on 13 islands**, exactly identical historical GSHHG polygon memberships, cited **exact species** GloBI original relationships and coordinate availability. After filtering for source-cited GloBI plant species records on the *same exact polygon where that plant is original-source native*:

| Original source-backed data unit | Strict observed count |
| --- | ---: |
| Georeferenced original GloBI distinct species×visitor group×source provider×citation×class×coordinate units among the 54 | 600 |
| Units geographically inside any of the exact 13 source-matched focal island polygons | **6** |
| Units also matching the ORIGINAL NATIVE PLANT SPECIES recorded on that same polygon | **6** |
| Unique matching native plant species | **3** |
| Distinct original islands containing those plant records | **3** |
| Same-island native plant species with 2+ documented GloBI guilds | **2** |
| Island native plant species with 3+ original distinct bibliography IDs | **0** |
| Island-polygon matches with explicit original GloBI `pollination` relation class | **3 rows** |
| Actual island population×mainland source population paired flower-visitation tests | **0** |
| Measured pollen transfer, seed set or verified evolutionary guild shifts | **0** |

**Five documented guilds at the same native plant on its exact original island:** flower-visiting birds **3 plant species/3 islands**; Lepidoptera (butterflies+moths together) **2 species/2 islands**; non-Bombus bees **1 species/1 island**; Bombus **0** and Diptera **0**. Records across visitor guilds can belong to the same plant and island; do not sum guild counts as independent species.

**Scientific significance without overclaim:** this is the first fd source-pinned result with **plant SPECIES NAME + originally native island flora ID + exact historical island polygon + original cited GloBI interaction type** simultaneously supported. It moves beyond assuming that a bird recorded anywhere on an island visits a particular flower. However, GloBI positions may reflect source-level georeferencing rather than standardized visit encounters; there is **no event timestamp in the original v0.9 TSV**, no original mainland counterpart visitor series, and no measured effective pollen transfer or fruit set. Thus a locally georeferenced `pollination` label must not be promoted to biologically verified effective pollination, and **none of the 3 species is a demonstrated evolutionary pollinator shift**. This 3-island sample is insufficient for a global macroecological regression.

[Anonymized 5-channel locality counts](native_54_exact_same_island_GloBI_v09_geo_5guild_counts.csv) · [verified original-source locality receipt](native_54_exact_same_island_GloBI_v09_locality_receipt.json).

## Source concentration and next macro-scale gate

Original exact-species GLOBAL data are strongly concentrated in one coarse archive-first label (~95% of 622 species×guild reference/class units); after removing that label the original 54-species globally cited coverage falls **27→7 species**, globally multiple-guild species **11→4**, and globally bird-documentation species **5→1**. Whether the **six island-local** source-cited records retain more than one genuinely independent original study, unique field site or archive group remains to be evaluated and must not be assumed from raw row count.

Because a study with only **3 species on 3 islands** cannot identify global island pollinator shifts, we additionally set up a **source-independent screening cohort of the full 928 Direct flora islands**, with native versus unresolved plant origin distinctly coded. The candidate macro screen matches the original Direct plant SPECIES name AND an original GSHHG polygon ID AND a GloBI original citation-bearing relationship; nonreproducible historical polygons are marked *not evaluable*, not pollinator absences. No result from that new full-928 screening is claimed before its source-verified run finishes.

The right manuscript architecture remains: **(1) robust floristic form/colour structure, (2) source-relative flora sorting, (3) an explicitly sparse but real source-cited plant-specific island-visitation bridge**, and (4) actual historical guild shifts/functional reproductive assurance as *untested hypotheses*, not conclusions.
