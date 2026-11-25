# Native plant–island species-local GloBI evidence: first verified source-georeferenced bridge

**12 November 2026, fd Draft PR #14 and island Draft PR #323.**

## What actually passed

Original [GloBI v0.9 SHA-verified workflow #43163702599](https://github.com/zuizui0223/island/actions/runs/43163702599) **succeeded**. Its [anonymous original-data receipt #13919492018](https://github.com/zuizui0223/island/actions/runs/43163702599/artifacts/13919492018) is consistent with the source side's [exact original GSHHG candidate polygon recovery #43162949729](https://github.com/zuizui0223/island/actions/runs/43162949729) and the original frozen 21 Bombus-qualified nondetection/bird-presence candidate islands.

Input source identity gates:
- GloBI June 4, 2026, v0.9 original `interactions.tsv.gz`, SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`.
- GSHHG 2.3.7 high-resolution original shoreline ZIP, SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`.
- Exact **21 historical candidate island polygon ID** matching; GSHHG global full 8,265 IDs do **not** fully replicate (8,260 rebuilt, 9 original missing, 4 new), and the reconstructed polygon subset must NOT be advertised as all-island identical.
- Candidate GIFT geo_k5 mainland *source-list name* overlap: 54 distinct strictly island-native source-matched accepted plant species in 49 genera / 197 island×species source opportunities / 13 focal islands, derived before inspecting GloBI interaction outcomes.

For these 54 strict native botanical species, the original source has **600** unique geo-cited plant×visitor guild×archive source×original reference×relation class×coordinate units after citation/provider qualifications. **6** are geographically inside one of the 13 exactly matched island polygons and match the **same native accepted plant species** in the original island flora. The 6 entries concern **3 accepted native species on 3 islands**; **2 native island×plant pairs have ≥2 documented visitor guilds**. There are **3 relations labeled `pollination`** and **none of the three original island×species pairs has 3+ distinct supporting bibliographic keys**.

| Documented original GloBI visitor channel | Distinct same-island native flower species | Original candidate islands |
| --- | ---: | ---: |
| Bombus | 0 | 0 |
| Other bees | 1 | 1 |
| Lepidoptera (moths and butterflies together) | 2 | 2 |
| **Flower-visiting birds** | **3** | **3** |
| Diptera | 0 | 0 |

These original GloBI relationships are **not** automatically independent observations of pollination on that island; GloBI location metadata may be inherited from contributing source records, the original v0.9 event-time column is absent, and `pollination` as a relation type does not imply measured pollen deposition/seed production. An island animal occurrence by itself never proves interaction with a named flower.

**Biological interpretation: the positive source-backed result is narrowly that there exist original source-cited georeferenced associations of named native plants and visitor guilds on verified island polygons**. A mainland GIFT list is not an ancestral plant population; no mainland population flower visitor or pollination-effectiveness comparison has been paired. Thus **zero certified evolutionary source-to-island pollinator shifts in the input** is evidence-unavailable, not absence of biological switching.

## Original provider/site replication is the next credibility check

The previous genus documentation analysis found a single coarse `sourceArchiveURI`-first group accounts for 95.96% of plant genus×guild reference units. For these 54 EXACT plant species, the dominant group accounts for **95.016%** of 622 citation/provider/class units; deleting it reduces global source-cited species from **27 to 7** and globally bird-documentation species from **5 to 1**. This is a **source-archive label sensitivity**, not clean independent-paper replication.

We have extended source scripts in island PR #323 to report, **among the actual 6 same-island source-cited relations only**: original distinct bibliography keys, source archive labels, unique coordinate sites, whether two guilds came from repeated or distinct local coordinates, and how many **remain when the dominant archive label was selected from all 54 original candidate species before local filtering**. Those numbers are **not known from the earlier output**; the renewed SHA-pinned GloBI workflow is currently executing. No favorable independence or study count is imputed.

## Macro-scale extension, not 21 hand-picked anecdotes

The full original 928-island Direct-floral cohort can be screened for evidence availability without using the 21 Bombus-undetected/bird-occurred islands as the entire target population. A new source-side module `src/island_v2/fd_globi_macro_928_locality_v09.py` and tests define the screening population **before the original GloBI outcome scan** and require:
- exactly 928 islands with ≥30 complete original species-direct high/medium flower colour/form/symmetry botanical taxa;
- exact accepted plant species + exact historical GSHHG polygon ID + valid original GloBI original provider/ref and site coordinates;
- **island native records separate from unresolved origin status**, with no shift or native status imputed from geographic presence;
- original ID-mismatched polygons left `unassessable`, not forced positive/negative;
- no ecological claim of independent plant-specific effective pollination or historical pollinator switching from plain citations.

The original full-928 preflight is triggered by source workflow [#43198808095](https://github.com/zuizui0223/island/actions/runs/43198808095). Its code/tests, original source/cohort checks, and original island geometry retrieval passed. **The original GloBI TSV re-scan has not produced a validated macro numeric receipt yet** at this report version. Do not back-fill fabricated 928-island visitor proportions.

## Decision

Retain **floral trait Rao diversity / nonrandom colour×architecture coordination and floristic mainland source-pool filtering** as the fd positive quantitative evidence. Treat a real ecological shift as a distinct question needing independent same-species source/island flower visitor observations, historical comparability and ideally pollen-effective visit/seed set data. Present the 6 exact-polygon interactions as a tiny source-backed potential follow-up, not a macroscopic bird-for-bee adaptation study.

Public original identity-protected grouped files: [5-guild 3-island counts](native_54_exact_same_island_GloBI_v09_geo_5guild_counts.csv); [original source/locality receipt](native_54_exact_same_island_GloBI_v09_locality_receipt.json). Original botanical species/island/latitude-longitude/source provider/reference IDs were never committed to fd. Both PRs remain Draft.


## Confirmed original local citation independence audit

An independent later SHA-pinned original GloBI source replay **SUCCEEDED**: [island Actions #43197981268](https://github.com/zuizui0223/island/actions/runs/43197981268). The previously reported 6 matched native plant×island original-GloBI geo-citation units (3 plant species / 3 islands) can be decomposed by site and original evidence identifiers, without disclosing any original identities:

| Original local source evidence | Verified number |
| --- | ---: |
| Distinct original bibliography keys across all six rows | **4** |
| Distinct original georeferenced site coordinates | **4** |
| Distinct `sourceArchiveURI`-first labels | **2** |
| Native island×plant pairs with two or more original archive labels | **0** |
| Native island×plant pairs with two or more distinct coordinate sites | **1** |
| Native island×plant pairs with ≥3 original bibliography keys | **0** |
| Original relation units left after deleting the archive label largest in the *whole 54-species* source cohort | **3/6** |
| Plant species / islands with locality-linked evidence after that omission | **1 species / 1 island** |

The post-deletion single surviving plant species retains documented **other-bee, Lepidoptera and bird guild** references on the same island polygon, but they are **not evidence of simultaneous visits, pollen-transfer effectiveness or independent historic local pollinator switching**. The dominant archive was determined from the **whole source-verified 54 native source-matched species corpus**, not selected to make those six locality records favorable or unfavorable.

Most importantly, **none of the 3 native species’ island pairs is independently corroborated by ≥2 archive-first source labels**. Across the source-wide original GloBI v0.9 archive, provider independence remains limited even where same-island GloBI coordinates are available. The six rows are 4 citation keys and 4 source coordinate sites, NOT six independent plant population experiments.

[Evidence-resolution waterfall](native54_evidence_waterfall_source_study_and_locality_gates.csv) is a nonidentifying aggregation with denominators for plant species, islands, citation units and island×plant observations kept separate.

## Full-island macro next step is genuinely registered

The new `src/island_v2/fd_globi_macro_928_locality_v09.py` and synthetic `tests/test_fd_globi_macro_928_locality_v09.py` in [island PR #323](https://github.com/zuizui0223/island/pull/323) extend beyond the conditionally selected 21 candidate islands to the **whole original 928 Direct-trait-complete island flora**. Only the original same accepted plant species, source-cited full GloBI relationship, and a matching original WKB-hashed exact GSHHG island polygon are accepted. Strict `native` is kept separate from unresolved flora origin, and the 9 unreplicated original island polygon IDs can never be silently counted as guild absence. The new [source-pinned workflow #43198808095](https://github.com/zuizui0223/island/actions/runs/43198808095) has **passed all synthetic checks and original source recovery/geometry verification**, then entered the 2.7GB original GloBI v0.9 scan. Macro ecological outcome counts are **pending** until the original run succeeds; no numbers are guessed from the source-only genus proxy.

### Inference ceiling

The original fd **928-island colour×flower-architecture coordination** is a source-reconstructed **flower trait result** independent of these source-geocoded animal rows. The 3-island local links can serve as a descriptive data-quality bridge but **do not establish an island bee→bird pollinator shift, compensation, new pollination syndrome, selective adaptation or causal reproductive service**. A 928-island evidence-coverage preflight may say whether a geographically broad documented plant–visitor association analysis is even supportable; only if enough independent island species, matched sampling opportunity, bibliographies and providers are present should a real macro pollinator-functional model be considered.
