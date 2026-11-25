# Partner E2 source eligibility — 25 November 2026

**Real source audit completed; NO new ecological coefficients.** Protected [Actions 45529579281](https://github.com/zuizui0223/island/actions/runs/45529579281) SUCCESS, anonymous artifact 14984725127, provider-reported ZIP SHA256 `0d16ef3cf068f92c7ca5491f715046200e9459b00591ee90f5bd87ea90be4262`. Exact logged aggregate JSON reproduced in [source_eligibility.json](source_eligibility.json); archive digest is from Actions metadata, not a second downloaded ZIP verification. Source original 4,505 botanical / 8,264 geography frame, frozen SHA-checked trait and flora inputs, exact 234 bee geometries, exact original Marshall archive MD5 and model/species CSV SHA256. Source script [island commit 91b369a](https://github.com/zuizui0223/island/blob/91b369a79506ab21e5f98589aa8f530424448005/scripts/fd_partner_source_eligibility.py).

## Fruit positive control and growth-form diagnostic
Current original species-Direct high/medium ledger has ten scored flower/mating trait endpoints. No fruit morphology, dispersal, growth or life form. Actual public GIFT 3.2 trait metadata successfully retrieved (109 rows; SHA256 ee20bafd237be38b4cb97d65219f65d1f02ecc5fca09867c9cdfc1640fb30af4):
| Trait | ID | Published metadata count |
|---|---|---:|
| Fruit dryness (dry, fleshy) | 3.18.1 | 2578 |
| Fruit type | 3.16.1 | 6567 |
| Dispersal syndrome (coarse) | 3.3.1 | 9602 |
| Dispersal syndrome (detailed) | 3.3.2 | 6518 |
| Growth form (coarse) | 1.2.1 | 222338 |
| Life form (coarse) | 2.3.1 | 100390 |

Counts are metadata, NOT Direct raw, unrestricted, accepted matched species or eligible islands. Fix **3.18.1** as the direct dry/fleshy outcome for first acquisition, **1.2.1** growth form as environmental/confounding diagnostic only. Need versioned raw records, non-derived/unrestricted reference screening, conflict/grade audit, accepted-name matching and protected island support census. Do not fill absent species as dry. If direct dryness source fails, future expanded unambiguous fruit-type mapping needs its own preregistered acquisition/reconciliation, no effect-driven substitution. The existing 38 provisional tropical bird/bee candidate islands cannot meet the 200/30 primary gate; no BirdLife polygons received. Positive control not yet executable on primary frame. Pollen or zoochory labels are not direct fruit morphology.

## Bombus source support
Marshall archive 305 model islands and 11,424 occurrence rows; 68 distinct corrected name strings with certain native Bombus records (not an independently taxonomically accepted count), 64 model islands with positive certain native Bombus records, 25 with uncertain native Bombus entries. 304 model richness lists reconcile. Cov 303 nonmissing, all in [0,1]. 87 model islands pass the frozen source checklist quality gate before geometry/trait/region join. 234 unique mutual-overlap island geometry joins: tropical 133, northern midlatitude 92, south extratropical 6, northern high latitude 3.

| Cohort | Islands | Groups | Positive | Nonrecord | Informative groups (presence varies) |
|---|---:|---:|---:|---:|---:|
| Northern midlatitude geometry + checklist | 92 | 28 | 47 | 45 | 7 |
| Northern source-qualified | 26 | 12 | 18 | 8 | 4 |
| Northern source-qualified + ≥5 true same-species symmetry/depth plants | **25** | **12** | **18** | **7** | **3** |
| All regions, source-qualified | 69 | 36 | 20 | 49 | 4 |

Nonrecord is not confirmed ecological absence, especially in first row. Source-qualified positives/nonrecords still not true local presence/absence validation. Even maximal northern source overlap is <200 islands; stricter source/trait frame fails the original 30-group gate by a large margin. No E1, E2, trait-by-exposure summary or p fit, no threshold relaxation, no Baker's law test. Need independent adequate range/checklist data before a primary Bombus test. IUCN coverage is taxon/assessment-specific, and absence of an assessment/map must not become zero. Existing 305-island archive is valuable for validating mapped positive records where geometry and origin agree.

## Common code
PR #66's original nectar-bird contract unchanged. New generic E1/E2, binary fruit-label normalization, plant-only null for genus-defined Bombus, two-sided conditional null and fixed Holm for diet-defined groups, clade-specific model with other focal diet taxa as separate control. [CI 45529725186](https://github.com/zuizui0223/fd/actions/runs/45529725186) SUCCESS, 28 tests (synthetic only). Fruit and Bombus source audit above is real data eligibility, not model evidence. IUCN source-adapter config and source rights still needed before animal range data execution.
