# 14 November 2026 — Original GSHHG centroid-near photo source audit for tropical island endemics

**Original protected source CI completed: [island Actions #43581171912](https://github.com/zuizui0223/island/actions/runs/43581171912) SUCCESS**, anonymous aggregate result artifact #14108556768. First attempt #43580858311 computed its original response but failed in the optional printed receipt due to a key typo; the corrected run alone is accepted as source-verified. No botanical taxon/island/photo IDs or precise observation coordinates are exported.

## Source and inference guard

- Independent Roeble 2024/Zizka 2022 geological labels: 434 geologically continental/oceanic island records in the original GSHHG source. LGM non-connection is NOT considered proof of oceanic geologic origin.
- Frozen original GIFT native + separately verified endemic status, high/medium Direct flower-colour category, five closest GIFT native candidate mainland source pools with independent accepted same-genus Direct colour evidence.
- **13 shared tropical endemic genera** represented by 26 different source-verified island endemic accepted species, selected BEFORE using photo API response, 13 oceanic and 13 continental. Sample spans nine oceanic and five continental original 10-degree source geographic blocks, and 11 verified source archipelago labels.
- Current iNaturalist public observations API v1, species-name **exact** unambiguous active species query; photos must be attached to a research-grade taxon-identical record, georeferenced with PUBLIC unobscured latitude/longitude and original positional accuracy ≤5 km, with supported CC0 / attribution / share-alike or selected NC open reuse terms. No photos are treated as flowers based solely on observation name.
- Candidate locality filter is a circle around the verified original source GSHHG island **CENTROID**, at a radius clipped to 5–75 km from the equal-area island circular radius. This circle can cross sea, neighboring islands, and mainland. It is explicitly NOT actual island polygon point-in-polygon validation or an exact island-local population occurrence. Obscured/private locations are never accessed or inferred.

## Actual original public provider observations

| Original botanical source condition | Oceanic endemic | Continental endemic |
| --- | ---: | ---: |
| Pre-selected true accepted endemic plants, one per genus | 13 | 13 |
| Unambiguously exact named public iNat species | 11 | 11 |
| Global-public licensed geo/grade photo candidate on queried original first page | 4 species | 3 species |
| **Candidate observed near original source island centroid in independent radius query** | **4 species** | **1 species** |
| **Qualified publicly open geophoto observations near centroid (first page only)** | **18 records** | **1 record** |
| Provider errors | 0 | 0 |

**ZERO** genus pairs have a near-centroid-qualified public photo for BOTH the geologically oceanic endemic species and the geologically continental endemic species. Note that the preceding GLOBAL source-only stage had **one** botanical genus with licensed photo records for both endemic species; that does NOT survive the more restrictive original island CENTROID proximity test.

There is no analysis of floral ROI, photogenic pigmentation, effective visitors, within-plant/within-species variation, ancestry or actual establishment. Photo availability is geographically and taxonomically biased; **4 vs 1 is not evidence of a real geological difference in floral colour, local presence or reproductive adaptation**. The source n=13 is too small and nonrandom for a geological effect.

## What this changes in the research strategy

The earlier idea to require a geographically local photograph of two different island endemic species (one in each geological class) from the SAME genus before even comparing phenotype is demonstrably very restrictive. Do **not** interpret failure to find an island-photographed same-genus pair as absence of suitable biological comparisons.

A stronger two-process design is:

1. **Assembly/colonization proxy**: source-region- and native-flora-adjusted floral trait frequency composition, genus/family turnover and categorical functional diversity across oceanic/continental islands, with independent archipelago and spatial-block weighting. Original severe support limit: only three fully independently geography-disjoint geological-class matched island pairs per tropical or northern-midlatitude broad region, so no promoted geological-causal effect.
2. **Post-arrival phenotype divergence**: sample multiple original source-verified endemic island populations together with candidate mainland close relatives/conspecific populations. Independently verify original island polygon point membership, continuous/palette flower ROI with genuine photo bytes and consistent methods, photo metadata/license/observer, native status, and phylogenetic proximity. Compare phenotype divergence island-minus-mainland across geological origins using phylogenetic and geographic clustering. A GIFT potential source is NOT necessarily the evolutionary ancestor.
3. **Functional mechanism**: independently observed actual floral visitor guilds, flower pollen deposition, autonomous seed set/pollen addition and plant fitness, which are not contained in current public photo metadata. A photographic source is NOT evidence of syndrome-driven selection or pollinator replacement.

A practical next data gate is target-first **four oceanic-endemic taxa with presently eligible source-island CENTROID-near public photography**, verify original polygon and blooming flower pixels and then find independent mainland same-genus source photo cohorts. Expand across genuine independently named archipelagos, not merely more photos per already photographed plant.

## Reproducibility

[fd source module](../../src/fd/oceanic_endemic_near_island_photo.py), [original prior global photo selection](../../src/fd/oceanic_endemic_inat_photo_pilot.py), [synthetic geoprivacy/radius tests](../../tests/test_oceanic_endemic_near_island_photo.py), [island original rights-safe runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_oc_endemic_near_island_inat_photo.py), [original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-oc-endemic-near-island-photo.yml).

The user-facing numbers are only source-verified *availability and observation effort*, not a biological effect. Keep the source and manuscript pull requests **DRAFT**.
