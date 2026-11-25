# All-928 original Direct floral island × original GloBI visitor-record coverage

**12 November 2026. Empirical preflight completed; these are source-cited locality matches, NOT measured effective pollination, historical visitor switching or pollination-syndrome selection.**

## Exact original data and reproducibility

[Original GloBI v0.9 source workflow #43198808095](https://github.com/zuizui0223/island/actions/runs/43198808095) **SUCCEEDED**. It verified frozen GloBI original v0.9 `interactions.tsv.gz` SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17` (2026-06-04, DOI `10.5281/zenodo.20546682`), the original Direct species floral-trait checkpoint, original published GBIF genus-family taxonomy and GSHHG 2.3.7 high-resolution islands source ZIP SHA256 `8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41`.

The originally locked **928** island floras with ≥30 directly measured, high/medium evidence, strictly single-state colour+form+symmetry species were sampled **before looking at GloBI visitor outcomes**. **927/928** island geometry IDs match the reconstructed original GSHHG WKB-based IDs; the unmatched original geometry is explicitly **UNASSESSABLE**, not a missing pollinator guild.

## New source-backed coverage: actual original plant species and exact historic island polygon

| Unit | Verified observed aggregate |
| --- | ---: |
| Strict Direct botanical species globally included in the 928-flora universe | **2,135** |
| Distinct strict original Direct island×plant species rows for 927 spatially evaluable islands | **81,593** |
| GloBI named species×guild original source-geo/citation candidate relation units across the 2,135 species | **23,693** |
| Original source-geocited relation units inside one of the 927 eligible polygons | 3,286 |
| Original relation units matching the SAME accepted plant species in original island flora and polygon | **3,182** |
| Island×plant species combinations with such original GloBI locality-relation evidence | **552** |
| Distinct botanical species with such locality-relation evidence | **249** |
| Distinct islands with such locality-relation evidence | **109** |
| Island×plant species pairs having ≥2 different GloBI documented visitor guilds | **213** |
| Island×plant species pairs with ≥3 distinct GloBI bibliography keys | **149** |
| Distinct bibliography keys across all matched original rows | 2,155 |
| Distinct original source-archive-first labels across those rows | 12 |
| Strict ORIGINAL NATIVE island×plant species pairs among 552 with documented original locality relations | **157 pairs across 40 islands** |
| Origin-UNRESOLVED island×plant species pairs with locality relations | **265** |
| GloBI original field `eventDateUnixEpoch` present | **NO** |
| Original mainland+island plant-specific effective visitor comparisons | **0 in input** |
| Measured pollen delivery or seed set | **0 in input** |

**Interpretation:** The first three large numbers are different units! 3,182 source-defined plant×visitor×reference×site relations are NOT 3,182 independent pollination field studies, and 552 island×plant documented source pairs are NOT 552 distinct plant species or independent sampling sites. The 213 multi-guild island×species associations are **globally source-cited, exact-name/geo-matched documentary co-occurrence**. Without common-visit sampling, independently recorded pollen deposition and study-source controls they are NOT confirmed ecological pollination generalism or compensatory guild replacement.

The source data coverage is HIGHLY selective: source-specific plant visitation is documented on only **109/927 spatially evaluable Direct floras (~11.8%)**, and strictly native plant pairs appear on **40/927 (~4.3%)**. Most of the 927 original island floras lack suitable cited flower-visitor evidence. Do not fill unsupported island×species cells with a zero or assume visitor absence. This is a conditional *evidence opportunity* diagnostic, not an unbiased all-island proportion of multiguild pollination.

### Original five guilds (all origin statuses retained with native explicitly distinguished)

| Documented original GloBI plant×visitor group | Islands with an exact same-flora-plant locality citation (any origin) | Species with such citations | Islands with a strict native plant in the original flora and this citation |
| --- | ---: | ---: | ---: |
| Bombus | 39 | 128 | 7 |
| Non-Bombus bees | 51 | 143 | 24 |
| Lepidoptera — moths+butterflies pooled | 62 | 127 | 21 |
| Flower-visiting birds | 16 | 18 | 11 |
| Diptera | 41 | 108 | 15 |

Records from the five groups overlap, and a source-cited `pollination` predicate is not a quantitative pollen transfer or seed set observation. Original GloBI v0.9 lacks event time, so historical visitor replacement cannot be reconstructed.

## Stronger original species-level GloBI universe, before geography restriction

Among **2,135** Direct measured plant species in the 928 original island floras, original GloBI global source documentation yields:
- **677** species with ≥1 exact plant-species GloBI cited flower interaction;
- **471** with ≥3 distinct bibliography keys;
- **475** with ≥2 documented visitor guilds (potential documentation breadth, not functional breadth);
- **314** with at least one original `pollination` label (not independent seed set);
- **588** species with any georeferenced relation;
- **44** archive-first source-provider labels, with the largest comprising **84.7026%** of 23,468 exact source species×guild×reference×evidence-class units. Removing that one label leaves **422** species with ≥1 documented original cited guild, **243** with ≥3 bibliography keys, and **295** with 2+ documented guilds.

This species-level original source backbone is substantially broader than the previous 54/13 candidate-only study and supports a meaningful new **SOURCE-DOCUMENTATION** prediction study if whole botanical families and source omission are respected, but cannot be used alone as direct evidence of local ecological pollinator shifts or an independently controlled study count.

## Predefined next independent-source / flower syndrome gates

The original source program [island Draft PR #323](https://github.com/zuizui0223/island/pull/323) now adds two strict paths:

1. **Original locality-source independence:** compare the SAME native and all-origin island×plant citation pairs, multi-guild and bird+other-guild counts before and after deleting the archive label selected as dominant across the full original 2,135-plant source set. Also count distinct provider labels and coordinate sites per island×plant pair. This directly tests whether the documentary multi-guild result remains without one source archive. This does not test an evolutionary mechanism.
2. **Exact botanical SPECIES source-label composition transfer** (not genus imputation): fit models to the global five-channel GloBI *documentation shares*, comparing reference-volume controls, open flower form, symmetry, colour, and pairwise flower-state interaction indicators. Hold out the ENTIRE botanical FAMILY, and keep an IDENTICAL botanical species set with ≥3 independent bibliography keys both before and after removing the globally dominant archive. Train-only scaling, ridge regularization, five-share simplex projection, equal-family bootstrap. Test whether floral morphology improves source-conditional out-of-family prediction, and whether an interaction term genuinely outperforms independent traits. If this fails, maintain only the Direct 928-island nonrandom morphological coordination result as the supported positive macro-ecological finding.

Both extensions have source-side code and synthetic tests, but the full verified source rerun **has not completed at the moment of this report**. No source-independent local multiguild or flora–pollinator-syndrome association estimate is invented.

## What is maintained in fd

The **separate 928-island Direct colour/architecture trait coordination** and its distinction from floral Rao diversity are still real reproducible flower assemblage patterns. This new 928-island source GloBI result adds that locally georeferenced **documented flower visitor relationships exist on 109 island floras, including native plant records on 40**. It does NOT demonstrate that morphological flower syndromes *cause* visits, that open/white flowers evolved after a historical bee→bird shift, or that 213 apparent multiguild local citation pairs represent effective pollination compensation.

No original plant species/island/geocoordinate/source publication/provider identities were uploaded. All released results are anonymous aggregate evidence gates with original source SHA and unit denominators preserved.
