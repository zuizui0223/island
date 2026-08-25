# Chapter 1 when / where frozen result — 2026-08-25

## Canonical observation-robust run

This checkpoint freezes the Chapter 1 result after explicitly separating the ecological trait-composition surface from the observation process.

- workflow run: `32845980788`
- artifact: `chapter1-when-where-main-32845980788`
- artifact id: `9562465057`
- artifact digest: `sha256:d4cd8ad0eea9cd2bda879af14f54c549e130de365faf1b56bebddad38b24f1c2`
- direct-trait input run: `32702160934`
- floristic-status input run: `32559322028`
- geographic/context input run: `29228212586`
- M3 genus-fixed null draws: `1000`

No pollinator variable enters the fitted Chapter 1 design.

## Data-generating view

GBIF occurrence data are **not treated as a census of island floras and are not assumed to be a simple random sample**. They are treated as an opportunistic, incompletely observed sample of realised island floras.

Raw record multiplicity is collapsed to unique `island × species` incidence before trait composition is calculated. An island with no flora record is retained in the observation-process layer but is never coded as trait absence.

The ecological estimand is therefore:

```text
P(trait state | directly trait-resolved observed realised flora, geography)
```

This is an assemblage/filtering probability conditional on the observed flora. It is not:

- the probability that a trait evolved on an island;
- the probability that an unrecorded species is absent;
- a source-pool-standardised colonisation probability.

The mainland-distance variable is interpreted as a **composite geographic gradient** that can represent both dispersal limitation and changing accessibility to mainland/source species pools. It is not labelled a pure causal isolation treatment.

## Primary question

> **When and where do floral/reproductive trait probabilities change along the mainland-distance/source-pool-accessibility gradient, and where do multivariate response vectors differ?**

## Full-island observation layer

All **8,265 islands** in the frozen geographic universe remain represented in the observation-process analysis.

- no realised-flora record: **3,760 islands**
- at least one realised-flora record: **4,505 islands**
- at least one direct Chapter 1 trait domain represented: **425 islands**
- direct evidence for all three core domains somewhere in the island flora: **405 islands**

Observation itself is geographically structured. For example, conditional on having flora records and controlling for area, climate, observed flora support and context, the direct-trait-evidence endpoint has a northern-midlatitude reference distance-gradient odds ratio of approximately **4.32 per 1 SD** (`p ≈ 1.2e-4`). This is why observation-process robustness is mandatory rather than assumed away.

## Trait-centric probability layer

For each directly resolved atomic trait state, models estimate the probability of encountering that state in the trait-resolved observed flora as the geographic/source-pool gradient changes.

The example motivating this design, self-compatibility, does **not** show a clear positive distance-gradient signal in the current direct-evidence sample:

- northern mid-latitude, all native: OR ≈ **1.018 per distance SD**, `p ≈ 0.305`
- tropical, all native: OR ≈ **1.049**, `p ≈ 0.120`

Thus the Chapter 1 result must not be reduced to “farther islands have more SC species.” Different trait components carry different parts of the regional response vector.

## Primary WHERE result

The response-vector omnibus tests ask, within each biogeographic context:

```text
H0: all supported atomic trait-state geographic slopes = 0
```

At the confirmatory threshold (each retained response represented on >=50 islands), the canonical species-count-weighted result supports geographic trait filtering in both:

- **northern mid-latitude island floras**;
- **tropical island floras**.

The result is present in both `all_native` and `native_nonendemic` strata.

## BETWEEN-WHERE result

For responses meeting the same support threshold in both contexts, the formal test is:

```text
H0: northern-midlatitude response vector = tropical response vector
```

The northern-midlatitude and tropical vectors differ confirmatorily in both status strata.

Canonical species-count weighting:

- all native: q = **2.35e-08**
- native non-endemic: q = **7.13e-07**

A significant slope in one region and a nonsignificant slope in another is not used as evidence of a regional difference; the vector contrast itself is tested.

## Observation-bias robustness

### Information-weight sensitivity

The canonical grouped-binomial model gives more information weight to islands with more trait-resolved species. Therefore the same trait shares were refit under:

1. canonical species-count weighting;
2. cap at 100 effective trials;
3. cap at 50;
4. cap at 20;
5. equal-island weighting.

Across 5 weighting modes × 2 status strata, the full headline was reproduced **10/10 times**.

Even under equal-island weighting:

- all-native northern WHERE: q = **5.25e-130**
- all-native tropical WHERE: q = **1.77e-32**
- all-native north-vs-tropical: q = **1.75e-11**
- native-nonendemic north-vs-tropical: q = **7.83e-14**

Thus the headline is not dependent on a few trait-rich islands receiving large species-count weight.

### Trait-resolution coverage adjustment

Direct-trait coverage varies strongly among regions and domains. For example, in all-native observed floras, median direct evidence fractions are approximately:

| trait domain | northern mid-latitude | tropical | southern extratropical | northern high-latitude |
| --- | ---: | ---: | ---: | ---: |
| flower colour | 0.769 | 0.614 | 0.727 | 0.815 |
| floral form | 0.254 | 0.168 | 0.141 | 0.097 |
| self-incompatibility | 0.357 | 0.182 | 0.171 | 0.473 |

A response-specific smoothed-logit trait-resolution fraction was therefore added to every atomic response model.

After this adjustment:

### All native

- northern WHERE: joint χ² = **418.26**, df = 21, q = **3.06e-75**
- tropical WHERE: joint χ² = **177.24**, df = 17, q = **1.02e-28**
- north vs tropical: joint χ² = **71.73**, df = 17, q = **2.17e-08**

### Native non-endemic

- northern WHERE: joint χ² = **447.11**, df = 21, q = **3.13e-81**
- tropical WHERE: joint χ² = **236.21**, df = 17, q = **1.35e-40**
- north vs tropical: joint χ² = **69.08**, df = 17, q = **3.11e-08**

The observation-bias checkpoint therefore records:

```text
observation_robust_headline = true
```

This does not prove missing species are missing at random. It shows that the headline cannot readily be explained by the measured differences in trait-resolution coverage or by the information dominance of trait-rich islands.

## Geographic-functional-form robustness

All distance sensitivities retain the full island universe. No 50-km, oceanic-island, positive-distance, or other geographic threshold is used.

The headline is reproduced for:

- `log1p(distance)`;
- `sqrt(distance)`;
- raw distance.

Across 3 distance forms × 2 status strata, the headline is reproduced **6/6 times**.

## Spatial-block influence

The northern WHERE, tropical WHERE and north-vs-tropical comparison were recomputed after deleting one spatial block at a time.

- all native: **84/84** block deletions preserve the full headline
- native non-endemic: **84/84** preserve the full headline

No single spatial block is necessary for the result.

## WHEN result

The confirmatory status-persistence classification is:

| context | all native | native non-endemic | endemic | current classification |
| --- | --- | --- | --- | --- |
| northern mid-latitude | supported | supported | not confirmatorily testable | **persists in native non-endemic** |
| tropical | supported | supported | not confirmatorily testable | **persists in native non-endemic** |
| northern high-latitude | not confirmatorily testable | not confirmatorily testable | not confirmatorily testable | **unresolved** |
| southern extratropical | not confirmatorily testable | not confirmatorily testable | not confirmatorily testable | **pilot signal / confirmatory unresolved** |

Therefore the northern and tropical responses are not patterns confined to endemic taxa.

## Domain decomposition after the WHEN/WHERE freeze

The multivariate result is not one classic island syndrome.

- **floral form:** strong northern WHERE; confirmatory northern-vs-tropical difference
- **flower colour:** WHERE signal in both northern and tropical contexts; weaker/inconsistent evidence for between-context difference after domain-wide correction
- **self-incompatibility:** much weaker contribution; SC itself has no clear positive distance slope in the current direct sample

Atomic/domain results explain what makes the already-supported regional vectors different; they do not replace the primary vector-level WHEN/WHERE test.

## Current Chapter 1 conclusion

> **Within opportunistically observed island floras, floral and reproductive trait probabilities change systematically along a mainland-distance/source-pool-accessibility gradient in both northern mid-latitude and tropical regions. The multivariate response vectors differ between those regions, persist in native non-endemic assemblages, and remain supported after equal-island weighting, explicit trait-resolution coverage adjustment, alternative full-universe distance transformations, and deletion of every single spatial block in turn. Current data do not support equivalent confirmatory conclusions for northern high-latitude or southern-extratropical floras.**

This is a WHEN/WHERE conclusion. It does not identify the causal ecological mechanism and it does not distinguish dispersal limitation from changing source-pool accessibility.

## What additional trait acquisition is now for

Additional trait acquisition must **not** be used to rescue the already-supported northern/tropical result.

Its primary purpose is now to expand the WHERE map into currently underpowered regions, using an outcome-blind acquisition priority based on evidence gaps rather than observed effect size.

Highest-value gaps are currently:

1. **southern extratropical floral form and SI/SC** — pilot regional signal exists but confirmatory response support is too thin;
2. **northern high-latitude floral form** — very low current direct coverage and too few trait-resolved islands;
3. **northern high-latitude SI/SC** — island count, rather than only within-island fraction, is the limiting factor;
4. colour is lower priority because direct coverage is already comparatively high.

The correct acquisition target is therefore `region × trait-domain testability`, not “species that would strengthen the expected effect.”

## Stronger source-pool extension

Where a defensible candidate source pool can be declared independently of island trait outcomes, a stronger secondary analysis can ask:

```text
P(species occurs on island | trait, candidate source pool, geography)
```

That is distinct from the present conditional trait-composition estimand and should be attempted only where the source-pool denominator is defensible.

## Chapter 1 -> Chapter 2 handoff

Chapter 1 ends with:

> **Why does the same mainland-distance/source-pool-accessibility gradient produce detectable floral/reproductive filtering in both northern mid-latitude and tropical island floras, yet generate different multivariate response vectors?**

Mechanistic discrimination belongs to Chapter 2 (`izu-core`).
