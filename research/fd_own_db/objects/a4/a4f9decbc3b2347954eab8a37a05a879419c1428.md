# Range-map nectarivorous birds × island floral form — frozen contract (2026-11-25)

**Status: frozen before any BirdLife range-map data has been read. No result exists yet.**
Code: `src/fd/range_map_nectarivore.py` (incidence), `src/fd/nectarivore_flora_test.py` (test),
`scripts/range_map_nectarivore_incidence.py` (runner), `scripts/nectarivore_identifiability_sim.py` (design check).

## Why this test

The 24 Nov exploratory result (`results/20261124_tropical_pollinator_source_intake/SPECIALIST_PARTIAL_FLORA_VERIFIED.md`)
found that islands with more **specialist** nectarivorous birds (EltonTraits `Diet-Nect` ≥ 50) carry a larger share of
tubular flowers, while **partial** nectarivores (0 < `Diet-Nect` < 50) have a negative coefficient. The contrast was positive in
11/11 nested panels, but it rests on 18–32 islands in 4–6 archipelagos, disappears under genus-constrained nulls on both the
bird and the plant side, and gives 0/160 Holm-supported tests on the original GIFT cohort.

With six archipelagos the result cannot distinguish two explanations:

- **Lineage co-occurrence.** Specialist bird clades (hummingbirds, honeyeaters, sunbirds, Hawaiian honeycreepers) and plant
  genera that happen to be tubular arrived at, or persisted on, the same archipelagos. Within each plant genus the share of
  tubular species does not track the birds.
- **Trait matching.** Where specialist nectarivores are present, the plant species present are more often tubular *than
  expected from the plant genera present*.

The bottleneck is the bird checklist crosswalk (38 tropical islands). Expert range maps can supply specialist and partial
nectarivore counts for most of the 4,505 GIFT islands, which gives enough archipelagos to separate the two explanations.

## Data (to be obtained by the author; none is redistributed)

- **Bird ranges.** BirdLife International & Handbook of the Birds of the World, *Bird species distribution maps of the world*
  (BOTW), the most recent release available on request from datazone.birdlife.org. Record the release year, file SHA256
  and the taxonomy version in the run output. The file is licensed for research use and must not be committed or exported.
- **Islands.** The GSHHG island polygons already crosswalked to the 4,505 GIFT floras in the island repository.
- **Diet.** EltonTraits 1.0 bird table (Wilman et al. 2014, doi:10.1890/13-1917.1), `Diet-Nect` (% nectar, pollen, plant
  exudates). Sensitivity: AVONET `Trophic.Niche == "Nectarivore"` (Tobias et al. 2022) on the BirdLife taxonomy.
- **Flora.** The protected GIFT Direct floral records, run only inside the island repository's protected CI. Only anonymous
  aggregates leave it.

## Frozen incidence rules

1. Keep range polygons with `PRESENCE ∈ {1 extant, 2 probably extant}`, `ORIGIN ∈ {1 native, 2 reintroduced}` and
   `SEASONAL ∈ {1 resident, 2 breeding}`. Sensitivity A adds `SEASONAL = 3` (non-breeding) and `4` (passage).
2. Compute areas in an equal-area projection (EPSG:6933). A species is **present on an island** when its retained range
   covers **≥ 50 %** of the island's area (primary). Sensitivities: ≥ 10 %, any overlap (> 0), and ≥ 50 % after a 1 km island
   buffer. Islands below 1 km² are excluded from the primary analysis and reported separately, because range maps are not
   drawn at that resolution.
3. Range maps give potential occurrence, not observed presence. A zero means "no retained range covers the island", never a
   confirmed absence.
4. **Name crosswalk.** Exact match of normalised binomials between BOTW and EltonTraits, then an explicit, versioned synonym
   table (`from,to,source`) and nothing else; no fuzzy matching. The run **fails closed** if more than 5 % of the species that
   occur on at least one qualifying island are unmatched. Unmatched species are reported by count, never given a diet of 0.
5. **Diet classes.** specialist `Diet-Nect ≥ 50`; partial `0 < Diet-Nect < 50`; other `Diet-Nect = 0`. Sensitivity B uses
   20 as the specialist cut, sensitivity C the AVONET trophic niche.

## Population and support gates

- Primary region: tropics, |latitude| ≤ 23.44° of the island centroid. Secondary (reported, not part of the primary family):
  northern midlatitude, southern extratropical, using the existing fd macroregion definitions.
- Island eligibility: ≥ 5 Direct-scored species for the endpoint, all-recorded status, island area ≥ 1 km².
- Grouping unit: archipelago (the Weigelt archipelago label already used in the island repository); islands without an
  archipelago form their own 10° physical block. Sensitivity: 10° physical blocks for everything.
- A group is **informative** when it has ≥ 2 eligible islands and its specialist count varies.
- **Gate:** ≥ 30 informative groups and ≥ 200 eligible islands in the region. Below the gate the status is
  `insufficient_support` and no estimate is interpreted. Thresholds may not be changed after the first run.

## Estimands

For island *i* with *n_i* scored species, *y_i* tubular species, plant pool defined per region:

- `share_i = y_i / n_i` (observed).
- `expected_i` = Σ over the island's scored species of the mean tubular indicator of all region-pool species in the same
  **plant genus × evidence-grade block** (the analytic expectation of a within-block label exchange), divided by *n_i*.
- `residual_i = share_i − expected_i`.

Within-group weighted least squares, each group's islands demeaned and each group given equal total weight:

`outcome_i ~ log1p(specialist_i) + log1p(partial_i) + log1p(other_i) + log(area_i) + log2(1 + distance_to_continent_km_i) + log2(n_i) + |lat_i|`

- **E1 (association)** = β_specialist − β_partial on `share`.
- **E2 (beyond plant lineage)** = β_specialist − β_partial on `residual`. This is the primary estimand.

## Inference

- Group bootstrap: resample groups with replacement, 1,999 draws, percentile 95 % interval. Draws that lose full rank are
  dropped and counted.
- Bird-side null: permute diet labels among species **within bird genus** (primary) and within family, recompute every
  island's counts from the fixed island × species incidence, refit. 1,999 permutations, one-sided p for E > 0.
- Plant-side null is built into E2.
- The bird-genus null only uses diet variation *within* bird genera. Where every genus is uniformly specialist (e.g.
  hummingbirds), it has no power and returns p ≈ 1 by construction. Report the number of bird genera with mixed diet
  classes among species present on qualifying islands, and read the genus null together with the clade replication.
- An outcome with no within-group variance (e.g. a trait state absent from the region) is `not_estimable`, not a null result.

## Decision rule (primary family = one test)

Endpoint: tubular form (`floral_form == tubular`), tropics, all-recorded, archipelago grouping.

| Outcome | Reading |
|---|---|
| Gate fails | `insufficient_support`; nothing is claimed |
| E1 > 0 (interval excludes 0) and E2 interval includes 0 | lineage co-occurrence: bird and flower lineages are co-distributed; no matching beyond plant genus composition |
| E2 > 0 (interval excludes 0) **and** bird-genus permutation p < 0.05 | association beyond both plant and bird lineage composition (still observational; no visitation or pollen transfer) |
| E2 > 0 but bird-genus p ≥ 0.05 | beyond plant lineage, but not separable from which bird genera are present |
| E1 interval includes 0 | the 24 Nov exploratory pattern does not generalise to range-map incidence |

## Secondary analyses (Holm across the listed set, reported whatever the result)

1. Endpoints: deep tube class, extended access (tubular, salverform, funnel/trumpet, spurred), red/pink colour × extended access.
2. Other regions (northern midlatitude, southern extratropical).
3. **Clade replication.** Re-estimate E2 with the specialist count restricted to one nectarivore clade at a time
   (Trochilidae; Meliphagidae; Nectariniidae; Drepanidini; others pooled) and the remaining specialists moved into a control
   term. Trait matching predicts a positive E2 in each clade's realm; lineage co-occurrence predicts a positive E1 that is
   carried by single clades.
4. Sensitivities A–C, the overlap thresholds, the island area floor and the 10° grouping.

## Claim limits

Range-map co-occurrence is not visitation, pollen delivery, or a historical pollinator shift. Floral traits are species-wide
records, not island populations. A positive E2 is consistent with matching but also with unmeasured island-level factors that
favour both specialist birds and tubular species within plant genera (habitat, elevation, age). Native-status restriction is a
sensitivity, not the primary population, because 83.65 % of island × plant records lack resolved origin.

## Design check before data

`scripts/nectarivore_identifiability_sim.py` simulates both explanations. E1 is positive under both. E2 keeps a 2–5 %
false-positive rate under co-occurrence at every archipelago count and detects matching; the archipelago count sets its
power (64 % at 6 vs 98 % at 15 archipelagos for the weaker effect), not whether the two can be told apart. Its output is in
`results/20261125_range_map_nectarivore_design/README.md`.
