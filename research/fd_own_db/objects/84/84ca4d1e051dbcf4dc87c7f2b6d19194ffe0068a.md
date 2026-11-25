# Equal absolute-distance floral source reference: N1d 1,000/2,000 km

**Verified empirical execution, 2026-11-09.** [Public-source full CI #42587303351](https://github.com/zuizui0223/fd/actions/runs/42587303351) completed successfully with the original pinned 496 floral colour/form state-set profile dissimilarity matrix, independently sourced island centroids and aggregate island-profile incidence. The workflow artifact fd-origin-strata-analysis (ID 13660371567) contains the full grouped CSVs, support table, and input SHA256 manifest. [Protocol and exact restrictions](../../docs/FIXED_RADIUS_DISTANCE_SES.md) · [source code](../../src/fd/fixed_radius_ses.py).

## Why a fixed geographic radius matters

The previous N1c reference draws from ten closest OTHER eligible islands, but the median tenth-neighbor distance differs by origin: approximately 387 km for continental and 1,149 km for oceanic island floras (All). This analysis instead fixes **absolute geographic reach** to the same 1,000 km or 2,000 km for both origins, taking *all other eligible islands in the focal island's same broad regime* within that radius. Donor profile weights are **island incidence counts**, never taxon carrier counts. At least three other recorded donor islands and at least n distinct reference profiles are required to sample n distinct types without replacement. N1b (all other eligible regional floras) is recalculated and compared **on the same focal island subset** for each radius.

The reference choices are deliberately **post-hoc**, informed by N1c's geography imbalance; neither radius was chosen before inspecting the previous result. All values below are **descriptive**. No new geological p-value or claim of adaptive floral contraction was generated.

## Source-pinned All-scope results: equal 30 distinct recorded floral profiles

Endpoint: within-island mean distance to the nearest OTHER recorded floral profile, on original equal-weight colour/form set-Jaccard D (range 0–1). Observed 30-profile score is averaged over 99 uniform thinnings of the same island; each null has 999 weighted draws. SES standardizes the difference by the null's own sample SD.

| Matched reference population / endpoint | Continental | Oceanic |
|---|---:|---:|
| **1,000 km** n focal islands | **127** | **115** |
| 1,000 km N1b regional SES median, same islands | +0.183425 | -0.526999 |
| **1,000 km N1d SES median** | **-0.036852** | **-0.226675** |
| 1,000 km median within-island N1d-minus-N1b SES change | -0.100419 | **+0.294707** |
| 1,000 km fraction negative N1d SES | 51.97% | **73.04%** |
| 1,000 km median observed-minus-N1d expected original D | -0.000864 | **-0.005024** |
| **2,000 km** n focal islands | **129** | **124** |
| 2,000 km N1b regional SES median, same islands | +0.183425 | -0.548493 |
| **2,000 km N1d SES median** | **+0.001062** | **-0.370289** |
| 2,000 km median within-island N1d-minus-N1b SES change | -0.079859 | **+0.227078** |
| 2,000 km fraction negative N1d SES | 49.61% | **75.81%** |
| 2,000 km median observed-minus-N1d expected original D | +0.000025 | **-0.008790** |

**Independent of reference reach selection:** The same-island oceanic SES becomes less negative under both fixed radii than its own broad-regional N1b comparator. Yet the negative median remains, including in unstandardized original D: **-0.0050 at 1,000 km and -0.0088 at 2,000 km**, respectively. These magnitudes are small relative to D's full 0–1 scale and are NOT causal biological effect sizes or evidence for independent geology.

**Direct nested evidence:** Equal30 nearest-D N1d median SES is -0.0374 continental (118 islands) and -0.1893 oceanic (102) for 1,000 km; -0.0794 continental (121) and -0.3667 oceanic (113) for 2,000 km. Direct is an overlapping source-evidence sensitivity, **not another experiment**.

An alternative within-island mean distinct-*pairwise* D gives All oceanic median N1d SES -0.159 at 1,000 km and -0.248 at 2,000 km, likewise negative but weaker; continental medians +0.025 and +0.085. These metrics share the same profiles.

## Support imbalance is NOT resolved by equal geographic reach

| Scope/origin | Radius | Eligible geography-typed islands | Equal30 reference possible | Median OTHER donor islands | Median donor profile types |
|---|---:|---:|---:|---:|---:|
| All continental | 1,000 km | 130 | 127 | **41** | 249 |
| All oceanic | 1,000 km | 130 | 115 | **9** | 197 |
| All continental | 2,000 km | 130 | 129 | **123** | 298 |
| All oceanic | 2,000 km | 130 | 124 | **24.5** | 226.5 |
| Direct continental | 1,000 km | 123 | 118 | 38 | 239.5 |
| Direct oceanic | 1,000 km | 120 | 102 | 8 | 192 |
| Direct continental | 2,000 km | 123 | 121 | 115 | 292 |
| Direct oceanic | 2,000 km | 120 | 113 | 22 | 211.5 |

Even within the SAME radius, oceanic target islands have only about **one fifth** as many eligible neighboring source-island observations as continental targets. Consequently the N1d comparison is **not** an equal-opportunity or source-lineage-conditioned test: spatial sampling density and island geography remain strongly imbalanced.

By broad region in All equal30, N1d nearest-D oceanic median SES is -0.175 (47 northern-midlatitude islands) and -0.248 (66 tropical islands) at 1,000 km; -0.317 (48) and -0.351 (73) at 2,000 km. The southern-extratropical oceanic cells have fewer than five scoreable islands and are suppressed. Treat these as spatially dependent descriptive regions, not external replications.

## Ecological decision gate

1. **Falsification of a strong claim:** Universal oceanic flower-packing/selection is NOT demonstrated by KDE, by raw signed 37-pair tests, by geography-disjoint 10-pair categorical tests, or by the current distance-SES spatial nulls.
2. **What survives:** A relatively small negative nearest-profile *recorded-category* departure for several oceanic island assemblages when conditioned on geography and donor-island profile frequency. It depends on reference radius and highly unequal nearby-island recording support.
3. **Next necessary identification:** Genus/trait donor conditional N2 with independently heldout evaluated taxon names and mainland-native source opportunities, predictive log-loss/calibration and expected-richness audits. Only its public aggregated summary may be released; current licensed input/checkpoints are insufficient to run it within this public CI. No new N1 variant should substitute for this lineage test.

## Reproduce without restricted private inputs

    python -m fd.fixed_radius_ses \
      --input work/island_public \
      --g0 work/fd_g0 \
      --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
      --output results/fixed_radius_ses \
      --draws 999 --thinnings 99

Only group-level effect summaries with small-cell suppression and anonymized donor support counts are in the published artifact. No island IDs, source-archipelago identifiers, taxon mappings, original restricted GIFT native-region records, or raw exact coordinates are included.
