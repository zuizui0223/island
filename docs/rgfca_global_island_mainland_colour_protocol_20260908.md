# Prospective RGFCA × island global island–mainland flower-colour test

Frozen: 2026-09-08 (JST), before any island-versus-mainland colour contrast was computed.

## Question

Do photographs of the same flowering-plant species taken on islands show more white and less chromatic floral signal than photographs taken on continental mainlands?

This lane deliberately drops the RGFCA global-common-boundary estimand. The treatment contrast is island versus continental mainland, and the primary unit of replication is species.

## Frozen data sources

### RGFCA reserve photographs

Use the already completed FCP reserve measurement artifact only:

- repository: `zuizui0223/fcp`
- workflow run: `34091091640`
- artifact: `rgfca-reserve-measured-photos-v1`
- artifact id: `10025692225`
- artifact digest: `sha256:fbe4e8ab93e665115b5560cd05ddd31d02b9c9138b8a30ce460b61b6c6775ef5`
- measurement head: `1f80af7f5db81d61f28ae2818c130ef058a9f850`

No image reacquisition or colour remeasurement is authorized for this first test.

### Island polygons

Use the frozen GSHHG v2 island-universe artifact already used by the island project:

- GSHHG version: 2.3.7
- shoreline level: L1
- resolution: high (`h`)
- primary island area threshold: >= 5 km2
- continental exclusion threshold used to create the island universe: > 7,000,000 km2
- workflow run: `28659417688`
- artifact id: `8066083419`
- artifact: `gbif-three-block-pilot-28659417688`
- artifact digest: `sha256:bee33e14672ec7ff4ed1f7acaea36b32cc6170a45c3d14b2caeb5f43bfaa623b`
- exact island polygon count: 8,265

The primary island label is an exact point-in-polygon join against `islands_v2.gpkg`. A conservative island-size sensitivity retains only polygons >= 20 km2.

### Continental-mainland reconstruction

The island artifact intentionally excludes continental landmasses, so this first bounded analysis reconstructs the continental mask using the locally bundled Basemap 2.0.0 GSHHG-derived coastline data (GSHHG 2.3.6) and applies the same >= 7,000,000 km2 component rule. This yields five large land components in the available global mask. The proxy is used only to classify continental points; exact GSHHG 2.3.7 island polygons remain authoritative for island membership.

This is a declared source-version limitation. A publication-grade replication must reclassify continental membership from the exact GSHHG 2.3.7 L1 source and require agreement with this result; no effect threshold may be changed between versions.

Points falling in neither an eligible island polygon nor a continental component are `other_or_unclassified` and are excluded rather than forced into either class.

## Eligibility fixed before colour contrast

A photograph enters the continuous-colour frame only when:

1. `roi_status == automated_colour_state_admitted`;
2. flower 12-anchor palette mass is > 0;
3. matched background 12-anchor palette mass is > 0;
4. geographic class is exactly `island` or `mainland`;
5. positional accuracy is <= 1,000 m for the primary analysis.

The already-selected reserve observer rules are retained unchanged. No species is selected using its colour values.

Species enter the primary crossover frame only if they have >= 5 eligible island photographs and >= 5 eligible mainland photographs. Metadata-only pre-outcome census gives 149 such species under the <=1,000 m rule. The <=5,000 m geographic-accuracy sensitivity has 157 species.

## Frozen colour outcomes

Let the twelve palette states be:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze, green, brown, black`.

For each photograph calculate from saved pixel counts only:

- `flower_white = white / total12`;
- `background_white = background_white / background_total12`;
- `white_differential = flower_white - background_white`;
- `flower_chromatic = (yellow + orange + red + pink + magenta + purple + blue) / total12`;
- `background_chromatic` analogously;
- `chromatic_differential = flower_chromatic - background_chromatic`.

Biological direction predicted a priori:

- island `white_differential - mainland white_differential > 0` (whitening);
- island `chromatic_differential - mainland chromatic_differential < 0` (dulling / reduced chromatic signal).

Raw flower-only white and chromatic contrasts are descriptive secondary outputs. Matched flower-minus-background differentials are the primary flower-specific outcomes.

## Repeated species-equal estimator

Primary repetition count: `R = 1000`.

For each realization and each primary crossover species:

1. sample exactly 5 eligible island photographs without replacement;
2. sample exactly 5 eligible mainland photographs without replacement;
3. calculate the island-minus-mainland mean difference for each outcome;
4. aggregate species differences with equal species weight.

Thus every species contributes one value per realization regardless of how many photographs it has.

Random seed: `20260908`.

Report the 2.5%, median and 97.5% quantiles of the 1,000 realization estimates and the fraction of realizations in the predicted direction. These are resampling-stability summaries, not independent-sample confidence intervals for all angiosperms.

## Null test

For each eligible species, keep the 10 sampled photographs per realization fixed and randomly permute the island/mainland labels while preserving 5 versus 5 labels. For each realization use one label permutation per species and compute the equal-species global statistic. Use 1,000 null realizations with seed `20260909`.

For whitening use the upper tail; for chromatic dulling use the lower tail. Report `(1 + extreme_null) / 1001`.

No threshold or outcome is selected according to these p-values.

## Prespecified sensitivities

1. geographic accuracy <= 5,000 m, same >=5 photographs per side;
2. island area >= 20 km2, otherwise identical;
3. stricter sampling depth: species with >=10 photographs per side and 10 sampled per side, positional accuracy <=1,000 m.

All are sensitivities; none can replace the primary result post hoc.

## Claim boundary

A supported matched differential permits the bounded statement that, within crossover species represented in this opportunistic global photo sample, floral colour differs between island and mainland photographs beyond the matched local-background colour contrast.

It does not by itself establish evolutionary white-flower syndrome, adaptation, pollinator loss, island colonization history, or causal convergence. Species-turnover / assemblage-level island effects are a separate estimand and must be analysed separately.
