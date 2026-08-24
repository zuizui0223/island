# Flora-status real-data checkpoint — 2026-08-22

## Scope

This checkpoint separates what is already supported by the first real GIFT
materialization from what is still unresolved.

Frozen inputs used for the diagnostic:

- GIFT status run `32555972362`, artifact `gift-floristic-status-32555972362`;
- latest reviewed direct-trait run on PR #132: `32551742699`;
- locked island distance/regime artifact: `purpose-shortest-distance-regime-29228212586`.

## GIFT materialization

GIFT v3.2 exported:

- 1,934 island-list metadata rows;
- 689,073 checklist rows;
- 1,412 exact-Island geometries with shapes.

Conservative GIFT-to-GSHHG geometry matching produced:

- 443 accepted GIFT entity matches;
- 426 frozen islands with at least one status-ledger row;
- 304,914 source-level ledger rows;
- 55,137 species represented in those source-level rows.

The first ledger used GIFT list-level endemism as if it were focal-island
endemism. That interpretation is **invalidated**. GIFT explicitly warns that a
species reported as endemic in one source/list may occur in other checklists.
Examples in the real artifact include widespread taxa carrying contradictory
`endemic` / `nonendemic` source-list states. Therefore GIFT endemism is now only
a source claim and cannot enter primary endemic analyses by itself.

## Origin status remains usable independently

Collapsing only native / introduced status, while ignoring endemism conflicts,
produces 198,559 unique `island x species` origin rows:

- native: 150,531;
- introduced: 19,442;
- unresolved: 28,586;
- native/introduced conflicts: 3,670 rows, retained unresolved.

Endemism disagreement never erases an otherwise consistent native origin in the
new resolver.

## Direct trait support among confirmed-native rows

Primary model support is now `>=1` direct trait trial on an island. The actual
number of direct species remains the grouped-binomial trial count, and stricter
`>=3`, `>=5`, `>=10` subsets are retained as robustness diagnostics. `30` and
`50` refer to numbers of usable islands, not numbers of species per island.

### Northern mid-latitude

| outcome | >=1 direct | >=3 | >=5 | >=10 |
|---|---:|---:|---:|---:|
| flower colour | 241 | 237 | 233 | 229 |
| floral form | 235 | 230 | 225 | 218 |
| SI/SC | 239 | 233 | 228 | 220 |

### Tropical

| outcome | >=1 direct | >=3 | >=5 | >=10 |
|---|---:|---:|---:|---:|
| flower colour | 138 | 137 | 134 | 124 |
| floral form | 134 | 124 | 113 | 98 |
| SI/SC | 135 | 129 | 116 | 99 |

### Southern extratropical

| outcome | >=1 direct | >=3 | >=5 | >=10 |
|---|---:|---:|---:|---:|
| flower colour | 33 | 31 | 28 | 25 |
| floral form | 31 | 25 | 21 | 17 |
| SI/SC | 32 | 28 | 25 | 17 |

Interpretation:

- northern and tropical confirmed-native composition is already well above the
  predeclared 50-island confirmatory-count target for all three focal outcomes;
- this remains true under the stricter >=10-direct-species robustness subset;
- southern extratropical support is limited mainly by the number of islands with
  source-backed floristic status, not by a global shortage of trait records.

The isolation range among supported native islands is broad rather than confined
to near-continent islands. Under the >=1 direct-trial definition, the median
continent distance is about 183 km in northern mid-latitudes and about 883–894 km
in the tropical subset, with long upper tails into remote oceanic islands.

## Endemism correction now in progress

The repository now uses a separate WCVP corroboration layer:

1. GIFT provides the source-reported endemic claim;
2. WCVP/POWO provides the accepted species' global native TDWG level-3 range;
3. the frozen island is independently mapped to a TDWG level-3 region;
4. `endemic` is accepted only when the GIFT claim is corroborated by exactly one
   WCVP native TDWG-L3 region matching the focal island region;
5. WCVP native range in more than one TDWG-L3 region is sufficient to classify
   the taxon as non-endemic at this regional scale;
6. all other cases remain unresolved.

This is deliberately conservative. The resulting `endemic` category is a
TDWG-L3-corroborated regional endemism class; it must not be described as
single-island endemism unless the TDWG region itself corresponds to that island.

## Current decision

Do **not** acquire more flower/reproductive traits for the northern or tropical
all-native analyses. Those models are already data-supported.

The next acquisition decision waits for the corrected endemism artifact. If an
endemic-stratum outcome has fewer than 50 usable islands but at least 50 islands
with corroborated endemic flora, target trait acquisition to unsupported endemic
islands. If the total number of islands with corroborated endemic status is below
50, the bottleneck is floristic-status/geographic resolution rather than trait
acquisition.
