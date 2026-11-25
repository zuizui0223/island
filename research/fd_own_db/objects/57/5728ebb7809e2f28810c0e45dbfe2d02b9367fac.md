# N1c: geographically nearest island profile frequency reference

**2026-11-09.** Direct companion to the validated N1b regional
profile-frequency null. This diagnostic asks whether **spatially nearby
recorded floras** already explain apparently unusually close island
flower-colour/form categories. It does not require or reconstruct any
restricted species, genus or GIFT mainland donor data.

## Prespecified source and statistic

- Input: the frozen 496-profile original colour/form set-Jaccard D,
  public island×profile support aggregates (All/Direct) and
  independently audited island-centroid geographic crosswalk.
- Focal population: the same geological origin-labelled eligible
  ≥50 complete-floral name islands as the public 17D and N1b analyses,
  with nested All/Direct scope flags.
- For each focal island, choose **the 10 geographically closest OTHER
  eligible islands within the SAME broad analysis_regime**, regardless
  of geological class, using haversine distances between audited
  centroids. Ties use sorted island keys for deterministic selection.
  Exclude the focal island completely.
- Build the union of all distinct profiles recorded by these ten
  donor islands. Give each profile a sampling weight equal to the
  **number of the ten OTHER islands recording it**, not the number
  of individual taxon-name carriers.
- Sample n distinct categorical profiles WITHOUT replacement per
  draw, where n is the observed focal distinct profile count.
  The local source pool may lack n profiles; mark such a comparison
  **unestimable**, never sample with replacement.
- Sensitivity: if the focal island has ≥30 recorded profile types,
  average 99 independently seeded uniform without-replacement
  30-profile observed thinnings; compare with n=30 neighbor draws.
  The 99 thinnings are dependent views of one observed flora.
- Each null uses 999 (default) reference draws; compute
  mean distinct off-diagonal D and mean within-island nearest
  OTHER profile D. Record SES and raw distance difference.
  If reference standard deviation is zero, SES is undefined;
  do not add a pseudocount.
- Recompute N1b alongside N1c on *the same island focal/support
  subset* for a direct **descriptive** comparison. N1b pools all
  other eligible islands in the same regime; N1c restricts to the
  fixed ten nearest.
- Publish scope/origin/regime aggregates with all effect entries
  suppressed when cell n<5; report the fraction of target
  islands with 10 eligible geographic donors and sufficient
  support plus median tenth-neighbor radius (only cells ≥5).
  **No individual island IDs, coordinates or neighbor membership.**
- As N1c will sometimes be unestimable, a separate
  `local_geography_matched_reference_comparison.csv` measures N1b and N1c
  on **the same physical islands**, separately by scope, origin, broad
  regime and endpoint, plus pooled regime summaries. It compares the
  median SES, fraction of negative SES and **raw original-D observed-minus-
  null differences**, not just unrelated group median values. Small cells
  (n<5) suppress effect entries. These remain descriptive contrasts,
  not an inference from geographically independent species/archipelagos.

## Why this is not N2

N1c is another **observational profile-incidence reference**.
It does not condition on observed genus composition, independent
mainland lineage opportunity, island shared species, observation
probability, plant ancestry or the true source flora. Indeed,
nearest neighbors may share numerous exact taxa with the focal
island. This is a deliberately conservative *geography* sanity
check, not pollination filtering or source-mediated trait evolution.

If N1b yielded a strongly negative SES on oceanic islands, but
N1c moves the same measured distance close to zero, local
spatial/available-category composition is a plausible sufficient
description. If N1c remains negative, this is **still not
independent source/genus filtering**: proceed to the private
taxon-heldout N2 proposed in [N2_SOURCE_DISTANCE.md](N2_SOURCE_DISTANCE.md).

No new geological-origin p values are calculated. The existing
32 N0/N1 tests and the four planned N2 local-genus tests remain
separate, explicitly post-hoc multiplicity families; avoid adding
N1c by searching for an attractive threshold or nearest-neighbor
count after seeing outcomes.

## Runnable public-source replay

```bash
python -m fd.local_geography_ses \
  --input work/island_public \
  --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/local_geography_ses \
  --neighbors 10 --draws 999 --thinnings 99
```

Only grouped CSVs (local_geography_ses_summary.csv,
local_geography_support.csv, local_geography_matched_reference_comparison.csv) and a schema/hash/limits manifest
are released in the public CI artifact. No required restricted
input can or should be downloaded as part of this step.

**Scientific interpretation ceiling:** local predictability of
*recorded floral phenotype category neighborhoods* conditional
on empirical neighboring recorded assemblages. Never claim
selection, adaptation, ecological niche volume or true
phylogenetic source independently.
