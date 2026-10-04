# Original seven-trait poster inference + WCVP

User-directed restoration on 4 October 2026. Broad all-analysis primary estimates, standard errors and p values are copied exactly from the original corrected seven-trait tables. Individual significance is nominal, unadjusted, two-sided normal P < .05; intervals are estimate +/- 1.96 SE. No aggregate score, t substitution or Holm adjustment is used. WCVP is fitted with the original seven-trait stacked cluster covariance and the same inference rules. Direct-only remains an evidence sensitivity. H2–H4 are unchanged.

The earlier `h1_traitwise_20261004` t/Holm replay is superseded as the primary analysis; its files are preserved for audit only. Original joint-Wald results have not been promoted into new evidence by this restoration.

## WCVP support

Before trait/covariate filtering: 513,320 island-species records, 2,372 islands and 90,693 distinct species. This includes existing source-native rows plus WCVP-compatible upgrades of unresolved rows. It is regional native compatibility, not verified nativity on each island.

| Region | Islands per trait (All) | Spatial blocks per trait |
|---|---:|---:|
| Northern mid-latitudes | 861–961 | 77–82 |
| Northern high latitudes | 127–261 | 41–58 |
| Tropics | 562–827 | 85–104 |
| Southern extratropics | 110–156 | 31–38 |

All seven traits in every region exceed the existing 50-island fit threshold, so the sensitivity is estimable. This is not a power guarantee. Some uncertainty is concentrated in a few blocks; e.g. the preceding influence diagnostic attributed about 77% of northern-high actinomorphy variance to one block. Species records are not independent replicates; regional coverage and intervals matter more than total records.

## Restoration verification

Broad estimates, SEs and p values match original source rows exactly (including the original Direct-only northern-high shallow/open-tube optimizer failure flag). That failed historical Direct-only fit must not be claimed as validated support; it is retained transparently rather than silently altered. The primary All fits and all WCVP fits converged. H2–H4 inputs are unchanged.

In primary All: original nominal P supports 17/28; the superseded finite-cluster replay supports 17/28 before Holm and 3/28 after Holm. These are different inference policies, not changing data or disappearing effects. Unadjusted individual P values do not control family-wise false positives.

## Replay

Use the pinned input downloads in `docs/REPLAY_H1_TRAITWISE_20261004.md`. Arrange the downloaded source artifact under INPUTS/sources/progressive-input, the status artifact under INPUTS/traitwise-inputs/status, and the WCVP artifact under INPUTS/traitwise-inputs/wcvp. Then run `python -m island_v2.chapter1_h1_poster_original --inputs INPUTS` from this repository with PYTHONPATH=src. The manifest records exact source hashes.
