# Target species/genus repetition sensitivity for GIFT floral prediction (20 November 2026)

The complete GIFT botanical source study retains 4,505 island floras and 13 high/medium genuinely species-Direct floral/reproductive endpoints. Source-regional Brier and log-loss improvements have been resampled by geographic 10-degree blocks and first candidate mainland-source groups. **Neither geography grouping makes a species-level flower trait on the same accepted botanical name independent across separate islands.** One species can be recorded on many island checklists, but the same published taxon-level trait state is reused.

## Protected-source original follow-up

The original [island PR #342](https://github.com/zuizui0223/island/pull/342) introduces two additional, explicitly different weighting estimands on identical previously held-out original Direct target plant carriers:

1. **Unique accepted species-equal:** for each distinct accepted plant name, first average its island-specific source prediction proper Brier/log-loss gains across ALL GIFT target islands on which it appears. Bootstrap that accepted species once as one outcome unit, not once per island (minimum 30 distinct scored species).
2. **Unique genus-equal:** average the already species-equal values within each true botanical genus, then bootstrap genera as units (minimum 10 distinct source-supported target genera). This helps diagnose congeneric dependence, but it is not a phylogenetic independent-contrasts model.

Both estimands evaluate geographic predicted-genus flower-category models relative to the same original THREE proper categorical source references: target-region genus-BLIND, other-region genus-BLIND and simple pooled regional other-island GENUS. Both losses, both source target units, 499 exploratory bootstrap replicates. Values of each accepted plant name/genera stay in memory of the protected botanical source; only anonymous scores, CI bounds, counts and mean accepted-species repeat ratios can be exported.

Per original focal island, the mean individual botanical accepted-species gain is checked to reproduce the identical original island-level proper-score difference to 1e-10, to ensure the original prediction has not changed. Source original whole botanical accepted-name embargo for the focal entire geographic block/first GIFT candidate source entity, pinned original GIFT Direct traits/statuses, 8264 GSHHG geo parent, 4505 GIFT parent, all 13 endpoints, 3 geographic regions and full/native nested strata remain frozen. **No 13-trait complete-case filter is applied.**

## Conservative inferential implication

If tropical floral symmetry or Southern extratropical primary flower colour are still positive against the three source references on BOTH Brier and log loss after accepted species/genus reweighting, their source-category prediction would be more robust to trait-identity pseudoreplication. If they do not survive, the earlier island-equal predictive result remains a real result for ISLAND checklists but cannot be interpreted as repeated independent taxon-level floral response. These are different estimands; no single best weighting is chosen post hoc.

No result of original PR #342 is claimed before its protected Actions completes and the original anonymous artifact is checked. All intervals are exploratory and multiple trait categories are dependent; genus units are not fully phylogenetically independent.

**Nonclaims:** not independent mainland ancestor measurements, within-species island floral adaptation, observed pollinator guild switching/fitness, actual propagule establishment, morphological 17D KDE or a corrected phylogenetic causal selection inference.

The independent fd auditor is scripts/verify_gift_target_taxon_source.py, validated by tests/test_verify_gift_target_taxon_source.py. It checks full original 78-panel family, privacy, accepted-name repetition counts, both taxon source units × all proper-loss signed gates, and exact preservation of the earlier original per-island model score and cohort counts.
