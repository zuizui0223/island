# Functional traits: matched-taxon ledger stability of mainland-source predictive gains

**2026-11-10 | exploratory real-data source/definition audit, fd PR #13.** All taxon, island, mainland GIFT source and spatial-block identities remain private. This report contains only aggregated quantities. No adaptive or independently replicated geological/pollinator selection effect is inferred.

## Why another test was needed

The prior functional-turnover tests used two related but **nonidentical** trait ledgers from the SAME original protected `island` progressive ZIP:
- `snapshot/chapter1_trait_ledger_{all_analysis,direct_only}.csv.gz` (the revised snapshot in the upstream study);
- `fixed/source/current/evidence/{all_analysis_eligible,direct_only}_trait_ledger.csv.gz` (the source-adjusted historical frozen layer used in the preceding donor analysis).

Earlier flower-symmetry genus predictions used the snapshot, but the recent source-placebo exercise used the source-adjusted ledger. The proper diagnostic is **not** to compare source gains on different evaluated species. We directly assessed strict-state consistency and reran **the same species × island outcomes** on both training ledgers, holding the five-fold name partition and nearest mainland source choices fixed.

## Strict binary evidence agreement

Only explicit single-state `SC` versus `SI` or `actinomorphic` versus `zygomorphic` is admissible. Multi-state descriptions are missing for these strict comparisons, never treated as zeros.

| Evidence scope | Binary trait | Snapshot known names | Source-adjusted known names | Shared known names | Conflicting recorded states |
| --- | --- | ---: | ---: | ---: | ---: |
| All | SC / SI | 21,885 | 17,885 | 17,009 | **71 (0.42%)** |
| Direct | SC / SI | 6,671 | 6,726 | 6,124 | **64 (1.05%)** |
| All | Actinomorphic / zygomorphic | 50,621 | 56,202 | 50,189 | **34 (0.068%)** |
| Direct | Actinomorphic / zygomorphic | 32,039 | 32,431 | 31,894 | **25 (0.078%)** |

The states are highly concordant *where jointly observed*, but eligibility differs. This matters especially for All SC, where an additional 4,876 strict-state species occur only in the snapshot ledger and 876 only in the source-adjusted one. The source-ledger SC denominator cannot be called the 21,885-name snapshot denominator.

## Same-name heldout genus prediction still succeeds

Global names (each species assigned once via SHA256 first byte modulo five) were held out from their genus predictors. When scoring ONLY the shared names with **identical recorded states**, genus Brier errors at exploratory tau=0.05 were:

| Trait / scope | Snapshot-train genus Brier | Source-train genus Brier | Identical target species |
| --- | ---: | ---: | ---: |
| SC / SI All | 0.05436 | 0.05531 | 16,938 |
| SC / SI Direct | 0.15073 | 0.15199 | 6,060 |
| Floral symmetry All | 0.01051 | 0.01047 | 50,155 |
| Floral symmetry Direct | 0.01662 | 0.01640 | 31,869 |

These refer to the **whole island-trait species catalog**, not independent native-source lineage histories. The high predictive information supplied by known genus is much more stable than any further geographic source gain.

## Direct GIFT source-versus-global-genus sensitivity on identical targets

**Population and controls.** Use the same 1,419 fd G0/G1 sampled islands, 7,618 overlapping original island-pair edges, and the same five geographic GIFT-native source entities (`geo_k5`) for each island. 748 native GIFT source regions contain 1,666,186 native-list × species records. The status flora contributes 879,010 distinct island × species records to these 1,419 islands. Focal species-name folds are globally excluded from both global and source-specific GIFT donor counts. GIFT source-list duplicates of the same species within the five-region union are collapsed.

For the **target**, retain ONLY island × species records that are strict-known **and state-identical in BOTH ledgers** (no favorable sample selection by source prediction). This yields 339,918 All SC, 456,612 All symmetry, 289,070 Direct SC, and 405,416 Direct symmetry *island × taxon occurrences*. Each case uses the same physical islands, target names, observed states, source entities and broad region regardless of the donor-training ledger.

For reproducibility of the smoothing contrast, global genus probability is shrunk by **0.05 pseudo names** to the fold-excluded global native GIFT frequency, and five-nearest-source genus probability by **8 pseudo names** toward that global-genus prediction. Both parameters are exploratory sensitivity choices; the arbitrary tau=8 global-genus prior inflated source-versus-global gains dramatically and is **not** used below.

The target-name Brier gain is `(y-p_global_genus)^2-(y-p_local_source)^2`. First average each island's recorded evaluated names, then give every spatial block equal weight (16 north-high, 66 north-mid, 22 southern, 65 tropical blocks), with 4,999 whole-block bootstrap resamples. The following values are **source-versus-global-genus**, NOT the more stringent donor-coverage-matched alternative source placebo result.

| Trait, All | Region | Snapshot model gain | Source-adjusted model gain | Paired source-minus-snapshot gain, 95% conditional block interval |
| --- | --- | ---: | ---: | ---: |
| SC/SI | North high latitude | +0.002564 | +0.001317 | **−0.001248 [−0.002285,−0.000512]** |
| SC/SI | North midlatitude | +0.000728 | +0.000266 | **−0.000462 [−0.000886,−0.000040]** |
| SC/SI | Southern extratropical | +0.000439 | +0.000718 | +0.000278 [−0.000063,+0.000628] |
| SC/SI | Tropical | +0.000278 | +0.000029 | **−0.000248 [−0.000406,−0.000085]** |
| Symmetry | North high latitude | +0.000325 | +0.000653 | **+0.000328 [+0.000141,+0.000483]** |
| Symmetry | North midlatitude | +0.000443 | +0.000298 | **−0.000145 [−0.000278,−0.000029]** |
| Symmetry | Southern extratropical | +0.000172 | +0.000159 | −0.000013 [−0.000061,+0.000038] |
| Symmetry | Tropical | +0.000037 | +0.000034 | −0.000004 [−0.000027,+0.000021] |

Direct evaluates the **same island/target evidence constraint independently** (but it overlaps All, so is NOT an independent external replication). Northern high SC changes by **−0.001742** [−0.002957,−0.000835], northern mid SC by **−0.000887** [−0.001442,−0.000345]; northern high symmetry by **+0.000382** [+0.000161,+0.000592]. Other Direct results and exact denominators are in the local, anonymous 16-row aggregate.

**Interpretation:** Even after the evaluated species and their observed states are identical, changing the *donor-training trait catalog* can materially change the **very small extra mainland-source predictive increment**. A narrow conditional block-bootstrap interval is not a preregistered significance claim; multiple regions, scopes and traits were examined post-hoc, and donor histories/training species overlap across blocks. Do not call these predictive score differences the ecological effect of geology or mainland ancestry.

## What these results do NOT test

- They do not repeat the donor-coverage- and pool-size-matched **real source vs alternative source assignment** of [SOURCE_CONDITIONED_FUNCTIONAL_TURNOVER.md](SOURCE_CONDITIONED_FUNCTIONAL_TURNOVER.md). The earlier finding of no robust symmetry source-identity advantage therefore remains unaltered.
- They do not compare independent taxon/clade evolution or model spatial/phylogenetic source dependence beyond block-resampled diagnostic uncertainty.
- They do not separate source-provider reporting, true species historical origin, native versus introduced status, field performance or selection.
- Exact trait-state matching on common names removes target outcome composition shifts but does **not** remove changes in training donors, which are the sensitivity being measured.

## Current robust conclusion

**Genus-associated predictability of recorded floral symmetry is strong and stable across both original source ledgers. Mainland source-geography increments are much smaller, depend on trait evidence availability and donor training, and do not yet establish an independent biogeographic assembly mechanism.** SC may retain geographically conditional source information, but the previous matched-placebo and strict-native caveats remain.

## Reproducibility and rights

Inputs are original licensed/restricted ZIPs `island` GitHub Actions artifacts `10058653212` and `10042377593`, plus public fd G0/G1 artifact `13741845313` (original SHA256s recorded in the parent PR). The pipeline used a private local source-replay script and preserved ONLY anonymized summaries:
- `fd_trait_definition_harmonization_20261110.csv` (four aggregate rows)
- `fd_functional_ledger_oof_20261110.csv` (48 grouped summary rows including tau sensitivity)
- `fd_functional_harmonized_source_common_v2_20261110.csv` (32 region/scoped source-versus-genus summary rows)
- `fd_trait_ledger_paired_source_delta_20261110.csv` (16 paired conditional block audits)

Public reusable checks are in `src/fd/functional_ledger_harmonization.py` with synthetic tests; they validate strict states, paired name/target identity, crossfold-name consistency and identity-free group aggregation. **The GitHub public CI can test logic but cannot independently replay licensed records.**
