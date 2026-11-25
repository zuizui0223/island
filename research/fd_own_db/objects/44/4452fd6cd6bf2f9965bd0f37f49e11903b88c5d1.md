# Fixed-target, fixed-placebo trait-ledger sensitivity

**fd PR #13 · 10 November 2026 · rights-protected original-data replay**. This is an exploratory *predictive* comparison, NOT evidence of mainland ancestry, adaptation, functional pollination equivalence, or phylogenetically independent replication. Source/island/species IDs and row-level predictions are not released.

## Test that previous analyses did not yet perform

The previous [ledger harmonization](LEDGER_HARMONIZATION.md) held the **island and evaluated species** constant but compared a nearby-source model against its own global-genus baseline. The previous [source-placebo test](SOURCE_CONDITIONED_FUNCTIONAL_TURNOVER.md) matched alternative source assignments but used only one of two archived trait ledgers.

**This replay uses both safeguards simultaneously:** the evaluated island × accepted species × observed binary state is identical, the actual five-nearest native-GIFT source entities are identical, and the *same two outcome-blind alternative source sets* are used when fitting **both** original trait ledgers.

The two measured binary traits are explicit recorded SC versus SI and actinomorphy versus zygomorphy, not flower-visitor effectiveness. Do not equate these strict original-source eligibility counts with the upstream Chapter 1 seven-atomic-trait contract.

## Audit of cohort and controls

- Original rights-protected `island` source archive SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`; separate trait-source ZIP SHA256 `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`; original public fd G0/G1 ZIP SHA256 `2efa7d2bb40adc60cc7d22b3c0c1043eb63c7c608e492fb0a422a14b3a563801`.
- Original fd cohort: **1,419 islands**, **7,618 overlapping pair rows**. The conditional source validation uses **1,068 geographically matched focal islands** across **105 distinct spatial blocks**. Local matching was selected for 85 northern-high-latitude, 876 northern-midlatitude, 71 tropical, and 36 southern-extratropical islands. The southern group fails the fixed >=40-island and >=10-block release gate, so no scored effect is reported.
- Source selection universe: 296 distinct five-nearest-GIFT-source combinations in the original 1,419-island cohort, from the original GIFT-native 748-region, 1,666,186 region × species archive. No reconstructive claim of historical source ancestry.
- Source controls: two *different* alternate five-region sets from other islands in the same broad latitude regime but a different spatial block. All sources use the same original island taxa and known genera. Select **before looking at observed 0/1 trait states**, matching within **±2 percentage points** the focal same-genus donor availability for **both traits in BOTH All archived training ledgers**, native GIFT species-pool size ratio **0.5–2**, and source-set Jaccard **below 0.8**; at most 128 candidates sampled per focal; two fixed seeds for this post-hoc design.
- Among matched groups, the reported average maximal known-genus donor coverage difference was <1 percentage point per target. Two source-assignment seeds must not count as independent replications: both produced the **same focal islands**, 85 northern-high and 71 tropical island predictions were identical in both seeds, and about **87%** of northern-midlatitude scores were identical.
- A species name has a single SHA256 first-byte modulo-five fold, globally excluded from all native-GIFT source genus donors when evaluating that fold. For each trait and All/Direct scope, the target is the **intersection of species known in both trait ledgers with the SAME binary state**. Every island/trait keeps identical evaluated names and statuses between ledger scorings.
- Predictive priors: global genus shrunk by **0.05** pseudo-donor names to the global heldout trait state; local five-region genus shrunk by **8** to the fold-heldout global genus predictor. These are exploratory choices rather than hyperparameters selected in nested training.
- Metric: per-name proper Brier error, averaged once **per focal island**, then once **per spatial block**. Report `mean(placebo Brier) - real-source Brier` (positive favours actual source), plus **4,999 whole-block paired bootstrap** diagnostic 95% intervals. The paired source-versus-snapshot ledger delta uses *the same block bootstrap indices*. This is conditional uncertainty, not p-values from independent evolutionary histories.

## Main measured results

The table uses the first control-selection seed, giving exactly the same focal islands for both trait ledgers. The second seed produces near-identical directions.

| Strict trait, region | Scope | Focal islands / blocks | Snapshot-source gain | Source-adjusted-ledger gain | Paired new-minus-old gain (conditional 95% block interval) |
| --- | --- | --- | ---: | ---: | ---: |
| SC/SI, northern high latitude | All | 85 / 12 | **+0.000697** | **+0.000537** | −0.000160 [−0.000663,+0.000243] |
| SC/SI, northern high latitude | Direct | 85 / 12 | +0.000933 | +0.000799 | −0.000134 [−0.000707,+0.000327] |
| SC/SI, northern midlatitude | All | 876 / 59 | +0.000135 | **−0.000036** | **−0.000170 [−0.000374,−0.000002]** |
| SC/SI, northern midlatitude | Direct | 876 / 59 | +0.000160 | −0.000018 | −0.000178 [−0.000405,+0.000019] |
| SC/SI, tropics | All | 71 / 33 | −0.000115 | −0.000151 | −0.000036 [−0.000101,+0.000029] |
| SC/SI, tropics | Direct | 71 / 33 | −0.000143 | −0.000179 | −0.000036 [−0.000122,+0.000042] |
| Symmetry, northern high latitude | All | 85 / 12 | −0.000051 | −0.000004 | +0.000047 [−0.000056,+0.000207] |
| Symmetry, northern midlatitude | All | 876 / 59 | −0.000031 | −0.000028 | +0.000003 [−0.000027,+0.000034] |
| Symmetry, tropics | All | 71 / 33 | +0.000016 | +0.0000004 | −0.000016 [−0.000035,+0.000002] |

The corresponding Direct symmetry residuals are also small and not consistently source-positive. Full 24 aggregate rows and the shared-focal two-seed sensitivity are in the locally retained **identity-suppressed** result CSVs, not source taxon ledgers.

### Source-placebo intervals and robustness

- Northern-high SC All: snapshot gain **+0.000697** [0.000200,0.001331], source-adjusted gain **+0.000537** [0.000296,0.000795]. Direct is also positive under both ledgers. Dropping any one of the **12** blocks retains a positive mean. These 12 blocks are not 12 demonstrated independent source histories; both control seeds yield identical island-level scores.
- Northern-mid SC All: snapshot +0.000135 [−0.000094,+0.000374], source −0.000036 [−0.000268,+0.000200]. **Neither ledger individually establishes a source advantage.** The *paired* training-ledger difference is negative under the conditional 59-block bootstrap in All, but Direct's 95% interval includes zero and the test is one of many post-hoc sensitivities. Across the second control seed, the All paired delta is −0.000191 [−0.000386,−0.000020].
- Tropical SC: both training-ledger source-placebo estimates are negative but intervals cross zero.
- Florally actinomorphic/zygomorphic symmetry: **no region** has a reproducible positive source-placebo gain with a 95% block interval excluding zero in both ledgers. The strong global-genus trait predictability from prior fd studies remains a different, surviving endpoint.
- The southern-extratropical region has only 36 eligible matched focal islands across 9 region blocks; it is **insufficiently supported**, not biologically confirmed zero.

## Scientific conclusion and limits

The same-source placebo vs two-training-ledger test **does not rescue a geographically universal mainland-source functional-trait signal**. Known genus remains strongly informative, especially for recorded floral symmetry. For strict SC/SI, a **small northern-high-latitude source-predictive advantage** survives the two archived ledgers and matched alternatives, but with only twelve spatial blocks and no external botanical/phylogenetic replication. The northern-midlatitude signal is strongly training-ledger-sensitive.

Observed flora status and raw trait evidence completeness, source-genus donor opportunity, repeated accepted species across islands, and historic dispersal origin all remain confounded. An out-of-name prediction plus whole-block bootstrap does **not** demonstrate within-lineage flower evolution, reproductive fitness, true native colonization or pollinator-mediated selection.

**Next identification gate:** new phylogenetically/geographically disjoint source histories and source-independent trait ascertainment, or direct trait/fitness observations, rather than more island-pair pseudo-replications. A preregistered independent geographic holdout is necessary to claim ecological generality.

## Reproducibility/rights

Reusable input validators and equal-block paired result scorer: `src/fd/functional_matched_ledger_audit.py` and `tests/test_functional_matched_ledger_audit.py`; the public CI uses *synthetic data* and public GBIF frequency aggregates only. Full private ZIP replay was performed in the private session, but protected raw species, genus, island, source and per-name model predictions are never uploaded to GitHub or a public workflow. Public outputs are 3 aggregate CSVs plus a SHA/provenance manifest; this report contains only broad-region group information.
