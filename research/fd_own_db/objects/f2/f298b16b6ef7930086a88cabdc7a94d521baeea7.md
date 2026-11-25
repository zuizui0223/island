# Functional mainland-source **independence audit** — overlapping donor regions

**2026-11-10 | fd PR #13 | real-source, exploratory sensitivity.** This audit addresses whether the previously observed, small northern-high-latitude SC/SI predictive advantage can be treated as replicated across *independent mainland source opportunities*. It cannot: the largest donor-overlap network collapses its candidate source sets to one connected component. This is a conservative label-overlap graph, **not** proof that all those islands share one evolutionary history.

## Frozen source and analysis unit

The same rights-authorized `island` progressive ZIP (SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`), source-scale trait ZIP, and fd G0/G1 pinned cohort were used as in [MATCHED_LEDGER_PLACEBO.md](MATCHED_LEDGER_PLACEBO.md).

The target is a recorded binary trait, strict **SC versus SI** or **actinomorphy versus zygomorphy**, scored on identical held-out species × island outcomes under both archived training ledgers (`snapshot` and `source`). The first predeclared source-control seed 20261110 is used; the second highly overlapping control seed is NOT an additional independent sample. Reported gains are proper Brier loss of two coverage-/pool-size-matched alternative native-mainland GIFT donor assignments minus the actual geographic five-region assignment (positive favours actual source).

Each island has exactly one five-region `geo_k5` GIFT donor *signature* from the original source assignments. The full source ZIP's native checklists and private island/source signature IDs are never released. Codes are validated against five distinct rank-1..5 source regions per included island.

## Independent-unit diagnostic — the number of island pairs is misleading

Two signatures are connected if they share at least one of the five assigned mainland-region entities. Transitive closure defines a **connected source-donor group**. This is intentionally strict and different from a phylogenetic model or an experimentally verified ancestral relationship. A further threshold of at least three shared regions is reported only as an exploratory sensitivity.

| Regime | Matched focal islands | Geographic blocks | Exact five-region signatures | Components, any overlap | Components, ≥3 shared regions |
| --- | ---: | ---: | ---: | ---: | ---: |
| Northern high latitude | **85** | **12** | **11** | **1** | 3 |
| Northern midlatitude | 876 | 59 | 146 | 6 | 14 |
| Tropics | 71 | 33 | 40 | 7 | 12 |
| Southern extratropical | 36 | 9 | 12 | 1 | not reported as a scored inference |

**Northern-high critical finding:** the two most frequently used source signatures cover **63/85 islands (74.1%)**. Under the any-shared-source graph, all 11 signatures are connected into one donor network. Even under the arbitrary >=3-shared-regions edge threshold, the largest network contains **8/11 signatures**. Thus the high-latitude SC source gain has NO source-*disjoint* replication in this matched sample. Individual source signatures, islands or spatial blocks may share GIFT donors and botanical taxa. Source graph components are not automatically independent histories either.

Midlatitude's largest any-overlap component includes **594/876** islands; tropical's includes **50/71**. Therefore even the regions with 6–7 linked components have highly imbalanced source sampling. The 7,618 original island-pair rows and two overlapping All/Direct scopes must not be presented as independent evolutionary replications.

## Robustness to the largest source signatures

For high-latitude strict SC/SI, the observed **equal-geographic-block** Brier advantage, and its changed-weight sensitivities, are:

| High-latitude SC/SI | All snapshot | All source-adjusted |
| --- | ---: | ---: |
| All 85 islands, equal geographic blocks | +0.000697 | **+0.000537** |
| All 11 distinct signatures, equal signature weights | +0.000679 | **+0.000450** |
| Exclude two most frequent signatures (only 22 islands left), equal remaining blocks | +0.000653 | **+0.000446** |

The sign remains positive, so it is not merely the two most frequent exact five-region signature groups creating the mean. But 22 residual islands have inadequate independence support for a confirmatory inference, and the original 11 signatures are **connected by donor sharing**. A bootstrap that treats the 11 exact signatures as independent units yields an apparently positive interval, but it is *not* a valid independent-source confidence interval. This inference must **not** be promoted to multiple independent source histories.

## Control matching strictness versus geographic support

The originally accepted alternative donor sets were chosen with known-genus donor-coverage discrepancy no larger than two percentage points (across both archived All ledgers and both traits). Without changing the chosen alternatives or focal outcomes, restrict the maximum discrepancy across both alternatives:

| Threshold | High-latitude focal islands | Geographic blocks | Source-adjusted All SC/SI block-equal Brier gain |
| --- | ---: | ---: | ---: |
| ≤2 percentage points | 85 | 12 | +0.000537 |
| ≤1 percentage point | 73 | **9** | +0.000739 |
| ≤0.5 percentage point | 49 | **7** | +0.000701 |

Positive signs remain, but both stricter cohorts lose the original **10-block minimum support gate**, so their point estimates are **diagnostic only**. The very small numeric gains should not be extrapolated to source-disjoint biogeographical histories.

## Changing the unit of inference changes conclusions elsewhere

For northern-midlatitude All SC/SI, the source-adjusted ledger yields **-0.000036** with equal spatial-block weighting, **-0.000085** with equal five-region signature weighting, and **+0.000137** when the **six donor-overlap components** are equally weighted (four components positive). The sign changes solely with the unit and weights, not the evaluated species or their trait states. This underscores that the data do NOT identify an unambiguous general positive mainland source effect.

For tropical All SC/SI, the corresponding values are -0.000151 (spatial block), -0.000232 (signature) and -0.000620 (seven donor-overlap components; only one positive). This is a contrast in predictive aggregation, not evidence for selection against self-compatibility in the tropics.

## Implication and stop rule

The **strongest positive finding remains the ability of known genus to predict held-out recorded floral symmetry and other categorical trait states**. The source geography residual is much smaller and sensitive to donor evidence. Northern-high SC has a reproducible *within the currently observed source-network* positive predictive contrast under both trait ledgers, but **none of the matched northern-high source sets represents a disconnected GIFT donor-region network**. It therefore cannot presently support a general macroevolutionary source-specific mechanism.

To make a new independently generalized finding, acquire or qualify **geographically/source-disjoint high-latitude donor networks and independent trait evidence**, or explicitly shift the estimand to a *single connected mainland-source system* and report only within-network predictability. The five-nearest mainland donor labels are geographic opportunities, not confirmed historical ancestors.

## Verification and release constraints

Public implementation `src/fd/source_signature_dependency.py` reconstructs the five-region source-signature overlap graph **from the licensed local source input**, checks missing/repeated donors, preselects one control seed, suppresses groups with <40 islands / <10 spatial blocks / <5 exact source signatures, and returns only region×trait×scope×ledger coarse aggregates. `tests/test_source_signature_dependency.py` contains synthetic shared-source graph, missing/duplicate source, repeated-seed and privacy checks. This module was exercised against the original locally retained private island-level score table, with its 24 supported aggregate rows independently reconciled against an earlier implementation (including region counts, source signatures, graph components and weighted score means). No protected island, genus, taxon, donor source ID, or per-island score is committed to GitHub; public GitHub CI tests logic only.


## New outcome-blind test: what if every original fd island is included?

The matched source/placebo subset is not the full FD island universe. To see if the lack of northern-high independent donor networks was created by matching/trait coverage, repeat the source graph construction on ALL **1,419 distinct G0 paired-cohort islands** *before using any floral or reproductive outcome*; use only the frozen G0 island IDs/regime and the same five `geo_k5` mainland donor entities.

| Regime | Original G0 islands | Exact donor signatures | Any-shared-donor components | Largest component by islands | >=3-shared-donor components |
| --- | ---: | ---: | ---: | ---: | ---: |
| Northern high latitude | **90** | **15** | **1** | **90** | 4 (largest 83 islands) |
| Northern midlatitude | 1,051 | 171 | 8 | 640 | 15 (largest 589) |
| Tropics | 218 | 90 | 8 | 127 | 22 (largest 60) |
| Southern extratropical | 60 | 29 | 5 | 44 | 8 (largest 43) |

All 1,419 original islands have their five geographical native-GIFT source labels. Northern-high **still has one connected geographic donor network even at the full eligible-island scale**, prior to self-compatibility measurement or placebo matching. Therefore increasing the number of pair comparisons using this unchanged data source cannot supply a source-*disjoint* northern-high test. This does **not** establish that the underlying plant lineages have one shared historical ancestry: the overlap graph follows the pragmatic GIFT source-region labels, not a resolved colonization phylogeny.

The independent, nonidentifying implementation is `src/fd/source_network_universe.py`, tested with `tests/test_source_network_universe.py`. It checks the original island's unique broad regime and exactly five distinct `geo_k5` donor entities, and publishes **only** region-level counts and connectivity; no species, island or donor IDs. This is an outcome-blind *feasibility/positivity* result, NOT another functional trait p-value or selection inference.
