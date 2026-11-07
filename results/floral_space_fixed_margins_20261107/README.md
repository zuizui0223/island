# Fixed-margin floral-profile sharing sensitivity

This **exploratory complementary null** uses the same island pairs as PR #321 but retains all taxon
names in the whole-island incidence. It holds each island's recorded number of floral profiles
and each profile's number of sampled islands fixed, then makes symmetric Curveball trades.
It **does not** condition on genus/source pool or remove shared taxa, so it is not a stronger
version of the original conditional test. A significant result would not establish evolution or
pollinator replacement. Markov-chain tail frequencies are approximate, not exact p-values.

| Scope | Region | Pairs | Islands | Observed J | Null median | ΔJ | Approx. BH q | Mixing |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| all | northern_midlatitude | 2000 | 1051 | 0.3406 | 0.3689 | -0.0283 | 0.0073 | pass |
| all | northern_high_latitude | 2000 | 90 | 0.4793 | 0.5116 | -0.0323 | 0.0073 | pass |
| all | tropical | 2000 | 218 | 0.3468 | 0.3587 | -0.0119 | 0.0073 | pass |
| all | southern_extratropical | 1618 | 60 | 0.3217 | 0.3293 | -0.0075 | 0.0133 | pass |
| direct | northern_midlatitude | 2000 | 1051 | 0.3368 | 0.3670 | -0.0302 | 0.0073 | pass |
| direct | northern_high_latitude | 2000 | 90 | 0.4777 | 0.5000 | -0.0223 | 0.0073 | pass |
| direct | tropical | 2000 | 218 | 0.3428 | 0.3529 | -0.0102 | 0.0073 | pass |
| direct | southern_extratropical | 1618 | 60 | 0.3108 | 0.3179 | -0.0071 | 0.0184 | pass |
| all_single_state | northern_midlatitude | 2000 | 668 | 0.5200 | 0.5455 | -0.0255 | 0.0073 | pass |
| all_single_state | northern_high_latitude | 171 | 40 | 0.6500 | 0.7353 | -0.0853 | 0.0073 | pass |
| all_single_state | tropical | 2000 | 104 | 0.6216 | 0.6296 | -0.0080 | 0.0073 | pass |
| all_single_state | southern_extratropical | 391 | 30 | 0.5778 | 0.5789 | -0.0012 | 0.9476 | pass |
| direct_single_state | northern_midlatitude | 2000 | 668 | 0.5227 | 0.5455 | -0.0227 | 0.0073 | pass |
| direct_single_state | northern_high_latitude | 171 | 40 | 0.6400 | 0.7273 | -0.0873 | 0.0073 | pass |
| direct_single_state | tropical | 2000 | 104 | 0.6071 | 0.6111 | -0.0040 | 0.4096 | pass |
| direct_single_state | southern_extratropical | 391 | 30 | 0.5556 | 0.5641 | -0.0085 | 0.1026 | pass |

The observed profile set of any island is reused across its pairs; pairs are not independent.
Rows with zero shared taxon names are enumerated only as a descriptive sensitivity,
because their numbers vary and the permuted incidence cannot preserve taxon identity.
Null 2.5–97.5% quantiles are randomization reference ranges, not effect confidence intervals.
See [summary](permutation_summary.csv), [draws](permutation_draws.csv.gz) and [manifest](manifest.json).
