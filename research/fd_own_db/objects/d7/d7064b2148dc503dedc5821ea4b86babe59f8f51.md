# Original categorical-distance packing: source-pinned public SES run

**9 November 2026; validated execution, not causal discovery.**
[CI #42579672053 (N1b-complete)](https://github.com/zuizui0223/fd/actions/runs/42579672053)
successfully ran the pinned public island aggregates, all existing tests, verified
geology crosswalk, 17D KDE diagnostic, spatial-independence audits and the new
[distance-native fixed-support SES analysis](../../docs/DISTANCE_SES.md).
The CI artifact named fd-origin-strata-analysis contains the *complete*
aggregate-only SES summary CSV (216 rows), the 32 paired spatial sensitivity
slots and SHA256 manifest. The repository stores this reduced report without
any potentially identifying island-level SES or full crosswalk.

## Main result: reference-relative proximity across BOTH island classes

Two endpoints use the ORIGINAL 496-profile colour/form set-Jaccard D:
mean distinct-pair distance (pairwise D), and mean distance to each
profile's nearest OTHER observed profile (nearest D). Nulls sample exactly
the same number of *distinct profiles*; SES=(observed-null mean)/null SD.

**N1 = leave-one-island-out pool of recorded profiles from other eligible
islands in the SAME broad analysis regime. This is NOT a genus-conditioned
source pool.** Equal30 = mean of 99 observed uniform thinnings to 30 profile
types, compared with 999 fixed-n=30 reference draws per relevant pool.
All and Direct datasets overlap and are NOT independent replications.

| Endpoint, N1, All | Continental (n=130) | Oceanic (n=130) |
|---|---:|---:|
| **Equal30 nearest D**, observed median | 0.292 | 0.276 |
| Equal30 nearest D, reference median | 0.334 | 0.337 |
| **Equal30 nearest D, median SES** | **−1.452** | **−2.205** |
| Fraction of islands with negative nearest-D SES | 100% | 100% |
| Full-support nearest D, median SES | −1.924 | −3.063 |
| Equal30 pairwise D, median SES | −0.008 | −0.715 |
| Full-support pairwise D, median SES | +0.018 | −1.280 |

Direct equal30 nearest-D N1 medians reproduce the negative reference
relative direction (continental **−1.500**, 123 islands; oceanic
**−2.154**, 119 islands; 100% below zero in both groups). Full support
Direct N1 nearest-D medians are **−2.094** (123 continental) and **−3.027**
(120 oceanic).

The unrestricted **N0 global496** reference likewise gives negative nearest-D
SES across both origins. For All equal30, N0 medians are **−1.833**
continental and **−2.390** oceanic (both 130/130 islands negative).
N0 and N1 show an *observational* concentration of recorded within-island
profile neighborhoods compared with uniform random profile combinations.
This does not by itself show environmental filtering: genus composition,
regional availability, trait reporting, and categorical repertoire structure
make observed states nonexchangeable with an arbitrary random reference.

## Critical falsification: regional frequency weighting removes the apparent universal packing

**N1b is NOT a new geology p-value.** N1b uses exactly the same leave-one-
island-out regional pool as N1 but weights each candidate profile by
the **number of OTHER eligible islands recording it**. Drawing n distinct
types without replacement this way allows common regional profiles to
be selected more often, so rare isolated categories cannot dominate
the null artificially. We retain N0 and N1 as the frozen 32 paired
exploratory comparisons and treat N1b solely as a post-hoc
description/negative-control sensitivity.

| Original D endpoint, equal30, All | Continental (n=130) | Oceanic (n=130) |
|---|---:|---:|
| **N1b nearest-D observed median** | 0.292 | 0.276 |
| **N1b frequency-weighted reference median** | 0.286 | 0.291 |
| **N1b nearest-D median SES** | **+0.181** | **−0.525** |
| Fraction of island SES < 0 | **40.8%** | **85.4%** |
| N1b pairwise-D median SES | −0.090 | −0.447 |
| N1b full-support nearest-D median SES | +0.151 | −0.973 |

The Direct result is closely aligned: equal30 nearest-D median SES
**+0.095 continental** (123 islands) versus **−0.537 oceanic**
(119 islands), with 42.3% versus 86.6% negative.

This **overturns the simplistic interpretation** that all sampled island
floras exhibit unusual packing relative to a well-formed regional
reference. Under regionally common-profile weighting, continental
floras are broadly in line with their comparison pool; oceanic floras
retain a *descriptive* excess of local proximity. The frequency
weighting is not genus or ancestor-source conditioning, and
island samples are spatially dependent. We did not apply
post-hoc pair p tests to N1b, so **no verified geological-origin effect**
is established by this pattern.

By broad region, the equal30 N1b All nearest-D median SES is
**+0.268 vs −0.287** (continental n=85, oceanic n=48) for northern
midlatitudes, and **+0.165 vs −0.667** (31 and 76) for tropical islands.
These are ecological/biogeographic descriptive contrasts, not
replicate-independent formal confirmations.

## Falsification: no geographically robust geological-origin separation

We retested **32 explicitly post-hoc scalar SES contrasts**: two nested
All/Direct scopes × two D endpoints × two nulls × two support sizes × two
geographical pairing definitions. The original 37 matched pairs have
only **5** independent sign groups under the conservative transitive
spatial-block graph; the alternative **both_verified** source/block-disjoint
rematching has **10** disjoint pairs and exact 2^10 sign configurations.
These are sensitivities, not a priori independent geological interventions.

- **32/32** scalar tests had uncorrected exact p >= **0.087890625**.
- Smallest Benjamini–Hochberg q = **0.446429**, and all BY q = **1.0**.
- Smallest raw p was **0.087891**, for Direct/N0/equal30 *pairwise D*
  in the 10 geographically verified matched pairs (mean
  oceanic−continental SES difference **−0.335**).
- The old 37 original pairs' conservatively grouped sign-test p values
  likewise do not support a unique oceanic category-packing effect.

Not rejecting here is not biological equivalence: the 10-pair analysis
has low independent geographic support, and N1 remains observational.
Unmatched islands' negative SES cannot be turned into independent 130-v-130
geological-origin significance by counting individual islands as
independent ecological histories.

## Interpretation and next identification gate

**Defensible observation:** Relative to *uniform* broad-regional profile
sampling, local recorded island profile sets are unusually close by nearest
D in both geological classes. Relative to **frequency-weighted** regional
sampling, the apparent universal pattern disappears in continental floras;
residual oceanic proximity is descriptive, and not lineage-, source- or
independent-origin-confirmed.

**Not demonstrated:** A universal selective shrinkage of island floral
niches; source/lineage-independent oceanic compression; pollinator
functional adaptation; or colonization-mediated category exclusion.

A stronger mechanism-oriented analysis needs an explicitly donor-held-out
genus × plausible mainland-source null (N2), independent detection and
flower trait provenance, plus better independent spatial replication.
The high-dimensional Gaussian KDE HDR volume remains an appendix
diagnostic rather than the primary inferential statistic.

## Reproduction and privacy

Run: python -m fd.distance_ses --input work/island_public
--g0 work/fd_g0 --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv
--output results/distance_ses --draws 999 --thinnings 99

Artifact: [complete public N1b-verified CI zip with SES aggregate outputs](https://github.com/zuizui0223/fd/actions/runs/42579672053).
Only grouped summaries with n>=5 and aggregate matched test statistics
are published. Original island identifiers, matched endpoint IDs,
source archipelago labels, GIFT donor identities and individual trait
assignments are not redistributed.
