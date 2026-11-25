# Native-mainland source geography versus genus-associated functional-trait prediction

**Exploratory empirical analysis, 10 November 2026.** This report adds a genuinely source-conditioned functional-trait comparison to fd PR #13. The original taxon/island/source rows remain private. This is **not** a causal pollination-service, native colonization or independent evolutionary adaptation test.

## Question and measurement

Do the five geographically closest native GIFT mainland source floras predict recorded **self-compatibility** (explicit SC versus SI) and **floral symmetry** (actinomorphic versus zygomorphic) beyond the known plant genus, and beyond other outcome-blind GIFT source assignments with comparable donor evidence?

Two category-specific binary indicators are measured, not experimentally realized functions. The analysis is not the earlier Chapter 1 H1 seven-trait combined coding: the explicit strict original-source ledger contains **17,885 All / 6,726 Direct** known SC/SI taxon names and **56,202 All / 32,431 Direct** known actinomorphic/zygomorphic taxon names before island-cohort filtering. The previously archived fd PR #13 traitwise denominators use a different eligibility/coding surface, so they are **not silently equated with these strict-source denominators**. Both archived ZIPs were compared and source-scale `direct_species_trait_ledger.csv.gz` exactly matches the corresponding progressive-checkpoint member.

## Pinned evidence and correct unit

- Original rights-protected `island` progressive ZIP: GitHub Actions artifact `10058653212`, SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Original rights-protected source-scale trait ZIP: artifact `10042377593`, SHA256 `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`.
- fd public G0/G1 source-pinned cohort: artifact `13741845313`, **7,618 sampled pairs, 1,419 distinct islands**. Pairs are heavily overlapping; All and Direct scopes are nested, not independent replications.
- The GIFT native mainland checklist has **1,666,186 region × taxon records across 748 regions**, from the original ZIP. Island source choices are the recorded `geo_k5` five nearest mainland entities; they are NOT known ancestral source populations.

A deterministic species-name SHA256 fold split (first digest byte modulo five) holds every evaluated species out from **all** native-mainland trait-donor training regions in its test fold, including when the same name occurs on many islands. Global-genus probability is estimated from geographically pooled heldout-species-excluded native GIFT donors, shrunk toward the fold-excluded global prior. Local-source genus probability is `(k_nearest_5 + tau * p_global_genus) / (n_nearest_5 + tau)`, with `tau=8` reported and `tau=2,20` sensitivity. No unknown trait states are set to zero.

Because one plant species is repeated across many islands, loss differences are first averaged within each island, then spatial blocks are given equal weight. 1,999 block-bootstrap resamples give **conditional descriptive intervals**; shared botanical histories and training donor pools still violate strict independent-lineage replication.

## First diagnostic: local GIFT source improves floral symmetry over a generic genus model

For the **original 1,419-island population** and `tau=8`, before any alternative-source placebo matching, equal-spatial-block Brier improvement for actinomorphy is positive in every region and both nested evidence scopes:

| Region | All improvement | Direct improvement |
|---|---:|---:|
| Northern high latitude | +0.000735 | +0.001364 |
| Northern midlatitude | +0.000676 | +0.001096 |
| Southern extratropical | +0.000539 | +0.000792 |
| Tropical | +0.000412 | +0.000725 |

For SC/SI this first diagnostic is positive in the two northern regions but negative in the tropical and southern-extratropical populations. **This is not yet source-geography specificity:** any pool with more suitable same-genus training examples might improve over the global prior.

## Decisive negative control: actual GIFT source versus donor-matched pseudo-sources

For each focal island, choose exactly **two alternative source assignments** from other islands in the same broad regime, different spatial block, with native source-pool size 0.5–2×, different source entity sets (five-entity Jaccard <0.8), and All-scope **known same-genus trait-donor coverage matched within 2 percentage points for BOTH traits**. Selection uses taxon/genus identities and observed trait *availability* but NOT the observed binary trait states. Cap the outcome-blind candidate search at 128 per focal; choose the nearest caliper matches. Repeat using two fixed random seeds.

The first matching yielded **1,157 focal islands** (88 northern high, 906 northern mid, 44 southern extratropical, 119 tropical); alternative seed **1,163** (88/913/44/118). Matching retained 1,129 focal islands in common. Median alternative/actual native checklist pool-size ratio = **1.02**, median geographic source-set Jaccard = **0.429**. Post-fold placebo-minus-real known-genus donor coverage is only about 0–0.7 percentage points across the reported trait-region scopes, although perfect evidence matching is not achieved. These are alternative geographically assigned donor sets, NOT independently reconstructed historical colonist sources.

**Main outcome:** actual source minus two matched alternatives, mean equal-block proper **Brier loss advantage**, positive favors actual source. Parentheses show conditional 95% spatial-block bootstrap intervals from **first** control seed:

| Trait | Region | All advantage | Direct advantage |
|---|---|---:|---:|
| SC/SI | Northern high latitude | **+0.000712** (+0.000254,+0.001356) | **+0.000930** (+0.000308,+0.001718) |
| SC/SI | Northern midlatitude | +0.000364 (+0.000033,+0.000690) | +0.000439 (+0.000055,+0.000885) |
| SC/SI | Southern extratropical | −0.000126 (−0.000296,+0.000058) | −0.000174 (−0.000475,+0.000096) |
| SC/SI | Tropical | −0.000089 (−0.000339,+0.000147) | −0.000118 (−0.000430,+0.000140) |
| Symmetry | Northern high latitude | −0.000001 (−0.000040,+0.000038) | +0.000015 (−0.000059,+0.000090) |
| Symmetry | Northern midlatitude | −0.000035 (−0.000110,+0.000034) | −0.000056 (−0.000146,+0.000034) |
| Symmetry | Southern extratropical | +0.000025 (−0.000152,+0.000212) | +0.000018 (−0.000196,+0.000222) |
| Symmetry | Tropical | −0.000019 (−0.000090,+0.000048) | −0.000010 (−0.000103,+0.000075) |

**The symmetry source signal disappears versus donor-matched source controls** across all four regions in the primary seed, with nested Direct confirming the lack of a stable positive advantage. An alternative matched control seed yields no robust positive symmetry advantage either.

The positive **northern-midlatitude SC** result is fragile: with a second source-control seed its intervals cross zero; on the **same 879 focal islands and 60 spatial blocks** common to both seeds, the SC All block-equal gains are +0.000345 / +0.000314, and Direct +0.000421 / +0.000381, **both seeds' bootstrap intervals include zero** on this common cohort. In contrast, the northern-high-latitude SC signal is +0.000712 All/+0.000930 Direct under both seeds—but **the selected control pairs are identical across seeds**, so this is one finding across only **14 blocks**, not an independent second replication.

## Native-status sensitivity and support loss

Strict island-native status from the original `chapter1_status_flora` is much less complete. Restricting to >=20 *recorded binary trait carriers* per island leaves the matched north-midlatitude native population at **142 All / 141 Direct focal islands** for SC, with 30/29 spatial blocks. On the original outcome-blind alternative-source assignment:

- All SC advantage **+0.000754** (block 95% interval +0.000069,+0.001516).
- Direct SC advantage **+0.000939** (−0.000032,+0.002010), unresolved.
- Symmetry on the native-only north-midlatitude target yields All +0.000053 (−0.000075,+0.000171) and Direct +0.000021 (−0.000188,+0.000195), both unresolved.

Other regions fall below the **40-island release/inference-support minimum** under these controls; no regional native extrapolation is licensed. The matching calipers were computed on the broader All known-trait universe, not reoptimized on native-only support. A WCVP regional-native-compatible label is not exact island nativeness.

## Direct source/genus expected pairwise functional difference

Using the original 7,618 pair cohort and excluding genera shared in each pair's **full** recorded flora (including unmeasured species), calculate observed absolute prevalence gaps and independent-Bernoulli predictions from the **out-of-name genus** and **out-of-name local-source genus** probability models, requiring 20 measured carriers per island. For the All source model, mean observed-minus-source expectation is:

| Region | SC/SI excess gap | Symmetry excess gap |
|---|---:|---:|
| Northern high latitude | +0.03455 | +0.00595 |
| Northern midlatitude | +0.02734 | +0.00510 |
| Southern extratropical | +0.02370 | +0.01689 |
| Tropical | +0.02267 | +0.01415 |

Direct SC/SI excess is also positive in the four regions (+0.02790,+0.01665,+0.02613,+0.02265 in that order). These **are descriptive model residuals** and are not a valid source/genus-conditioned evolutionary permutation test. The independent-Bernoulli model does not reproduce species phylogenetic/detection correlations, source-arrival probabilities or trait covariance; the same taxa and islands appear across many pairs. Random nonoverlapping-block pair re-selections are correlated sensitivity diagnostics, NOT p-values.

A 700-draw Monte Carlo check of the normal `E|P_a-P_b|` approximation on 20 sampled eligible pairs per region/scope/trait (320 sampled pairs total) gives pooled mean normal-minus-direct-simulation biases of about −0.00045 to +0.00008 and mean absolute per-pair differences ≈0.0010–0.0016 depending on trait/scope. The model residuals are not dominated by this *numerical approximation*, but could absolutely arise from missing ecological, taxonomic or observation processes.

## Interpretation and nonclaims

The result discriminates between **predictability from known genus** and **additional information specifically associated with the nearest mainland GIFT source assignment**. Recorded floral symmetry is highly predictable from genus and has no robust additional nearest-source assignment gain once same-genus donor evidence is matched. Recorded SC/SI, a more variable trait, has a small northern geographic source signal with strong qualifications and does NOT show a common global signal.

This does NOT demonstrate geographic isolation *causes* self-compatibility, source-specific selection, convergent within-lineage evolution, historical mainland ancestry, equivalent floral function, pollinator-mediated trait evolution, or higher seed set. A model-residual difference is not evidence of natural selection. A strict ecological mechanism would require observed floristic status, independent source-archipelago histories, phylogenetic residual controls, trait-ascertainment modeling, and biological pollen-delivery/seed-set endpoints.

**Methods/artifact safety:** Only this broad-region report, aggregate CSV/manifest, and synthetic validation code should enter fd. All original GIFT source names, accepted species, genus/island IDs, source-assignment rows and per-species/per-pair predicted scores remain in protected local inputs. The public GitHub CI cannot replay the licensed source, and CI success alone must not be called independent validation of these empirical estimates.
