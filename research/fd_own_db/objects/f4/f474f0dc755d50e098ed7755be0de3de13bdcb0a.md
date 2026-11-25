# fd: Functional diversity consequences of pollinator shortage and putative guild switching

**Research design v1 · 13 November 2026 · source-backed exploratory test, not a confirmed causal effect.**
Independent analysis: [fd PR #15](https://github.com/zuizui0223/fd/pull/15), original restricted-data runner: [island PR #325](https://github.com/zuizui0223/island/pull/325).

## Question and non-overlapping hypotheses

**What changes in island floral functional proxy diversity are associated with an independently supported shortage of potential pollinator channels, and when is change in the floral trait repertoire distinct from functional contraction?**

- **H1, contraction:** for an island with independently supported pollinator *availability shortfall*, standardized floral-trait dispersion and joint-profile diversity decrease relative to a valid comparator.
- **H2, assembly reorganization:** floral joint trait states, colour–architecture coordination, and open-flower prevalence change *even when marginal trait diversity remains stable*. The same marginal Q can accompany distinct joint Q; use a fixed-profile null instead of assuming they are interchangeable.
- **H3, visitor/guild switch:** a directly documented, temporally or island–mainland *same-plant* change in actual visitor guild composition predicts a shift in the floral trait distribution, controlling source/lineage and visit effort. An occurrence of a flower-visiting bird on an island with adequately nondetected Bombus is a **candidate opportunity**, not a measured bee-to-bird shift.

The ecological estimand is the **effect of observed partner-community change on flora composition and/or reproductive output**, not how well colour predicts GloBI bibliographic visitor categories. Morphological traits alone are imperfect functional proxies.

## Current measurable source-backed endpoints

1. Direct species colour, floral form, floral symmetry — exactly 30 accepted species sampled without replacement per local flora; sample equal 30 from its independent candidate **five-nearest GIFT native mainland region lists**, if that pool has ≥30 species with valid Direct trait evidence. **Not** a reconstructed ancestor.
2. Marginal categorical Rao Q, \(Q_m=\frac13\sum_{j}(1-\sum_s p_{js}^2)\). With equal trait weights, this is **blind** to recombination of traits on individual species.
3. Joint categorical Rao Q, \(Q_J=1-\sum_{c,f,s}p(c,f,s)^2\), and Hill number \(N_{2,J}=1/\sum_{c,f,s}p(c,f,s)^2\), both sensitive to *actual recorded* colour×form×symmetry profiles. Multistate and genus-consensus trait records are **not** expanded into hypothetical flowers.
4. Colour × intact architecture \(I(C;[F,S]) - E_{\rm perm}[I(C;[F,S])]\) using the **analytic hypergeometric expectation under the exact fixed-marginal colour shuffle**, not noisy 49-shuffle Monte Carlo.
5. Open-radial form prevalence and actinomorphic symmetry prevalence, **structural visitor-access proxies**, not validated species-specific guild matching.
6. A subsequent *separate* reproductive-assurance endpoint will use original directly supported SC/SI, autonomous selfing and mating-system states with trait-specific coverage and frozen Chapter 1 coding. Do NOT demand complete cases for all seven traits or infer unknown reproductive traits from floral form.

## Pollinator exposure and loss identification

Freeze *source available* + *observed qualified guild occurrence* = **retained**; *source available* + *zero qualified guild occurrence after the prospectively fixed within-guild effort gate* = **adequately nondetected/disrupted**. Treat all other channels as **unresolved**, never absent. The five guilds are Bombus, other bees, Lepidoptera (moths and butterflies pooled), Diptera and flower-visiting birds.

For each guild × broad geographic regime, an associational comparison is **suppressed unless** both observed states have at least 10 islands and each ≥10 spatial blocks **and** at least 3 nearest main-source signatures contain both states. This is a *minimum positivity/support test*, not a license for historical or causal claims. Richness of observed channels is highly correlated with evidence opportunity and must never be used as an uncorrected ecological exposure.

Strict native only is reported separately and on the **identical eligible island identities** in a within-island paired sensitivity; all-observed and native-only are nested observations, not external confirmations.

## Original source rights / independence

- Original progressive Direct flora + GIFT native source and geography: upstream original source artifact ID **10058653212**, SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Frozen source-positive GBIF island guild evidence: original source run **34750153702**, artifact **10314984302**, SHA256 `c643b8432e04de0392867028810606f0f60f5a1d45b94bca1e28a1696ffc292d`.
- Analysis is performed inside `island` Actions with separately checked-out `fd` branch. Only anonymous region and source/block-level aggregate JSON is released; no private taxon/species, island/source region, provider/reference, original geography or interaction IDs.
- Geographic spatial-block and nearest GIFT-source signature bootstrap uncertainty are **two complementary, overlapping, conditional descriptive levels**, not phylogenetic histories or independent evolutionary replication. Use both for exploratory consistency but do not treat either as a causal selection null.

## What would identify actual pollinator decline or shift?

Build standardized **dated island × grid/site × guild × survey-dataset** *effort-controlled* repeat-detection histories, with at least two independent windows, for example 2000–2010 and 2015–2025. Define change in effort-corrected detection/occupancy, **not** GBIF record-count change. Connect temporally matched plant-island floras or site observations so that change in floral assemblage diversity can be estimated. Existing static global plant species trait assignments allow **species sorting/community reassembly**, not within-species phenotype evolution.

A real shift needs the **same accepted plant species or independently justified island–mainland sister lineage**, actual flower-linked visitors in each population, sampling years/site effort, named guilds with adequate nondetection, and preferably pollen deposition or seed set. GloBI v0.9 has no event-time field and disproportionately lacks same-island native species-linked independent references; its bibliographic frequency is not pollinator abundance or measured service. Separate Lepidoptera moths from butterflies only when independent visitor records support it.

For any newly acquired longitudinal cohort, freeze:
- treated and control sampling-frame dates and environmental/observation covariates **before** reading floral outcomes;
- within-island paired change of trait-composition responses and guild-specific exposure; island and year/campaign effects, spatial/source blocks, modern land-use and climate nuisance controls;
- isolation/source-proximity and comparable surrounding mainland/islands as negative controls;
- guild effort/opportunity and native-status strict sensitivity;
- exclusion of same taxon donor leakage and source/provider-only effects;
- report unidentifiable contrasts explicitly, instead of treating low counts and null non-rejection as confirmation.

**Status: current source-qualified cross-sectional guild data can test floral assembly differences but do NOT identify the impact of historical pollinator decline or an evolutionary shift.**
