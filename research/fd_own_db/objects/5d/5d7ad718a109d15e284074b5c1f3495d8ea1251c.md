# Floral-category repertoires are locally predictable, but not independently archipelago-specific

**9 November 2026 — complete exploratory fd public-aggregate analysis.** [Methods and predictor audit](../../docs/ARCHIPELAGO_FLORAL_HOLDOUT.md); [source pipeline](../../src/fd/archipelago_cv.py); [all 12 comparisons](archipelago_holdout_tests.csv). The scientific outcome is **prediction of recorded floral colour × floral-form state-set distributions**, not historical colonization, phylogenetic independence, or plant-pollinator functionality.

## 1. Why the independent archipelago analysis was necessary

Previous fd PR4/PR8 showed exact categorical differences in 37 continental–oceanic matched islands when pairs were swapped independently. But PR9 found those 37 original pairs belonged to only **five connected geographic groups**; whole-group exact swaps did **not** support the origin effect. PR10 re-matched completely geographically disjoint origin pairs (10 fully labelled pairs) and found **zero BH/BY-supported geological contrasts**. The 17D projected Gaussian KDE volume difference also disappeared after profile-richness adjustment. Increasing the count of similar islands is therefore not sufficient to identify a geological mechanism.

This project instead asks: **Does an archipelago's observed flower-colour/form repertoire help predict the distribution of floral categories on a different island, after accounting for generic same-region information and physical island proximity?**

Earlier island functional studies already covered island plant trait space and trait equilibrium: Barajas Barbosa et al. (2023; Nature, DOI 10.1038/s41586-023-06305-z) documented trait assembly in Tenerife, and Schrader et al. (2023; JBI, DOI 10.1111/jbi.14526) established an Equilibrium Theory of Island Biogeography for Traits. This cross-sectional fd analysis is **not the first island trait-space study** and does not test temporal trait equilibria. Its distinct diagnostic is **heldout island floral repertoire prediction with an explicit geographic-only negative control and source-group-level inference**.

## 2. Data and spatially honest units

The same source-pinned public island×profile aggregate table as fd PR321 and a Weigelt archived named-archipelago crosswalk. Of the **1,442 islands with ≥50 complete colour/form names under both All and Direct**, **711 islands carry a nonmissing Weigelt archipelago label** (112 named groups). There are **45 region × archipelago subgroups with ≥4 eligible islands**, but **44 distinct named archipelagos** supply evaluated targets. Keep All and Direct on the same physical island and donor cohorts; do not treat them as two independent validations.

Each island is encoded as a normalized **equal-distinct-profile** vector across the same 496 original recorded flower colour-state-set × flower-form-state-set profiles. No multistate description is expanded into hypothetical flowers. Only the region and archived source group determine eligibility; target flower profiles do **not** determine neighbor selection.

Predict each held-out target using either an island-equal **outside-target-archipelago regional average**, or **three geographically nearest other islands of its own archipelago**, or **three nearest islands regardless of archived archipelago membership**. In each local model, use three donor islands mixed 50:50 with the identical outside-group regional average. The target itself is excluded from every source. A stricter control uses three **different-archipelago** islands selected by matching source-to-target distances within a 3-fold tolerance, with within-group third donor ≤1500 km.

Evaluate error as squared distance in exact 496-profile identity space (L2) and in the full original colour/form set-Jaccard similarity kernel. Neither is a measurement of pollinator performance, flower-colour frequencies within real plant individuals, or 2D/17D KDE volume.

## 3. Results: local prediction succeeds, geographic source identity adds no clear signal

**Main heldout sample:** **606 islands** representing **44 named Weigelt archipelagos** (45 regime-specific subgroups). Treat the **archipelago**, not the individual island, as the uncertainty/sign-flip/bootstrap unit; one average gain per archipelago. There are 9,999 directional group-sign permutations and 4,999 group-bootstrap draws. All 12 correlated post-hoc tests (three comparisons × two endpoints × All/Direct scopes) share one BH/BY family. The closest-island negative control was added **after an initial 8-test pilot** revealed possible distance confounding; this is fully transparent post-hoc work, not an independent preregistered test.

### Mean predictive errors (smaller is better)

| Evidence scope / original-profile metric | Same-regime outside group | Same archipelago, three donors | Geographic-nearest three (group label ignored) |
|---|---:|---:|---:|
| All, exact-profile L2 | **0.006861** | **0.005077** | **0.005078** |
| All, original categorical similarity-kernel | **0.010609** | **0.006203** | **0.006209** |
| Direct, exact-profile L2 | **0.007351** | **0.005473** | **0.005473** |
| Direct, original categorical similarity-kernel | **0.011137** | **0.006574** | **0.006573** |

The same-archipelago model improves mean island-level exact-profile error by **26.0% All and 25.5% Direct** compared with outside-group regional average, and original categorical kernel error by **41.5% All and 41.0% Direct**. At equal source-group weighting, mean improvements remain positive; all four broad contrasts have group-based directional p=0.0001, **BH q=0.0003 / BY q≈0.00093** over 12 declared tests. This is evidence for **geographically local predictability of recorded flower-profile composition**, not an independently identified archipelago-identity mechanism.

### The decisive negative control: geography alone predicts equally well

The three geographically closest eligible islands were **exactly the three within-archipelago donors for 581 of the 606 targets (95.9%)**. As a result, the source-identity model has practically NO incremental predictive skill beyond geography-only nearest neighbors:

| Source group vs geography-only nearest three | All mean archipelago gain | All p | Direct mean archipelago gain | Direct p |
|---|---:|---:|---:|---:|
| Exact-profile L2 | **−0.00000564** | **0.5913** | **−0.00000891** | **0.7033** |
| Original categorical similarity kernel | **−0.00002162** | **0.8565** | **−0.00003251** | **0.9296** |

Positive gain would indicate an advantage of selecting donors using explicit archipelago identity over selecting neighbors from distance alone. The point estimates are instead tiny and negative. All confidence intervals for the gain straddle zero, and none has BH/BY support.

A second, more restrictive geographical negative control replaces same-group donors with **different-archipelago donors selected to match target distances**. The original sample with all three matching calipers satisfied is only **51 targets across 18 named archipelagos** (33 north-midlatitude, 18 tropical; other regimes lack source overlap). The same-group exact-profile advantage is not supported: All **p=0.0891**, Direct **p=0.0813**; original-kernel All **p=0.1564**, Direct **p=0.3515**. Their group-bootstrap 95% improvement intervals all include zero.

**Therefore:** floral trait profiles are more predictable from **nearby island floras** than from a broad regional average, but the current evidence does **not** show a distinct Weigelt-archipelago effect after accounting for physical proximity. Lack of support in the strict 51-target subset is not proof that all genuine archipelagos have identical floras.

## 4. Effective geographic replication and novel interpretation limits

- Regional group support (all scopes use the same 606 targets): northern midlatitude **450 targets/31 named groups**, tropical **123/11**, southern extratropical **25/2**, northern high latitude **8/1**. Named source groups may span multiple geographic regimes, leaving **44** distinct archipelagos rather than 45 independent units.
- Median distance to third nearest same-group donor: **47.5 km** in northern midlatitude, **201.7 km** in tropical. A nearby donor can share plants, introduced species, climate, geology or horticultural history with the test island.
- A Weigelt 'archipelago' is a published geographic label, **not** a recovered historical evolutionary source pool. Even the 44 named units may share plant species and source lineages; group-sign flips assume more independence than the original botanical records can demonstrate.
- Source-specific species and genus identity are **not available in the public island×profile aggregates**. A leave-one-island-out comparison may use donor flora containing *the exact same species* present on the target island. This means its successful prediction may reflect taxonomic overlap and shared environmental conditions; it is **not out-of-species evolutionary transfer** or a demonstration of selection.
- The two evidence scopes overlap; the three source baselines share regional data; all 12 statistical tests and donor choices are **post-hoc** and correlated. Conservative BY correction does not identify a causal spatial/lineage mechanism.

## 5. Paper positioning and next discriminating step

The results support this empirical proposition:

> **Recorded island floral trait repertoires are strongly locally predictable, but a source-archipelago label adds no demonstrable information beyond the actual nearest neighboring island floras.**

They do **not** establish an intrinsic floral hypervolume collapse on geologically oceanic islands, adaptation to depleted pollinator services, an archipelago-specific causal dispersal filter, or species colonization/failure probabilities.

**The next scientifically decisive separation is species turnover vs floral-profile turnover**: if two islands share few plant species yet preserve the same colour×form repertoire, that would suggest assemblage-level phenotypic redundancy. If floral similarity disappears once shared species and genera are removed, the current neighborhood predictability may simply reflect taxonomic overlap. Testing this requires the separately rights-checked identity-level island flora/trait checkpoints, and must release aggregate group summaries only. Independently resurveyed floras could further distinguish arrival from historical observation; the current single-snapshot data cannot.

## 6. Reproducibility and rights

- Full implementation: [src/fd/archipelago_cv.py](../../src/fd/archipelago_cv.py), with **66 successful synthetic and integration tests** (two preexisting nonblocking warnings) and replay from the pinned original public aggregate workflow.
- [Methods and multiple-comparison protocol](../../docs/ARCHIPELAGO_FLORAL_HOLDOUT.md), including the documented post-hoc geography-only negative control.
- [Complete twelve statistical comparisons](archipelago_holdout_tests.csv), [predictive error summaries](predictive_errors.csv), [regional and geographic donor-overlap audit](archipelago_support_by_region.csv).
- GitHub Actions artifact from the source-pinned full run contains the exact floating-point original results and input SHA256 provenance. Only public aggregates/source checksums are retained; no named source group labels, island IDs, matched donor identities, species names, genus assignments or restricted source evidence are published.
