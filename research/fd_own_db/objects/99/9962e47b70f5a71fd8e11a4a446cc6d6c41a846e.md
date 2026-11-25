# Floral donor geography improves similarity-aware probability scores, not exact-category log loss

**2026-11-10 JST — privately reconstructed, source-verified, exploratory macroecological prediction audit.**

[Executable five-geographic-fold predictor](../../src/fd/n2_crossfit_macro.py) · [proper scoring implementation](../../src/fd/proper_scoring.py) · [source/donor feasibility audit](../20261109_n2_donor_support/README.md). No GIFT raw native floras, island/name mappings, species trait ledgers or geographic source assignments are republished. This is **not** a new geographic-origin selection test, fossil dispersal reconstruction or pollinator function measurement.

## Scientific question and data unit

**Does a candidate nearby native-mainland source repertoire help predict floral *similarity* beyond the globally known genus of each recorded island taxon?**

Frozen upstream [island PR #321](https://github.com/zuizui0223/island/pull/321): 496 full colour-state-set × flower-form-state-set profiles (467 observed in Direct scope), original equal-weight state-set Jaccard D. The original traits and exact profile frequencies were reconstructed from the SHA-pinned source and progressive checkpoints. A 32-bit independent FNV audit against the pinned public profile reference verified **all 496 profile IDs, intact colour/form state sets and per-profile All carrier-name counts**. No unobserved colour/form Cartesian combinations are generated.

All includes **1,533 eligible islands / 328,933 island×recorded taxon-name carriers / 16,386 distinct complete-profile names**. Direct includes **1,443 / 299,235 / 8,244**, and is an overlapping evidence subset, **not independent replication**.

## Holdout scheme and exact predictive estimand

- Entire source geographic spatial blocks are allocated to **five outcome-blind, disjoint folds**, balanced by All eligible island numbers (307/307/307/306/306); Direct inherits the same spatial-block fold assignment.
- Each fold removes **EVERY accepted taxon name observed in its entire test-island fold from EVERY GIFT 3.2 native mainland source region**, *before* global and local genus trait distributions are constructed. Repeating a name across 20 GIFT entities cannot leak it into training.
- Known actual island genus is supplied to both global-genus and local-genus models; neither predicts genus arrivals. The local model updates the independently held-out global-genus prior with intact colour/form donor profiles listed as native in a candidate island's five geographically nearest GIFT mainland source entities (geo_k5). Both use hierarchical shrinkage 8 to global genus.
- Per test name, evaluate a full 496-category predicted probability vector, with **equal weight across islands** for the main result. Evaluate the same names under *global no-genus prior*, *global genus* and *geographically selected local genus*. No unknown or unrecorded island taxon is coded as biological absence.

The existing expected untransformed D loss is **not a proper probability score**: its Bayes action can collapse to one modal category. We instead add two proper scores:

1. **Brier**: sum over category types of (predicted probability minus observed one-hot indicator) squared. **Strictly proper** for the full categorical distribution.
2. **K-score**: \\(p^T Kp -2(Kp)_y +K_{yy}\\), where \\(K=1-D\\). The frozen 496×496 source kernel is numerically positive semidefinite (minimum eigenvalue approx **−2.85×10⁻¹⁴**, a floating-point zero). Its rank at eigenvalue threshold 1e-8 is **76**, so this score is **proper but not strictly proper**; it tests similarity-aware predictions, **not** uniquely the probabilities of all 496 exact categories.

All improvements below are *loss(global genus) minus loss(local genus)*, positive indicating a better geographically assigned native-source predictor. Intervals were computed from 4,999 complete **spatial-block** bootstrap resamples (170 All blocks; 159 Direct blocks), retaining all islands within each sampled block. All estimands and scoring choices are exploratory/post-hoc; these intervals are not nominal preregistered hypothesis confirmations.

## Entire eligible island population: spatial-block inference sensitivity

| Increment when including local source (vs global genus) | All (1,533 islands) | Direct (1,443 islands) |
|---|---:|---:|
| Exact-profile log loss, proper | **+0.00011** [−0.00451,+0.00379] | **+0.00019** [−0.00606,+0.00517] |
| Brier, strictly proper | **+0.00221** [+0.00095,+0.00335] | **+0.00178** [−0.00022,+0.00349] |
| Original D K-score, proper non-strict | **+0.00829** [+0.00633,+0.01004] | **+0.00948** [+0.00655,+0.01191] |
| Expected original D, *not* proper | +0.00952 [+0.00769,+0.01129] | +0.00898 [+0.00686,+0.01077] |

**Global genus identity itself** adds a much larger heldout exact log-loss improvement relative to the non-genus baseline: All **+0.64073** [0.58738,0.69103], Direct **+0.36006** [0.29703,0.41254]. Nearest-source identity thus does **not** meaningfully improve *exact* heldout category likelihood after genus, even though local donor geography still shifts probability mass toward **similar** flower colour/form profiles.

Source support remains sparse: local independent genus donor support covers only **24.26% All / 19.56% Direct** of evaluated island-name carriers, versus globally available genus donor support **69.07% / 60.36%**. For most names, local model falls back to global genus or reference priors.

**Fold stability:** K-score source improvement is positive in **5/5 disjoint spatial folds**, All range **+0.00382 to +0.01153**, Direct **+0.00442 to +0.01349**. Exact-category log-loss improvement is positive in only 3/5 All and 2/5 Direct folds. K-score improvements are also positive as descriptive region-level means in all four broad biogeographic regimes, but these regions/folds share an underlying species profile map and are **not independent experiments**.

[Machine-readable island-weighted source/geography blocked summary](island_block_bootstrap.csv) · [all five fold signed score audit](five_fold_sign_audit.csv).

## Falsify repeated-species pseudoreplication: one taxon name = one vote

A single recorded name appears in many islands. **91.6% All / 95.5% Direct** of the evaluated island×name carrier rows come from names recorded on **at least 10 islands**. Spatial-block bootstrap alone cannot eliminate cross-block correlation caused by identical taxa.

We therefore (a) average an accepted name's heldout score across its island occurrences **within each heldout fold**, (b) average equally across folds in which that name is tested, and (c) count each globally unique taxon name once. Confidence intervals now resample entire **genera**, rather than island records; this is a stronger taxonomic sensitivity, **not** fully independent phylogeny or family-level sampling.

| Name-balanced local-source gain (genus bootstrap 95% interval) | All (16,386 names / 2,262 genera) | Direct (8,244 / 2,214) |
|---|---:|---:|
| K-score | **+0.01354** [+0.01027,+0.01716] | **+0.00905** [+0.00696,+0.01109] |
| Brier | **+0.00649** [+0.00097,+0.01190] | +0.00095 [−0.00060,+0.00254] |
| Exact-profile log loss | +0.01872 [+0.00443,+0.03390] | −0.00388 [−0.01423,+0.00478] |

The K-score gain persists even when common species no longer receive greater weight merely because they occupy more islands. **The log-loss verdict is weighting-sensitive**, so no all-scale exact-category predictive improvement is claimed. This is an in-sample-to-source-admissible heldout taxon-name prediction benchmark, not a reconstruction of evolutionary trait gains.

[Machine-readable taxon-equal, genus-bootstrap summary](taxon_name_balanced_genus_bootstrap.csv).

## Interpretation and the next decisive falsification

**Supported observational prediction result:** candidate local mainland genus source lists carry positive **similarity-aware** information over and above global genus traits in geographically blocked out-of-taxon-name tests. They do **not** robustly improve exact-category log loss across island- and taxon-weighted definitions.

**Still unverified:** that the *specific* geographically nearest assignment is uniquely informative. Sparse local donor coverage, source-list sampling and hierarchical smoothing could allow generic local subpool estimates to improve a non-strict K-score. The next test must **replace each island's five nearest mainland-source entities with geographically shuffled, same-regime source entity sets under the same training, shrinkage and candidate-pool size**, and compare the *proper K-score gain against this negative control*. This tests whether native mainland **geographic identity** is informative, not merely whether modeling local subpools changes soft-trait error.

Even a positive real-versus-shuffled source result would concern recorded flower categories and taxonomic source opportunity; it would not demonstrate evolutionary selection, successful colonization, real mainland ancestor sources, functional pollinator replacement or confirmed continental–oceanic geological filtering. No geological origin effect is estimated or declared here.

## Reproduction and rights

The code paths in this PR, `src/fd/n2_crossfit_macro.py` and `src/fd/proper_scoring.py`, implement full categorical scoring and same-fold name/genus reweighting. Real-source results were computed locally using the **rights-restricted original checkpoint inputs**; the public CI only executes synthetic unit tests and source-approved public aggregates. The published files contain grouped aggregates, score intervals and fold counts only—never accepted names, species/profile mappings, GIFT source region IDs, island IDs or private geographic crosswalks.

Input provenance: archived source-batch and progressive checkpoints specified by the frozen [island PR #321 reproducibility runner](https://github.com/zuizui0223/island/blob/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/scripts/run_floral_space_exploration.py). The 496-profile reference was hash-checked against the public source tuple and profile-frequency maps; training donor removal was validated per fold.
