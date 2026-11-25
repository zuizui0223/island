# fd — category-appropriate functional **metric volume**, original macro study (2026-11-13)

**Research question:** How much *realized original floral categorical state space* do island floras cover relative to geographically plausible mainland source floras, **beyond differences in the recorded number of flower profiles**?

**Status:** Original-source result is separately executed in `island` protected GitHub CI; read the run artifact and the source validation receipt before making numerical biological claims. This is a distinct data/estimand from the 17D Gaussian smoothed KDE and from G1 fixed-margin island/profile-incidence tests.

## Why the previous analysis was not a full functional hypervolume

- The 2026-11-08 2D PCoA 90%-HDR Gaussian KDE area retains only **15.35%** of positive reference inertia and is not a full 496-profile space volume.
- The 17D Gaussian KDE uses 17 of the positive PCoA axes (80.49% inertia), two bandwidths **h=0.55, 0.85**, two repeated 12,000-draw MC HDR volume integrals. Both are **legitimate numerical calculations in an artificial Euclidean embedding** of unordered original colour/form sets, but continuous smoothing adds impossible phenotype combinations; volume is near an increasing multiple of distinct recorded profiles when most kernel centers barely overlap. Pairwise profile-richness-adjusted 37-pair origin contrasts all include zero, even before the stricter five-geographic-component dependence correction.
- G1 fixed binary island×profile margins actually tests **equal-profile Rao Q and between-island nearest-neighbor category distance**, not KDE volume itself. No G1 BY-adjusted comparisons survived the original 32-slot family.

**Avoid making a positive 17D volume statement.** Its absolute numbers change multiplicatively with bandwidth and its observed origin effect is essentially a recorded-profile-number effect.

## A defensible finite category-space alternative

Freeze the original **Direct species-high/medium, single-state, unequivocal** trio of flower-primary-colour, floral form, floral symmetry. The original global Direct ledger has **2,234** three-trait complete plant species. They form a finite universe of **actually observed** three-category phenotype profiles `R`. The 928 eligible island floras each have ≥30 such species. This reference is fixed globally before reviewing any island results; it is not a set of hypothetical Cartesian combinations.

Define the original categorical Hamming dissimilarity as
[
d(x,y) = ig[{f 1}(C_x\ne C_y)+{f 1}(F_x\ne F_y)+{f 1}(S_x\ne S_y)ig]/3.
]

Let `A_i` be one equally species-rarefied set of 30 original source-trait-complete island records. The **finite support hypervolume** is the measured uniform reference mass of the union of radius-`r` categorical metric neighborhoods around those ACTUAL profiles:
[
V_r(A_i)= {1\over |R|}\sum_{x\in R}
 \mathbf 1\{\min_{a\in A_i}d(x,a)\le r\}.
]

- **r=0** equals the number of exact occupied original reference profiles divided by |R|. It is **profile richness**, not an independent geometric result; report only as a baseline.
- **r=1/3** covers all originally observed source profiles reachable by **at most one** original categorical trait change. It can differ among sets with identical original exact-profile richness; this is the primary *metric support volume*, not a continuous density estimate. No unrecorded Cartesian combination is reported as an observed source profile.
- **r=2/3** is solely a saturation diagnostic. A regional population with ≥80% of islands near the 98% reference coverage ceiling fails the meaningful neighborhood-resolution gate; do not promote volume differences in a saturated metric.
- One global **species-frequency-weighted** `V_1/3` is secondary; original global Direct species frequencies, not each island's original n, determine this fixed reference measure.

This is an operational **reference-measure metric-ball volume in a categorical space**, not *convex-hull Euclidean FRic*, not a 17D Lebesgue volume and not effective pollination functional space. Colour and form/symmetry are partially ontologically dependent; changing the metric is a legitimate sensitivity, not evidence that the raw equal Hamming weights represent the strength of pollinator matching.

## Three prespecified controls that isolate the mathematical artifact

1. **Exactly 30 directly characterized species per island/source** (25 stochastic thinnings; 19 frozen colour permutations each; seed 20261113). Original GIFT 5-nearest mainland source flora are separately thinned to 30 Direct species. They are **candidate native source floras, not demonstrated ancestors**.
2. **Same distinct-profile count k**: For each rarefaction, compute the difference between observed `V_1/3` and mean `V_1/3` of **199 null subsets with the same k distinct reference profiles**, selected without replacement from the global original reference. This explicitly tests whether metric-ball coverage adds anything beyond counting observed original flower profiles. This geometry null is **not** a genus-conditional, geography-conditional or trait-margin-conditional biological null.
3. **Exact original trait marginals**: Shuffle only flower-colour labels over the 30 sampled Direct species while keeping full actual (floral form, floral symmetry) *pairs intact*. The difference between original and mean permuted `V_1/3` reflects colour-architecture state-pairing contribution to one-trait coverage. Null centers can be counterfactual trait combinations; only the ORIGINAL GLOBAL REFERENCE profiles carry mass.

For any source-relative comparison, report the island−source metric volume, the island−source profile-richness difference, and the difference in **k-conditioned geometric residual**. Also evaluate an exclusion of **every taxon species that is shared between the observed island and its five-nearest source flora**; if no ≥30 species remain in the reference, mark unassessable rather than relaxing support.

## Biological support and nonclaims

- Analyze *all-recorded* 928 Direct islands and *strict-native-recorded* subset independently; then perform a within-exact-island all-versus-native diagnostic, because the two cohorts are nested and source-native missingness selects special islands.
- Regional grouping: northern midlat, northern high, tropics, southern extratropics. Report equal weighting of independent 10° **spatial blocks**, and separately equal weighting of the **first assigned GIFT mainland-source signature**. Conditional 1,999-level bootstrap intervals do not reconstruct independent colonization or phylogenetic histories.
- The same species has one immutable label in this original database at every location. Thus any island–source metric result measures **documented floristic sorting and representation**, not an evolved within-species floral phenotype, ecological plasticity, actual pollinator depletion or guild switching.
- Genus-conditioned and plant-source-status-conditioned nulls, spatial range dependence, under-recorded islands, morphological observation error and real pollination service all remain unresolved.

## Stop/go decisions

- If `V_1/3` saturates in many islands, the metric **cannot detect additional ecological composition change**. Do not tune the radius post-hoc to force significance: use original exact trait-pairing MI, three-trait joint Rao Q, and taxonomically informed sources instead, or obtain numerical flower geometry and colour/nectar/UV traits to define a genuine continuous morphology volume.
- If island−source volume is directional but the same k-conditioned geometry residual is not, **the apparent volume difference is primarily distinct original profile richness/ascertainment**, not extra functional space spread.
- If native-status and source-excluding-shared-plant controls disagree with the broad all-recorded result, report the selection dependence; do not conclude universal island syndrome or pollinator-induced floral volume reduction.
- Only a source/genus/phylogenetic and observation-effort controlled geometric separation across independent geographic blocks can advance this macro ecology beyond descriptive floristic association.

**Audit source code:** `src/fd/discrete_functional_volume.py` and original source-protected `island/scripts/run_fd_discrete_functional_volume.py`. No original taxon×island/source identity or map coordinates may be made public.
