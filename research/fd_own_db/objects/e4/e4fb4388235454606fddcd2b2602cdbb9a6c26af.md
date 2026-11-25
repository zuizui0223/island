# G1 functional-distance null and projected hypervolume — verified 2026-11-08

**Research home:** [zuizui0223/fd](https://github.com/zuizui0223/fd)  
**Original source:** [island PR #321](https://github.com/zuizui0223/island/pull/321), public aggregate source commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`.  
**Successful execution:** [GitHub Actions #42438861743](https://github.com/zuizui0223/fd/actions/runs/42438861743), [results artifact #13600563834](https://github.com/zuizui0223/fd/actions/runs/42438861743/artifacts/13600563834).  
**Validation:** 11 tests passed (2 warnings); all G0, G1 and projected KDE steps completed. No restricted species-to-trait ledger was imported or published.

## 1. What was actually computed

G0 measured original colour/form category-set distances for 496 profiles. G1 tested **two distinct equal-profile-weight statistics**, with 2 independent Curveball chains × 200 draws each, preserving island profile richness and each profile's prevalence over the selected island population. Scored subsets contain up to 128 islands and 128 matched pairs per regime/scope; the null incidence matrix includes **all islands in the matched-pair population**, not only the scored subset. Pair and local-island subsets are chosen by identifiers, without using trait distances. All/Direct and their single-state subsets are overlapping sensitivities.

All **32/32 G1 comparisons** passed mixing checks (Rhat 0.9975–1.0098; maximum absolute lag-1 autocorrelation 0.196). BH supports 13/32, including 9 pair-NN comparisons and 4 local-Q comparisons. **BY supports 0/32** (minimum BY q = 0.05889). The 400-draw two-sided p minimum is 0.00499 and contributes to limited BY resolution. Both BH and BY belong to one fixed 32-comparison family. The G1 results are **exploratory MCMC-approximate conditional tests**, not independent exact permutation tests or causal proofs.

## 2. Strongest positive directions (effect = observed − fixed-margin null median)

| Regional scope | G1 endpoint | Observed | Null median | Delta | BH q | BY q |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Northern midlatitude, All | Local equal-profile Rao Q | 0.77019 | 0.77375 | **−0.00356** | 0.0145 | 0.0589 |
| Northern midlatitude, Direct | Local equal-profile Rao Q | 0.76915 | 0.77247 | **−0.00332** | 0.0145 | 0.0589 |
| Northern midlatitude, All | Inter-island equal-profile nearest-distance | 0.11067 | 0.09955 | **+0.01112** | 0.0145 | 0.0589 |
| Northern midlatitude, Direct | Inter-island equal-profile nearest-distance | 0.11135 | 0.10017 | **+0.01118** | 0.0145 | 0.0589 |
| Northern high latitude, All | Inter-island nearest-distance | 0.08297 | 0.07223 | **+0.01074** | 0.0145 | 0.0589 |
| Northern high latitude, Direct | Inter-island nearest-distance | 0.08311 | 0.07496 | **+0.00815** | 0.0145 | 0.0589 |
| Northern high latitude, All single-state | Inter-island nearest-distance | 0.10232 | 0.07236 | **+0.02996** | 0.0145 | 0.0589 |
| Northern high latitude, Direct single-state | Inter-island nearest-distance | 0.10628 | 0.07548 | **+0.03080** | 0.0145 | 0.0589 |

Additional BH-supported results include local Q deficits in southern extratropical All/Direct and pair NN increases in northern-midlatitude single-state All/Direct and southern-extratropical All. Tropical pair NN effects are **not** BH-supported. No inference about a globally universal syndrome or tropical functional replacement follows from these comparisons. Full 32-row estimates, null reference ranges, diagnostic values and BH/BY adjustments are in the linked artifact's `results/run_fd_g1/g1_summary.csv` and all per-draw statistics in `g1_draws.csv`.

**Interpretation:** Northern midlatitudes show *both* slightly reduced within-island functional dispersion and increased between-island functional separation under the binary-margin null. The local Q difference is small (approximately 0.4–0.5% of the null magnitude), whereas observed nearest-distance differences are around 11% of the null. These effects are not independent, as pair IDs reuse islands and evidence scopes overlap. The northern high-latitude single-state cohort contains just 40 sampled islands and a restricted set of 171 eligible original pair slots; the current G1 samples 128.

## 3. Hypervolume measurements: **computed**, but limited to two-dimensional PCoA

All 4,722 eligible scope×island records received: (a) finite categorical occupancy of the 496 recorded reference profiles and (b) 90% highest-density **Gaussian KDE area in the common two-dimensional PCoA projection**, with standardized global reference axes and `h=0.35,0.55`, once equally weighting recorded profiles and once weighting name carriers. These are projected smoothing-dependent areas, not full biological floral functional hypervolumes.

| Region | All islands | Median occupied cells | Median projected KDE 90% area h=0.35 | Area h=0.55 |
| --- | ---: | ---: | ---: | ---: |
| Northern high latitude | 95 | 67 | 12.873 | 16.521 |
| Northern midlatitude | 1,111 | 80 | 13.503 | 17.071 |
| Southern extratropical | 70 | 50.5 | 14.011 | 17.808 |
| Tropical | 256 | 65.5 | 14.487 | 17.964 |

North-midlat h=0.35 area 13.503; tropical 14.487 is a **7.3% larger descriptive projected median**. There is **no conditional null, regional effect test, area-adjustment, or independence claim** for this regional difference. Increasing the bandwidth notably raises all projected areas. The smoothed density invents mathematical locations between originally discrete categories.

The PCoA 2-axis reference retains **15.35% of positive inertia**; reaching 80% requires **17 axes**. KDE areas therefore must be named **projected 2D KDE HDR area** in titles, figures and manuscript. They cannot be called total floral niche volume. Finite recorded-category occupancy is simply profile richness in the 496-profile reference and is not an independent hypervolume result. All generated values are preserved in `results/run_fd_hypervolume/*_hypervolumes.csv`.

## 4. What is not established and the next scientific gate

G1 preserves fixed binary island/profile margins **but not** genus composition, source pools, floristic status, taxon identities or shared-name exclusion. A positive effect here therefore does not overturn the distinct PR #321 genus-conditioned shared-taxon-excluded result (no BH/BY-supported island-pair difference). Pollination service, selection, historical trait evolution, species-independent convergence and functional replacement remain unresolved. The main next test is a **source/genus-conditioned null for the same distance endpoints** using rights-cleared identity-level data, plus dependence-aware and geography/coverage sensitivity checks. Increasing only MCMC iterations without revising the biological null would not identify the mechanism.

## 5. Provenance and checksums

- Source `source_manifest.json` SHA256: `a4e8a36eed3dff65d1d5ba11f71caabb3478e87be49f9cda17e2d4f6e834bf5e`
- G1 `g1_summary.csv` SHA256: `83ee82b2fce177485825cd5c007852c628d042d9283000e0db41053fc8644b1f`
- KDE `regional_hypervolume_summary.csv` SHA256: `0853144a4e1d0886bf67cf26ca2c151c37110b0b3dbd6fd054a28ec43d1191ce`
- Artifact contents: CSV outputs, 32-slot test table, null draws, cohort audit, fixed-parameter manifests and checksums; source-rights audit prohibits copying private taxon→trait rows.

Analysis definitions, caveats and replay: [G1/HV methodology](../../docs/G1_HYPERVOLUME.md).
