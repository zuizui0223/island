# fd functional-trait turnover: empirical count-matched diagnostic

**10 November 2026. Status: exploratory OBSERVED association, not a geological-origin, phylogenetic, or pollination-effect hypothesis test.** Source-level species, trait and island records were processed only in a private container. This note is restricted to aggregate results with >=40 pair comparisons, >=40 distinct islands and >=10 spatial blocks per reported analysis cell. It contains no source taxon or island identifiers.

## Question

Does removing shared plant taxa or shared genera change the **recorded prevalence** of binary reproductive/access traits between island assemblages by more than one would expect simply because each island then has fewer known trait carriers?

This tests an elementary measurement/sampling rival explanation for functional differentiation, **not** a biological functional-retention null. Neither an identical trait state nor a difference in its prevalence establishes equivalence or differentiation of actual pollination service.

## Source and population

- Rights-authorized original `island` frozen checkpoint: ZIP SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999` (private input not redistributed).
- Validated fd public aggregate / pair artifact: ZIP SHA256 `2efa7d2bb40adc60cc7d22b3c0c1043eb63c7c608e492fb0a422a14b3a563801`.
- **7,618 original fd sampled island-pair comparisons**, sharing **1,419 distinct islands**: 2,000 northern-midlatitude, 2,000 northern-high-latitude, 2,000 tropical and 1,618 southern-extratropical pairs. These pairs overlap and are NOT 7,618 independent evolutionary or spatial units.
- Source trait scope: **All** plus nested **Direct**; evaluated separately, not independent replications.
- Seven binary trait indicators from strictly **single-valued** source records: self compatibility (SC vs SI), predominantly/obligate selfing vs predominantly outcrossing, autonomous/delayed vs absent, generalized/open floral form vs prelisted restricted form, actinomorphy vs zygomorphy, open/shallow vs deep flower tube, and broad plain colour vs other stated colours. Mixed or multistate records and all other categories become missing, NOT inferred zero. Therefore this strict diagnostic is **not identical** to Chapter 1's seven H1 model outcomes. Each pair/trait requires >=20 valid carriers **per island** after each exclusion.

The original eligible group is compared sequentially to species-exclusive floras, then species-exclusive to genus-exclusive floras; the two steps are distinct. For each pair/trait, a deterministic-seed 999-draw hypergeometric reference independently thins each baseline island to the **same number of known trait carriers** as the retained subset. This removes the elementary carrier-number artifact, not the correlation of genus membership and traits. It does not reconstruct historical mainland donor pools.

## Real data: incremental genus exclusion beyond shared species (primary descriptive comparison)

Entries are the **mean of per-pair observed-minus-equal-carrier-reference absolute trait prevalence gaps**, in units of proportion (positive means the genus-exclusive assemblages are MORE different than carrier-number matching alone predicts). This is NOT a p-value, causal effect, or difference of the two marginal medians.

| Region | Trait | All excess | Direct excess | All/Direct eligible pairs |
| --- | --- | ---: | ---: | ---: |
| Northern midlatitude | Self compatibility | +0.0027 | +0.0002 | 1,769 / 1,720 |
| Northern high latitude | Self compatibility | -0.0014 | -0.0012 | 1,037 / 847 |
| Tropics | Self compatibility | +0.0187 | +0.0088 | 1,736 / 1,523 |
| Southern extratropical | Self compatibility | +0.0063 | -0.0024 | 1,448 / 1,310 |
| Northern midlatitude | Generalized form | +0.0203 | +0.0122 | 1,671 / 1,635 |
| Northern high latitude | Generalized form | -0.0021 | +0.0204 | 802 / 710 |
| Tropics | Generalized form | +0.0053 | +0.0029 | 1,304 / 1,091 |
| Southern extratropical | Generalized form | +0.0352 | +0.0217 | 1,331 / 1,285 |
| Northern midlatitude | Actinomorphic symmetry | +0.0117 | +0.0111 | 1,807 / 1,802 |
| Northern high latitude | Actinomorphic symmetry | +0.0202 | +0.0214 | 1,148 / 1,069 |
| Tropics | Actinomorphic symmetry | +0.0186 | +0.0170 | 1,967 / 1,955 |
| Southern extratropical | Actinomorphic symmetry | +0.0236 | +0.0239 | 1,530 / 1,509 |

Across all four broad regions, the **sign of the actinomorphy excess** was positive in both nested source scopes. In contrast, self compatibility changed direction or largely vanished across regions and evidence scopes. This is a measured contrast in recorded trait-state composition, **not** proof of stronger adaptive selection on floral symmetry.

A separate all-to-genus-exclusion comparison without stagewise matching gives larger differences for several traits; its confounding by shared-species removal explains why the *sequential* diagnostic above is necessary. Plain colour, selfing mating system, autonomous selfing and tube-depth effects were also computed under the same rule but are not promoted as separate findings. Support for tube depth is particularly sparse (no northern-high-latitude genus-excluded pairs pass its >=20-per-island gate); no universal seven-trait composite is defined.

## What remains open

1. **Spatial and shared-species dependence:** Many sampled pairs contain the same islands, genera and species. The 10–66 sampled geographic blocks per regime are NOT a permutation design. No p or confidence interval based on treating the 7,618 pair comparisons as independent is valid.
2. **Source pool and phylogeny:** The thinning reference conditions on carrier counts, **not** actual source genus distributions, introduction status, or plant phylogeny. Genus removal is itself trait-associated; excess does not identify adaptation.
3. **Trait ascertainment:** Only explicitly single-state raw categories are retained; many mixed records disappear nonrandomly. A strict-native and within-archipelago, source-opportunity-matched comparison is not completed.
4. **Functional meaning:** Self compatibility and floral structure are **proxies**; they do not measure pollen receipt, reproductive output or visitor effectiveness. Upstream Chapter 1 GloPL associations are context, not independent fd confirmations.

## Conclusion at this stage

**Different taxonomic exclusions reveal heterogeneous recorded floral/reproductive trait turnover.** Floral symmetry shows a recurring descriptive excess beyond random loss of known carriers, whereas self compatibility's excess is regional and evidence-sensitive. These results justify a source/lineage-conditioned analysis **but do not yet establish functional differentiation or functional preservation**.

Code for the count-based, rights-protected aggregate reference: `src/fd/functional_rarefaction.py`, with tests `tests/test_functional_rarefaction.py`. The original private taxon-level checkpoint and intermediate pair-level sufficient statistics remain outside the public repository. The source ZIP replay was not verified end-to-end in this PR; no claim of completed public recomputation is made.
