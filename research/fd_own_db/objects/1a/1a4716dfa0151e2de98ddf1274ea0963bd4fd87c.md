# Floral marginal trait diversity and syndrome-like pairing across island floras

**11 November 2026 | exploratory, Direct-only original-data audit | NOT independent pollinator-syndrome validation.**

## Central question

**Do island floras retain diversity in recorded floral trait states while exhibiting nonrandom associations among colour, form and symmetry?** These are TWO estimands. The first measures equal-weight categorical morphological disparity; the second measures *cross-trait combinations beyond observed marginal trait frequencies*. Neither automatically measures actual pollination function or named pollinator guilds.

The prior fd work measured colour/form 496-category profiles and genus-associated predictive transfer. The original equal-colour+form additive kernel cannot identify full trait pairing. This analysis explicitly tests pairing, using the original RIGHTS-RESTRICTED direct-species floral data, while keeping the biological inferential ceiling.

## Frozen inputs and observed support

- `island` progressive source archive: GitHub Actions artifact `10058653212`, ZIP SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Source tables used **locally only**: `fixed/source/current/evidence/direct_only_trait_ledger.csv.gz`, `fixed/canonical/input/chapter1_status_flora.csv.gz`, `fixed/source/isolation/results/purpose_shortest_island_data.csv`; optional historical template audit `fixed/source/current/syndrome/direct_only/species_syndrome_concordance.csv.gz`.
- Exactly one recorded accepted species x trait state; **no family/genus imputation**, no mixed-state values (`A|B`) assigned a single state. The selected three morphological traits are original flower primary colour, floral form and floral symmetry. This is a **flower morphology proxy space**, NOT realized ecological performance or the full Chapter 1 seven-trait reproductive system.
- **2,234 distinct species** have complete single-state Direct records for all three traits, among 62,360 species with at least one trait in the Direct source ledger.
- These species occur on 3,085 recorded islands. **928 islands** have at least 30 such recorded species, allowing the strict fixed-support analysis (743 northern midlatitude, 106 tropical, 53 northern high latitude and 26 southern extratropical). Because the same species recur on islands, 928 is NOT the number of phylogenetically independent plant radiations. Source-backed native-only records yield **87 islands** reaching the same criterion; in this selected status subset 67 are northern midlatitude, with sparse other regions.

## Two distinct estimands

**A. Marginal floral diversity.** On the same trait-complete species in each island, after drawing exactly 30 distinct species (five independent rarefaction samples per island):

`Q_morph = (1/3) sum_{trait j} (1 - sum_{state s} p_(j,s)^2)`

This is the categorical-equality-distance version of Rao's quadratic diversity, with equal trait weights and species occurrence weights, under the WITH-replacement Rao convention. It is **exactly invariant** to reshuffling which individual species carries which cross-trait combination if the three trait marginals stay fixed. It does not identify colour×form or colour×symmetry pairing.

**B. Recorded floral colour/architecture association.** Define `C`=colour, `F`=form, `S`=symmetry. Measure mutual information `I(C;[F,S])` in natural-log units, holding actual form/symmetry PAIRS intact. Independently shuffle colour labels among the 30 species **49 times per rarefaction** and subtract the mean chance-biased information. This null preserves all observed colour states and the full `(F,S)` joint frequency table, unlike an overly generous shuffle that also destroys ontologically overlapping form/symmetry. Average the five observed-minus-null differences per island. A positive score shows *recorded colour–architecture association beyond island marginal opportunity and finite-sample bias*, not causal pollinator selection.

A companion `I(F;S)` statistic scrambles symmetry but preserves its own margin. A third `I(C;F)` omits symmetry. These tests separate correspondence between colour and architecture from the partly definitional form–symmetry overlap.

**C. Identifiability of beyond-genus pairing.** Also shuffle colour only within **congeneric species present on the SAME island**. This keeps each island–genus state distribution fixed and removes only within-genus pairing. Such permutations have little scope when most genera are singletons or colour-monomorphic, so a zero residual does NOT prove that genus composition fully explains the association.

This is a **POST-HOC, DESCRIPTIVE** study, not a confirmatory randomization experiment. Across-island confidence intervals below are conditional whole-spatial-block bootstrap diagnostics. Spatial blocks still share species, genera, floristic sources and trait-observation providers. Replicating 49 random shuffles or five rarefactions does not create independent islands.

## Measured results: 30-species-standardized Direct flora

| Broad region | Eligible islands | Spatial blocks | Median marginal Rao Q | Median excess I(colour; form+symmetry), nats | Share positive colour/architecture association |
| --- | ---: | ---: | ---: | ---: | ---: |
| Northern high latitude | **53** | 7 | **0.507** | **+0.193** | 100.0% |
| Northern midlatitude | **743** | 53 | **0.553** | **+0.182** | 99.5% |
| Southern extratropical | 26 | 13 | 0.623 | +0.083 | 96.2% |
| Tropical | 106 | 42 | 0.606 | +0.066 | 86.8% |

Across the **928 islands**, 909 (98.0%) show a positive mean five-subsample colour×architecture excess; 874/928 show positive colour×form excess; 920/928 show positive form×symmetry excess (these fractions are descriptive and do NOT constitute independent evolutionary replicates). The full three-way total correlation also exceeds a complete-independence within-island marginal shuffle in 924/928 fixed-30 islands, but the source of this association could be taxonomic and categorical. The median Q is quite high even where colour/architecture association is low: these axes are distinct, not two names for the same KDE or profile-richness estimator.

**Equal spatial-block-weighted colour/architecture excess and conditional 95% block resampling interval:**
- Northern high +0.173 [0.132, 0.208], **only seven blocks**;
- Northern mid +0.159 [0.140, 0.176], 53 blocks;
- Southern extratropical +0.105 [0.082, 0.134], 26 islands / 13 blocks;
- Tropical +0.073 [0.059, 0.086], 42 blocks.

These conditional bootstrap intervals are **NOT** tests of distinct independently selected pollination syndromes, genus-independent trait convergence, or between-region causal adaptation. Trait-complete ascertainment is nonrandom, and area/climate/source-genus baselines are not matched. In particular, observed latitudinal medians cannot be promoted as formal between-region differences.

### Source-verified native-only sensitivity

Only **87 native-record** islands retain 30 complete Direct species; **83/87** show positive excess `I(C;[F,S])` under 30-species rarefaction. Northern midlatitude has 67 such islands across 16 spatial blocks, median colour/architecture excess +0.156. In the tropical native subset (17 islands/10 blocks), the positive excess is less stable with equal-block conditional interval spanning zero; southern extratropical has only three eligible islands. Do not claim four-region independent native-flora replication.

### Important failure of genus-conditional identification

The median fraction of complete-trait species in each island for which **colour can meaningfully be shuffled among congeneric cooccurring species** is only 7.8% in northern high, 10.8% in northern mid, 14.0% in southern extratropical, and 5.9% in the tropics. The median extra colour/architecture association relative to these *within-genus* permutations is approximately **zero** in three regions and −0.0038 in northern high. This is mostly a **low-permutation-opportunity failure** (the shuffles often make no change), NOT independent evidence that phylogeny or source sorting fully accounts for the association. A more relevant future test needs multiple biologically independent genera/colonizations and source-conditioned whole-trait-vector references; simply repeating current shuffles will not solve this.

## Why the existing 'guild-like' syndrome templates cannot be interpreted as pollinator identity

The original preregional `island` score tables label **large-bee-like, butterfly-like and bird-like** *morphological template concordance*, not measured pollinator identities. Restricting to Direct species with **>=3 informative source traits for EACH of the THREE templates**, only **205** comparable species remain. At the previous fixed high-concordance score threshold `soft_membership>=0.75`:
- **195/205 (95.1%)** have high concordance with **at least two** named templates;
- **168/205 (82.0%)** have high concordance with **all three**.

These candidate templates overlap heavily in their preferred architecture/color/size features and the source score is not an independently calibrated classifier. These numbers are an explicit **ambiguity audit** on a highly selected subset, NOT evidence that the same flowers are effectively pollinated by birds, bees and butterflies simultaneously. Older positive tropical/temperate template findings must not be turned into observed guild shifts. Reproductive assurance/selfing is a separate trait axis, not a pollinator identity template.

## Interpretation and manuscript decision

**Supported observational claim:** across hundreds of recorded island floras, the marginal diversity of directly recorded colour/form/symmetry traits does not describe the full arrangement of these traits into flowers. Colour/architecture co-association is frequently greater than island-level marginal-shuffle expectations, even when the form/symmetry relationship is left intact.

**Not established:** ecological functional diversity in pollination efficiency, pollinator syndrome classification, bird/bee/butterfly turnover, adaptive source-independent convergence, within-lineage evolution, or whether morphologically diverse island floras maintain equal pollen-transfer services. Direct trait completeness, species/genus overlap and correlated source evidence remain dominant limits.

A defensible synthesis is: **floral marginal diversity and multivariate trait coordination are distinct assemblage properties.** A causal island syndrome hypothesis remains a separate future test against independent observed pollinator interaction/effectiveness data. The pollination-syndrome literature supports considering multiple correlated traits *with measured pollinator function*, not calling a high template score a pollinator observation (Dellinger 2020, doi:10.1111/nph.16793; Rosas-Guerrero et al. 2014, doi:10.1111/ele.12224).

## Reproduction and rights

- `src/fd/floral_syndrome_coupling.py` is a private-input / aggregate-output implementation with strict species-direct and mixed-state guards.
- `tests/test_floral_syndrome_coupling.py` uses fully synthetic botanical profiles, validates the mathematical RaoQ-vs-joint-association distinction and confirms the grouped reports contain no island, species, genus or spatial-block IDs.
- The anonymous region-level CSV omits sparse native-only regional subsets. Original species × trait records, original island IDs, and precomputed syndrome species scores are **not** copied to the fd public repository.
- Empirical computations: five rarefactions to 30 species × 49 independent label shuffles per rarefaction for colour/architecture and within-genus checks; companion 99-shuffle unrarefied/genus check and 4,999 block-resampling diagnostic. Full private-source ZIP bytes and SHA remain under original source permissions.
- These are real-data exploratory checks but **not** independent external biological replication. Public GitHub CI runs synthetic logic tests only.
