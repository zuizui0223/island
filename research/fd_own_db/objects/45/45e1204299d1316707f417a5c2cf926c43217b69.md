# Spatial independence, strict native status, and genus-predictive controls for fd functional traits

**10 November 2026 | empirical exploratory audit | PR #13.** This replaces any strong claim that repeated rarefaction-positive island-pair statistics alone establish source-independent functional diversification. All figures are restricted-source-derived, pooled aggregates; no original species, island, genus, geographic-block IDs or private source rows are released.

## Identifiable question

How much of the observed turnover in binary floral/reproductive trait prevalence can be attributed to (i) losing shared species, (ii) losing shared genera and changing observed carrier number, and (iii) predicting the traits of the remaining species from other species **in the same genus**?

The observed island-flora data are NOT a census, and the traits are **proxies** rather than experimentally measured pollinator performance.

## Data and independence

- Rights-authorized `island` original frozen progressive species/trait snapshot; independent public fd G0 cohort: **7,618 sampled pairs, 1,419 distinct islands**. All/Direct are nested evidence, not replicates.
- Four broad geographical regions. 999 independent *within-pair carrier draws* per comparison are a **sampling-support null only**, not independent biogeographic replicates.
- The spatial-block graph of the sampled pairs is completely connected **within each region**, so the 7,618 edges cannot be reinterpreted as independent geographic histories.
- 16 high-northern, 66 middle-northern, 22 southern extratropical and 65 tropical observed blocks under this cohort. At most 8/33/11/32 **mutually block-disjoint island pairs** can be selected within these regional pair graphs. Maximum cardinality is a structural bound; 200 outcome-blind random GREEDY MAXIMAL, not necessarily maximum-cardinality, block-disjoint re-selections were used as diagnostic perturbations. They are correlated with one another, so positive fractions are NOT p-values or confidence levels.
- Leave-one-spatial-block-out drops all pairs using that block as **either endpoint**, and recomputes the descriptive residual mean.

## Stagewise genus exclusion: positive floral-symmetry direction but no independent biological null

Statistic: mean over pairs of `observed |trait-frequency difference| after excluding shared genera` **minus** the expectation under carrier-count-matched random thinning from the **species-exclusive** flora. A positive number means more divergence than random carrier loss alone predicts, without implying selection or realized functionality.

| Region | All symmetry extra gap | Direct symmetry extra gap | All/Direct eligible pairs | Max disjoint block pairs |
| --- | ---: | ---: | ---: | ---: |
| Northern high latitude | +0.02018 | +0.02138 | 1,148 / 1,069 | 8 |
| Northern midlatitude | +0.01167 | +0.01105 | 1,807 / 1,802 | 33 |
| Southern extratropical | +0.02362 | +0.02396 | 1,530 / 1,509 | 11 |
| Tropics | +0.01856 | +0.01697 | 1,967 / 1,955 | 32 |

Random block-disjoint re-selections were positive in All 94%/75%/87.5%/100% and Direct 83.5%/80%/93.5%/99%, in the regional table order. A substantial portion of these group means survives single-block deletion, but this does not remove shared botanical species/lineages across blocks or infer independent histories.

### Sensitivity using ONLY source-verified island-native records

Strict source-backed native status retains 150,531 of the original 1,039,757 island×species records; unresolved origin is **not** promoted to native. Shared-genus exclusion is conservatively defined using **the full recorded flora** (including origin-unresolved taxa) before retaining the known natives. The result is severe eligibility loss:

- Tropical symmetry: **All 219 eligible pairs / 82 islands / 34 spatial blocks, extra gap +0.02534**; **Direct 200 / 81 / 33, +0.01790**.
- Nonoverlapping-block greedy reselections are positive only **56% All / 63.5% Direct** for this native-only sample.
- Other cells frequently fail the required 40 pairs, 40 islands and 10 spatial blocks reporting gates. Do not generalize the native sensitivity to all regions or treat its sample as random.

This means that the original broad positive tropical direction persists on the selected known-native subset, but its **geographically distinct reproducibility is weak**.

## Crucial lineage check: floral symmetry is overwhelmingly genus-conserved in recorded data

Strict single-state source traits, distinct accepted species names:
- **All**: 50,621 known binary symmetry names, 2,214 genera with at least two known species; **only 81 (3.66%) have both recorded symmetry states**.
- **Direct**: 32,039 known binary names; 2,200 multi-known-species genera, **81 (3.68%) variable**.
- By comparison, in All, within-genus variation was observed in 348/999 multi-known self-compatibility genera (34.8%), and 171/948 multi-known generalized-form genera (18.0%).

These are *recorded* trait states under strict coding, not a phylogenetic evolutionary-rate estimate: limited traits, correlated evidence sources and within-genus underreporting may reduce observed polymorphism.

### Strict species-name holdout

Assign each species name to one of five deterministic SHA-256 folds. Its trait state is excluded from **every** donor estimate when that fold is scored. The predictor uses the other folds' genus frequency shrunk toward their fold-excluded global frequency. This is out-of-*taxon-name* prediction, **not** phylogenetic-/source-archipelago-disjoint validation.

For recorded actinomorphic symmetry, using the small smoothing factor `tau=0.05` (chosen *post hoc* based on species-level proper score sensitivity):
- All global Brier **0.2202** → genus Brier **0.01113** (50,621 names; 96.3% have a same-genus training donor).
- Direct global Brier **0.2381** → genus Brier **0.01768** (32,039 names; 94.2% donor coverage).
- The extremely strong genus predictive skill is **conditional on known genus and recorded binary state**, not on successful island arrival, trait-associated fitness or pollinator effectiveness.

### Genus-based expectation of geographic trait distance

For each sampled genus-exclusive island pair, score ONLY taxon names whose strict binary trait was recorded. Use their out-of-name genus probabilities and an explicitly approximate independent-Bernoulli expectation for the **absolute difference between island trait prevalences**. Probability variance is handled with the normal absolute-moment formula. A species has one global train/test fold across all its island occurrences. No novel trait state or ecological absence is imputed.

| Region | All observed mean gap | All genus-predicted gap | Direct observed | Direct genus-predicted |
| --- | ---: | ---: | ---: | ---: |
| Northern high latitude | 0.12276 | 0.12266 | 0.12300 | 0.12548 |
| Northern midlatitude | 0.08630 | 0.08303 | 0.08516 | 0.08270 |
| Southern extratropical | 0.15266 | 0.13551 | 0.15983 | 0.14238 |
| Tropics | **0.09960** | **0.09566** | **0.10799** | **0.10306** |

In the tropics the observed-minus-predicted residual is only **+0.00394 All / +0.00493 Direct**, corresponding to 0.39/0.49 percentage points of trait prevalence, after providing genus information and sampling variance. The *numerical expected-to-observed mean ratio* is about **96.0% / 95.4%**, NOT a causal fraction of variation explained. Under 200 block-disjoint greedy re-selections, **68.5%/69%** are positive: the residual sign is not geographically stable enough for a claim of novel functional sorting. Southern-extratropical differences are larger but have only **11 maximum possible disjoint block pairs** and the same source/observation limitations.

### Numerical approximation cross-check

To test whether the *small* tropical genus-model residual is an artifact of the normal absolute-moment approximation, we drew an outcome-blind sample of 350 candidate island pairs per evidence scope and compared the normal formula with **1,500 direct independent-Bernoulli trait-state draws per eligible pair** (tau=0.05, the same fixed genus predictions). After the >=20 known-carrier-per-side gate, **341 All and 340 Direct pairs** remained. Mean `normal estimate - simulated expectation` was **-0.0000365 All / -0.0000146 Direct**; the 95th percentile of *absolute per-pair* discrepancy was **0.00141 / 0.00159**. This does not account for botanical shared-source correlation or model uncertainty, but indicates that the 0.00394/0.00493 **aggregate mean** tropical residual is not driven primarily by the numerical normal approximation.

**Smoothing sensitivity is crucial:** with tau 0.01 to 0.20, tropical residual is **+0.00363 to +0.00503 All / +0.00452 to +0.00637 Direct**. An earlier arbitrary stronger shrinkage tau=8 inflates the tropical residual to +0.02320/+0.02804. A genus model's apparent failure must not be interpreted as selection if its donor-strength parameter was chosen poorly.

## Main scientific interpretation (evidence-limited)

The macroscopic turnover in **recorded floral symmetry** is strongly linked to which genera are represented by each island flora. Once the model predicts held-out species trait states from genus, nearly all the observed tropical *mean absolute symmetry difference* is predictable, and the small residual is not robust to geographically disjoint pair reselection or strict-native support attrition.

This supports **lineage-associated assembly and strongly genus-conserved floral architecture**, not demonstrated morphological evolution or realized pollination functionality. Genuine fitness-neutral versus adaptive selection models remain unidentified. Self-compatibility is more variable among congeners and is a promising contrasting response axis, but its native-only and block-disjoint evidence is unstable.

**Next identification gate:** source-/genus-availability-conditioned whole-trait-vector permutations using independent floristic opportunities and spatial components or new source archipelagos; also validate the approximate absolute-gap expectation against direct Bernoulli simulations and prospectively held-out independent botanical sources. Do not present the present post-hoc smoothing choice or five folds as a confirmatory test.

## Computational/rights reproducibility

- New shared pure functions `src/fd/functional_genus_reference.py` with tests: SHA-heldout genus probability estimation, proper Brier/log scores and expected absolute gap approximation.
- Existing `src/fd/functional_rarefaction.py` and `src/fd/functional_support_null.py`: count-matched sequential attrition checks.
- Exact protected source checkpoint and all identity-level matrices remain outside git. This public report stores aggregate denominators and scientific limitations only. GitHub CI runs synthetic checks and public-source aggregate verification; **it does NOT rerun private species/trait results**, and therefore CI success alone does not independently verify the private replay.
