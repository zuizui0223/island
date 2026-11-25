# GloBI five-channel documentation composition and direct floral morphology

**11 Nov 2026 — exploratory, source-pinned fd PR #14. This analyzes GloBI DOCUMENTED visitor-channel reference shares, not observed focal-island pollinator assemblages, local effective pollination, or evolutionary adaptation.**

## Why a new analysis?

Earlier fd analyses found that Direct floral colour×architecture pairing is common across island floras, and that globally documented GloBI effective visitor-channel *breadth* covaries with island floral coordination. But adding morphology did NOT improve prediction of breadth in an entirely unseen genus (486-genus heldout, RMSE 0.705 source-evidence only versus 0.729 with marginal traits, 0.784 with joint profiles).

Breadth and the **relative composition of GloBI-documented visitor channels** are different ecological *recording* endpoints. We now test the latter, directly at one genus per prediction test unit.

## Frozen original inputs and population

- Rights-protected `island` species-level evidence checkpoint, GitHub artifact `10058653212`, ZIP SHA-256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- Outcome-independent GloBI genus-breadth v2 source artifact `10328696801`, ZIP SHA-256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`. The GloBI source original was compiled before the fd morphology relation; plant-side morphological data did not produce GloBI interaction identities.
- Original trait ledger path `fixed/source/current/evidence/direct_only_trait_ledger.csv.gz`. Only `evidence_scope=species_direct`, High/Medium, and **explicit single-state** flower primary colour, floral form, and floral symmetry records are accepted. No `genus_consensus` labels, mixed states or unmeasured plant traits are imputed.
- From **2,234 strict complete Direct flower-morphology species** representing 893 botanical genera, **1,597 species / 486 genera** have a qualifying GloBI source profile with >=3 independent flower-interaction references. The parent GloBI genus catalog has 2,228 genera (1,534 pass >=3 references across the full catalog). This matched 486-genus subset is **not all islands, all plants, or an independent sample of plant clades**.

Each genus' GloBI outcome is the proportion of its reported **five-channel independent reference-channel units** allocated to `bombus`, `non_bombus_bees`, `lepidoptera` (unseparated butterflies and moths), `flower_visiting_birds`, and `diptera`. The five documented fractions sum to 1, conditional on the original reference window. A **zero count means no documented link in the compiled source**, never true pollinator absence; references are not abundance-weighted visits or pollen delivery.

## Exact predictive comparison

All models use identical 486 GloBI genera and a shared fivefold *whole-genus* heldout evaluation for each of three fixed seeded fold assignments. Source-evidence-only base predictors are `log1p(number GloBI references)`, `log1p(number matched GloBI plant species)` and `log1p(number directly measured complete-flower species per genus)`. Two candidate floral designs are:

- **Simple morphology (3 additional predictors):** fraction actinomorphic, fraction open-radial floral form, and fraction green/brown inconspicuous colour, among complete Direct species in a genus.
- **All trait marginals (22 additional state proportions):** the seven colour categories, twelve floral-form categories and three symmetry categories.
- **Marginals plus full joint profiles:** add all observed colour×form×symmetry combination proportions (115 total feature columns), a deliberately high-dimensional stress test, NOT a prevalidated pollination-syndrome classifier.

Multivariate ridge alpha is fixed to 10 or 100 in this posthoc sensitivity; the split is by genus, and all feature standardization is learned **inside training folds**. Out-of-fold predicted five-channel fractions are Euclidean-projected to the 5-simplex: nonnegative and summing to one. The primary reported point values use alpha=100, seed=20261111. MSE is the five-channel mean squared documentation-share error (lower is better).

| Model | Out-of-genus five-channel MSE | Gain versus evidence-only |
| --- | ---: | ---: |
| Evidence-only | 0.041396 | baseline |
| + 3 simple floral morphology rates | **0.039785** | **+0.001611 / 3.89% lower** |
| + 22 floral category marginal rates | **0.038272** | **+0.003124 / 7.55% lower** |
| + complete joint morphology rates | 0.043707 | **−0.002311 / 5.58% worse** |

This ranking survives both fixed ridge strengths and all three random genus folds. The source-only coefficient and the morphology source are distinct. **The earlier breadth-prediction failure remains fully valid**: the target has changed from number/diversity of documented channels to their documented relative frequencies.

With the genus pointwise OOF squared-error difference retained locally (never published), a conditional 4,999-genus bootstrap yields:
- simple morphology gain **[+0.000936,+0.002288]**;
- all marginal traits gain **[+0.001517,+0.004625]**;
- full joint trait model gain **[−0.004599,−0.000045]**.

These are **not family- or provider-disjoint intervals**; genera may share botanical family, source publications and survey preferences. They do not prove transfer beyond family or a causal syndrome.

### The five channels are not equally precise

With 3 simple traits, alpha100/first split, per-channel MSE improvements are Bombus **+0.001491**, other bees **+0.001403**, Lepidoptera **+0.001599**, flower-visiting birds **+0.000974**, Diptera **+0.002589**. The marginal genus-bootstrap interval crosses zero for **other bees and Lepidoptera**, so do not promote five individually proven guild-specific effects. The five proportion components are also compositional and statistically dependent. For a pooled ecological interpretation, the **equal-weight five-channel MSE** is a better primary reporting unit than choosing a single favorable insect guild.

## Coarse source-effort-matched randomization

To ask whether any three-trait forecast advantage is simply a GloBI source-documentation artifact, permute each genus' *whole triplet of morphology features* within 5 bins of log independent GloBI reference counts × 3 bins of GloBI represented plant-species support × 3 bins of Direct complete-morphology species count. Forty observed strata, 34 with >=2 genus members; **480/486 genera** can be shuffled. Refit the same 5fold alpha100 model to each of **999** independently shuffled morphology labels. Outcome GloBI channel shares, source evidence features, and 486 test genus folds stay fixed.

- Observed pooled gain **+0.001611**.
- Permuted pooled gain median **−0.0001945**, middle 95% approximately **[−0.0005026,+0.0001573]**.
- Posthoc one-sided tail fraction **0.001** (minimum possible with 999 shuffles and +1 correction).
- All five individual documented-channel gains exceed their respective shuffled 97.5th percentiles; individual tail fractions 0.001–0.004, **not independent confirmatory guild tests**.

This falsifies the *narrow* rival that three-trait prediction benefits arise solely from broad GloBI reference effort and number of studied plant taxa. It does NOT condition on botanical family or phylogeny, island native history, correlated trait reporting sources, local pollinator sampling, or actual effective pollen transfer. The morph predictor family and randomization are posthoc discoveries.

## Scientific interpretation and failures

**Supported as a predictive source-documentation result:** coarse Direct floral morphology carries a **small, transferable signal for the relative GloBI-documented visitor-channel references of an unseen botanical genus**, after simple measures of reference/documentation opportunity are supplied.

**Not demonstrated:** (1) the composition of real pollinator assemblages on any of the 928 focal islands; (2) selection toward actual bee-, moth-, bird- or fly-pollination syndromes; (3) successful pollen deposition/seed production; (4) new evolutionary convergence, family/phylogeny-independent syndromes, or nocturnal moth-vs-diurnal butterfly channel identity. There are no independent local visitor effectiveness records in this GloBI genus reference matrix. All large-bee/butterfly/bird source templates remain heavily overlapping and do not independently assign realized pollinators.

Do not present the successful three-feature model as evidence that full colour×form×symmetry *trait coupling* determines pollinators: adding high-dimensional joint features **worsens** heldout performance. Conversely, this penalty does not biologically reject all syndrome hypotheses: 115 predictors for 486 observed genera are susceptible to limited-data variance and could be improved by different preregistered modeling/covariates.

**Next identification gates:** family-/phylogenetic-clade-heldout predictions, source/provider-stratified external visitor records, and matched island-specific effort-normalized local visitor assemblages. A legitimate causal claim of functional flower–pollinator fit further requires visitation effectiveness or pollen deposition. The prior exact-island five-channel dataset has just 33/928 mixed retained/disrupted islands and zero of 20 region×guild cells passing the stringent two-class support gate; it cannot establish real local pollinator loss.

## Implementation and release limits

Reusable code and synthetic provenance/holdout tests:
- `src/fd/globi_channel_composition_cv.py`;
- `tests/test_globi_channel_composition_cv.py`.

Nonidentifying public aggregates:
- `fivechannel_documented_prediction_aggregate.csv` (four-model main comparison);
- `fivechannel_effort_stratified_label_null.csv` (six aggregated posthoc shuffled test rows).

Local, **identity-suppressed** output: all 3 seeds × 2 ridge penalties × 4 models × 5 channels, including pooled five-channel scores, and stricter 999-reassignment sensitivity. Private genus names, plant identities and citations are not included in the public files or GitHub CI. Source SHA256s were checked against prior reported originals, and the new synthetic tests do not independently replay protected empirical inputs.
