# N2: private held-out genus × native-mainland source floral distance reference

**2026-11-09; implementation and synthetic CI only unless a separate
source-authorized real-data run is attached.** Do not claim that N2 has
produced real p/SES values from the public PR321 aggregates. Their
island × profile carrier table does not contain the island × taxon × genus
join needed for this conditional analysis.

## Question, units and why N1b is not sufficient

The N1b reference weights each distinct colour/form profile by the count
of *other regional islands* bearing it. On original D and uniform
30-profile rarefaction, All N1b mean nearest D median SES was +0.181
for continental and -0.524 for oceanic islands (130 each), while the
independent 10-pair geology test and old 5-component grouping did NOT
provide multiple-testing-supported differences. These are descriptive,
not adaptation.

N2 asks: **Given a test island's recorded genera and taxon carrier counts,
does its set of distinct recorded flower-colour/form profiles have smaller
internal nearest-neighbour D than a model trained on OTHER taxon names
from native GIFT mainland candidate floras would predict?** The known
island genus mixture is intentionally an *oracle conditional reference*,
not source-only colonization prediction. The source regions are
geographically plausible opportunities, not demonstrated ancestors.

## Donor leakage, data rights and full provenance

Run only locally with authorized original PR321 per-scope
`species_profiles.csv.gz` tables, private `island×accepted_species`
flora, archived GIFT 3.2 native mainland entity×species names,
island×mainland source assignments, and the existing pinned
public G0 and geological crosswalk. Do not add any private inputs to Git.

1. Reconstruct the same original outcome-blind 37 origin-matched pairs
   and the *post-hoc distinct* 10 `both_verified` geographic pairs
   using public covariates. Take the union of their physical islands,
   per All/Direct scope.
2. Join island name membership with its **single intact, observed
   colour/form state-set profile**. Validate *each island's per-profile
   taxon-carrier counts*, not just totals, against its pinned public
   G0/carrier aggregates; abort on any mismatch.
3. Collect ALL test taxon names across BOTH matched cohorts in each
   evidence scope. Remove each such name from EVERY native-mainland
   GIFT source entity BEFORE training, then join remaining native
   species to independently held-out flower-profile donors. This
   disallows exact name reuse across any test island, potential
   source feature or trait donor. Genus sharing is deliberately retained.
4. Use the existing `DonorModel` with prior pseudo-count 0.02/profile
   and hierarchical smoothing (default 8 taxon-pseudocounts).
   `N2_global_genus` models traits from globally available
   non-test native GIFT names in each observed genus.
   `N2_local_genus` updates them with those in the assigned mainland
   entity set (default `geo_k5`). For genus with no non-test donors,
   the model shrinks/falls back to a global prior rather than
   inventing a false biological absence. REPORT donor coverage.

Source donor names are originally mapped to traits by the PR321
island-observed profile map; this is out-of-tested-taxon, **not** a
separate primary mainland trait survey or phylogenetic independence
demonstration.

## N2 fixed-distinct-profile working reference

Let `c_ig` be the number of recorded focal-island taxon-name
carriers of observed genus g, and `p_ig(a)` the heldout donor
probability for complete categorical profile a under global/local
genus model. Define approximate inclusion propensity:

    u_ia = 1 - product_g [1 - p_ig(a)] ** c_ig

This is the exact *marginal presence probability* of a category
under independent donor-assigned per-name profiles, conditional on
the genus/name counts. The sum of all u is that **generative model's
expected distinct-profile richness**, which is reported and compared
with observed n.

The script then samples **exactly n distinct profile types**
without replacement using probabilities proportional to u; n is
the island's own observed distinct-profile count. This additional
weighted sampling is ONLY AN APPROXIMATION TO the exact joint
distribution of occupancy sets conditioned on n: inclusion
indicators are correlated, and sequential without-replacement
sampling is not necessarily the same conditional distribution.
Thus N2 is **a predictive fixed-support benchmark, NOT an exact
genus-preserving randomization test**. This distinction is
essential even if its SES is striking.

For each of 999 (default) simulated n-profile supports, calculate
`mean off-diagonal D` and `mean nearest-other-profile D` on
the ORIGINAL 496×496 colour/form state-set Jaccard matrix.
Report `SES=(T_obs-mean(T_null))/sd(T_null)` with SD ddof=1,
and observed-minus-expected D, alongside N1b on the same target
islands (sample size fixed and focal island excluded). The separate
equal-30 sensitivity averages 99 uniform thinnings of observed
island support and uses n=30 in all three null references. All/Direct
are nested, not independent biological replications.

## Predictive calibration is a separate gate, not a discovery

Before interpreting a local-genus N2 SES, score each tested island's
**held-out name-level** actual recorded profile under three models:
the non-genus `global_prior`, `global_genus` and `local_genus`.
Use exact categorical **negative log probability** as the proper
primary predictive score, and the expected original state-set
Jaccard distance to each held-out profile as an ecological
distance comparator. Calculate per-island mean losses, then
aggregate *equally by island* within All/Direct and geological
classes. No individual target name, genus or island is released.
The output `n2_heldout_prediction_calibration.csv` reports
local-minus-global-genus improvement and the fraction of
tested islands on which local genus is better. Positive local
increment is a necessary adequacy check, not a guarantee of
correct source ancestry, trait predictions or selection.
Evaluate the fixed shrink=8 model first and disclose sensitivity
for shrink=2/20; do not tune shrinkage on the same held-out
outcomes and subsequently present the optimized loss as an
independent validation.

## A second generative check against artificially forcing observed n

The primary fixed-n N2 uses predicted inclusion propensities with
weighted sampling **without replacement**. This is NOT equivalent to
the exact joint genus-by-name model conditional on observed profile richness.

Therefore a separate **N2g** model check samples one complete profile
for EACH observed genus-name carrier from its held-out local-genus donor
distribution, independently within the model, without forcing the number
of distinct types to match observation. Across 199 simulations (default),
report the generated distinct-richness median and 95% interval, the
fraction of focal islands whose observed richness falls outside this
interval, and the fraction of simulated communities with >=30 distinct
types. Only simulations meeting that last criterion are rarefied uniformly
to 30 types and compared with the observed equal30 metric as an additional
**descriptive** distance check. If fewer than 20 such simulations are
usable, no generative distance SES is reported for the island.

The N2g model check must not enter the four geographical p tests.
If N2g and fixed-n N2 disagree, or richness support is very poor,
interpret fixed-n N2 as a misleading reference rather than cherry-pick
the smaller p-value. N2g remains a simulated generative model, not proof
of independent taxa, source ancestry or ecological filtering.

## Inference and sampling safeguards

There is **one fixed post-hoc four-slot hypothesis family**, only
`N2_local_genus` / equal30 / nearest D × All and Direct ×
(1) original 37 pairs with all five spatial connected components
sign-flipped together and (2) 10 distinct `both_verified` pairs
with exact pair sign flips. Exact two-sided p and BH/BY q
are conditional on these observational geographic pairing assumptions.
If <8 matched pairs have computable SES, retain the test slot
as `insufficient_complete_matched_pairs` (p=1 in the four-slot
family), not as a result. All other N2 comparisons are
**descriptive only**; do not add unreported post-hoc significance tests.

Out-of-target taxon holdout and source conditioning do not eliminate
shared genus phylogeny, source-area misclassification, undocumented
taxon detection/missing flower trait evidence, hybrid assignments,
native/introduced uncertainty or true missing species.

A meaningful interpretation gate needs ALL of the following:
- Real private inputs audited and source/run digest matched, no
  hidden taxon-name leakage
- Coverage reported separately: proportion of focal name carriers
  with independent global/local genus donor support (previous
  mainland analyses found local donor support often sparse)
- Compare observed distinct-profile n with the generative model's
  expected distinct-profile count to detect major model mismatch
- Local-genus baseline substantially improves prediction over
  global-genus or global reference in external taxon holdout;
  otherwise apparent negative SES may only diagnose misspecification
- Geography-replicated comparisons and explicit detection/trait
  coverage cautions, no claim of evolution/selection unless
  supported by independent causal evidence.

## Local authorized replay (never public GitHub Actions)

The default output is ignored by `.gitignore`
(`results/run_*/`); use only aggregate outputs in a separately
reviewed public report.

```bash
python -m fd.n2_source_distance \
  --private-all /secure/island/floral_space_20261106/all/species_profiles.csv.gz \
  --private-direct /secure/island/floral_space_20261106/direct/species_profiles.csv.gz \
  --flora /secure/island/private_flora.csv.gz \
  --mainland-native /secure/gift/native_mainland.csv.gz \
  --source-assignments /secure/island/mainland_assignments.csv \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --g0 work/fd_g0 \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --input work/island_public \
  --output results/run_private_n2 \
  --source-mode geo_k5 --shrink 8 --draws 999 --thinnings 99 --generative-draws 199
```

The path basenames above are illustrative; they are not claims that
the source archives have those exact names. File readers require
the documented column schema. Also run explicit shrinkage
sensitivities (e.g. 2 and 20) to check robustness;
they are not independent hypothesis tests.

Aggregate outputs only:
- `n2_heldout_prediction_calibration.csv`: held-out per-name categorical log loss and expected original D, island-averaged by group, plus local-vs-global increment
- `n2_distance_ses_summary.csv`: per-scope/origin/broad-region
  observed, modeled, delta, SES; suppress all scores for n<5
- `n2_source_donor_coverage.csv`: donor support and expected vs
  recorded distinct richness by origin/scope
- `n2_paired_spatial_sensitivity.csv`: four exact coarse or
  independently geography-disjoint exploratory test slots
- `manifest.json`: source hashes, parameters, run status and
  interpretation limitations, **no taxon/island/source identifiers**

**Stop rule:** Without authorized private checkpoint + GIFT source
mapping, release NO empirical N2 SES or p values. Synthetic tests
only establish the code's basic behavior, not data adequacy.
