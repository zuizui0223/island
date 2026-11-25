# Island lineage opportunity, predicted floral repertoires and 17D projected hypervolume

Date: 2026-11-08 JST. Independent `fd` research line. Source-checkpoint-coupled to [island PR #321](https://github.com/zuizui0223/island/pull/321) and its controlled data-rights audit.

## A. Distinguish three biological quantities

**1. Mainland source opportunity:** a lineage is listed as native in a potential mainland region, with recorded distance from the island. This does **not** mean it dispersed across the ocean.

**2. Recorded island genus:** at least one accepted taxon-name unit of that genus is listed from an island. **No record ≠ biological absence**, so a source-only ranking of recorded versus unrecorded genera cannot be interpreted as estimated successful colonization probability.

**3. Floral phenotype space:** the set of recorded flower-colour and form state sets of the sampled island flora. This is downstream of genus assembly and trait assignment, not a direct test of pollinator performance or evolutionary change.

The key predictive chain is therefore **possible source mainland → ranking of candidate genera → GIFT genus-specific floral profiles → observed island floral profile composition**. Evaluation may be done without supplying the **actual observed genus identity** to the source-only predictor, with a separate post-hoc known-genus upper benchmark.

## B. Source-only genus opportunity ranking (private source inputs; aggregate outputs)

Input: original [GIFT 3.2](https://gift.uni-goettingen.de/) native mainland lists (748 source entities; 1,666,186 entity×taxon rows), frozen outcome-blind `geo_k5` source assignments, original `island` taxon presence/status table, and the independently audited Roeble/Zizka geological-origin crosswalk. **Strict independent prediction MUST first remove every name observed on any of the 434 geology-classified islands from all native GIFT source lists**; the original non-held-out source overlap is only a descriptive sensitivity, not independent prediction. The strict source donor sample contains **923,786 entity×name records** (742,400 originally overlapping records removed). No island labels are fitted as rank-model coefficients.

For each island, the **candidate set** is the union of genera with ≥1 name in five assigned mainland native regions. A test genus is **recorded** if observed in the island source list, otherwise *unrecorded/unlabelled*. No unrecorded genus is declared absent. Genus outside the candidate union is handled in a separately reported **source coverage fraction**, not coded as a false negative.

Compare four fixed source-only rank scores within the **same candidate genus universe**:

- `global`: fraction of all 748 native GIFT mainland entities with the genus.
- `local_count`: number of the five nearby mainland entities listing that genus.
- `local_distance`: presence weighted by `1/sqrt(1+distance_km/500)`.
- `local_density`: weighted `(native name units of genus / all native name units in that source entity)`, controlling source-list richness.

Scores evaluated with tie-aware **average precision (AP)** and pairwise rank **ROC AUC** for recorded-versus-unrecorded candidate genera. AP must be interpreted relative to each island candidate prevalence; metrics are **not** calibrated arrival probabilities.

**Support**: island ≥10 recorded genera, ≥5 among potential source candidates and at least one unrecorded source genus. Under strict taxon-name holdout this yields **399 physically distinct islands** for all-recorded status (**174 continental / 225 oceanic**) and **301** for explicitly native status (**141 / 160**). The older overlapping-source exploratory values were 406 (174/232) and 306 (141/165), and must not be described as independent unseen-species prediction. The native status subset is not representative of all flora because most island occurrence statuses are unresolved.

Individual-island bootstrap intervals for source gains condition on the selected island population; **37 original matched geological origin pairs** (31 native-eligible pairs for the ranking score) are evaluated as a separate geographical comparison. The **12 paired tests** (three local models × AP and AUC × two status classes) use one BH/BY family. Even a zero-bias population with species detection probability depending on genus can generate an opportunity rank effect; detection still needs independent modelling.

## C. Two-stage source → genus → floral repertoires

Use the **same original matched 37 island pairs / 74 islands** as fd PR #4. All scoped trait input has 13,141 island×complete-name records /3,316 unique target name units, Direct 11,664/2,605. Train on the GIFT mainland flora with the original PR321 private **complete** taxon→colour/form profile map. Hold out every one of the **74 islands' target taxon names globally from BOTH the mainland candidate genus opportunity source features and every flower-trait donor record**, including every mainland source entity in which that taxon is listed. An initial pilot excluded targets only from donor traits; its larger apparent source benefit is not the strict primary result.

Create four model predictions of the **taxon-carrier-weighted flower-profile frequency vector** on each island:

- `global_prior`: native GIFT profile frequency across all non-target mainland donor names (0.02 pseudo-count per all 496 reference categories).
- `source_genus_prevalence`: mixture of non-target genus-specific profiles weighted by each candidate genus's **worldwide GIFT entity prevalence**, normalized **within the island's source-candidate genera**.
- `source_genus_density`: same genus-specific floral profiles, weighted by `local_density` from section B.
- `known_genus_oracle`: observed island genus carrier counts, but still **non-target mainland genus→trait profiles**. This is an *upper benchmark*, NOT a prospective source-only model.

Genus-specific distributions are shrunk to the global profile prior using eight prior donor pseudo-counts. If a source candidate genus has no independent floral donor, the source-only predictor falls back to the global prior; report coverage rather than silently treating missing donors as absence.

Evaluation: proper name-carrier profile **categorical log loss**, expected full original colour/form **set-Jaccard distance loss**, squared frequency-vector error and similarity-kernel frequency error (`K=1−D`, full original 496×496 dissimilarity). All/Direct are overlapping sensitivity definitions. Compute 4,999 cluster-bootstrap comparisons over 37 matched island-pairs; non-overlapping target names are required for training leakage control, but multiple test islands may share species. **Source-only prediction cannot estimate island genus richness** here; all scoring conditions on the evaluated island's recorded taxon count.

These models do not estimate true establishment, phylogenetic colonist ancestry, plant persistence, or pollinator-driven floral selection. They test whether a source-only **genus-weighted flower repertoire** improves prediction beyond the mainland-global profile frequencies.

## D. Two complementary hypervolume measurements

**Existing 2D projection:** globally standardized first two PCoA axes, 15.353% of positive reference inertia, Gaussian 90%-mass highest-density region, fixed bandwidth 0.35/0.55. Descriptive only.

**New 17D projection:** build the exact 496-profile `D = .5*(colour state-set Jaccard distance + floral-form state-set Jaccard distance)`. Embed `sqrt(D)` by positive classical PCoA; 17 axes retain **80.4899%** of positive inertia. Use **one common scalar standardization** across the 17 axis coordinates, not per-axis whitening. Two isotropic KDE bandwidths `h=0.55` and `h=0.85`. For each island, simulate 12,000 draws from the equal-recorded-profile-weight Gaussian mixture, estimate density `f(x)`, select its 90%-mass HDR via the empirical tenth density percentile, then estimate volume by `E_f[1_HDR(x)/f(x)]`. Repeat **twice independently** and report the relative Monte Carlo integration error and replicate log-volume discrepancy.

Run on **all independently geological-typed, ≥50 complete name islands** in public All/Direct scopes (260 All, 243 Direct), and the same 37 pairs for geographically matched contrasts. BH/BY four-family applies to 2 bandwidths × All/Direct raw 37-pair comparisons. Separately regress per-pair difference in log volume on **pairwise log difference in number of distinct recorded floral profiles**, using a pair-bootstrap for the adjusted intercept; the latter tests whether projected volume adds more than profile richness.

**Neither KDE is the entire categorical flower trait-space volume**: 17 axes still omit ~19.51% of positive inertia; Gaussian smoothing fills impossible/unobserved category combinations; continuous Lebesgue hypervolume is not intrinsically defined for unordered category-set combinations. Raw volume magnitude depends exponentially on bandwidth and dimension. Do not compare numeric 2D and 17D volumes directly. The exact original-distance Rao Q and exact-profile repertoire incidence remain primary.

## E. Privacy and method gate

The repo can execute 17D KDE on source-approved public **island×profile aggregates** and profile IDs; its CI reruns the independent geology reference join and emits only group summaries, matched-test aggregates, numerical error diagnostics and digests. The genus opportunity and two-stage modules take **private input paths** at execution time; public GitHub CI only runs synthetic tests, not raw GIFT or name-level source archives.

Never commit or attach private original taxon×trait/name rows, exact matched island/source IDs, historic source lineages or GIFT mainland entity checklists. Release only code, aggregate source/model outcomes and source-hash manifests.

To claim an actual **colonization mechanism**, independently estimated detection-adjusted island absence (surveys with repeat opportunities), observed historical access/propagule supply, source-region ancestry and phylogenetically informed entry survival will be needed. AUC on GBIF/name checklists does NOT provide this on its own.

## Reproducibility

Public geometric stage (the [origin-analysis workflow](../.github/workflows/origin-analysis.yml) executes it automatically):

```bash
python -m fd.hypervolume17 --input work/island_public --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/hypervolume17 --scopes all,direct --samples 12000
```

Private source-only stages (local, authorized input paths; do not upload to GitHub Actions):

```bash
python -m fd.source_genus_opportunity \
  --flora /secure/island/chapter1_status_flora.csv.gz \
  --mainland-native /secure/GIFT/gift_native_mainland_flora.csv.gz \
  --source-assignments /secure/GIFT/island_source_assignments.csv.gz \
  --geology /secure/fd/island_origin_crosswalk.csv \
  --matched-pairs /secure/fd/matched_pairs_private.csv \
  --output work/private_genus_opportunity
python -m fd.two_stage_assembly \
  --profile-all /secure/fd/all_profiles_private.csv.gz \
  --profile-direct /secure/fd/direct_profiles_private.csv.gz \
  --membership-all /secure/fd/all_matched_membership_private.csv.gz \
  --membership-direct /secure/fd/direct_matched_membership_private.csv.gz \
  --mainland-native /secure/GIFT/gift_native_mainland_flora.csv.gz \
  --assignments /secure/GIFT/island_source_assignments.csv.gz \
  --matched-pairs /secure/fd/matched_pairs_private.csv \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --output work/private_two_stage_repertoire
```
