# Supplementary Information draft

## Article

**A recurrent global floral island syndrome extends beyond the selfing syndrome**

This Supplementary Information draft is tied to the corrected 8,264-unit Chapter 1 submission surface selected by `config/chapter1_submission_current.json`. It does not use superseded v13/v14 geography or historical WHEN/WHERE branches.

## Appendix S1 — Geographic exposure and geography correction

### S1.1 Corrected island universe

The plant analysis uses 8,264 island units. One GSHHG split component corresponding to continental Eurasia was excluded after the geography audit.

Primary machine-readable sources:

- `results/geography_20260924/corrected_universe_exclusion.csv`
- `results/geography_20260924/confirmed_distance_and_universe_issues.json`
- `results/geography_20260924/corrected_geometry_gate.json`

### S1.2 Source-matched coastline distance

Isolation is measured as minimum minor-great-circle arc separation between island and continental coastlines from the same GSHHG 2.3.7 archive on a mean-radius sphere.

Sources:

- `results/geography_20260924/corrected_geography_covariates.csv`
- `results/geography_20260924/gshhg_spherical_distances_all.csv`
- `results/geography_20260924/distance_same_coastline_audit.json`
- `results/geography_20260924/spherical_geometry_validation.json`

### S1.3 Repaired and retained zero distances

The correction repaired 1,113 spurious island zero distances produced by the previous mixed-resolution geometry. GloPL retains 996 site coordinates that truly lie on continental land.

Sources:

- `results/geography_20260924/formerly_zero_islands_recalculated.csv`
- `results/geography_20260924/glopl_corrected_site_distances.csv`

**Figure S1:** old versus corrected island distance for the 1,113 formerly zero units.

---

## Appendix S2 — Island flora and trait evidence

### S2.1 Flora and taxonomic universe

The trait database contains 106,295 accepted angiosperm species and three raw evidence axes: flower colour, floral structural complexity and reproductive assurance.

### S2.2 Trait evidence and missingness

The frozen scientific ledger contains 222,688 resolved cells of 318,885 possible species × axis cells (69.83%). Missing values remain missing. All-analysis-eligible and Direct-only evidence scopes are analysed separately.

### S2.3 Provenance and redistribution rights

The complete scientific database is identified by:

`sha256:a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a`

Current redistribution audit:

- resolved scientific cells: 222,688;
- redistribution-authorized cells: 46,274;
- blocked/review-required cells: 176,414.

The public subset is not the canonical analysis database.

Supporting provenance:

- `docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md`
- `docs/CHAPTER1_DATABASE_PUBLIC_SUBSET.md`
- `config/chapter1_database_versions/v1.0.0.yml`

**Table S1:** trait coverage and evidence scope by axis.

**Figure S2:** resolved-cell coverage and Direct-only versus all-analysis support.

---

## Appendix S3 — H1 recurrent multivariate island response

### S3.1 Primary all-analysis models

Atomic isolation slopes:

- `results/geography_20260924/all/beta_binomial_within_slopes.csv`

Joint seven-response tests:

- `results/geography_20260924/all/beta_binomial_within_omnibus.csv`

Original-versus-corrected comparisons:

- `results/geography_20260924/all/comparison_beta_binomial_within_slopes.csv`
- `results/geography_20260924/all/comparison_beta_binomial_within_omnibus.csv`

### S3.2 Direct-only sensitivity

- `results/geography_20260924/direct/beta_binomial_within_slopes.csv`
- `results/geography_20260924/direct/beta_binomial_within_omnibus.csv`
- `results/geography_20260924/direct/comparison_beta_binomial_within_slopes.csv`
- `results/geography_20260924/direct/comparison_beta_binomial_within_omnibus.csv`

The main claim is multivariate recurrence across four geographic strata. Individual atomic coefficients are retained to show that recurrence does not imply uniform positivity.

### S3.3 Northern-high-latitude Direct-only optimizer audit

The frozen corrected Direct-only northern-high-latitude table marked the shallow/open-tube component as `optimizer_success=false`, while the seven-response joint test remained FDR-supported. Because `vector_supported` in the historical runner was based on q-value and did not itself require every optimizer flag to be true, we performed a dedicated numerical audit rather than treating the flag as cosmetic.

Audit source:

- `results/geography_20260924/h1_direct_northern_high_convergence_audit.json`
- `scripts/geography_correction/audit_h1_direct_convergence.py`

Results:

- frozen shallow/open-tube estimate: 0.3544814;
- enhanced re-fit estimate: 0.3544784;
- absolute change: 2.93 × 10^-6;
- enhanced re-fit: converged;
- fully converged seven-response replay: q = 3.793 × 10^-8;
- fully converged six-response sensitivity excluding shallow/open tube: q = 1.430 × 10^-8.

Thus the Direct-only northern-high-latitude H1 conclusion is numerically stable and does not depend on the warned component.

**Table S2:** complete H1 coefficient, SE, P and FDR-adjusted q table.

**Figure S3:** Direct-only H1 coefficient forest plot beside the primary all-analysis result.

---

## Appendix S4 — H2 conditional decomposition and raw colour/architecture

### S4.1 Reproductive-assurance-adjusted models

Primary:

- `results/geography_20260924/all/h2_decomposition_models.csv`
- `results/geography_20260924/all/h2_summary.json`

Direct-only:

- `results/geography_20260924/direct/h2_decomposition_models.csv`
- `results/geography_20260924/direct/h2_summary.json`

Corrected-versus-parent comparisons:

- `results/geography_20260924/all/comparison_h2_decomposition_models.csv`
- `results/geography_20260924/direct/comparison_h2_decomposition_models.csv`

The tropical Direct-only accessibility estimate is nominally positive but not FDR-supported after correction (q = 0.1196).

### S4.2 Raw colour states

All-analysis:

- `results/geography_20260924/all/raw_patterns/raw_colour_model_results.csv`
- `results/geography_20260924/all/raw_patterns/raw_colour_joint_omnibus.csv`

Direct-only:

- `results/geography_20260924/direct/raw_patterns/raw_colour_model_results.csv`
- `results/geography_20260924/direct/raw_patterns/raw_colour_joint_omnibus.csv`

### S4.3 Colour × architecture

Joint colour–architecture prevalence:

- `results/geography_20260924/all/raw_patterns/raw_colour_architecture_model_results.csv`
- `results/geography_20260924/direct/raw_patterns/raw_colour_architecture_model_results.csv`

Architecture conditional on colour:

- `results/geography_20260924/all/raw_patterns/raw_colour_conditioned_architecture_model_results.csv`
- `results/geography_20260924/direct/raw_patterns/raw_colour_conditioned_architecture_model_results.csv`

Comparison files in the same directories preserve every original-versus-corrected row.

**Table S3:** complete H2 decomposition results.

**Table S4:** all raw colour and colour × architecture results, including null and opposite-sign rows.

**Figure S4:** regional raw-colour response matrix.

**Figure S5:** colour-conditioned architecture heatmap.

---

## Appendix S5 — H3 experimental pollen limitation

### S5.1 Corrected site exposure

- `results/geography_20260924/glopl_corrected_site_distances.csv`

### S5.2 Model-ready GloPL data

- `results/geography_20260924/h3_corrected_effect_rows.csv.gz`
- `results/geography_20260924/h3_corrected_measurement_cells.csv.gz`

### S5.3 Corrected-versus-parent comparison

- `results/geography_20260924/h3_original_corrected_comparison.json`

Primary corrected estimate:

- beta = 0.09191;
- SE = 0.03806;
- two-sided P = 0.01575.

The no-zero-constant sensitivity remains supported; the supplemental-only sensitivity remains positive but unsupported.

**Table S5:** H3 primary and sensitivity models.

**Figure S6:** H3 primary, no-zero and supplemental-only coefficient comparison.

---

## Appendix S6 — H4 exact-species functional bridge

### S6.1 Literal H2 score bridge

- `results/geography_20260924/h4_exact_corrected.csv`
- `results/geography_20260924/h4_exact_original_reproduced.csv`
- `results/geography_20260924/h4_exact_comparison.csv`

Model-ready effect rows:

- reproductive assurance: `results/geography_20260924/h4_ra_corrected_effect_rows.csv.gz`
- accessibility: `results/geography_20260924/h4_arch_corrected_effect_rows.csv.gz`

Primary corrected estimates:

- selfing_core: beta = -0.29830, P = 0.00396;
- generalized_accessible: beta = -0.29566, P = 0.02187.

### S6.2 Atomic-trait sensitivities

- `results/geography_20260924/h4_atomic_corrected.csv`
- `results/geography_20260924/h4_atomic_original_reproduced.csv`
- `results/geography_20260924/h4_atomic_comparison.csv`
- `results/geography_20260924/h4_atomic_refit_manifest.json`

Autonomous selfing is the strongest atomic functional association; actinomorphy is also negative and supported. Self-compatibility is negative but imprecise.

**Table S6:** H4 score models and sensitivities.

**Figure S7:** atomic-trait H4 coefficient forest plot.

---

## Appendix S7 — Inferential boundaries and reproducibility

### S7.1 Claim hierarchy

The current paper permits:

1. recurrent multivariate floral/reproductive response to isolation;
2. partial separation of reproductive assurance and floral accessibility;
3. an independent isolation-associated pollen-limitation gradient;
4. post-hoc functional compatibility of response trait states with lower current pollen limitation.

It does not establish:

- historical mediation;
- global pollinator abundance or visitation decline;
- realized pollinator identity from floral traits;
- uniform direction of every atomic trait;
- species sorting versus within-lineage evolution;
- prospective confirmation from the geography repair.

### S7.2 Repository replay

- `results/geography_20260924/repository_replay_verification.json`
- `results/geography_20260924/refit_completion.json`
- `results/geography_20260924/raw_pattern_refit_completion.json`

### S7.3 Submission data-policy boundary

Ecology Letters data-policy status:

- `submission/chapter1_current/ECOLOGY_LETTERS_DATA_GATE.md`

Rights-aware data statement:

- `submission/chapter1_current/DATA_ACCESSIBILITY_DRAFT.md`

The SI is not a substitute for permanent data/code archiving.

---

## Supplementary table assembly map

| Supplement item | Source |
| --- | --- |
| Table S1 | trait database coverage/provenance manifests |
| Table S2 | all/direct H1 within-slope + omnibus CSVs |
| Table S3 | all/direct H2 decomposition CSVs |
| Table S4 | all/direct raw-pattern CSVs |
| Table S5 | H3 comparison JSON + corrected effect/cell files |
| Table S6 | H4 exact + atomic corrected/comparison CSVs |
| Table S7 | geography, replay and redistribution-rights receipts |

## Supplementary figure assembly map

| Supplement item | Source |
| --- | --- |
| Figure S1 | formerly-zero distance audit |
| Figure S2 | trait coverage/evidence-scope audit |
| Figure S3 | H1 all-analysis vs Direct-only coefficient tables |
| Figure S4 | raw-colour model results |
| Figure S5 | colour-conditioned architecture results |
| Figure S6 | H3 primary/sensitivity comparison |
| Figure S7 | H4 atomic corrected results |
