# Own-database FD: species-name held-out genus prediction

25 November 2026. Exploratory H1/H2 extension, not a new submission lock.

**Execution:** [run 45546558379](https://github.com/zuizui0223/island/actions/runs/45546558379) SUCCESS, [artifact 14991284126](https://github.com/zuizui0223/island/actions/runs/45546558379/artifacts/14991284126), ZIP SHA256 `sha256:b88fdc9383750723416c3e302bf57add256acb9e4b132c2c20958abfc4412ace`. Twelve tests passed. All 54 estimable parent triple slopes reproduced for each of two fixed hash-fold seeds (108 fits; maximum |delta| 4.86e-17). Full 324 aggregate model rows, 84 proper-score rows, 16 support rows preserved from verified job log in `genus_crossfit_full_from_log.json`. Original aggregate artifact is retained separately. No original taxon/island-level rows exported.

## First gate: actual held-out name prediction

Entire accepted-name test fold is excluded from regional training, including blind prior frequencies. Fixed profile vocabulary, known genus and observed island membership; unseen genera fall back to blind. In the primary tau=1 model, seed 20261125:
- North-mid 1,435 target names, 1,012 with a training congener (70.5%); Brier .938026 -> .647935 and log loss 3.598406 -> 2.502635. Genus-equal gain intervals positive: Brier [.08412,.13042], log loss [.37193,.57098].
- Tropical 902 names, 500 with a congener (55.4%); Brier .959705 -> .835638 and log loss 3.873703 -> 3.266175. Genus-equal gain intervals positive: Brier [.02410,.06567], log loss [.15909,.32020].
- Both proper scores improve for both fixed seeds and both shrinkage values (tau=1,8) in these two regions. The same-fold label-change synthetic test verifies target-fold label exclusion. These are known-genus prediction results on documented species, not independent phylogenetic or geographic holdouts.

## Second gate: island functional structure

Primary tau=1, same coverage-adjusted corrected geography cohorts. MSE is island-equal squared error of predictive means. Block-equal gain intervals resample fixed predictions, not model refits; not formal familywise tests.

| Region (focal metric) | Seed | Islands | Observed slope | Blind forecast | Genus forecast | Island MSE reduction | Forecast/observed slope ratio | Block-equal MSE gain 95% interval |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| northern_midlatitude | 20261125 | 742 | 0.02003 | -0.00023 | 0.00979 | 54.6% | 48.9% | [0.000906, 0.003125] |
| northern_midlatitude | 20261126 | 742 | 0.02003 | 0.00031 | 0.00862 | 51.2% | 43.0% | [0.000914, 0.002923] |
| tropical | 20261125 | 105 | -0.01108 | -0.00014 | -0.00444 | 53.0% | 40.0% | [0.000225, 0.000577] |
| tropical | 20261126 | 105 | -0.01108 | -0.00038 | -0.00399 | 49.5% | 36.0% | [0.000205, 0.000508] |

North-mid focal metric is excess colour–architecture MI; tropical metric is marginal Rao Q. The model predicts the observed direction and lowers island-level MSE for both metrics in both seeds, including tau=8 sensitivity. It does not reproduce the full slopes. Genus forecast/observed slope ratios are numerical, not causal explained proportions. Predictive draw ranges include zero for both focal tau=1 forecasts; positive conditional skill does not establish a nonzero predicted gradient. Residual pointwise intervals exclude zero but are post-hoc conditional on the fitted predictive reference; no discovery claim. Compared with in-sample genus-exchange reference (63.3% and 78.3% numerical ratios), honest name-held-out forecasts are attenuated; the estimands differ and this is not a formal attenuation test. Mean island-level carrier donor coverage is about 72–74% north-mid and 45–47% tropical. Incomplete same-genus donors and shrinkage, limited trait evidence, annotation/lineage processes and regional composition may create forecast gaps.

Across all 216 estimable nonblind model rows (nested variants/seeds/metrics), 208 have positive island MSE gain and 194 have positive lower conditional block interval. These counts are descriptive, not independent confirmations. Literal-native north-mid 67-island MI prediction improves (tau1 MSE reduction 39.3% and 25.5% by seed); its observed slope +.00484 is weak and residual intervals cross zero. Tropical native remains below the parent model support gate.

## Existing H2 matched-carrier adjustment, descriptive sensitivity

| Region | Seed | n | Observed beta | Genus-predicted beta | Residual beta |
|---|---:|---:|---:|---:|---:|
| northern_midlatitude | 20261125 | 739 | 0.01770 | 0.00925 | 0.00845 |
| northern_midlatitude | 20261126 | 739 | 0.01770 | 0.00822 | 0.00948 |
| tropical | 20261125 | 68 | -0.01087 | -0.00544 | -0.00543 |
| tropical | 20261126 | 68 | -0.01087 | -0.00512 | -0.00575 |

The actual selfing covariate is held fixed, not predicted jointly with flower morphology. No new pathway-independence claim.

## H1–H4 decision

The new evidence supports **genus-linked predictive information about floral functional assemblage structure**, distinct from a pure in-sample genus null. It does not show that genus composition explains all isolation-linked differences or that residuals are adaptive evolution. The parent H1 SC recurrence and H2 conditional accessibility interpretation remain unchanged; H3/H4 functional GloPL links are separate. The next identification gate is to compare regions on a common distance/covariate-support cohort and to evaluate predictions trained outside target spatial blocks and target-name folds, with identical donor opportunity and no substitution of a different cohort. Prior GEB freeze and all negative results retained.

See [design](../../docs/fd_self_db_genus_crossfit_design_20261125.md).
