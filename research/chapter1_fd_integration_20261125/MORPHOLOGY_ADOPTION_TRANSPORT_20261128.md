# Final morphology adoption review: from local composition to geographic transport (28 November 2026)

Author-side update to the previously rendered PR354 review; not a new submission freeze. Original H1–H4 and six main figures are unchanged. PR362 and PR363 are same-database, outcome-motivated diagnostics and not independent replication.

## Evidence chain

| Analysis | Main result | Interpretation |
|---|---|---|
| PR354 tropical paired symmetry × depth | 257 islands/718 species, beta +.027620, Holm .0002206 | Exploratory composition contrast, high-only unsupported, literal native n44 insufficient |
| PR362 local genus margin | observed and expected +.027620 exactly, zero exchangeable genus-grade strata | Within-genus association not identifiable; arithmetic reconstruction is not independent prediction |
| PR363 external held-block transport | conditional contrast RMSE .251445→.159921; Brier gain +.050872 | Known-genus conditional composition prediction improves |
| PR363 held-out isolation slope | observed +.027620; external +.003975; other tropical +.003707 | Good composition prediction does not transport the isolation response |
| PR363 coverage partition | supported and unsupported genus components both positive in the broad primary cohort | Donor coverage alone does not account arithmetically for the residual; limited congener sample sizes and shrinkage remain alternatives |

## Coverage decomposition (all original panel fields replayed)

[Run 46113724276](https://github.com/zuizui0223/island/actions/runs/46113724276), job 140882736641, artifact 15250753067, executable head da09942e282bf17ae6630d324550fd315a4f6439. The [full eight-panel receipt](FLORAL_GENUS_TRANSPORT_COVERAGE_RESULT_20261128.json) was recovered from the successful log and checked against the original receipt, with SHA256 metadata (ZIP bytes not independently checked). All original fields match within 1e-10. Anonymous-only outputs.

| Cohort | Training region | Residual isolation slope | Donor-supported genus contribution | Donor-unsupported genus contribution | Unknown genus carrier mass |
|---|---|---:|---:|---:|---:|
| Broad n257 | External | .02364552 | .01298261 | .01066291 | .534763 |
| Broad n257 | Other tropical | .02391297 | .01383389 | .01007908 | .470515 |
| Broad n122 | External | .02082714 | .00873588 | .01209127 | .548624 |
| Broad n122 | Other tropical | .02240573 | .00858028 | .01382544 | .478205 |
| High-only n62 | External | .00162123 | .00186814 | -.00024691 | .701163 |
| High-only n62 | Other tropical | not estimated | not estimated | not estimated | 6 failed folds |
| Literal native n44 | External/other tropical | not estimated | not estimated | not estimated | insufficient parent support |

For each island, split observed and predicted four-cell proportions into disjoint contributions from carriers with or without an eligible donor congener. Evaluate the conditional contrast on two possible orders of replacing predicted component vectors by observed component vectors; average each component's two marginal differences (two-factor Shapley accounting). Both replacement vectors use target outcomes strictly after prediction. The two terms sum exactly to the residual at each island and hence on the frozen linear slope design. This is neither a new predictive model nor a causal partition. The 499 block resampling ranges of the original forecasts are descriptive, and no new test is assigned to these residual components.

In the primary external panel, unsupported-genera error accounts for .033759 of .047038 total mean squared category error by signed projection (71.77% arithmetic ratio), but only .010663 of .023646 slope residual (45.09%). Supported genera contribute .012983 (54.91% of slope residual). Average carrier weights across 0, 1, 2–4 and ≥5 eligible donor species are 53.48%, 14.55%, 13.26%, 18.71%. Thus most carriers have little or no congener information, and nonzero donor support is not dense sampling. In n122 the unsupported-genus slope component is larger; no universally dominant component is identified. No causal inference that filling traits would recover the slope, no adaptive local divergence inferred, no claim of independent evidence or confirmed cross-region effect.

## Adopt in a single paper

**Main conclusion:** recurrent self-compatibility across four regions; contingent floral composition; current reproductive constraint and functional compatibility. Morphology contrasts are useful secondary examples, with the exact genus-margin and cross-validation limitation adjacent. No general FD contraction, no identifiable within-genus coupling, no named pollinator replacement.

The new text-only [author-review manuscript](geb_review/MANUSCRIPT_BLINDED_REVIEW_TRANSPORT_20261128.md) adds the diagnostic to the prior rendered candidate, not the active GEB submission. The present addendum is intended as supplemental author-side evidence, not a finished blinded SI file. Figures 1–6 and S8 unchanged. Previous rendered artifact does not certify the new manuscript. Re-render and proof after author acceptance; source rights, editor access and authorship gates remain.

## Stop rule and next material action

Do not add more small residual partitions, change regularization or screen new endpoints merely to recover the original tropical slope. The existing positive H1–H4 spine does not require this slope to transport. If a stronger floral reassembly claim is desired later, prioritize new paired native-flora high-quality trait observations and eligible independent congener coverage, then evaluate a frozen cohort. The current same-roster analysis cannot identify that claim. No new empirical run is required for the present editorial disposition.
