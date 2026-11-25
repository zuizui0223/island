# Own-database FD: coherent regional genus null — verified result

25 November 2026. Follow-up to the exploratory corrected-isolation model, not a GEB submission replacement.

**Completed:** [run 45544038800](https://github.com/zuizui0223/island/actions/runs/45544038800) SUCCESS; [artifact 14991150352](https://github.com/zuizui0223/island/actions/runs/45544038800/artifacts/14991150352), ZIP digest `sha256:72f52097456fadba4400db5b6a2cffef266b2f2e1ee623a81ad199803f093368`. Nine tests passed. All 144 aggregate rows and eight support rows reproduced from verified job log and archived in `genus_null_full_144_from_log.json`. All 54 parent estimable triple slopes were reproduced, maximum absolute beta difference 2.78e-17. 54 residual fits; 90 nonfit slots preserved. No raw island or species identities committed. Original 399 whole-profile within-region, within-genus exchanges reuse each species assignment across islands, with fixed source carrier and rarefaction masks. The native status subset has its own regional species pool. See the [locked design](../../docs/fd_self_db_genus_null_design_20261125.md).

## Main decision: attenuated but not zero; no multiplicity-robust coherent-null departure

| Focal broad coverage-adjusted result | Islands/blocks | Observed beta | Genus-null mean beta | Residual beta (conditional CR1 95% interval) | Coherent permutation p / focal Holm2 | Residual BY q |
|---|---:|---:|---:|---:|---:|---:|
| Northern midlatitude excess MI | 742/52 | +0.020026 | +0.012683 | +0.007343 [+0.002634,+0.012053] | 0.050 / 0.100 | 0.05378 |
| Tropical marginal Q | 105/41 | −0.011084 | −0.008676 | −0.002408 [−0.003845,−0.000971] | 0.090 / 0.100 | 0.05378 |

Expected-to-observed slope ratios are 63.3% and 78.3% respectively. These are numerical comparisons under a specified in-sample genus-exchange null; **not causal explained variance or mediation proportions**. The signed residuals are not zero, but neither passes the two-focal-test Holm gate or the 12-slot BY-adjusted permutation gate. No permutation test passes BY across any variant/scope family. Conditional residual t inference treats the estimated reference as fixed and is not equivalent to the coherent label-exchange distribution.

## Conditional residual sensitivities retained, not promoted

- Northern midlatitude unadjusted-coverage baseline MI residual +0.008008 [0.003431,0.012585], conditional BY q=.03500, permutation p=.040, permutation BY q=.74477. Parent raw baseline was not BY-supported.
- On parent same-morphology-roster selfing-adjusted sample, north-mid MI residual +0.006585 [0.002411,0.010760], conditional BY q=.04837, coherent permutation p=.055 (BY=.68271); tropical marginal Q residual −0.002353 [−0.003717,−0.000989], conditional BY q=.04782, coherent p=.055 (BY=.68271). Those two sensitivity intervals are not independent confirmation; same-roster selfing is scored on a subset, not the entire morphology sample.
- Whole-profile null does not preserve morphology–reproductive trait coupling within each species; the observed selfing covariate is held fixed. Therefore H2 residual analyses are diagnostics, not isolation effects independent of all reproductive and lineage pathways.
- Literal-native north-mid 67 eligible islands; tropical 17, southern 3, high north zero. 18 native residual models fitted; 0 native models pass residual or permutation BY. Do not generalize this to native-only global absence or equivalence.

## Exchange opportunity

All-recorded regional pool contains 1,472 Direct triple species/605 genera in north-mid and 980/541 in the tropics (including islands below the 30-carrier target threshold). Median fraction of target-island carriers with a variable genus pool is 48.9% and 38.1%, respectively; high north 23.7% and southern extra-tropical 42.2%. This is nonzero but limited exchange support, not phylogenetic independence. Accepted-name genus token is a proxy, not a reviewed phylogenetic backbone. There were zero malformed genus singleton tokens in the analyzed pools. The same-name taxa are repeatedly represented among islands; within-region coherent assignments preserve their dependence but separate regional pools do not give a formal region interaction. The difference of expected/observed slopes is not a proof of historical colonization or selection.

## Relation to H1–H4

H1 remains the frozen traitwise result (SC recurrence and context-dependent floral responses). FD adds the distinction between marginal diversity and colour–architecture coupling, and the new result supports compatibility with regional genus composition as an important explanation, not proof that all residual responses are explained by genus. H2 retains its conditional accessibility decomposition. The FD evidence does not establish a new genus-independent or selfing-independent direct pathway. H3/H4 remain independent experimental-pollen-limitation context and exact-species compatibility; no new pollinator mechanism or causal chain. A more stringent next diagnostic needs held-out accepted-species genus donors, donor-support matching, and exact common-island/coverage controls, or an explicit regional common-distance-support contrast. None of the present post-hoc tests warrants promotion to frozen GEB conclusions.
