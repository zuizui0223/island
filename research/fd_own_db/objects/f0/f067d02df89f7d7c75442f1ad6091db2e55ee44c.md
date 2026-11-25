# GIFT 4,505-island Direct flower prediction after removing repeated target-species and genus pseudoreplication

**Completed-source-artifact verification is now independently GREEN.** The original 78-panel numerical source fit completed and was rechecked by [independent source-artifact Actions #44757114610](https://github.com/zuizui0223/island/actions/runs/44757114610): **SUCCESS**, exact frozen JSON SHA256 `1370f468cb8c4fb559c1a23070f685cf2d12c64931f434107adba4ad9acba5b8` matched alongside prior original GIFT category JSON `7cbdbc34d20e1d5902c464cc1d93786004f3718848317db3e83038d6063e7201`. The dedicated source-side CI invoked the already passing independent fd auditor on BOTH exact downloaded original artifacts, checked 4,505 parent/78 panel original identity and previous proper-score invariance, and reproduced **0 genus-equal** and **2 species-equal all-recorded** strict triple-reference positives. This is a **new successful verification of a previously completed numerical artifact**, not a claim that the original full fit Actions #44669409016 (failed on a stale metadata check) succeeded or that a fresh second fit has completed. Protected original [island PR #342](https://github.com/zuizui0223/island/pull/342), initial source-original [Actions #44669409016](https://github.com/zuizui0223/island/actions/runs/44669409016): 17 source tests PASS, full original 4,505 GIFT botanical flora / 8,264 GSHHG parent restored with source SHA checks, original 13 true species-Direct flower/reproduction endpoints × 3 disjoint macroregions × nested all-recorded/native botanical status = exactly **78 original panels**, model fit finished, full anonymous artifact **#14603722607**, downloaded ZIP **SHA256 2828a3549cff10fb5cf3843466bfdb93e8501c783cbcf52e90b564a6cf85993e**, one 1,139,332-byte anonymous JSON. The run's subsequent output-integrity step failed solely because its workflow copied obsolete **local-vs-external blind flag assertions** from a different source experiment: these keys are not produced by PR342's otherwise completed original source model. [Corrected source workflow](https://github.com/zuizui0223/island/pull/342), [corrected Actions #44750154943](https://github.com/zuizui0223/island/actions/runs/44750154943) rerunning original data. Thus the **original long-running fit #44669409016 still correctly has a failed overall GitHub conclusion**, but the exact source-output artifact has passed a separate, SHA-pinned, original-source independent verification workflow #44757114610. The corrected fresh full-fit source run #44750154943 is still being checked separately; disagreement would require a reconciliation, not post-hoc selection. Numerical findings in this report are supported by the GREEN source-artifact verification, with that transparent provenance caveat.

An independent downloaded-file audit confirms source status, 78 unique exact panels, original full-flora holdout declarations, no exported species/genus/island/source identifiers, original 4,505 parent, exact original scored source island counts and source group support flags, and complete **original geographic-group Brier and log-loss point estimates plus confidence intervals for three baseline comparators** compared to prior verified GIFT original 13-endpoint Actions #44519589815 (ZIP SHA256 015b4910119d2ba76266ecb2af15addc30ed3c43fc6fd50d976bf00c72995052), all identical under 2×10^−12. All 78 taxon-reweighting stored positive flags were independently recomputed from their paired Brier/log-loss bootstrap source intervals. The fd source-safe independent auditor is scripts/verify_gift_target_taxon_source.py with CI-tested null-cohort support guards, corrected in [fd PR #58](https://github.com/zuizui0223/fd/pull/58) after source evidence showed 17 genuinely unsupported panels omit meaningless repeat counts.

[Full 78 × 2 distinct-taxon bootstrap interval and reference benchmark audit table](GIFT4505_SOURCE_TARGET_SPECIES_AND_GENUS_EQUAL_78.csv).

## Different estimands, not refitting the source model

All sources use the previously verified strict no-actual-focal-genus-oracle geography→predicted botanical genus→other-island congeneric Direct category model. Source-fold botanical species names excluded from training before source genus and Direct flower stage are never used to predict held-out botanical species. However, species-wide Direct traits are reused on repeated island floras. Block/group source resampling treats geographical groups as units but does not remove dependence caused by the **same species** appearing on multiple target islands.

New sensitivity weights are calculated after prediction on the same exact held-out source plant carriers, relative to three already verified source category references: target-region genus blind, externally sourced genus blind and simple pooled source-region GENUS mixture.

- **Unique accepted target species-equal**: mean the same distinct botanical species' paired Brier/log contributions across all its scored island occurrences, then bootstrap each distinct name as one unit (minimum 30 source names).
- **Unique true target genus-equal**: average each species' already species-equal scores within its genus and bootstrap genera as units (minimum 10 genera). This is a coarse similarity sensitivity; biological genera are NOT independent phylogenetic lineages.

The 499-resample 95% intervals are exploratory and must be positive for **BOTH** Brier and log loss for a baseline-positive verdict. A strict signal must beat all THREE source references under each weighting estimand. These are **different source-target estimands**, not interchangeable source-region island average estimates. Source scopes and traits overlap; the 78 cells and thousands of taxa do not provide 78 independent island evolutions.

## Source-original findings independently verified on SHA-pinned completed artifacts

| Source original all-recorded 13 regional endpoints | Source three-reference gains using unique species units | Gains using unique genus units |
| --- | ---: | ---: |
| Northern midlatitudes (13 endpoints) | **0** | **0** |
| Tropics (13 endpoints) | **0** | **0** |
| Southern extratropical (13 endpoints) | **2** | **0** |
| **Total, 39 source region×trait cells** | **2** | **0** |

The two species-equal only supported Southern endpoint×region cases are **floral symmetry** (original target n178 island floras) and **TRUE same accepted-species primary flower colour×symmetry** (n162), but neither clears all three source references under genus-equal intervals. The nested literal native-only source sensitivity has one species-equal three-reference positive (tropical symmetry), but **zero genus-equal** three-reference positives across 39 native-scoped cells and zero Southern-native source-qualified cases.

### Original two previously strongest island-equal category + nominal Q candidates

| Property | Tropical floral symmetry | Southern extratropical primary flower colour |
| --- | --- | --- |
| Original source target island floras | **925** | **184** |
| Distinct accepted source-target species across them | **22,914** | **7,964** |
| Distinct target genera | **3,016** | **2,249** |
| Mean source target island occurrences per distinct botanical species | **5.12** | **3.78** |
| **Species-equal** gain vs LOCAL blind Brier [95%] | +0.000726 **[−0.000226,+0.001546]**, ambiguous | +0.006819 **[+0.004620,+0.009168]**, positive |
| **Genus-equal** gain vs LOCAL blind Brier [95%] | +0.003319 **[+0.000921,+0.005756]**, positive | +0.011031 **[+0.007554,+0.014760]**, positive |
| **Species-equal** gain vs EXTERNAL blind Brier [95%] | +0.000764 **[+0.000628,+0.000911]**, positive | −0.000147 **[−0.000547,+0.000282]**, ambiguous |
| **Genus-equal** gain vs EXTERNAL blind Brier [95%] | −0.000158 **[−0.000474,+0.000149]**, ambiguous | +0.000932 **[+0.000217,+0.001714]**, positive |
| **Species-equal** gain vs pooled source-region GENUS Brier [95%] | +0.001380 **[+0.001181,+0.001604]**, positive | −0.000214 **[−0.000433,+0.000023]**, ambiguous |
| **Genus-equal** gain vs pooled source-region GENUS Brier [95%] | −0.000219 **[−0.000674,+0.000196]**, ambiguous | +0.000351 **[−0.000064,+0.000744]**, ambiguous |
| **All three source-reference proper losses positive?** | **NO for species and genus units** | **NO for species and genus units** |

Proper log loss intervals were assessed simultaneously and produce the stricter overall verdict; the displayed Brier intervals alone do not constitute complete support. The ORIGINAL island/group-equal source prediction and corrected nominal categorical Rao Q positives for these two cases remain real source outcomes for botanical island checklist weights. They cannot be advertised as independent across species/genus or demonstrated island evolutionary flower transitions.

## Strong conclusion and ceiling

**Result changes with the denominator.** Geography- and target island-equal GIFT floral category source prediction improves in the two original special cases. Species-equal and genus-equal weighting can still support subsets of their benchmark comparisons, including Southern colour versus local blind under BOTH taxon equalizations. But the more stringent universal **all-three-reference** information gain disappears at genus level in ALL 39 original all-recorded region×trait panels. This is a direct empirical warning against claiming broadly taxonomically replicated floral adaptation or universal channel-gated pollinator-syndrome loss from species-wide island checklist traits alone.

These are original botanical source numerical results, not simulated island communities, and the **exact completed artifact** has now passed independent source Actions #44757114610. This supports including the source-weighting sensitivity in the fd evidence ledger, while preserving a strict distinction: original full-refit Actions #44669409016 failed on a stale output-key assertion; corrected independent full-refit #44750154943 has not yet been confirmed complete. If that fresh fit later disagrees, reconcile rather than cherry-pick these results.

**Nonclaims:** no verified pollinator guild transitions, actual mainland progenitor phenotypes, within-species island flower evolution, independent genus phylogenetic histories, causal sources of selection, measured dispersal/establishment, general mixed-trait Gower-Rao morphology or 17D floral hypervolume. Missing data/Direct trait documentation and nested native cohort shortfalls remain, and 95% intervals are exploratory/non-multiplicity adjusted.
