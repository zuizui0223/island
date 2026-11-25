# GloBI v0.9 original provider-label omission: matched botanical-family prediction

**11 November 2026 — confirmed real-data result, not actual pollination effectiveness.**

Source pipeline: [island PR #323](https://github.com/zuizui0223/island/pull/323), **successful [Actions #43029060199](https://github.com/zuizui0223/island/actions/runs/43029060199)**, anonymous original-source receipt [artifact #13860344048](https://github.com/zuizui0223/island/actions/runs/43029060199/artifacts/13860344048). Exact GloBI June 4 v0.9 `interactions.tsv.gz` verified SHA256 `095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17`, combined with the previously source-pinned original `island` Direct plant-trait ledger and frozen GBIF family crosswalk. The pre-specified genus universe was the original **2,228 GloBI source genera**, selected before island floral outcomes.

## The real new comparison uses identical genera and families

One source-provider label, selected because it had the largest share of the **archive-first** cited reference units, was removed. GloBI cited-reference channel fractions were reconstructed from the original rows on both sides. For the matched genera that retained **>=3** cited interaction references after the removal, a fresh documentation-volume-only baseline was fitted separately before and after removal; adding the **Direct source-backed open-radial flower fraction** was scored with an entire botanical family absent from training.

- Matched target sample: **344 botanical genera / 77 unambiguously assigned families**, not the earlier unfiltered 484 genera / 88 families.
- Same botanical test genera and same 77 heldout-family training exclusions were used before/after. The response is a *different GloBI documented reference composition after removing the largest archive label*, so the two target measurements are correlated rather than independent experiments.
- Source representation: 62 **archive-first labels**; one label accounts for **95.9584%** of 369,659 genus×channel×provider×citation units. The label's **within-record source metadata** contains **one** distinct `sourceNamespace`, **one** `sourceCitation` and **zero** `sourceDOI` values. This does **not** establish that the label equals one independent human study; it could represent an umbrella archive.
- Source genus support with >=3 citation-backed units: **1,674 before** versus **846 after**, within the original 2,228-genus universe; matched before/after **846**. Median within-genus documented five-channel composition change after removal is **0.260904** total variation.

### Prediction: open flower form added to source reference-effort baseline

| Same 344-genera/77-families | Before deleting largest archive label | After deleting it |
| --- | ---: | ---: |
| Genus-equal source-effort-only five-channel MSE | **0.034927** | **0.034739** |
| Genus-equal + Direct open-form model MSE | 0.033868 | 0.034274 |
| Genus-equal open-form incremental MSE gain | **+0.001059** | **+0.000465** |
| Family-equal open-form incremental gain | **+0.001769** | **+0.000687** |
| Conditional 4,999-draw one-family-one-vote 95% interval | **[+0.000838,+0.002735]** | **[−0.000196,+0.001535]** |
| Family means with positive gain | **52/77** | **45/77** |

**Mean family-equal gain falls by 61.2%**, and the after-omission conditional interval includes zero. The genus-equal gain falls by ~56.1%. These are *same-cohort descriptive attenuation estimates*. The current independently resampled before and after intervals **do NOT yet provide a paired CI on the change itself**; a further same-77-family paired resampling/guild-level extension is implemented in [island PR #323](https://github.com/zuizui0223/island/pull/323) and a fresh original-source rerun was started. Do not present the 61.2% reduction as a statistically confirmed difference before its paired interval is available.

**This is an important negative-control / sensitivity result:** a previously positive, weak morphology-to-GloBI-documentation relation depends strongly on the original source archive. The reduced mean gain remains positive, so it is not evidence that the effect is biologically zero. GloBI source documentation is not independently measured local pollinator functional composition. The huge removed fraction makes archive-level provenance definition paramount.

## Independent source ID / biological interpretation

The source grouping favored `sourceArchiveURI`, then `sourceDOI`, then `sourceNamespace + sourceCitation`. Even though the dominant label has only one recorded namespace and citation string, provider labels do not establish independent biological observations or distinct publications. A large archive may aggregate many studies, and the `reference_id` identities may derive from shared metadata rather than independent sampling. The original scanner excluded 32 reference keys occurring across provider labels, but cannot guarantee original-field DOI or citation completeness. This limitation is explicit.

The model's five channels (Bombus, other bees, Lepidoptera, Diptera, flower-visiting birds) are **global GloBI documentation fractions**, not visitation frequencies, actual island pollinator communities, effective pollen delivery or adaptation. Lepidoptera does not distinguish butterflies and moths. Source provider dependence cannot prove or disprove a pollination syndrome.

## Results and next gated test

**Supported:** A same-genus, same-family comparison of before/after source archive deletion reveals substantial attenuation of a documented visitor-channel composition *prediction improvement* attributed to coarse open flower morphology.

**Not established:** independent multi-provider replication, after-omission model gain distinctly positive on equal-family uncertainty, paired attenuation confidence interval, true island guild presence or pollination service.

The next source-pinned rerun is [island Actions #43034251795](https://github.com/zuizui0223/island/actions/runs/43034251795). It adds a **proper paired family resampling of the after-minus-before incremental gain** using shared 77 family units, and evaluates the individual five documented channels in the same heldout cohort. If source support is insufficient, the exact code must report `NOT_IDENTIFIABLE`. The original data and taxon/provider reference IDs remain in restricted input; fd receives anonymous aggregate summary only. Draft PRs are kept unmerged pending source integrity review.


## 12 Nov 2026: the paired 77-family result has completed

The new original-GloBI provider/whole-family workflow [#43034251795](https://github.com/zuizui0223/island/actions/runs/43034251795) **SUCCEEDED**, with anonymous [source receipt artifact #13864221183](https://github.com/zuizui0223/island/actions/runs/43034251795/artifacts/13864221183). Both source conditions use exactly **344 genera and 77 botanical families**. Unlike the previous separate confidence intervals, the completed analysis resamples the **same 77 families as PAIRS**:

- Source-archive omission change in mean family-equal flower-open-form MSE gain: **−0.001082**, paired conditional 95% interval **[−0.002177,−0.000107]**; negative throughout. This establishes an exploratory, **conditional source-documentation attenuation** under this coarse archive-label omission, not absence of actual pollination service.
- Non-Bombus bees: family-equal incremental gain **+0.004348 before / +0.001087 after**, *paired change* **−0.003262 [−0.006630,−0.000129]**.
- Lepidoptera: **+0.001609 before / −0.000069 after**, paired change **−0.001679 [−0.003375,−0.000039]**.
- Bird flower-interaction documentation: **+0.002243 before / +0.002621 after**, *paired change* **+0.000378 [−0.000508,+0.001350]**, includes zero. Importantly, the **after-source-omission bird documentation prediction gain** alone has a positive conditional family-equal 95% interval **[+0.000879,+0.004586]**. It is NOT a verified bird-pollination syndrome.
- Bombus and Diptera paired changes have conditional intervals containing zero.

A caveat on the artifact metadata: the original source-coverage receipt retains a legacy top-level string `status=PROVIDER_COVERAGE_ONLY_NOT_FAMILY_PREDICTION`, even though its explicitly nested `provider_omitted_family_transfer.status=COMPLETED_EXPLORATORY_PROVIDER_OMITTED_FAMILY_TRANSFER` and boolean `provider_family_prediction_rerun_completed=true` show **actual completed family transfer**. Interpret the more specific nested fields, which are supported by matching runner logs and GitHub job success. This is a metadata legacy labeling issue, not a failed biological fit.

This source-label omission is NOT guaranteed independent-publication-source filtering: the dominant `sourceArchiveURI` may aggregate several human research providers. The five GloBI outcomes are compositional and overlap in botanical/documentary origins, and all results are exploratory after prior model selection.

**Updated interpretation:** the original GloBI documentation-composition pattern is partly source-label dependent; a narrower **documented bird-reference composition** prediction signal remains after this broad source-label omission. The result supplies an evidence filter for future true pollinator shift work, not evidence that birds replace bees on islands.

See the [paired six-guild anonymous aggregation](paired_provider_omission_77families_fiveguilds.csv) and the separate [strict local island-pollinator evidence gate](../20261112_pollinator_shift/POLLINATOR_SHIFT_EVIDENCE_AUDIT.md).
