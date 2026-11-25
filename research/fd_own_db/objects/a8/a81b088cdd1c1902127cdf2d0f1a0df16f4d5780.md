# Original GIFT stage-2 baseline reversal: regional versus external genus-blind priors (19 Nov 2026)

This is an **independently re-derived diagnostic of an EXISTING original dataset**, not a new ecological source sample or a new model fit. It reopens the exact original source-anonymous ZIP from [island Actions #44513833517](https://github.com/zuizui0223/island/actions/runs/44513833517), artifact #14525321320, ZIP SHA256 **9cc97d7e572c043d586bda4a661ba69573fc2049bdd35e2c9f38de35dd7bf725**, reads all 12 original GIFT panels, and evaluates the reported proper-loss identities. Original 4,505 botanical island parent, two separately prespecified TRUE Direct endpoints, and all-name-disjoint source model remain unchanged.

## Identity and estimand

For the EXACT SAME original Direct focal plant carrier set and geographical source grouping, let the Brier/log loss of target-region genus-BLIND floral-category reference be L_local, other-region genus-BLIND reference L_external, and two-stage genus→floral forecast L_model. The original source anonymized report already stores (L_local−L_model) and (L_external−L_model). Their DIFFERENCE is therefore **exactly** L_local−L_external, an estimator of **baseline portability/calibration** rather than genus-composition effects.

This identity was separately checked against all FOUR original source genus models (region-pooled genus source, geographically nearest-40, environment-nearest-40, environment+ranked candidate-source nearest-40), each of the TWO proper losses and both physical 10° block-equal and first candidate source-group-equal means, across the six original all-recorded region×endpoint panels. **Maximum algebraic deviation 1.11e−16** in floating point. This proves arithmetic consistency, **not** statistical significance of the direct difference nor independent biological replication.

## Observed exact group-equal mean source loss differences

Positive L_local−L_external means **external-region categorical frequency baseline itself scores better** than the nominally local-region category reference. Negative means the target-region categorical reference scores better. This is a contrast between Brier/log loss of two blind reference models, not a contrast of observed botanical trait frequencies.

| Original GIFT macroregion | Genuine Direct endpoint | n evaluable islands | 10° block-equal local minus external Brier | first-source-group-equal Brier | 10° block-equal local minus external log loss |
| --- | --- | ---: | ---: | ---: | ---: |
| Northern midlatitudes | Floral form | 1,661 | **−0.103769** | −0.084899 | −0.391010 |
| Northern midlatitudes | Colour × form | 1,463 | **−0.049768** | −0.043667 | −0.455833 |
| Tropical | Floral form | 771 | −0.003242 | −0.008786 | −0.009717 |
| Tropical | Colour × form | 582 | −0.005398 | −0.006208 | −0.068525 |
| Southern extratropical | Floral form | 157 | **+0.057487** | **+0.048160** | +0.226373 |
| Southern extratropical | Colour × form | 140 | **+0.050341** | **+0.046233** | +0.171959 |

Numerically, the original southern stage-2 gains versus the *strongly mismatched local reference* are **+0.060389 Brier** for floral form and **+0.051092** for true same-species colour×form. The original model versus the *better-performing external blind reference* is only **+0.002902** and **+0.000751**, respectively. Thus the large southern comparison is predominantly a **local-versus-external blind baseline difference**, and the much smaller remaining model-over-external score must not be conflated with a large ecological taxonomic/genus source effect. This is just a signed score identity: do NOT interpret the ratio as a causal variance or adaptation component.

The direction of local-versus-external blind advantage differs between the northern midlatitudes and southern extratropics, under both independently defined source grouping averages. This is a striking source-regional predictive calibration difference. But there is **NO directly paired bootstrap interval for the local minus external baseline difference in the anonymously aggregated file**; subtracting the separately stored bootstrap marginal interval endpoints is invalid. Therefore the table contains numerical group-equal point means, **NOT a new confirmatory p-value or source-based interval**.

## Research implication

The verified southern locally benchmarked two-stage predictive positive remains numerically and source-proper-score valid. What changes is the **scientific explanation**: it cannot be used as standalone evidence that geographical source genus assembly caused functional floral diversification, because an outside-region genus-blind prior nearly matches its Brier score on those same floral outcomes. Regional trait category frequency source/annotation differences and heterogeneous source data coverage can drive the apparent gain over a poorly transferring local reference. Meanwhile northern results favour the local reference strongly, so global prediction is not a single-shape phenomenon.

Next original-data test should directly bootstrap the **paired local-blind minus external-blind loss** by physical block and first-source group on the same Direct carriers, and check category-ontology and annotation differences before attributing regional reversal to ecological guild dynamics. The ongoing independent [expanded external-congener sensitivity](https://github.com/zuizui0223/island/pull/332) and [all-13 Direct endpoints](https://github.com/zuizui0223/island/pull/333) may further distinguish model-information effects from local/external baseline mismatch. Neither has a new verified outcome in this note.

Related [original stage-2 source-verified report](../20261119_stage2_name_disjoint_floral_forecast/README.md).
