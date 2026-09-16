# Chapter 1 trait-first family response audit — 2026-09-15

Status: **completed; pre-model family averaging does not recover the atomic response difference**.

## Provenance

- workflow run: `34964919955`
- artifact: `chapter1-trait-first-family-34964919955`
- artifact ID: `10394986981`
- digest: `sha256:14d8f67d1105251d1932f1f2a880cb07545fd8b376b076f9eb85a7fb49c76cf8`

## Question

Previous diagnostics established that:

- the six atomic North--Tropical vector is supported;
- a post-estimation two-family linear contrast of the six logit-scale slopes is also supported;
- the historical species-first two-axis concordance score is not supported;
- relaxing its `minimum_informative_traits` gate from 2 to 1 does not restore support.

This audit asks whether the problem is specifically species-first averaging. It instead aggregates **trait first**: each island gets one proportion for each atomic trait, then those trait proportions are combined into the two historical response families using the predeclared weights.

Historical weights retained:

- accessibility/generalisation: generalized form = 1, actinomorphic symmetry = 0.75, shallow/open tube = 1;
- reproductive assurance: self compatibility = 1, selfing mating system = 1, autonomous selfing = 1.

The resulting two continuous family scores are then analysed with the same distance, area, climate and spatial-block robust linear design used for continuous response summaries.

## Result

### All-analysis evidence

Using every available family component:

- **3,626 islands**, 187 blocks;
- joint two-family p = **0.50765**.

Requiring all three component traits in each family:

- **2,982 islands**, 178 blocks;
- joint p = **0.11072**.

### High/Medium direct evidence

Using every available family component:

- **3,569 islands**, 187 blocks;
- joint p = **0.25947**.

Requiring all three components:

- **2,886 islands**, 178 blocks;
- joint p = **0.06374**.

For the all-analysis / available-components case, the family interaction estimates are small on the raw family-score scale:

- accessibility/generalisation: `+0.01070`, SE `0.01315`, p=`0.4158`;
- reproductive assurance: `+0.00855`, SE `0.01200`, p=`0.4760`.

## Interpretation

Species-first averaging is **not** the sole cause of the discrepancy. Even trait-first averaging on the raw response scale weakens the North--Tropical contrast.

Combined with the response-compression contrast audit, this localizes the issue to **when and on what scale the response is compressed**:

1. fit each atomic count response with its own beta-binomial/logit model, then contrast family slopes -> supported;
2. average atomic proportions into family scores before modelling -> unsupported;
3. average species-level concordances into family scores before modelling -> unsupported.

Therefore the defensible expanded analysis should keep the six atomic outcomes through model fitting and derive accessibility/generalisation and reproductive-assurance family summaries **after estimation as planned contrasts on the model scale**, rather than treating a pre-averaged two-axis score as the primary response.

This is not an argument for six unrelated traits. The two biological families remain useful, but statistically they should be contrasts of atomic estimands, not compressed island-level response variables.
