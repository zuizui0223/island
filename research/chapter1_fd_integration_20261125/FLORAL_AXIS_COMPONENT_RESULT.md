# Floral-axis decomposition: what drives the tropical response?

[Successful run 45590946804](https://github.com/zuizui0223/island/actions/runs/45590946804) · [full verified aggregate panels](FLORAL_AXIS_COMPONENT_RESULT.json) · [candidate figure](FLORAL_AXIS_COMPONENT_FIGURE.svg).

The same H1 and H2 designs and island masks were reused. All four previously fitted difference slopes were reproduced by the runner. The outcome difference equals the deep-group radial prevalence slope minus the shallow/open-group radial prevalence slope. Units: absolute proportion per unit of the frozen H1 log-distance transform; not standardized effects or temporal changes.

| Model | Islands | Radial among deep | Radial among shallow/open | Difference | Joint genus/geography 2.5–97.5% ranges, deep / shallow / difference |
|---|---:|---:|---:|---:|---|
| H1: no reproductive adjustment | 257 | +0.01988 | -0.00774 | +0.02762 | [+0.00532, +0.03771] / [-0.02484, +0.00498] / [+0.00741, +0.04957] |
| H2: whole-flora selfing_core | 244 | +0.01516 | -0.00451 | +0.01967 | [+0.00091, +0.03279] / [-0.02272, +0.00979] / [-0.00123, +0.04565] |
| H2: pair-carrier SC | 188 | +0.01837 | -0.00626 | +0.02463 | [+0.00431, +0.03470] / [-0.02165, +0.00852] / [+0.00273, +0.04588] |
| H2: pair-carrier autonomy | 106 | +0.01382 | -0.00707 | +0.02089 | [-0.00067, +0.03020] / [-0.02272, +0.00759] / [+0.00087, +0.04294] |

H1 (257 islands, 66 spatial blocks, 718 distinct paired-trait species): deep +0.019876 (95% spatial-cluster interval [0.009232,0.030521], Holm2 p=.000811); shallow/open −0.007744 ([−0.014161,−0.001327], Holm2 p=.01879). The exact contrast is +0.027620. The ratio of the deep-group coefficient to the difference is about 72%, purely arithmetic, not a causal explained fraction. The deep component alone retains positive joint-weight sensitivity [0.005322,0.037706], whereas shallow/open crosses zero [−0.024844,0.004982]. Thus the better supported ecological description is an isolation-associated **increase in radial-symmetry prevalence among recorded deep-tubed tropical species**; opposite shallow/open response is a pointwise result, not a robust second process.

After conditioning on the original whole-flora selfing_core, the deep coefficient remains +0.015161 with joint sensitivity [0.000913,0.032786], though the overall difference's joint range crosses zero. The same-pair SC adjustment retains the deep component's positive joint range [0.004314,0.034698]. The autonomy-adjusted deep range crosses zero; its overall difference remains positive. These component-conditioned results are descriptive, not another Holm-controlled family or independent confirmations. Joint ranges are sensitivity distributions from 499 positive genus weights and geographical block resampling, not calibrated CIs. All original scope/grade caveats apply.

**Editorial decision:** use this as the single focused exploratory H1–H2 addition. The prior joint-Q/forecast analyses remain supplementary; do not add another paper or change the H1–H4 definitions. No claim of named pollinator replacement, adaptive evolution, biological mechanism or robust region interaction. H4 four-cell gate failed (global 32 species, tropical 8), so do not fit, reinterpret, or replace the existing H4. High-only and native-only support remain unresolved. All frozen GEB files stay untouched.
