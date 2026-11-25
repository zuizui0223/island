# H2 paired floral-axis conditional results — source-backed

[Successful run 45584985080](https://github.com/zuizui0223/island/actions/runs/45584985080); [full 24 anonymous panels](FLORAL_AXIS_H2_RESULT.json). Original Direct audit and morphology/geometry pins passed; min2 score covers 2651 species and complete3 covers 564. Only 124 of the 1112 morphology-pair species have the min2 score. All results exploratory extensions of original H2.

| Score carriers / adjustment | Islands | Matched baseline | Conditional slope | Conditional 95% interval | Holm p |
|---|---:|---:|---:|---|---:|
| whole_flora / min2 | 244 | 0.027956 | 0.019675 | [0.008605, 0.030745] | 0.003613 |
| whole_flora / complete3 | 72 | 0.011934 | 0.012133 | [0.000545, 0.023721] | 0.1220 |
| whole_flora / SC | 257 | 0.027620 | 0.023562 | [0.011692, 0.035433] | 0.001115 |
| whole_flora / autonomy | 257 | 0.027620 | 0.020428 | [0.008935, 0.031921] | 0.003613 |
| pair_carriers / min2 | 62 | 0.006893 | 0.005608 | [-0.001118, 0.012335] | 0.1979 |
| pair_carriers / SC | 188 | 0.028888 | 0.024633 | [0.014811, 0.034455] | 0.00004177 |
| pair_carriers / autonomy | 106 | 0.022256 | 0.020894 | [0.012568, 0.029221] | 0.00005184 |

Whole-flora min2, SC and autonomy, and morphology-pair SC/autonomy subroster scores meet the 8-slot Holm gate. Pair-carrier means are computed on the known subset; the morphology response still uses all paired morphology carriers. No claim of a within-individual interaction or causal independence. High-only fits have no supported positive slope; every literal-native panel is below support. Complete3 same-pair support is two islands. The whole-flora min2 conditional coefficient is 70.4% of its own matched baseline; no mediation fraction is identified.

Next diagnostic: replay all 24 original panels on exact masks and verify slopes to numerical precision. In every estimable panel, assign one exponential weight per accepted-name genus, reused for every appearance of both morphological and reproductive traits; unparsed names each receive separate weights. Recalculate both morphology conditional difference and reproductive mean; keep original coverage and covariates fixed. Refit 499 geography-only, genus-only and joint genus/geography reweightings, preserving paired raw and conditioned estimates. These are sensitivity distributions, not calibrated intervals or new p values. No new endpoint search, no H1 promotion, no frozen submission changes. 
