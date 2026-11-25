# Shared-taxon-exclusive floral dissimilarity: turnover or nestedness?

This analysis partitions **observed** Jaccard dissimilarity of floral-profile
sets after removing shared taxon names in each sampled island pair.
Turnover is replacement of profile identities; nestedness-resultant dissimilarity
reflects unequal profile richness. No new p-values are calculated.

| Scope | Region | Pairs | Median turnover | Median nestedness | Median turnover share | Median relative Jaccard |
|---|---|---:|---:|---:|---:|---:|---:|
| all | northern_midlatitude | 2,000 | 0.610 | 0.127 | 0.830 | 0.445 |
| all | northern_high_latitude | 2,000 | 0.727 | 0.078 | 0.902 | 0.317 |
| all | tropical | 2,000 | 0.638 | 0.114 | 0.850 | 0.413 |
| all | southern_extratropical | 1,618 | 0.650 | 0.106 | 0.861 | 0.403 |
| direct | northern_midlatitude | 2,000 | 0.619 | 0.133 | 0.826 | 0.440 |
| direct | northern_high_latitude | 2,000 | 0.737 | 0.085 | 0.897 | 0.314 |
| direct | tropical | 2,000 | 0.667 | 0.111 | 0.856 | 0.395 |
| direct | southern_extratropical | 1,618 | 0.667 | 0.105 | 0.864 | 0.383 |
| all_single_state | northern_midlatitude | 2,000 | 0.444 | 0.117 | 0.793 | 0.592 |
| all_single_state | northern_high_latitude | 171 | 0.571 | 0.078 | 0.883 | 0.453 |
| all_single_state | tropical | 2,000 | 0.308 | 0.143 | 0.690 | 0.720 |
| all_single_state | southern_extratropical | 391 | 0.375 | 0.138 | 0.724 | 0.656 |
| direct_single_state | northern_midlatitude | 2,000 | 0.444 | 0.117 | 0.796 | 0.590 |
| direct_single_state | northern_high_latitude | 171 | 0.571 | 0.078 | 0.883 | 0.450 |
| direct_single_state | tropical | 2,000 | 0.333 | 0.137 | 0.720 | 0.692 |
| direct_single_state | southern_extratropical | 391 | 0.375 | 0.147 | 0.724 | 0.659 |

Per pair: a=shared, b=only A, c=only B. Jaccard beta=(b+c)/(a+b+c),
turnover=2min(b,c)/(a+2min(b,c)); nestedness-resultant=Jaccard beta−turnover.
Turnover share is omitted for identical profile sets (beta=0).
Relative Jaccard=observed J / [min(KA,KB)/max(KA,KB)] standardizes the
**attainable maximum only**; it does not fix profile prevalence or genus composition.
Component medians do **not** add to the median total; pair-level identities do.
Pairs sharing an island are dependent. Equal-island means are separately in CSV.
The [whole-island fixed-margin reference](../floral_space_fixed_margins_20261107/README.md)
includes shared taxa and is NOT the null distribution of this pair-exclusive result.
