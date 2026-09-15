# P1 paired spatial-block uncertainty result — 2026-09-15

## Status

Post-baseline robustness extension. This result does not reclassify the historical H3 gate.

Pre-result protocol commit: `39a238c996e395031d4ba94456f826fba489cf78`.
Implementation/run head: `2ad01d86d9803810f00187fcdef43b5f344356ff`.

Canonical run:

- workflow run `34939113182`;
- artifact ID `10383754545`;
- artifact `chapter1-p1ab-assembly-depth-34939113182`;
- digest `sha256:822c69e9a1be4eb2c340391c1fb84dfba91b45b8a765549bb8ad1ec6379ef480`.

Pinned biological input remained PR142 run `34232450884`, artifact `10058653212`, digest `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.

## P1a — exact paired-support reconstruction

The frozen taxonomic decomposition was independently reconstructed from the pinned artifact and tracked taxonomy.

- direct-only maximum absolute reconstruction difference: `1.11e-16`;
- all-analysis maximum absolute reconstruction difference: `1.11e-16`;
- the family and genus stages use common observed species;
- the same island × response rows and `n_species` information weights feed observed, family-adjusted and genus-adjusted stages;
- family/genus labels are grouping variables, not trait imputations.

Therefore stage-specific sample loss is not a viable explanation for the frozen attenuation within this decomposition.

## P1b — paired spatial-block uncertainty

The pre-result extension resampled the existing `spatial_block` unit 2,000 times with seed `20260915`. Within each draw the same resampled blocks were used for observed, family-adjusted and genus-adjusted stages, and the frozen point-estimate regression structure was refitted.

### Direct-only primary safeguard

Across the eight source-mode × floristic-stratum Palearctic profiles:

- observed total genus attenuation was large: `0.7881–0.8024` of the observed two-axis vector;
- the 95% lower bounds for total genus attenuation remained positive, ranging approximately `0.2765–0.4021`;
- observed family-to-genus *additional* attenuation was `0.5028–0.5467` of the observed vector;
- however its 95% lower bounds ranged approximately `-0.2095` to `-0.1193`;
- therefore **0/8** direct-only profiles had a 95% interval for the family-to-genus incremental attenuation wholly above zero.

The conditional genus attenuation relative to the family-adjusted remainder also had wide intervals that crossed zero in all eight direct-only profiles.

### All-analysis complementary scope

Only `1/8` all-analysis profiles had the family-to-genus incremental attenuation interval wholly above zero. This scope is complementary, not independent replication.

## Interpretation

Two statements now have to be separated.

**Supported:** the broad Palearctic response is strongly sensitive to taxonomic composition, and observed genus adjustment removes a large fraction of its vector magnitude. This is not produced by changing the analysed species between taxonomic stages.

**Not established by paired spatial uncertainty:** the *incremental* attenuation is precisely localized at the family-to-genus transition. The point-estimate pattern remains concentrated there, but its spatial-block-resampled difference from the family stage is too uncertain for a strong quantitative localization claim.

P1c asks a different question and remains informative: whether true genus boundaries attenuate more than arbitrary within-family partitions of matched complexity. A positive P1c result could establish genus-specific structural information even though the exact family-to-genus attenuation magnitude remains spatially imprecise.

## Manuscript consequence

Until P1c is known, avoid the strongest wording `the syndrome has a precisely measured family-to-genus assembly depth`.

Safe current wording:

> The strongest Palearctic floral-island response is strongly structured by taxonomic composition, with the observed point-estimate attenuation concentrated between family and genus; paired spatial-block resampling shows that the incremental family-to-genus attenuation itself is imprecisely localized.

Do not interpret this result as causal mediation, dispersal-only filtering, absence of within-lineage evolution, or evidence for/against pollinator causation.
