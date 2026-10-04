# Final H1 — seven traits, finite-cluster t, no Holm

### H1: seven separate trait responses

The primary analysis is broad contemporary island flora with All evidence. Each of seven binary traits is fitted separately in each of four regions by beta-binomial logit regression, adjusting for standardized corrected log isolation, log island area and climate PC1–4. Reproduction, colour and structure are labels, not aggregate scores. Spatial-block sandwich standard errors use a finite-cluster t reference with G−1 degrees of freedom and pointwise 95% intervals. All individual p values are two-sided and unadjusted; no Holm correction is applied. The final rule was user-approved after comparison of alternatives on 4 October 2026 and is retrospective, not preregistered. Individual significance does not establish family-wise significance or formal differences between regional slopes.

WCVP regional-native-compatible flora is the sole active origin sensitivity; Direct-only is a separate evidence-quality sensitivity. Both use the identical model and inference rule. WCVP compatibility retains existing source-native records and upgrades unresolved records only under accepted TDWG-L3 native-range compatibility; introduced records are not overwritten. This establishes regional compatibility rather than exact focal-island nativity. H2–H4 and their previously defined functional covariates are unchanged.

### H1: recurrent reproductive responses with regional floral differences

All 112 fits converged. Primary All evidence supports 17 of 28 individual associations at nominal two-sided P < 0.05, the same set of significant traits as the original poster's normal approximation. Self-compatibility increases with isolation in all four regions (P = 0.00569, 0.00671, 0.01812 and 0.01123 in northern mid-latitude, northern high-latitude, tropical and southern extratropical floras, respectively). Selfing mating system increases in the first three regions. Autonomous selfing increases in northern high latitudes and the tropics. Generalized form increases in northern mid-latitudes, the tropics and southern extratropics; actinomorphy increases in northern high latitudes and the tropics. Shallow/open tubes increase in northern high latitudes but decrease in southern extratropics. Plain colour increases in southern extratropics.

In WCVP All, 10 of 28 individual associations meet the same nominal threshold. Increasing self-compatibility remains supported in all four regions. Northern-high autonomous selfing, actinomorphy and shallow/open tubes remain positive and supported; tropical actinomorphy remains supported. Plain colour increases in both tropical and southern extratropical WCVP floras. Several broad-flora associations, including the southern shallow/open-tube decrease, no longer meet the threshold. This can reflect changed composition, precision and coverage, rather than proving an introduced-species mechanism. Direct-only yields 15/28 supported associations in broad flora and 11/28 in WCVP flora. Complete coefficients, intervals and unadjusted p values are reported in results/h1_final_traitwise_t_20261004/traitwise_results.csv.

## WCVP coverage

Before trait/covariate filtering: 513,320 island-species records, 2,372 islands, 90,693 species. All fits exceed the existing 50-island threshold, but this is not a power guarantee. Some variance contributions are concentrated in few blocks. WCVP is regional native compatibility, not exact island nativity.

## Supported All associations (nominal, unadjusted)

| Flora | Region | Trait | Estimate | 95% CI | P |
|---|---|---|---:|---|---:|
| broad | northern_midlatitude | self_compatibility | 0.05606 | [0.01678, 0.09533] | 0.00568533 |
| broad | northern_midlatitude | selfing_mating_system | 0.04634 | [0.01807, 0.07460] | 0.00161561 |
| broad | northern_midlatitude | generalized_form | 0.04613 | [0.01690, 0.07536] | 0.00235126 |
| broad | northern_high_latitude | self_compatibility | 0.11691 | [0.03360, 0.20023] | 0.00670944 |
| broad | northern_high_latitude | selfing_mating_system | 0.06963 | [0.01708, 0.12218] | 0.0102843 |
| broad | northern_high_latitude | autonomous_selfing | 0.06427 | [0.00057, 0.12798] | 0.0480593 |
| broad | northern_high_latitude | actinomorphic_symmetry | 0.24086 | [0.07494, 0.40677] | 0.00514165 |
| broad | northern_high_latitude | shallow_open_tube | 0.36877 | [0.00916, 0.72838] | 0.0446567 |
| broad | tropical | self_compatibility | 0.11111 | [0.01934, 0.20287] | 0.0181166 |
| broad | tropical | selfing_mating_system | 0.32823 | [0.08676, 0.56970] | 0.00822599 |
| broad | tropical | autonomous_selfing | 0.09778 | [0.00082, 0.19474] | 0.0481392 |
| broad | tropical | generalized_form | 0.11706 | [0.03949, 0.19463] | 0.00345866 |
| broad | tropical | actinomorphic_symmetry | 0.10990 | [0.03961, 0.18018] | 0.00247323 |
| broad | southern_extratropical | self_compatibility | 0.17606 | [0.04221, 0.30992] | 0.0112335 |
| broad | southern_extratropical | plain_colour | 0.07419 | [0.03763, 0.11074] | 0.000199803 |
| broad | southern_extratropical | generalized_form | 0.10792 | [0.02062, 0.19522] | 0.0167883 |
| broad | southern_extratropical | shallow_open_tube | -0.21660 | [-0.34046, -0.09275] | 0.00112082 |
| wcvp | northern_midlatitude | self_compatibility | 0.07865 | [0.02355, 0.13376] | 0.00571959 |
| wcvp | northern_high_latitude | self_compatibility | 0.14424 | [0.03698, 0.25150] | 0.00928864 |
| wcvp | northern_high_latitude | autonomous_selfing | 0.11263 | [0.00908, 0.21619] | 0.0335792 |
| wcvp | northern_high_latitude | actinomorphic_symmetry | 0.25616 | [0.02994, 0.48237] | 0.0271929 |
| wcvp | northern_high_latitude | shallow_open_tube | 0.44745 | [0.06866, 0.82624] | 0.0217817 |
| wcvp | tropical | self_compatibility | 0.13510 | [0.03330, 0.23690] | 0.00981011 |
| wcvp | tropical | plain_colour | 0.11908 | [0.04961, 0.18855] | 0.000961242 |
| wcvp | tropical | actinomorphic_symmetry | 0.11800 | [0.01915, 0.21684] | 0.019785 |
| wcvp | southern_extratropical | self_compatibility | 0.39945 | [0.17301, 0.62588] | 0.00101248 |
| wcvp | southern_extratropical | plain_colour | 0.16214 | [0.09373, 0.23055] | 2.59267e-05 |

## Reproducibility

All 112 fits converged. Exact input hashes and software versions: manifest.json. Full table: traitwise_results.csv. Figure: [PDF](traitwise_H1.pdf), [PNG](traitwise_H1.png). Replay: [instructions](../../docs/REPLAY_H1_FINAL_TRAITWISE_T_20261004.md). H2–H4 are unchanged; the poster PPTX itself has not yet been updated. Earlier Holm and restored-normal tables remain historical, not active.
