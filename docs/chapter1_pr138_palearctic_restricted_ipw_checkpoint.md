# PR138 within-Palearctic restricted IPW + Aegean stress checkpoint

Status: **robustness / claim-boundary checkpoint; not a new primary estimand**.

Run: `32980009101`  
Artifact: `pr138-palearctic-restricted-ipw-32980009101`  
Digest: `sha256:668cda94b88a46366b3235440c41d7373f1fdcfbce8099a521c7e625ec2f3c98`

This analysis combines three already-predeclared stressors without changing the response definition:

1. external GIFT mainland source-pool subtraction;
2. exact reproductive restrictions (`si_only`, `predominantly_outcrossing`, `si_and_predominantly_outcrossing`);
3. within-Palearctic outcome-blind observation-selection IPW, with and without deletion of the preidentified Aegean/eastern-Mediterranean block `lat12_lon20`.

The response is source-adjusted `attraction_shift = (-large_bee_like + generalized_accessible) / 2`. Propensity models use only score availability plus frozen measured distance, area and climate covariates. Syndrome-score values, fitted ecological effect sizes and p-values are not used in propensity fitting.

## 1. Full Palearctic realm

All 72 full-Palearctic fits have a positive distance slope.

Across four source definitions × two floristic strata:

| Restriction | Unweighted | IPW 0.05 | IPW 0.025 |
|---|---:|---:|---:|
| SI only | 8/8 FDR-supported | 8/8 | 8/8 |
| Predominantly outcrossing | 8/8 | 4/8 | 0/8 |
| SI AND predominantly outcrossing | 4/8 | 8/8 | 8/8 |

Effect ranges remain positive:

- SI only: `+0.064 to +0.095` across all selection specifications;
- predominantly outcrossing: `+0.034 to +0.061`;
- SI AND predominantly outcrossing: `+0.052 to +0.079`.

Thus, within the full Palearctic realm, the attraction/accessibility response is not confined to self-capable taxa and remains especially strong for exact SI taxa. The outcrossing-only subset is more sensitive to aggressive IPW.

## 2. Delete the preidentified Aegean/eastern-Mediterranean block

After deleting `lat12_lon20`, **all 72 fitted slopes remain positive**, but inferential support changes sharply.

Across four source definitions × two floristic strata:

| Restriction | Unweighted | IPW 0.05 | IPW 0.025 |
|---|---:|---:|---:|
| SI only | 8/8 FDR-supported | 0/8 | 0/8 |
| Predominantly outcrossing | 0/8 | 0/8 | 0/8 |
| SI AND predominantly outcrossing | 8/8 | 0/8 | 0/8 |

The weighted estimates attenuate rather than reverse:

- SI only: unweighted `+0.070 to +0.082`; IPW `+0.017 to +0.037`;
- predominantly outcrossing: unweighted `+0.042 to +0.055`; IPW `+0.007 to +0.028`;
- SI AND predominantly outcrossing: unweighted `+0.075 to +0.086`; IPW `+0.013 to +0.041`.

The combined stress also reduces support and effective sample size. Under IPW, minimum effective-sample fractions are roughly `0.27–0.51`, with realized maximum weights `3.38–6.38` after Aegean deletion.

## 3. Interpretation

This result changes the strength of the decomposition claim but not its direction.

Supported:

1. Source-adjusted Palearctic attraction/accessibility reorganization is positive in exact SI and outcrossing subsets.
2. In the full Palearctic realm, exact SI restriction remains robust under both IPW specifications, so measured selfing ability is not a sufficient explanation for the full-realm pattern.
3. The sign of the restricted response does not reverse when Aegean leverage and measured observation selection are stressed simultaneously.

Not supported:

1. The exact-restriction result is **not invariant to the simultaneous Aegean deletion + IPW stress**.
2. It is therefore too strong to describe the SI/outcrossing restriction as decisive proof of a geographically homogeneous, selection-robust attraction pathway across the entire Palearctic realm.
3. The analysis still does not identify Bombus loss, relaxed attraction selection, realized pollinator identity, or causal mediation.

## Updated claim ceiling

The strongest defensible statement is:

> The Palearctic attraction/accessibility shift cannot be reduced to a simple measured-selfing by-product: it is present in self-incompatible and predominantly outcrossing taxa and survives source-pool adjustment. However, the restricted signal is spatially and observation-selection sensitive; after simultaneously removing the influential Aegean/eastern-Mediterranean block and applying IPW, effect estimates remain positive but FDR support is lost.

This makes the main Chapter 1 conclusion **biogeographic contingency of floral/reproductive assemblage responses**, with a strong but leverage-sensitive Palearctic attraction/accessibility branch. It does not justify a uniform Palearctic causal mechanism.
