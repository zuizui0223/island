# Current H1 taxonomic representation-depth audit — 2026-10-03

## Purpose

This audit asks whether the **current seven-response H1** is explained by taxonomic
composition at family or genus level. It is not a causal assembly or evolutionary test.

The audit uses the corrected 24 September 2026 geography and the exact current H1
response definition. For each atomic outcome, trait-resolved species are retained only
when their family and genus each contain at least one other scored species. Family and
genus expectations are leave-one-species-out means, so taxonomy never fills missing traits.

Before residualization, the current beta-binomial H1 is re-fit on this exact common
species support. The common-support H1 remains FDR-supported in all four regions in all
four analysis surfaces (all-observed / WCVP regional-native-compatible × all-analysis /
Direct-only). Taxonomic attenuation is therefore not caused simply by dropping
taxonomically isolated species.

## Main result

The strongest reproducible representation-depth result is tropical.

The seven-response tropical vector remains FDR-supported after genus residualization in
all four analysis surfaces:

| Flora/evidence | post-genus q |
| --- | ---: |
| all-observed, all-analysis | 0.00733 |
| all-observed, Direct-only | 0.00559 |
| regional-native-compatible, all-analysis | 0.00144 |
| regional-native-compatible, Direct-only | 0.000311 |

Self-compatibility is individually positive and supported after genus residualization in
all four surfaces (P = 0.00293, 0.00355, 9.93e-5 and 7.09e-5). Plain colour and
generalized form become additional supported post-genus components in both
regional-native-compatible surfaces.

Paired 300-draw spatial-block bootstrap shows that genus structure still contributes to
the tropical vector in the all-analysis scopes: incremental family-to-genus attenuation
has median 0.458 (95% interval 0.212–0.662) in all-observed flora and 0.341
(0.150–0.552) in regional-native-compatible flora. The corresponding Direct-only
intervals cross zero. Thus taxonomy matters, but it does not exhaust the tropical H1
response.

Northern high latitudes show substantial family-level attenuation: the family attenuation
interval is positive in all-observed all-analysis (0.300–0.702), all-observed Direct-only
(0.181–0.573) and regional-native-compatible all-analysis (0.156–0.670), while the
regional-native-compatible Direct-only lower bound is -0.0019. The additional
family-to-genus attenuation is imprecise in every northern-high surface, so a sharp genus
breakpoint is not identified.

Northern mid-latitude depth is model/evidence sensitive: the equal-island observed stage
does not reproduce H1 in two of four surfaces. Southern extratropical depth is strongly
floristic-status sensitive: genus residuals remain supported in all-observed flora but
not in regional-native-compatible flora, while bootstrap attenuation intervals do not
support a precise native family/genus breakpoint.

## Interpretation boundary

The tropical result rules out a simple explanation in which the current H1 pattern is
only a consequence of which genera are present. It does **not** demonstrate repeated
within-lineage evolution. A post-genus residual can still reflect within-genus species
sorting, unmeasured source-pool structure, persistence filters, introductions not removed
by the broad flora surface, or genuine evolutionary change.

The defensible synthesis is therefore:

> A recurrent functional island response can be represented at different hierarchical
> depths. In the tropics, a substantial H1 residual persists below genus composition;
> in northern high latitudes, lineage composition contributes strongly, although the
> exact family-to-genus breakpoint is uncertain.

## Provenance

- validation workflow run: **37095612710**
- artifact: **11264117016**
- artifact digest: `sha256:d2b46479d6eb40457157490b667660bf9728361faae0ab7321ec1cc8d6d1db5b`
- WCVP-compatible flora replay: 362,789 unresolved rows upgraded; 513,320 regional-native rows
- bootstrap: 300 paired spatial-block draws per context and analysis surface

Compact committed outputs:
- `classification_summary.csv`
- `bootstrap_summary.csv`
- `tropical_post_genus_slopes.csv`

Full stage-level slopes, omnibus tests and bootstrap draws remain reproducible through
`src/island_v2/chapter1_h1_taxonomic_depth_current.py` and
`.github/workflows/validate-chapter1-h1-taxonomic-depth-current.yml`.
