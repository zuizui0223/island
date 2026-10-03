# Tropical taxonomic-depth status-provenance diagnostic — 2026-10-03

## Question

The current H1 taxonomic-depth audit showed a tropical isolation-response vector that
remained supported after genus residualization in the all-observed and WCVP
regional-native-compatible surfaces. This diagnostic asks where that below-genus residual
sits with respect to floristic-status provenance.

It is **not** a test of within-lineage evolution. Endemicity is used only as a status
stratification, and WCVP regional-native compatibility does not establish exact island
nativeness.

## Frozen status strata

| Stratum | rows | islands | species |
| --- | ---: | ---: | ---: |
| source-backed native non-endemic | 118,429 | 425 | 26,200 |
| source-backed endemic | 17,444 | 143 | 17,429 |
| source-backed native, endemism unresolved | 14,658 | 380 | 8,469 |
| WCVP-upgraded from origin-unresolved | 362,789 | 2,281 | 60,117 |

The WCVP-upgraded stratum is exactly the subset that was originally origin-unresolved
and was upgraded only because the focal island's TDWG-L3 unit occurs in the species'
WCVP native range. Known introduced rows are absent.

## Tropical result

| Stratum / evidence | common-support H1 q | post-genus q | interpretation |
| --- | ---: | ---: | --- |
| native non-endemic / all-analysis | 0.003758 | 0.117741 | compatible with genus structuring |
| native non-endemic / Direct-only | 1.85e-5 | 0.075393 | compatible with genus structuring |
| endemic / all-analysis | 0.724571 | 0.211405 | parent H1 not reproduced; no depth inference |
| endemic / Direct-only | not testable | not testable | one retained outcome only |
| source-native, endemism unresolved / all-analysis | 0.012655 | 0.000648 | below-genus residual retained |
| source-native, endemism unresolved / Direct-only | 0.360019 | 0.002023 | beta-binomial parent H1 gate fails; no depth inference |
| WCVP-upgraded origin-unresolved / all-analysis | 9.88e-6 | 0.002147 | below-genus residual retained |
| WCVP-upgraded origin-unresolved / Direct-only | 0.002415 | 0.000603 | below-genus residual retained |

The strict source-backed native-nonendemic tropical response is therefore strongly
present before genus adjustment but is not supported after genus residualization. The
source-backed endemic subset cannot adjudicate the question because the parent H1 is
unsupported or not testable. The most stable below-genus residual occurs in the large
WCVP-upgraded status-unresolved pool.

## Interpretation

The tropical below-genus result from the broader WCVP regional-native-compatible surface
is **floristic-status-provenance sensitive**. It must not be presented as evidence that
native island lineages repeatedly evolved the syndrome below genus.

The strongest defensible inference is:

> Tropical H1 contains a finer-than-genus residual in broad and regionally
> native-compatible assemblages, but the residual is not reproduced in the strict
> source-backed native-nonendemic subset and cannot be evaluated in source-backed
> endemics. The finer-grained component is therefore localized mainly to records whose
> exact island floristic status remains unresolved.

This pattern remains compatible with within-genus species sorting, regional source-pool
structure, imperfect island-level status assignment, persistence filters, or genuine
within-lineage evolution. The current data do not distinguish among them.

## Provenance

- validation workflow run: **37096771148**
- artifact: **11264951254**
- exact current H1 depth implementation:
  `src/island_v2/chapter1_h1_taxonomic_depth_current.py`
- corrected geography: 24 September 2026 baseline
- no-bootstrap point diagnostic; paired block-bootstrap uncertainty for the parent
  four-surface audit remains in
  `results/h1_taxonomic_depth_current_20261003/`
