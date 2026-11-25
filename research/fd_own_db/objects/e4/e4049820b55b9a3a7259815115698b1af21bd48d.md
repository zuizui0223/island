# Source-holdout feasibility audit for genus × mainland floral N2 (2026-11-09)

**Result type:** independently retrieved and SHA-verified rights-restricted source
checkpoint, non-identifying source-coverage summaries only. This is a
*methodological support audit*, not a confirmed ecological null rejection.
[Original frozen island PR #321](https://github.com/zuizui0223/island/pull/321)
is the reference. Source checkpoint names, island names, source entity
identities, per-taxon floral profiles, and exact geographic joins remain private.

## Verified original source-reconstruction gate

The archived PR321 raw Direct ledger, axis coverage, All/Direct trait
snapshots, and island-status flora each match the original upstream frozen
SHA256 digests documented in the source runner. Reconstructed intact
colour-state-set × form-state-set combinations reproduce exactly:

| Field | All | Direct |
|---|---:|---:|
| Complete and lexically eligible taxon-name units | **16,739** | **8,455** |
| Unique categorical colour/form profile types | **496** | **467** |
| Islands with at least 50 recorded complete name units | **1,533** | **1,443** |
| GIFT-native mainland *distinct donor taxa* with a matching floral profile | **9,017** | **6,562** |

The original GIFT native-mainland source archive contains **748 source
entities** and **1,666,186 entity × taxon rows** (not 1.66 million
distinct trait donors). Cross-island sharing of the same accepted
taxon names and data reporting is a core inferential limit.

## Why all-island full-name holdout is underidentified

If *every* taxon recorded on *any* eligible island is simultaneously held
out from *every* GIFT native-mainland entity, the remaining known-profile
taxon donors are only:

| Scope | Donors remaining | Fraction of original known-profile donors |
|---|---:|---:|
| All | **228 / 9,017** | **2.53%** |
| Direct | **172 / 6,562** | **2.62%** |

Therefore a single global train/test split with all eligible islands
reserved as tests cannot identify a useful native-source genus trait
distribution in this shared-taxon data. Merely producing large
absolute sample counts would overstate independent donor support.

## More defensible geography-block-disjoint five-fold donor audit

All 1,533 All-eligible islands are allocated to five **nonoverlapping
spatial-block sets** by an outcome-blind greedy count balancing rule
(approximately 306–307 target islands per fold). Direct inherits the
same heldout spatial-block assignments. **Every accepted name** recorded
on *any* heldout island in a fold is removed from EVERY mainland source
entity before training, including from global and local trait donors.

Observed non-identifying approximate ranges over the five folds:

| Coverage diagnostic | All | Direct |
|---|---:|---:|
| Known-profile donor taxa remaining in fold training | **3,697–5,272** | **2,103–3,527** |
| Focal name-carriers with independent **global genus** donor | **64.1–72.7%** | **53.8–63.9%** |
| Focal name-carriers with independent **five-nearest native-mainland source genus** donor | **16.6–31.2%** | **12.9–25.3%** |

The donor coverage is **much lower** near putative mainland sources
than in the global genus pool. A smoothed local genus posterior can
be calculated, but local source location is only weakly informed
by independent taxon traits on many islands.

[Reproducible audit code](../../src/fd/n2_source_coverage_audit.py)
implements this donor-grain and spatial-block exclusion with fail-closed
checks and aggregate outputs. It is private-input-only;
GitHub Actions runs synthetic tests but never downloads or redistributes
the rights-restricted plant trait/flora/GIFT source tables.

## Pilot predictive checks (exploratory, not full N2 origin inference)

On an outcome-blind geographically distributed pilot, **118 All islands**
and **116 Direct islands** were evaluated with per-geographic-fold
taxon-name exclusion. A global genus-informed donor model improved
exact-profile log loss relative to the global no-genus prior;
however, adding the geographically nearest native mainland genus
donor data produced almost no further *mean exact-profile* benefit:

| Mean heldout name-level log loss | All pilot (118 islands) | Direct pilot (116 islands) |
|---|---:|---:|
| Global non-genus trait prior | 5.4524 | 5.6191 |
| Heldout **global genus** predictor | 4.8083 | 5.2545 |
| Heldout **local native-mainland genus** predictor | 4.8078 | 5.2543 |
| Increment from local source over global genus | **+0.0005** | **+0.0002** |

The sampled **model-based** nearest-profile reference still yields
negative distances relative to its expected floral combinations:
median fixed-observed-richness SES about **−3.20 All / −3.42 Direct**,
and per-name trait-simulation equal-30 SES about **−1.53 / −1.88**.
These are **not** evidence of evolutionary selection or geological
origin dependence. Low independent *local* genus donor support,
observational taxonomy, reporting bias, the independent name-by-name
simulation assumption, and the fact that the local source predictor
hardly improves exact taxon prediction are serious model adequacy limits.

A pilot isolation-gradient regression on the N2 fixed-richness SES,
adjusting broad region and log island area with geographic block
cluster-robust errors, gave non-supported isolation slopes
(All n=118, p≈0.90; Direct n=116, p≈0.56).
These post-hoc pilot regressions are **not** confirmatory independent
biogeographic discoveries.

## Claims boundary and next gate

**Do not claim:** full-sample globally independent N2 validation;
adaptive island floral convergence; confirmed geological origin
selection; or independent mainland lineage ancestry.

The useful inference is a separation:
(1) *recorded genus composition* is informative about floristic
phenotypes; (2) *local mainland donor location* has not yet passed
a convincing independent heldout name-level predictive gate;
(3) reported narrow within-island profile neighborhoods can
also result from spatial species/trait composition and model
misspecification. These are different estimands.

A future source-supported ecological claim requires training/validation
on geographically or phylogenetically independent islands, richer
**independent** native-mainland flower-trait donors, detection-aware
record quality, and sufficient geographic replication. Strict original
37/10 matched geological N2 remains **unevaluated** here; this pilot
does not replace it. No new geological p-value has been issued.

All published entries are aggregate counts or summary scores.
Source archives, island IDs, source names, exact species profiles and
geographic crosswalks are not included in this repository.
