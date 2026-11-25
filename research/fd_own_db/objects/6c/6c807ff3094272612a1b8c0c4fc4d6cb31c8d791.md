# Supplementary figures and complete numerical audit tables

Source-verified aggregate supplement, 24 November 2026. All figures generated from `supplement_data.json` using `scripts/render_integrated_supplement_figures.mjs`. No new statistical fits or protected species/island identities. These are editable draft vector figures, not final journal plates.

| Item | Contents | Source/interpretation |
|---|---|---|
| [Fig. S1](supplement/figS01_geology_connectivity.svg) | Geological origin and LGM connectivity, including empty joint cell | Main verified evidence ledger; separate contrasts do not identify independent effects |
| [Fig. S2](supplement/figS02_all13_nominal_Q.svg) | All 78 category-vs-corrected nominal Rao-Q cells | Full original crossgate table, 13 traits × 3 regions × 2 statuses; unsupported preserved |
| [Fig. S3](supplement/figS03_animal_trait_sensitivity.svg) | Bird HWI; butterfly wingspan and documented host-family breadth across source masks and thresholds | Runs 45341369223 and 45358140791; recorded resource breadth is not actual local host use |
| [Fig. S4](supplement/figS04_original_null_families.svg) | 996 qualified comparisons across three original permutation families | 336 source-relative joint, 500 form/depth, 160 two-sided IUT; zero adjusted passes |
| [Fig. S5](supplement/figS05_origin_uncertainty.svg) | Native-status coefficient decomposition and unresolved-origin identification | Original runs 45385821139 and 45386194167; 0/460 feasible bounds identify sign |
| [Table S6 / provenance](DATA_SOURCES_AND_REPRODUCIBILITY.md) | Full source manifest, pinned runs and remaining publication gates | Source permissions and primary citations still require reconciliation; no need for a redundant figure |

## Machine-readable supplementary tables

- [Table S2 — all 78 corrected Q/category evaluations](supplement/table_S02_all13_crossgate.csv): unchanged archived source CSV.
- [Table S3 — full butterfly trait-filtering reported panels](supplement/table_S03_butterfly_traits.csv): all reported traits and models, including documented host accounts and unsupported fits; exact log values where available.
- [Table S4a — all 996 qualified original null comparisons](supplement/table_S04_original_null_tests.csv): coefficient, raw and adjusted p values, test family, and available null/leave-group data. BH and Holm are separate families; do not pool these tests as independent.
- [Table S4b — source panel/cohort status](supplement/table_S04_panel_support.csv): all source-relative, form/depth, dual-null and origin-bound panels, including insufficient-support cases.
- [Table S5 — all origin bound outputs](supplement/table_S05_origin_bounds.csv): 460 feasible endpoint-level intervals plus infeasible cap rows; intervals are conditional identification ranges and not sampling intervals.

Reproduction: `node scripts/render_integrated_supplement_figures.mjs --check`. The original Actions logs were parsed and the 996/460 totals verified. See main evidence ledger and protected run links in the source registry. Original source metadata for the native bounds inherited an outdated warning; actual status, source code and endpoint-wise bound outputs control interpretation.

**S1 denominator correction (25 November 2026 JST):** geology contrast is oceanic 124 versus continental 37; LGM is unconnected 347 versus connected 257. Earlier vector labels listed the two counts in reverse contrast order. Only the displayed order changed, not the effect estimates or joint support cells.
