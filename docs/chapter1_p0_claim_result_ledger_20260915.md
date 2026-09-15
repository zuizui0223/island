# P0 immutable claim-to-result ledger — 2026-09-15

Status: **P0 central-artifact reconciliation complete for the current headline chain.**

Purpose: identify the immutable empirical source for each manuscript-level claim before any P1 post-baseline robustness analysis is implemented. This file does not recalculate a biological result.

## Live artifact verification performed 2026-09-15

The following workflow artifacts were queried directly from GitHub Actions and were live (`expired=false`) with the expected IDs and digests:

| role | run | artifact ID | digest |
|---|---:|---:|---|
| final progressive H1/H2/H3 input | 34232450884 | 10058653212 | `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999` |
| effect-size / attenuation synthesis | 34796763611 | 10330230959 | `sha256:a2db683b6d46ca09d4eae6d9a3bfcaf85a4f172b5589227958bc388df8edba60` |
| historical formal V5 MNAR run | 33141250895 | 9674100317 | `sha256:273cf037611f041e540d8e4f845d7316999623d5da6ef5958c20ba50b9788663` |
| N1 independent pollinator-channel test | 34750351009 | 10315671249 | `sha256:01d938a1f42b7f64634232a2fcb76c68a60183824684c5f135f567107a4c5d85` |
| GloBI source-breadth V2 | 34790842133 | 10328696801 | `sha256:475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254` |
| V6 species-detection tipping | 34800498716 | 10331282464 | `sha256:095506a214143a2312bcb1d115644ab21d2d130b7813261019bc386206fd3849` |
| observed response geometry | 34765070227 | 10320630085 | `sha256:5ac34231b0cfa4acc1d536ab798a574a3cbcdd9423b6f0f2995195aab78a2f23` |
| H5c observed biotic-vs-wind specificity | 34803837463 | 10332871066 | `sha256:fb8927bcf2b6e1c762b3ecbe59c42360faac679cf37b1ad0ed31386631a0f405` |
| H5d distributed-threshold identifiability | 34803545574 | 10332671123 | `sha256:13ea1c3cd0066950d37fa3399f771d41ea2755c9746d02ac18594e224d895539` |

The figure artifacts are presentation surfaces and are not substituted for the empirical sources above.

## Claim ledger

### C1 — one universal floral island syndrome is not recovered

**Manuscript role:** first biological result; this is the failure of the classical one-direction global expectation, not proof that every region differs.

**Immutable empirical source:** final progressive analysis run `34232450884`, artifact `10058653212`.

**Primary tables:**

- `branching/all/global_branch_distance_slopes.csv`;
- `branching/direct/global_branch_distance_slopes.csv`;
- associated frozen direct between-context / multivariate branching outputs in the same artifact.

**Primary estimand:** isolation-associated two-axis response vector for `accessibility_generalization` and `reproductive_assurance` under the predeclared H1/H2 contract.

**Exposure:** `log_distance_to_continent_km`, interpreted as composite source separation/connectivity.

**Fixed universe:** 8,265 islands; Database 1.0 denominator 106,295 accepted angiosperms; 222,688 / 318,885 resolved raw-axis cells in the final snapshot.

**Scopes:** all-analysis-eligible and direct-only; all-native and native-nonendemic where support permits.

**Claim ceiling:** no temporal claim; no pollinator-loss mechanism; no inference from separate significance tests in different regions.

### C2 — the classical syndrome components can decouple by biogeographic context

**Immutable source:** same final progressive artifact `10058653212`.

**Headline frozen examples:**

- Palearctic all-analysis/all-native: accessibility/generalization `+0.0795`, q `0.00463`; reproductive assurance `+0.0534`, q `0.00463`;
- Palearctic direct-only/all-native: accessibility/generalization `+0.0616`, q `0.0389`; reproductive assurance `+0.0966`, q `8.48e-14`;
- Tropical direct-only/all-native: accessibility/generalization `-0.1005`, q `0.00275`; reproductive assurance `+0.1360`, q `0.00428`.

**Formal inference:** direct between-context vector tests remain the basis for context heterogeneity. The coefficient examples above are descriptive anchors, not the heterogeneity test themselves.

**Secondary decomposition:** the Palearctic floral-architecture response remains positive after conditioning on `selfing_core`; this is decomposition, not mediation.

### C3 — the strongest Palearctic island syndrome is localized near the family-to-genus transition

**Empirical source:** taxonomic-depth outputs inside final progressive artifact `10058653212`:

- `taxonomic-depth/all_analysis_eligible/slopes.csv`;
- `taxonomic-depth/direct_only/slopes.csv`.

**Descriptive effect-size synthesis:** run `34796763611`, artifact `10330230959`, especially:

- `taxonomic_effects_long.csv`;
- `taxonomic_vector_attenuation.csv`.

**Existing implementation identity:** `src/island_v2/chapter1_taxonomic_depth_decomposition.py` explicitly constructs family and genus residuals on the same observed island species after requiring both family and genus to have eligible source positions and source availability. Family/genus are grouping variables, not missing-trait imputers.

**Frozen result:** across 16 source-mode x scope x stratum profiles:

- family attenuation: `19.6–33.4%`, median `26.2%`;
- genus attenuation relative to observed: `78.8–85.9%`, median `81.3%`;
- conditional attenuation of the family-adjusted remainder at the genus stage: `70.6–79.1%`, median `75.7%`;
- canonical support ladder: `4/4 -> 4/4 -> 0/4`.

**Current uncertainty status:** the percentages are descriptive re-expressions of frozen slope estimates. A paired uncertainty interval for attenuation itself has **not** yet been produced. This is the main P1 gap.

**Claim ceiling:** compatible with genus-level assembly beyond family; not causal mediation; does not prove dispersal alone, absence of within-lineage change, or absence of evolution.

### C4 — the Palearctic core is not trivially explained by the tested observation processes

#### V5: trait-resolution MNAR

**Historical formal source:** run `33141250895`, artifact `9674100317`; canonical interpretation is recorded in `docs/chapter1_v5_mnar_tipping_point_checkpoint.md`.

**Role:** trait-state-dependent resolution sensitivity conditional on recorded flora.

**Claim:** the broad Palearctic primary vector survives the finite predeclared MNAR grid; some reproductive details and cross-context contrasts remain bounded rather than universally robust.

**Do not claim:** arbitrary MNAR is excluded.

#### V6: species-list detection

**Source:** run `34800498716`, artifact `10331282464`.

**Primary output:** `species_detection_tipping_surface.csv` plus manifest/summary.

**Frozen result:** Palearctic accessibility survives `99/100` baseline-supported surfaces; all `80/80` surfaces with distance-dependent under-survey in the concerning direction survive. Tropical accessibility survives `35/75`; the North–Tropical vector contrast survives `70/75`.

**Do not claim:** the grid estimates true flora completeness or a posterior distribution over detection.

### C5 — one common global nonlinear transition is not identified

**Source:** observed geometry run `34765070227`, artifact `10320630085`.

**Primary output:** `observed_geometry_cross_scope.csv`.

**Frozen result:** `0/12` broad atomic cells meet the cross-scope promotion rule; all are `monotonic_or_unresolved` under the calibrated contract.

**Interpretation:** no promoted common global assemblage breakpoint under the tested design.

**Do not claim:** local or lineage-specific thresholds are absent.

### C6 — independent global pollinator evidence does not identify the upstream mechanism

#### N1 channel heterogeneity

**Source:** run `34750351009`, artifact `10315671249`; lock `config/chapter1_nee_n1_result_lock.json`.

**Frozen result:** joint isolation x channel Wald `W=1.6187`, df `3`, p `0.65516`; N1 failed and the preregistered chain stopped before N2.

**Power boundary:** limited power for modest channel differences; this is non-identification, not equivalence.

#### effort-matched GloBI source breadth

**Source:** run `34790842133`, artifact `10328696801`; lock `config/chapter1_globi_source_breadth_v2_result_lock.json`.

**Frozen result:** `0/4` context x stratum cells promoted; N1 is not rescued and N2 is not reopened.

#### H5c independent biotic-vs-wind specificity

**Source:** run `34803837463`, artifact `10332871066`.

**Frozen result:** only one prospectively qualified observed cell; distance x biotic interaction `+0.06495`, 95% CI `[-0.09030, 0.22020]`, p `0.41221`; classification `no_pollination_mode_specificity_support`.

**Do not claim:** pollinators are irrelevant.

#### H5d distributed-threshold identifiability

**Source:** run `34803545574`, artifact `10332671123`.

**Frozen result:** `0/8` designs qualified; classification accuracy approximately `0.733–0.778`; smooth heterogeneous clines falsely selected as distributed-threshold generators `19.0–25.5%` of the time.

**Do not open:** observed genus-specific threshold distributions under this failed design.

## P0 discrepancies and boundaries

1. **No central artifact/digest mismatch was found** among the live artifacts listed above.
2. The user-supplied parent interpretation reference `e4a73796a` is not used as an empirical source in this ledger; the current main branch and immutable run receipts are the evidence anchors.
3. H3 attenuation percentages currently lack paired uncertainty for the attenuation estimand itself. This is an explicit P1 task, not a hidden completed result.
4. The existing H3 implementation already enforces common observed species for family/genus decomposition. P1 must therefore test **uncertainty and grouping specificity**, not claim novelty from re-establishing common support alone.
5. All-analysis and direct-only reuse much of the same underlying flora. They are robustness scopes, not independent datasets.

## P0 decision

**PASS for proceeding to P1.**

Reason: every central manuscript claim in the current island-first chain has a traceable immutable source and the major live artifacts match their recorded IDs/digests. The remaining central vulnerability is not provenance; it is whether the large family-to-genus attenuation is statistically precise and biologically specific to genus structure rather than generic fine-group flexibility.
