# Chapter 1 publication figures and tables — design

Date: 2026-09-17
Status: approved direction; implementation pending spec review
Base: `main` at `c9ebc7bec83ab49106a22ae2cefb432dfb09b8e3`
Target branch: `chore/chapter1-publication-figures`

## 1. Goal

Build a reproducible publication-figure pipeline for the current Chapter 1 H1–H4 paper directly from frozen result locks and pinned GitHub Actions artifacts.

Publication-facing artwork must not contain internal version labels such as `v13`, workflow IDs, artifact IDs, branch names, or repository implementation terminology. Provenance remains in code, manifests, source-data tables and documentation.

The figure set is designed for a Nature Ecology & Evolution Article: concise main-display items, exact sample sizes, explicit uncertainty, restrained colour, and editable vector output.

## 2. Journal-format constraints

Use Nature-branded research-journal artwork conventions:

- main figures prepared at 180 mm double-column width;
- target figure-body height <= 170 mm unless a specific panel requires less;
- Arial or Helvetica throughout;
- final-size text 5–7 pt;
- panel labels lower-case bold (`a`, `b`, `c`, `d`);
- RGB colour mode;
- colour-blind-safe palette and redundant shape/line encodings; do not rely on red-versus-green contrast;
- editable vector output as PDF and SVG, plus 300 dpi or higher PNG for review;
- line weights >= 0.25 pt at final size;
- figures must remain legible at final print size;
- legends explain error bars, sample sizes and statistical tests rather than embedding long prose inside panels;
- no figure titles such as `Main Figure 2` inside artwork; manuscript numbering belongs in the legend/caption;
- no internal labels such as `v13`, `H1`, `H2`, `H3`, `H4` inside the publication artwork unless scientifically necessary; use biological panel titles instead.

Nature Ecology & Evolution Article format permits up to six main display items and up to ten Extended Data figures/tables. The primary design therefore uses three main figures and keeps detailed audit/sensitivity material in Extended Data.

## 3. Statistical and provenance rules

Every plotted numeric value must come from one of:

1. `config/chapter1_v13_unified_island_syndrome_result_lock.json`;
2. `config/chapter1_v13_functional_bridge_result_lock.json`;
3. an artifact explicitly pinned by those locks or their frozen parents;
4. deterministic descriptive summaries of those frozen values.

No coefficient, standard error, confidence interval, sample size or P value may be copied manually into plotting code.

For every generated display item, emit a companion source-data CSV and a JSON provenance manifest recording:

- input lock paths and blob SHAs;
- GitHub Actions run IDs, artifact IDs, artifact names and artifact digests;
- source filenames within artifacts;
- row counts and uniqueness checks;
- output SHA-256 digests.

The plotting workflow must fail if pinned digests, expected row counts or required columns do not match.

## 4. Important sample-definition distinction

The global geography/covariate artifact contains 8,265 frozen islands. Of those, 4,453 have `n_trait_species > 0` and are **trait-informed plant-analysis inputs**.

That 4,453 count must not be described as the exact H1 multivariate model sample. The frozen H1 lock reports model-support counts of:

- primary all-analysis: 2,164 northern midlatitude + 409 northern high latitude + 1,462 tropical + 299 southern extratropical = 4,334 context-specific analysis islands;
- Direct-only sensitivity: 2,137 + 400 + 1,432 + 287 = 4,256 context-specific analysis islands.

Accordingly:

- the global map labels 4,453 as `trait-informed island inputs`;
- coefficient panels use their own frozen model-specific `n` values;
- no panel equates 4,453 with the complete-case H1 model sample.

The GloPL preflight contains 2,969 effect rows, 1,248 unique coordinate sites and 919 publications. Exactly 197 GloPL sites occur on 37 frozen islands, and all 37 are within the 4,453 trait-informed input set.

## 5. Main Figure 1 — Global study coverage and evidence design

Purpose: establish geographic scale and independent evidence layers before presenting effect estimates.

### Panel a — Global coverage map

Plot from frozen real data only:

- all 8,265 frozen islands as very light grey points;
- 4,453 trait-informed island inputs as the principal plant-side point layer;
- all 1,248 GloPL sites as small neutral crosses;
- the 37 frozen islands containing exact GloPL sites as larger outlined star/diamond markers;
- the four geographic replication strata differentiated with a colour-blind-safe palette plus redundant marker shapes if needed.

Do not draw schematic ellipses or region labels over the map.

### Panel b — Evidence design

A compact schematic, not a result panel:

`Geographic isolation -> recurrent plant response`

`Geographic isolation -> experimental pollen limitation`

with the plant response branching into:

- reproductive assurance;
- floral accessibility/generalization.

A dotted or visually distinct connector from frozen trait state to current GloPL response indicates post-hoc functional triangulation. A historical mediation/selection arrow is shown only if necessary and must be dashed and labelled `not identified`.

Keep this panel visually subordinate to the data map.

## 6. Main Figure 2 — Recurrent plant response across geographic replications

Purpose: show the core plant-side result with estimates and uncertainty rather than symbolic positive/negative icons.

### Panel a — Primary six-atomic coefficient forest

Display the six identically oriented atomic isolation slopes for each of the four replication strata in the primary evidence scope.

Requirements:

- x-axis is standardized isolation-effect estimate;
- point estimate + 95% CI for every atomic response;
- vertical zero line;
- traits grouped visually into reproductive-assurance and accessibility/generalization families;
- contexts encoded with colour-blind-safe colours and/or facets, not ranked;
- exact context-specific model `n` shown in strip labels or legend;
- joint multivariate P value shown once per context, not repeated per atomic coefficient;
- no significance stars.

### Panel b — Direct-only sensitivity

Same coefficient layout and axis limits as panel a so visual differences reflect evidence scope rather than rescaling.

### Panel c — Classic-island orientation

Show the descriptive mean of the six identically oriented slopes for each context in primary and Direct-only scopes.

These means are descriptive summaries. Do not draw inferential confidence intervals unless a frozen covariance-based uncertainty is available. Use paired points/lines or a compact dot plot and label them `descriptive mean direction` in the legend.

### Panel d — Two response families

Show the frozen descriptive family mean slopes for:

- reproductive assurance;
- floral accessibility/generalization;

for all four contexts in both evidence scopes. Use the same zero-centered axis and distinguish family from evidence scope with shape/line style.

Primary visual message: both families are positive in all four geographic replications, while atomic magnitudes vary.

## 7. Main Figure 3 — Experimental pollen limitation and functional compatibility

Purpose: connect the plant pattern to an independent experimental layer without implying historical mediation.

### Panel a — Global pollen-limitation distance effect

Prefer an estimate/CI display over a synthetic scatterplot unless the frozen artifact contains the exact model-level partial-effect data required to reproduce the fitted relationship.

Show:

- global standardized distance estimate `0.079368...` with 95% CI computed from the frozen SE `0.037734...`;
- exact support: 2,969 effect rows, 1,248 sites, 919 publications;
- two-sided P and one-sided positive P in the legend or compact annotation.

### Panel b — Frozen sensitivities and offshore shape diagnostic

Forest plot with separate rows for:

- primary global distance effect;
- supplemental-only sensitivity;
- no-zero-constant sensitivity;
- within-offshore gradient;
- offshore-only gradient;
- mainland-to-offshore step if retained as a diagnostic.

Visually distinguish confirmatory/frozen sensitivity rows from post-hoc shape diagnostics.

### Panel c — Exact-species trait-state associations with current pollen limitation

Forest plot of evaluable traits:

- autonomous selfing;
- self-compatibility;
- actinomorphic symmetry;
- generalized floral form.

Show estimate + 95% CI, exact P value in the source table, and sample support (`n_cells`, `n_sites`, `n_species`, `n_publications`) in a compact adjacent column or legend.

Group rows by the two plant-response families.

### Panel d — Robustness and negative-result boundary

Use a compact matrix/forest showing:

- autonomous-selfing supplemental-only, no-zero-constant, within-publication and within-publication-by-site checks;
- architecture sensitivity checks where evaluable;
- the frozen parent distance-by-trait moderation families as `not supported` and not reclassified.

The panel must visually separate:

- global association evidence;
- within-group robustness;
- non-evaluable tests;
- previously failed moderation tests.

No causal arrows appear in Figure 3.

## 8. Extended Data display items

Use no more than ten Extended Data items.

### Extended Data Figure 1 — Atomic support and sample-size audit

Per outcome × context × evidence scope, show analysis `n`, informative species/trials and missingness/support where available from pinned artifacts.

### Extended Data Figure 2 — Geographic sampling and GloPL overlap audit

Show the 8,265 -> 4,453 input attrition, 1,248 global GloPL sites, 197 island GloPL sites and 37 exact frozen islands, plus geographic-regime counts. This is an audit-oriented companion to Figure 1.

### Extended Data Figure 3 — Detailed GloPL sensitivity and shape diagnostics

Include model variants, offshore diagnostics and any influence diagnostics that are frozen and relevant to the promoted global distance claim.

### Extended Data Figure 4 — Functional-bridge sensitivity matrix

For all evaluable traits, show primary, supplemental-only, no-zero-constant, within-publication and within-publication-by-site estimates or explicit `not evaluable` cells.

### Extended Data Table 1 — Data layers and support

Columns: evidence layer, estimand, observations, sites/islands, publications, frozen source, inferential role.

### Extended Data Table 2 — Six-atomic plant coefficients

One row per context × evidence scope × atomic response, with estimate, SE, 95% CI, support counts and joint context P value.

### Extended Data Table 3 — Descriptive plant summaries

Context × scope rows for classic-orientation mean and the two family means, clearly labelled descriptive.

### Extended Data Table 4 — GloPL global and sensitivity estimates

All promoted/frozen distance and offshore diagnostic estimates with uncertainty, P values, sample support and evidence role.

### Extended Data Table 5 — Functional trait-state associations

All H4 evaluable/non-evaluable traits and every sensitivity, with estimate, SE, P values and support counts where available.

### Extended Data Table 6 — Exact GloPL-tested islands

The 37 frozen islands with island ID, coordinates, analysis regime, trait-informed status, distance to continent, number of GloPL sites, rows and studies.

## 9. Output structure

Create a dedicated publication surface:

```text
figures/chapter1/
  main/
    figure1_global_coverage.{pdf,svg,png}
    figure2_recurrent_plant_response.{pdf,svg,png}
    figure3_pollen_limitation_bridge.{pdf,svg,png}
  extended_data/
    extended_data_figure1_atomic_support.{pdf,svg,png}
    extended_data_figure2_sampling_overlap.{pdf,svg,png}
    extended_data_figure3_glopl_sensitivity.{pdf,svg,png}
    extended_data_figure4_functional_bridge_sensitivity.{pdf,svg,png}
  source_data/
    figure1_source_data.csv
    figure2_source_data.csv
    figure3_source_data.csv
    extended_data_table1_data_layers.csv
    extended_data_table2_atomic_coefficients.csv
    extended_data_table3_descriptive_summaries.csv
    extended_data_table4_glopl_estimates.csv
    extended_data_table5_functional_bridge.csv
    extended_data_table6_glopl_tested_islands.csv
    publication_figure_manifest.json
```

Generation code lives in focused modules under `src/island_v2/`, with a single CLI entry point that renders the complete set from materialized pinned inputs.

## 10. GitHub Actions workflow

Add one workflow that:

1. installs `.[dev,figures]`;
2. downloads every pinned artifact required by the publication figures;
3. verifies artifact names/digests and mandatory files;
4. materializes tidy publication source-data tables;
5. renders all main and Extended Data figures;
6. runs structural figure/table tests;
7. checks expected dimensions, file existence and non-zero size;
8. uploads the complete publication bundle as one artifact.

The workflow must not refit scientific models. It is a deterministic rendering/repackaging workflow over already frozen results.

## 11. Tests

Add tests before implementation for:

- 8,265 frozen geography rows;
- 4,453 trait-informed input islands;
- explicit distinction between 4,453 inputs and H1 model-specific sample counts;
- 1,248 global GloPL sites;
- 197 exact-island GloPL sites;
- 37 exact GloPL-tested frozen islands, all inside the 4,453 input set;
- Figure 2 source table contains all six atomic traits × four contexts × two scopes;
- Figure 3 primary coefficients equal the frozen locks;
- all publication-facing labels exclude `v13`, workflow IDs and artifact IDs;
- all main PDFs/SVGs use the expected panel count and output dimensions;
- source-data tables contain no undocumented hand-entered values.

## 12. Manuscript synchronization

Update `docs/chapter1_v13_submission_figure_sync_20260917.md` and the manuscript figure legends so the publication-facing captions describe the new data-first figures without internal version wording.

Figure legends should begin with a brief title sentence, describe panels and symbols, define uncertainty and exact `n`, and keep methods details in Methods rather than legends.

The manuscript scientific claims and frozen result locks do not change.

## 13. Claim ceiling

The figures may show:

- recurrence of the classic island-syndrome direction across all four replication strata;
- increasing experimental pollen limitation with isolation globally;
- positive reproductive-assurance and accessibility/generalization family directions;
- negative current-pollen-limitation associations for frozen trait states where supported;
- robustness and negative results.

The figures must not imply:

- historical causal mediation from pollen limitation to trait evolution;
- global decline in pollinator abundance or visitation;
- a named pollinator-loss mechanism;
- equal response vectors across regions;
- between-region ranking as a primary result;
- within-lineage evolutionary change.

## 14. Acceptance criteria

The publication bundle is ready when a fresh GitHub Actions run can regenerate every main and Extended Data figure/table from pinned sources and all of the following are true:

1. all numeric claims match frozen locks/artifacts exactly;
2. main figures contain no repository/version jargon;
3. Figure 1 correctly distinguishes 8,265 geographic islands, 4,453 trait-informed inputs and the 37 GloPL-tested frozen islands;
4. Figure 2 uses exact coefficient estimates and uncertainty rather than symbolic signs;
5. Figure 3 separates confirmatory, sensitivity, post-hoc and negative-result evidence visually;
6. PDF/SVG artwork remains editable and publication-size text is 5–7 pt;
7. every figure has a companion source-data CSV and provenance manifest;
8. manuscript legends and the figure-sync document match the generated artwork.
