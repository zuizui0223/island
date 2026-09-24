# Chapter 1 corrected submission baseline — 24 September 2026

This is the **primary submission baseline**, replacing the uncorrected v14 geographic exposure. The user explicitly selected this role after the geometry audit; the correction is not relegated to a sensitivity analysis. The frozen v14 locks remain unchanged as superseded provenance. This measurement repair does not create a new prospective confirmation and does not change the post-hoc role of H4.

## Geography and analysis population

The original GSHHG 2.3.7 high-resolution island polygons were measured against coarse Natural Earth 110m continental polygons in projected coordinates. In the broad H1 union, 1,113 of 4,380 units had zero distance because those geometries overlapped. Recomputing against continental coastlines from the **same original GSHHG archive** restores positive distances for all 1,113. Their median corrected distance is 0.927 km (range 0.001655–100.579 km). These are recomputed exposures, not recovered ground-truth distances; metre-level numerical output is not metre-level positional accuracy.

A second error admitted GSHHG `0-W`, a split component of Eurasian sibling `0`, as an island. Its geometry is exactly the locked unit `gshhg_2.3.7_h_8b13189234b949ee1ff6`. Excluding that component changes the universe from 8,265 to **8,264**, and the broad H1 union from 4,380 to **4,379**. It was not one of the 1,113 old-zero units. Future acquisition recombines `sibling_id` before the area filter. The current correction retains the other locked island identities instead of silently rebuilding trait joins.

Distance is the minimum separation between minor great-circle coastline arcs on a **mean-radius sphere, R = 6371.0088 km**, not a WGS84 ellipsoidal distance. Artificial dateline split edges are omitted. Continental siblings are selected with fixed interior seeds in Africa, Eurasia, North America, South America, Australia and Antarctica (L5 ice front). Segment indexing uses triangle-inequality bounds and exact arc minima, not a fixed nearest-vertex sample. Exhaustive comparisons, subdivision checks and 24 independent numerical minimizations validate the geometry implementation. Shoreline resolution and spherical approximation remain limitations.

GloPL exposure uses each site's coordinates, not an island-level proxy. Of 1,248 sites, **996 truly lie on seeded continental land** and correctly retain zero distance. This correction does not remove legitimate mainland zeros.

## Results replacing the uncorrected estimates

All original trait values, model families, controls, geographic strata, spatial/publication clusters, evidence scopes and multiplicity families are retained. Both all-analysis and direct-only results, native and native-nonendemic H1 strata, and full raw-colour/architecture tables are included. No sign or significance threshold was tuned after correction.

### H1 — recurrent, region-dependent floral responses

The seven-trait joint isolation response remains supported in all four broad strata: northern midlatitude q = 3.216e-10 (2,173 islands), northern high latitude q = 2.433e-5 (411), tropical q = 3.498e-7 (1,493), southern extratropical q = 2.504e-18 (302). These counts sum to the analytic union; individual-trait and native and native-nonendemic fits have their own denominators. Direct-only joint support also remains in all four regions.

Joint support does not establish a uniformly positive seven-trait syndrome. Twenty-six of 28 broad all-analysis coefficients are positive, but the southern shallow/open-tube coefficient is negative (beta = -0.2166, nominal p = 0.000370). Northern high-latitude generalized form weakens to p = 0.1018; southern selfing mating system weakens to p = 0.0540. Thus the submission claim is recurring functional components with region- and trait-dependent expression, not one universal floral response.

### H2 — partially separable responses, with limits

Selfing-core slopes are positive in all four regions but nominally supported only in tropical and southern regions. Selfing-adjusted accessibility has q = 0.1703, 0.000847, 0.01616 and 0.2873 in the four regions respectively. The tropical direct-only accessibility result weakens to q = **0.1196** (nominal p = 0.0448); it must not be described as FDR-supported. Southern adjusted plain colour remains supported (q = 0.000847).

The raw colour, joint colour, architecture and colour-conditioned architecture analyses were each reproduced on original inputs and then refit. Detailed patterns differ by region and conditioning model. For example, northern-midlatitude yellow/orange–large-bee-form coupling becomes negative and FDR-supported; northern-high-latitude blue/purple architecture associations remain negative; tropical direct yellow/orange–butterfly/deep-tube coupling remains positive (q = 0.03949). These descriptive trait labels do not establish the identity or historical loss of pollinators. All rows, including null and opposite-sign results, are in `all/raw_patterns` and `direct/raw_patterns`.

### H3 — independent pollen limitation

The primary standardized distance slope is **0.0919104 (SE 0.0380631; two-sided p = 0.0157489)**, compared with 0.07937 before correction. Denominators remain 2,969 effect rows, 1,408 measurement cells, 1,248 sites and 919 publications. The supplemental-only sensitivity remains unsupported (beta = 0.0441003; p = 0.310819). The no-zero-constant sensitivity gives p = 0.0186943. The positive primary association is therefore not evidence of equal robustness across measurement definitions, nor a direct reconstruction of historical pollinator limitation.

### H4 — post-hoc functional associations

The exact-species selfing score associates with lower pollen limitation (beta = -0.298301; SE = 0.103521; p = 0.003957; 455 species / 409 publications). Accessibility also associates negatively (beta = -0.295660; SE = 0.128965; p = 0.021873; 143 species / 143 publications). Supplemental-only selfing remains weak (p = 0.289619), whereas accessibility gives p = 0.039706. Atomic-trait reconstruction remains a sensitivity, with self-compatibility p = 0.1495, autonomous selfing p = 2.89e-8, generalized form p = 0.04434 and radial symmetry p = 1.197e-5.

These are explicitly post-hoc functional associations: no causal mediation, prospective validation, or proof of historical adaptation is claimed. Chapter 2 simulations are unchanged by this Chapter 1 geography correction and do not prove the historical causes of observed regional patterns.

## Reproducibility and current entry points

- Current machine-readable selection: [`../config/chapter1_submission_current.json`](../config/chapter1_submission_current.json).
- Complete corrected tables, comparisons and geometry receipts: [`../results/geography_20260924/`](../results/geography_20260924/).
- Reproduction instructions: [`../scripts/geography_correction/README.md`](../scripts/geography_correction/README.md).
- Parent model/method text: [`chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md`](chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md). Its uncorrected geographic estimates and result claims are superseded by this correction; it must not be submitted unchanged.

A full model replay is distinct from fast CI, which checks code, hashes, population invariants and current-surface consistency. Green CI alone is not a scientific replication claim.

The alpha1 database workflow remains a historical reproduction lane with its original 8,265-unit contract. It restores the exact hash-locked island artifact rather than rebuilding old identifiers with the repaired acquisition rule. This does not select alpha1 geography for the corrected submission.
