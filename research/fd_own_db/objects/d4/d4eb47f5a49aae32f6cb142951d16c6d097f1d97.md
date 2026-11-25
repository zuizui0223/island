# Original GIFT 4,505-island species-Direct flower evidence quality is associated with recorded category states, even within genera

**Real protected original botanical data, independently source-audited, 21 November 2026.** [Original island Actions #44838287115](https://github.com/zuizui0223/island/actions/runs/44838287115) **SUCCESS** (15 source/synthetic guard tests), original [island DRAFT PR #344](https://github.com/zuizui0223/island/pull/344); original source-anonymous artifact **#14678744384**, downloaded ZIP **SHA256 cb037c3ac1ba60240d9c769b23bc800aaebbe11593f2495ddbda7592eb44e5be**, exactly one JSON **SHA256 692e21ee2511305a1661e089c65b9b2f932eda05496786c8e290882c7d19627e** (120,397 bytes). Original SHA-locked GIFT accepted botanical island floras **4,505**, independent GSHHG geography island parent **8,264**, and source species-Direct trait ledger remained unchanged.

This study addresses a core limitation of the previous source-original geography→genus→floral-category model: on all 59 source-group qualified original region×trait×botanical-status panels, the largest remaining category forecast residual was externally learned **congeneric floral trait-profile mismatch**, not geography-predicted botanical genus composition. That prediction residual could arise from real congeneric biology, trait provenance, annotation uncertainty, missing species evidence and/or source-biogeographic domain shifts; it does **not** itself identify which one caused it.

**Critical source identifiability boundary:** The original Direct ledger has **one normalized category and one high/medium quality grade per accepted species×trait**, NOT replicate measurements of the same species on multiple islands. Any high-versus-medium source result is a difference between **different botanical species** graded differently, even when compared within the same genus. It is NOT proof of measurement error and CANNOT establish the proportion of the source residual caused by annotation versus biological heterogeneity.

## Prespecified study and source-safety checks

All 10 genuinely species-Direct single endpoints and 3 true same-accepted-species floral pairs are considered, in each of Northern midlatitude / Tropics / Southern extratropical original GSHHG macroregions, for both all-recorded and nested literal known-native floras: **13 × 3 × 2 = 78** complete source-anonymous scope panels. Each accepted species counts ONCE within a geographic region/status source collection even when its species-wide trait is assigned to numerous island checklists. The study retains the full GIFT original 4,505-island botanical universe and records all unsupported panels rather than requiring complete 13-trait cases. A paired trait receives HIGH grade only when **both trait values are high grade on the same accepted species**; otherwise it is MEDIUM if both were individually eligible Direct high/medium.

For each source endpoint×region×scope, compare the observed nominal category-frequency distribution of HIGH species with that of MEDIUM species by total variation distance (TV, 0–1). More importantly, restrict to original **certified botanical genera with at least 2 distinct HIGH and at least 2 distinct MEDIUM source Direct species in the SAME genus** and average within-genus TV **equally by qualifying genus**, not by island occurrence. This reduces, but does not eliminate, differences driven by which genera have records.

Under the conditional exchangeability hypothesis that trait category and source evidence-quality label are independent *WITHIN genus*, randomly reassign the high/medium labels among distinct species WITHIN each qualifying genus while retaining each genus's exact original high/medium totals, repeated **199** times with original frozen label-blind seeds. Require ≥10 matched original genera and ≥30 unique directly scored HIGH and MEDIUM accepted plant species to conduct inference. Report the observed equal-genus TV, null mean and 95th percentile, empirical one-sided p=(1+#null ≥ observed)/200; apply an exploratory Benjamini–Hochberg adjustment across ALL testable 78 nested region×status×endpoint tests, not only source-positive panels. Original species, genus, island and source-region identifiers, exact coordinates and categorical botanical trait labels never enter released JSON or CSV.

[All 78 original-anonymous source panels with evidence quality coverage, matched genus sample sizes, observed vs permuted category TV, exploratory p and BH q](GIFT4505_DIRECT_HIGH_MEDIUM_GENUS_NULL_78.csv).

## Results: source evidence quality and category state are often associated beyond genus mixture

Of 78 original region×trait×status panels, **55** had enough distinct HIGH/MEDIUM species and ≥10 genera for the genus-conditioned permutation. **43 of 55 testable panels** had observed within-genus category TV greater than expected from permuting evidence grade labels, with **exploratory BH q<0.05** across all 55 qualified tests. The 78 panels include overlapping source species/statuses and correlated flower/reproduction categories; neither 55 nor 43 is a count of independent evolutionary studies, and BH q does not cure nonrandom source trait grading.

| Source macroregion | All-recorded: testable | All-recorded: BH q<0.05 | Known-native: testable | Known-native: BH q<0.05 |
| --- | ---: | ---: | ---: | ---: |
| Northern midlatitude | 11/13 | **11** | 11/13 | **9** |
| Tropical | 11/13 | **9** | 8/13 | **3** |
| Southern extratropical | 11/13 | **9** | 3/13 | **2** |
| **Total** | **33/39** | **29** | **22/39** | **14** |

The original independent GIFT Southern known-native stage-2 *forecast* was unsupported because it lacked ≥30 independent geographical donor islands; the present source-quality association can still be computed for THREE endpoints using DISTINCT region-level botanical species. These are DIFFERENT estimands and support conditions, so do NOT use the source-quality subset to claim that stage-2 Southern native floral prediction has been replicated.

### Two prior strongest island-weighted flower-category source forecast endpoints

| Same actual species-Direct endpoint among all recorded original regional botanical species | n HIGH | n MEDIUM | Matched certified genera | Genus-equal TV observed | Mean TV under within-genus grade exchange | Observed minus exchange mean | Empirical p / exploratory BH q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| **Tropical floral symmetry** | **1,467** | **21,511** | **77** | **0.028290** | **0.015321** | **+0.012969** | **0.005 / 0.00764** |
| **Southern extratropical primary flower colour** | **3,747** | **4,222** | **252** | **0.423775** | **0.308269** | **+0.115506** | **0.005 / 0.00764** |

In tropical symmetry, only **6.38%** of directly documented distinct regional species have HIGH grade (1,467/(1,467+21,511)); this does not prove MEDIUM labels are wrong, but it makes strict high-only source trait coverage and positive-source-genera retention especially important. Southern primary flower colour has about **47.02%** HIGH-grade directly documented species; nevertheless HIGH/MEDIUM category distributions differ strongly even after genus conditioning. The quality-label permutation p-value has minimum 0.005 with 199 random draws; do not overstate the precision of BH-adjusted q or call a coarse source quality flag a direct repeated measurement quality calibration.

Several other endpoints show within-genus grade/category association. Example: original Northern all-recorded floral form has genus-equal observed TV 0.5808 vs permutation expectation 0.3076 across 118 matched genera (difference +0.2732). This suggests source botanical category annotation/quality, within-genus biological trait diversity, or a mixture may play a substantial role in source congener transfer. The direction and causal origin of these discrepancies **remain unidentifiable** from one species-level Direct state.

## Next original discriminating prediction check

The direct next test is a **genuine HIGH-only external source congener prior model compared against the original HIGH+MEDIUM source predictor on EXACTLY the SAME heldout GIFT target accepted species and island Direct categories**, retaining full source accepted-name exclusion, genus mixture weights, regional categories, geographical block/source-group inference, and explicit attrition of genera/species that lack sufficient HIGH grade source support. An improvement from HIGH-only donors would mean quality filtering changes source predictive signal; failure to improve could reflect fewer donors and greater estimator variance, not that MEDIUM is reliable or that biological source heterogeneity is absent. Matching source donor counts within genus, predeclaring shrinkage strength, source taxonomic support and label ontology are needed to distinguish source-quality effects from donor sample size. Original source models and native-only inferential support must not be silently redefined.

A FULL independent comparison of same-accepted-species repeat floral phenotypes measured separately on islands and putative mainland sources, independent source trait data records, cross-source coding audits and effective local pollinator visits/pollen/seed fitness are still necessary to attribute biological selective mechanisms or pollinator syndrome evolution.

**Nonclaims:** no direct observation of bee, butterfly, bird or moth guild switching; no measured mainland ancestral flower phenotype; no identified local island adaptation; no controlled phylogenetic evolution or independent biological genus replicates; no demonstrated functional hypervolume contraction, morphological 17D KDE validity or causal attribution of previous signed category residuals to annotation error alone.
