# Functional traits across unseen mainland-source regions: strict predictive transfer

**fd PR #13 · 2026-11-10 · real, rights-restricted original-source replay.** This is a new **conditional trait-prediction** analysis, NOT observed pollination effectiveness, independent phylogenetic replicates, true historical mainland ancestors or evidence that island selection caused the retained trait states. Only broad anonymous aggregate results may be released.

## Why this is different from previous island-pair tests

The previous `fd` SC/symmetry results showed strong *known-genus* predictability, but evaluated taxa or source-region data could recur in different source training comparisons. A former component-only source exclusion gave positive genus-vs-global Brier score gains for **72/72** scope × trait × ledger × large-source-component comparisons. However, different test components could share species that remained in each other's mainland training donors. The **global all-island-name holdout below is the stronger, controlling result**.

## Exact data source and fold/unit construction

- Original protected `island` progressive checkpoint artifact `10058653212` (SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`), and the separate original source-scale trait ZIP artifact `10042377593` (SHA256 `5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd`).
- Original public fd G0/G1 `all_pairs.csv` ZIP artifact `13741845313` (SHA256 `2efa7d2bb40adc60cc7d22b3c0c1043eb63c7c608e492fb0a422a14b3a563801`).
- **1,419 distinct original G0-paired-cohort islands**, containing 879,010 unique island × accepted-species records and **110,587 distinct accepted plant names**. The regional source-overlap graph defined by `geo_k5` donor-region sets yields nine component cohorts with >=40 islands: northern-high **90**, northern-mid **640/165/120/57/57**, southern-extratropical **44**, and tropical **127/59**.
- Each donor component is defined before inspecting traits. Every island contributes one five-nearest native GIFT mainland source signature. The component is a connected subgraph under *any shared donor source entity* within a broad region. Distinct components are not necessarily independent botanical lineages, and training/reference sets overlap across folds and regions.
- Original geographically plausible native mainland GIFT source lists have 1,666,186 region × accepted-name records across 748 regions. These represent potential checklist opportunities, **not proven arrival routes**.

## Leakage-controlled models and scores

For EACH test component:
1. Entirely exclude its set of five-nearest native GIFT mainland region entities from the training reference.
2. In the **strict `all_g0` test**, exclude **every one of the 110,587 accepted names appearing on ANY of the 1,419 fd island floras** from *all remaining GIFT entities*, not just the test component's name set. Hence a test species in a different component cannot become a donor, and all original test species remain excluded even if their trait is missing.
3. Use only strict recorded, single-state SC/SI and actinomorphic/zygomorphic species trait evidence. All/Direct and the two original training trait ledgers (`snapshot`, `source-adjusted`) are distinct **nested sensitivity definitions**, not external replications. Unknown and multistate records are missing.
4. Predict each **distinct target taxon name once within a heldout component** using the observed known genus. Baseline probability is the global remaining native-mainland trait frequency; genus probability is `(same_genus_positive_donors + tau * global_prior)/(same_genus_donors + tau)`. Report proper Brier score. The main exploration uses tau=**1**, with 0.05 and 8 as predeclared numerical sensitivities for this replay. This does **not** predict which genera colonize an island.
5. Independently evaluate `all_recorded` and **source-verified native-only** test names. The *training* exclusions do not change for native-only: ALL island taxa, including introduced and origin-unresolved taxa, remain globally excluded from donors.

There are no ecological null p values or confirmed spatial/phylogenetic confidence intervals in these tests. The nine source-connected component outcomes may share botanical species and mainland donors across geographic regions, so a `9/9` pattern is not nine independent evolutionary experiments.

## Main strict `all_g0` result: the functional-trait contrast becomes sharp

Brier gain = global native-mainland baseline error minus out-of-global-island-name-and-out-of-source genus-model error, **positive means better genus-based trait prediction**.

| Scope | Strict SC/SI (two ledger analyses × 9 components) | Strict floral symmetry (two ledger analyses × 9 components) |
| --- | ---: | ---: |
| **All** | **18/18** positive at tau=1 | **18/18** positive |
| **Direct** | **8/18** positive at tau=1 | **18/18** positive |
| Direct, tau=0.05 | **0/18** positive | **18/18** positive |
| Direct, tau=8 | **15/18** positive | **18/18** positive |

Thus **62/72** strict all-species-at-once heldout comparisons have a positive tau=1 genus prediction gain, as opposed to **72/72** under the weaker component-only name holdout. The difference lies entirely within the **Direct SC/SI** sensitivity; this was concealed by the easier name holdout.

For the strict `source-adjusted Direct` ledger, the following figures illustrate representative *proper-score* sizes:

| Region / source-connected component size | Trait | Global Brier | Genus Brier (tau=1) | Improvement |
| --- | --- | ---: | ---: | ---: |
| Northern high / 90 islands | SC/SI | 0.21719 | 0.21918 | **−0.00199** |
| Northern mid / 640 islands | SC/SI | 0.21634 | 0.21191 | **+0.00443** |
| Tropical / 127 islands | SC/SI | 0.22339 | 0.22097 | **+0.00242** |
| Northern high / 90 islands | Floral symmetry | 0.20966 | 0.17337 | **+0.03629** |
| Northern mid / 640 islands | Floral symmetry | 0.21896 | 0.16948 | **+0.04947** |
| Tropical / 127 islands | Floral symmetry | 0.24272 | 0.16869 | **+0.07403** |

All nine strict source-adjusted Direct symmetry gains are positive (minimum **+0.03236**). **All nine remain positive under tau=0.05 and tau=8**, as well as with the other Direct ledger. This demonstrates robust **genus-associated trait-state predictability out of evaluated taxa and geographic source-region identity**, conditional on known genus. It is *not* phylogenetically independent proof of selection or independent convergent evolution.

### Why reproductive assurance is much less stable

The strict globally excluded **source-adjusted Direct SC/SI** reference retains only **104–105 distinct trait-known mainland donor species per heldout source component**, versus **806–886** for Direct floral symmetry. Known same-genus donor support covers only approximately **4–13%** of SC/SI target species, compared with about **20–46%** for floral symmetry. Many SC genus predictions are thus pure/global-prior fallbacks; their failure to improve cannot distinguish weak biological genus conservation from severe data sparsity.

This is an **evidence-availability/identifiability failure**, not positive evidence that self-compatibility evolves independently of genus. Earlier broad All SC/SI source-name heldout gains are positive (18/18), but All includes less direct evidence and must not override Direct failure.

## Stringent source-native test

A strict source-backed island-native target subset was evaluated without allowing ANY original island names to return to mainland donor training. Release only groups with >=40 original component islands and >=50 distinct strict-known native targets:

| Test source | SC/SI | Floral symmetry |
| --- | ---: | ---: |
| All (both ledgers) | 16/16 positive | 15/15 positive |
| Direct (both ledgers) | **8/14 positive** | **14/14 positive** |

Therefore **53/59** eligible source-native strict comparisons are positive at tau=1. The **six nonpositive comparisons are all in Direct SC/SI** (8 positive of 14); All SC/SI and both evidence tiers of floral symmetry retain positive gains. One northern-high Direct SC/SI source-ledger cell has just **42** strict native species, so is **withheld** by the >=50 species gate, not treated as a negative or missing ecological effect.

## Source independent holdout remains conditional prediction, not an island syndrome

Important distinctions:
- The nine large *within-region* donor components are outcome-blind, but not globally independent. Under a single all-region GIFT-donor overlap graph the 1,419 islands fall into nine components, **one containing 1,333 islands**; consequently the 72 scored cells share taxa, genera and mainland donor opportunities.
- This experiment tests **unseen species trait states with their genus given**. It does not model genus dispersal, plant establishment, native history, pollinator performance, pollen limitation, or within-lineage evolution.
- The two trait ledgers and All/Direct scopes are overlapping, not independent replications. Exact strict trait states may share evidence providers and within-genus recording conventions.
- Smoothing sensitivity was assessed post hoc; a positive gain against the global mean is a predictive result, not a causal source-geography advantage. The separate actual-geographic-GIFT versus donor-matched placebo test in [MATCHED_LEDGER_PLACEBO.md](MATCHED_LEDGER_PLACEBO.md) shows that **the specific nearest mainland source labels add little robust extra information** beyond genus.
- These outcomes are aggregated by **distinct species within each component**, not weighted by island frequency. The same species can occur in multiple evaluated components, but no name ever leaks into the strict *training* donor pool.

## Reproducibility and publication boundary

The private replay was completed on the original source ZIPs and aggregate-only outputs were compared with a separate implementation: all **32** source-adjusted Direct (strict all-recorded plus retained native-only) component-level summary rows, including source-disjoint donor counts, matched evaluation species counts, Brier baselines, genus gains and tau sensitivities, agree to floating-point precision (max differences below 1.5×10⁻¹⁶). This detects implementation drift but does **not** make the original donor trait data independent.

Code: `src/fd/source_disjoint_trait_cv.py` with synthetic leakage/duplicate/unknown-trait tests in `tests/test_source_disjoint_trait_cv.py`. Inputs are rights-authorized private `island` and GIFT rows; the public GitHub CI validates logic on synthetic fixtures and does **not** replay licensed source-level calculations.

**Claim ceiling:** the reproducible biological signal is **transferable genus-linked recorded floral symmetry**, not source-specific island adaptation. For Direct SC/SI, insufficient source-disjoint known genus trait donors prevent a robust global result. Future work must improve truly independent reproductive trait measurements and independent source histories, not multiply overlapping island-pair comparisons.
