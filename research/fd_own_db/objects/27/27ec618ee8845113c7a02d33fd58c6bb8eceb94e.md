# Island mainland-source genus opportunity, two-stage floral prediction and 17D projected KDE

**2026-11-08 JST — post-hoc, source-audited, aggregate-only.** [Methods and privacy contract](../../docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md). Original trait universe: [island PR #321](https://github.com/zuizui0223/island/pull/321). **No biological arrivals, confirmed absences or selection processes are directly estimated.**

## 1. Crucial source-name leakage audit and corrected scientific interpretation

The initial opportunity-ranking pilot allowed **exact same species names** to appear in both the island occurrence outcome and its geographically nearby native-GIFT mainland source checklist. For **describing genuine mainland supply/overlap**, an identical mainland species is biologically informative. But for claims of **independent prediction of an unseen taxon-name assemblage**, this source identity overlap is a shortcut.

We therefore require a **stricter explicit holdout**: globally remove *every accepted taxon name recorded on any of the 434 geology-labelled islands* from the 748 GIFT native mainland source checklists **before constructing source-genus predictor features**. This removes 742,400 source-region×name records from the original 1,666,186, leaving **923,786**. It is conservative to the point that it removes true potential supply species, so the two analyses estimate different questions.

**Revision:** The previous apparently positive source-rank AUC gain on continental islands (**0.676 → 0.762**) does **not survive** this unseen-name holdout. It must **not** be cited as a validated predictive discovery. Retain it only as a source-overlap descriptive result, not an out-of-name prediction.

## 2. Genus source opportunity under strict taxon-name holdout

Candidate genera are the union of native genera in the island's five nearest retained GIFT mainland regions. Source-only models use: global number of mainland source regions with genus; count of five nearby regions; distance-weighted presence; and source-list-size-normalized native genus species density. Rank evaluation uses *recorded* vs *unrecorded* candidate genera (not biological presence vs confirmed absence).

**399 of 434** geology-labelled islands have adequate all-recorded genus outcomes (174 continental/225 oceanic); **301** have enough explicitly native source-status outcomes (141/160).

| Outcome | Geological type | n | Global mainland genus AUC | Nearest source-density genus AUC | Interpretation |
|---|---|---:|---:|---:|---|
| All recorded | Continental | 174 | **0.693** | 0.656 | Nearby source model is worse after strict holdout |
| All recorded | Oceanic | 225 | **0.750** | 0.637 | Nearby source model is worse after strict holdout |
| Explicit native | Continental | 141 | **0.682** | 0.675 | No clear incremental nearby-source gain |
| Explicit native | Oceanic | 160 | **0.753** | 0.662 | Nearby source model is worse |

The source-eligible share of **observed genus records** also falls under strict holdout (median all-recorded continental **41.9%**, oceanic **36.8%**), because the exclusion deliberately removes much of the actual geographically plausible mainland supply. Source genus coverage is a data-support statistic, not a biological lack of connectivity.

The same original **37 geographically matched continental/oceanic island pairs** (native evaluable 31 pairs) were compared for gains in AP and AUC from three local models over the global source baseline. Under one fixed **12-test BH/BY family**, there are **zero supported geological contrasts** (minimum BH q≈**0.487**). This rejects any robust inference, from this dataset, that source-rank predictive improvement differs by geological class.

**Biological limitation:** Non-recorded ≠ truly absent; ranked genus AUC is not arrival, establishment, survival, persistence or true occupancy probability. Observed source taxa are a legitimate supply description but cannot simultaneously serve as completely independent name-held-out predictor without changing the estimand.

## 3. Two-stage mainland source → genus → floral phenotype prediction

Using the **original 37 matched island pairs / 74 islands** and the PR321 complete-profile taxon-name mapping, the training candidate source lists and genus-specific floral donor counts both exclude **all** matched test-island names globally, within each All/Direct scope. Hence this stricter two-stage result does NOT rely on the same evaluated species as source predictors or as floral-profile donors.

The source-only model weights candidate genera by normalized nearby source-native genus density, predicts each genus's profile using other GIFT-native mainland taxa, then forms the **island taxon-carrier-weighted frequency vector** across 496 colour/form profiles. It does **not** consume observed island genus identities. A separately labeled `known_genus_oracle` benchmark is allowed to know the actual island genera and is NOT prospective.

| Scope/model | Mean held-out taxon exact-profile log loss ↓ | Mean full original categorical-kernel repertoire error ↓ |
|---|---:|---:|
| All / global mainland prior | 5.3563 | 0.020476 |
| All / **strict source-only genus-density** | 5.3505 | **0.018107** |
| All / known-genus oracle (post-hoc bound) | 5.2956 | 0.012699 |
| Direct / global mainland prior | 5.4687 | 0.024446 |
| Direct / **strict source-only genus-density** | 5.4721 | **0.021582** |
| Direct / known-genus oracle (post-hoc bound) | 5.4473 | 0.018572 |

**Strict two-stage improvement in the original categorical-distance repertoire error:**
- **All:** absolute +0.002369, relative **11.6%**, paired 37-cohort bootstrap 95% **[+0.001317,+0.003448]**.
- **Direct:** absolute +0.002864, relative **11.7%**, 95% **[+0.001616,+0.004249]**.

The improvement in **name-level exact-profile log loss is not clearly positive** (All +0.0058, interval [−0.0035,+0.0152]; Direct −0.0034, [−0.0099,+0.0032]). Thus a source-only prediction improves **assemblage-level categorical-similarity composition** more consistently than the exact flower type of an individual unknown taxon.

The source-only model is still not a colonization model: the evaluation uses the original island species/name counts; genera without independent floral donors fall back to a global prior. A fully independent island survey, source ancestry data, and phylogenetic generalization test would be required for causal assembly claims.

## 4. 17D hypervolume: measured, but no unique geological effect beyond profile richness

From 496 source-audited exact colour/form state-set profiles, classical PCoA of sqrt(original distance) has **80.4899%** positive inertia in 17 dimensions (versus **15.353%** in the previous 2D projection). Equal-profile-weight 90%-mass isotropic Gaussian KDE density with common scalar reference scaling, h=0.55 and h=0.85; **12,000 Monte Carlo draws ×2 repetitions/island**. Public aggregate + independent geological crosswalk CI calculated All **260 islands** (130 continental/130 oceanic) and Direct **243** (123/120), altogether **1,006 scope×island×bandwidth estimates**. Numerical relative integration SE median≈2%, independent log-volume replicate discrepancy median≈0.04.

| All-scope median 17D projected log hypervolume | Continental | Oceanic |
|---|---:|---:|
| h=0.55 | 19.537 | 19.368 |
| h=0.85 | 26.460 | 26.315 |

On the original 37 matched pairs raw All/Direct differences are negative (oceanic smaller), with four-test **BH q≈0.0306**, **BY q≈0.0638** (no BY support). **The key falsification is adjusting for each pair's number of recorded distinct profiles:**

| 37-pair adjusted oceanic−continental 17D logvolume intercept | Estimate | Pair-bootstrap 95% interval |
|---|---:|---:|
| All, h=0.55 | −0.0086 | −0.0240 to +0.0074 |
| All, h=0.85 | −0.0137 | −0.0340 to +0.0066 |
| Direct, h=0.55 | −0.0181 | −0.0360 to +0.0017 |
| Direct, h=0.85 | −0.0208 | −0.0432 to +0.0026 |

All adjusted intervals include zero. **Raw high-dimensional volume differences are substantially explained by numbers of recorded profiles**, rather than demonstrating extra geometrical floral simplification. Results are smoothed 17-axis *projected* Lebesgue volumes, not entire discrete trait functional space, nor a physical pollination niche volume. Bandwidth radically changes absolute volume.

## 5. Current paper conclusion and next inferential gate

**Two supported observations with strong limits:** (1) recorded island floral composition differs across geological classes (prior 37-pair exact-profile test); (2) a genuinely unseen-taxonomic-name **source-only genus-to-floral two-stage mixture** provides a modest improvement in predicting island **assemblage-level categorical similarity** compared with a genus-free global mainland reference. Strict name-held-out source genus **ranking individual island records** does *not* independently demonstrate nearby-source predictive success, and no origin-specific predictive gain is supported after multiple tests. The 17D projected KDE does not establish additional volume compression after recorded profile richness.

A defensible question is whether **source lineage assembly**, including genus-level source availability and which taxon names can reach a given island, can explain floral repertoire differences. The causal mechanism remains unresolved because historical dispersal, establishment, detection and true absence are not identified from these checklists.

**Reproduction:** [source opportunity ranker](../../src/fd/source_genus_opportunity.py), [source-only two-stage assembly](../../src/fd/two_stage_assembly.py), [17D hypervolume](../../src/fd/hypervolume17.py); [full method/rights contract](../../docs/SOURCE_OPPORTUNITY_HYPERVOLUME17.md). Public CI uses only source-approved aggregated flower profiles, coordinates and independently published geological references; GIFT/rights-restricted species-level runs are restricted to local authorized execution. All committed results are non-identifying summaries.
