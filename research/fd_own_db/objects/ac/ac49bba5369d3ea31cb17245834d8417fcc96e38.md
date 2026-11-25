# Independent pollinator evidence bridge: what floral trait coordination can and cannot predict

**11 November 2026. fd PR #14, exploratory macroecological analysis; Direct-only complete plant morphology.** This joins *three distinct evidence domains*, without treating one as the other: (1) recorded plant colour/form/symmetry, (2) GBIF **exact-island pollinator-channel detection / qualified nondetection**, and (3) global **GloBI genus-associated documented flower-interaction channel breadth**. Actual visitation effectiveness, abundance and island-specific pollen delivery are NOT measured.

## Pinned sources and join

- Original `island` protected progressive species/traits/floristics checkpoint, workflow artifact **10058653212**, original ZIP SHA256 `b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`.
- `island` prospective pollinator-channel qualification artifact **10314984302** from run **34750153702**, ZIP SHA256 `c643b8432e04de0392867028810606f0f60f5a1d45b94bca1e28a1696ffc292d`. The five-channel `nee_projected_channel_state.csv` differentiates `retained`, `disrupted`, `unresolved` using explicitly supported GBIF exact-island evidence. The companion channel count is *capped Search evidence*, NOT abundance.
- `island` GloBI genus-source interaction-breadth v2 artifact **10328696801** from run **34790842133**, ZIP SHA256 `475dda8b0eef4ce22ad1aa736fa82cc049042d801217c65d52641489288ca254`.
- Original prospective N1 fit artifact **10315671249** from run **34750351009**, ZIP SHA256 `01d938a1f42b7f64634232a2fcb76c68a60183824684c5f135f567107a4c5d85`, with frozen channel-isolation heterogeneity Wald `p=.65516` (unsupported).
- Focal flora: **2,234 Direct single-state accepted species** with flower colour, floral form and symmetry jointly known; **928 islands** with at least 30 complete species. The main plant response is the prior 5×30-species-rarefied, 49-colour-shuffle-centred **I(colour; form×symmetry)** in nats, not merely Rao's Q.

All original species, island IDs, source entities, pollinator taxon IDs and local row-level outputs remain in a protected working session; ONLY region×channel counts and macro summary diagnostics are released.

## 1. Direct local channel support: a falsification of naïve richness

Every one of the 928 flora islands has exactly five projected GBIF channel rows (4,640 rows), covering:
- Bombus: **514 retained / 22 adequate nondetection-disrupted / 392 unresolved**;
- non-Bombus bees: **545 / 7 / 376**;
- Lepidoptera (butterflies + moths, inseparable): **759 / 0 / 169**;
- Diptera: **659 / 3 / 266**;
- flower-visiting birds: **767 / 1 / 160**.

Only **409/928** islands have all five channels evaluable; **381** of those have all five retained. **147/928** islands have zero evaluable channels. **33** of 928 have both at least one retained and one properly disrupted channel (12 northern-midlatitude, 19 tropical, 1 each northern-high and southern-extratropical). Of the 928 islands, **759** have two or more evaluable channels, but this does not mean they contain independent positive/negative channel contrasts.

Critically, Spearman `rho=0.977904` between **number of retained** and **number of evaluable** channels; with total GBIF background record counts `rho=0.795224`. Thus observed richness of supported channels is overwhelmingly tied to **observation / evidence opportunity**. Treating unresolved channels as biologically absent would create a false pollinator assemblage gradient.

A strict outcome-blind channel×region qualification (>=10 independently observed retained AND >=10 properly non-detected source-available islands, each state represented in >=10 distinct spatial blocks) is satisfied by **zero of 20 region×guild cells**. Tropical Bombus has 16 disrupted but just 2 retained islands; northern-midlatitude Bombus has 480 retained but just 5 disrupted; Lepidoptera has **zero confirmed disrupted anywhere**. The published N1 earlier channel-heterogeneity gate also failed (W=1.6187, df=3, p=.65516). There is therefore **NO defensible current macro-scale effect estimate of local pollinator-guild loss on floral colour–architecture pairing** from these binary data.

This is an *identifiability and sampling* result, not evidence that pollinator communities have no influence. Search counts in the qualification workflow can stop at fixed caps (some channel maxima 300); they cannot be substituted for species abundance or flower-visit frequency. The 5 channels do not resolve butterfly from moth pollination.

## 2. GloBI source-genus breadth is associated with floral co-association, but NOT a local assemblage

Global GloBI source records are independent of the floral source ledger but linked to the **same plant genera**. At >=3 independent documented flower-interaction references per genus, **1,534 genera** have a source-level effective number of five documented guild channels. Among 928 Direct-floral islands, **927/104 spatial blocks** have complete ecological and source-evidence covariates.

Observed Spearman `rho = +0.5328` between island colour×architecture MI excess and the mean documented visitor-channel breadth of its plant genera. A spatial-block-clustered posthoc OLS adjustment for floral marginal Rao's Q, n recorded complete species, fraction of complete species whose genus has qualifying GloBI evidence, mean log independent GloBI references per genus, island area, isolation distance, climate PC1–4, and geographic regime yields:
- partial breadth slope **+0.258341 nats per 1 source-genus effective channel**, cluster SE **0.053829**, conditional normal 95% interval **[+0.152839,+0.363843]**;
- **northern-midlatitude** 743 islands/53 blocks: slope **+0.251718** [0.143917,0.359520];
- **tropical** 105 islands/41 blocks: slope **−0.045956** [−0.225657,+0.133744] (unsupported and directionally different);
- small northern-high and southern-extratropical regional regressions are **not** promoted because support/overfitting is weak.

The full-sample coefficient is positive under every original 104-spatial-block leave-one-block-out deletion (**+0.1959 to +0.2783**). Varying the GloBI reference threshold 2/3/5/10/20 yields full-sample partial slopes **+0.266/+0.258/+0.309/+0.306/+0.302**, but source-carrying trait coverage falls from median ~95.7% to ~74.5%. These are *multiple exploratory sensitivities*, not independent p-value replications.

An additional **1,999 outcome-blind GloBI genus-label permutations** within six bins of independent reference count crossed with 3 bins of represented plant-species count gives a posthoc two-sided permutation tail fraction **0.0045** (observed +0.2583, reference median +0.0225, reference 95% range −0.1775 to +0.2058). This shows that the observed link is not created *solely* by coarse global GloBI genus reference count and sampled-plant-count gradients.

**Unresolved confound:** whole-genus GloBI breadth and island morphology both carry plant taxonomy, phylogenetic family, and source/provider documentation patterns; the 18-bin shuffle is NOT a phylogenetically conditioned, mainland-source-conditioned or effective-service null. In particular, this is *NOT* evidence of a named bird-, moth-, butterfly- or bee-pollination syndrome on any focal island.

### Geographic-block prediction is directionally improved but not independently supported

Fivefold heldout **spatial blocks** predicting the same island-level MI excess (927 islands), using fixed ridge alpha 1:
- nuisance/effort/trait-marginal/area/climate/region model OOF MSE **0.006095** (R²=0.229);
- adding GloBI genus documented breadth OOF MSE **0.005396** (R²=0.317);
- island-equal gain **+0.000699**; equal-block bootstrap 95% interval **[−0.000041,+0.001416]** includes zero, positive in 61.5% of 104 blocks;
- in 105 tropical islands the OOF MSE **worsens** (gain −0.000169).

Spatial holdout still shares genera and GloBI reference providers across blocks. This adds a directional geographic prediction diagnostic, **not a confirmed pollinator function signal**.

## 3. Does this plant morphology predict visitor-channel breadth for an unseen genus?

To check the shared-lineage explanation, we projected the 2,234 complete Direct floral species to their **486 genera with >=3 independent GloBI references** (1,597 Direct species). Train/test split is by **WHOLE PLANT GENUS** (five folds, fixed ridge alpha 10 and training-only feature standardization). The same heldout genera and source interaction outcomes are scored under:

| Predictive model | Features | Out-of-genus RMSE (lower is better) | Out-of-genus R² |
| --- | ---: | ---: | ---: |
| GloBI reference/plant-support sampling controls only | 3 | **0.7050** | **0.246** |
| + directly measured colour/form/symmetry marginal state fractions | 25 | 0.7295 | 0.192 |
| + full observed joint floral-profile state fractions | 115 | 0.7836 | 0.068 |

A second seeded 5-fold genus split repeats the rank: RMSE **0.7033** (effort only), **0.7268** (+marginals), **0.8230** (+joints). The joint model's worse prediction partly reflects sparse high-dimensional morphology predictors; this fixed ridge comparison **does not prove that no other model can predict guild breadth**.

Nevertheless, the currently positive island-level correlation is **NOT an established generalizable morphology-to-pollinator functional rule across unseen genera**. The burden of proof for a true pollination syndrome still requires independent local pollinator interactions/effectiveness and stricter phylogenetic controls.

## Interpretation and decision

**Supported evidence-level finding:** recordable floral marginal diversity and colour–architecture coordination are separate dimensions, and colour–architecture coordination co-varies geographically with globally *documented* interaction breadth of the resident plant genera, especially in northern midlatitudes. The association persists after several sampling-depth diagnostics and within-source-effort genus-label placebo swaps.

**Falsified shortcuts:** 1) 'More observed channels on an island = more true pollinator guilds' (recording opportunity dominates); 2) 'a documented channel missing = biologically absent' (most are unresolved); 3) 'the named bee/butterfly/bird template identifies a real pollinator' (previous 95.1% 2+ template overlap); 4) 'the genus-linked morphology pattern predicts GloBI breadth for unseen genera' (fixed out-of-genus model comparisons do not improve).

**Not supported:** historical island pollinator loss, local moth-versus-butterfly shift, independent ecological functional diversity, floral adaptation to a new pollinator guild, visitor-to-plant phenotypic causal selection, or pollen-transfer success. H5 in upstream `island` also found general isolation-associated pollen limitation but did NOT identify its proposed North–Tropical pollinator mechanism. No new GloPL mediation is claimed here.

### Next genuinely identifiable macro-scale step

Collect independently sampled **local species-level pollinator assemblages** (not GloBI global source breadth), with a standardized protocol, known survey visits and explicit adequate nondetection across multiple channels, on islands already carrying >=30 Direct complete floral species. Require per-channel × region both observed positive and adequately sampled negative support and enough **spatial/source-disjoint** islands *before* testing colour×architecture associations. Independent direct pollen deposition or seed set is needed before claims about realized service. Do not chase a positive effect by treating unobserved GBIF insect taxa as absences or adding unrestricted high-dimensional named-guild scores.

## Public reproducibility and confidentiality

Reusable, input-gated modules:
- `src/fd/pollinator_assemblage_bridge.py`: true retained/disrupted/unresolved validation, strict channel×region support gates and spatial block cluster OLS.
- `src/fd/globi_morphology_falsification.py`: reference-effort-stratified GloBI genus reassignment and entire-genus heldout predictive comparison.
- Synthetic tests are under `tests/test_pollinator_assemblage_bridge.py` and `tests/test_globi_morphology_falsification.py`. They run in public CI without restricted source records.
- Source-pinned empirical original records and intermediate island x genus x visitor data remain in private local inputs; public release is aggregate-only: `exact_island_pollinator_channel_support.csv`, `globi_documented_genus_breadth_sensitivity.csv`, `outofgenus_globi_guild_breadth_prediction.csv`, `pollinator_bridge_aggregate_receipt.json`.

GloBI upstream interaction records require attribution to their original contributors, and its compiled documented interactions are a distinct evidence class from real, current island visitation. No raw licensed user source data, island IDs, original taxon IDs, local donor labels, or georeferenced observation rows are committed.
