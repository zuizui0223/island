# 15 November 2026 — Do floral genus trait categories transfer across biogeographic regions after strict taxon-name holdout?

**Original botanical source-protected numerical analysis SUCCESS:** [island Actions #43700248327](https://github.com/zuizui0223/island/actions/runs/43700248327), anonymous [artifact #14156231290](https://github.com/zuizui0223/island/actions/runs/43700248327/artifacts/14156231290); extracted original anonymous JSON SHA256 `c38d753e3d9242bc37f288f29996d168b2889ec70212a4be765b40c090a2dded` (741,887 bytes). Original GIFT 4,505 accepted botanical island floras and original high/medium single-state species-Direct trait archive/first GIFT potential source assignment were SHA256 audited; physical GSHHG 8,264-island regional/geographic crosswalk rebuilt from original independent sources; mathematical donor-name exclusion/equal genus source counts/negative proper scoring controls passed synthetic CI before original data use. The earlier source-only comparison run [#43700017352](https://github.com/zuizui0223/island/actions/runs/43700017352) succeeded, and the new full run adds a **stronger external-region genus versus original genus-blind proper predictive score** on exactly the same real directly documented plants.

## Global frame and distinct question

The **entire original 4,505** GIFT botanical island population is retained, including **2,205 northern-midlatitude**, **1,558 tropical**, **317 southern-extratropical**, **424 northern-high-latitude** and **one geography-unresolved original flora**. All **10** separately original Direct recorded flower/reproductive traits and **three** TRUE same-species two-trait endpoints are treated separately. Original direct colour×form ≥5 source flora opportunity remains **2,324 all-recorded**, **327 source-native**. No old n30 928 3-trait-complete-only filtering, artificial trait cross-products, or 17D continuous KDE as biological volume.

The previous [global original genus decomposition and single-island whole-flora-name holdout](GENUS_ASSEMBLY_EXACT_DECOMPOSITION_AND_NAME_HELDOUT.md) showed category frequency prediction improves with original plant GENUS identity, especially whole-flower floral symmetry, but name-withheld donor phenotypes were drawn primarily from other original island plant floras within the SAME geographic region. **Here we test true across-region transfer using botanical species names observed ONLY OUTSIDE the target original region**, and compare them against different, matched-count same-genus plant species observed ONLY INSIDE the target region.

### Source-set names and leakage gates

For a focal original geography R, source unique accepted species are partitioned **before** reading test flower outcomes:

1. **Local-exclusive source plants:** recorded on original plant floras in R and NEVER on any original botanical flora in another independently labelled GSHHG geography (the one geography-unresolved flora also excludes local exclusivity).
2. **External-exclusive source plants:** recorded in one of the OTHER THREE independently labelled GSHHG geographic regions, and on NO original botanical flora anywhere in R.
3. **Region-shared species:** any accepted species documented on original plant floras in R and some other geographic regime. Entirely **excluded from BOTH source donor sets**. They are NOT used to manufacture spurious regional cross-validation on an identical plant species.
4. Additionally remove **EVERY full botanical plant NAME** of the tested island (including species with no Direct flower annotation or with unverified nativeness) from each model's candidate source before predictions. Original source-native is a separate nested subset, NOT independently replicated by another flora survey.
5. For every target genus g represented among original directly recorded target plants with at least one valid source species in each donor universe, choose **exactly the same number of original DISTINCT source species** from local-exclusive and external-exclusive, capped at 20/genus, ranked by frozen accepted-name SHA256 hash independent of measured flower states or model performance. If one side lacks a congener, BOTH models use the SAME external-source category prior; no zero-state or missing-visitor classification is fabricated.
6. Both models use the **same outside-original-region-exclusive genus-blind trait-category prior** with Jeffreys 0.5 pseudocount/state, and both genus models shrink 5 source-label prior counts toward that single prior. Thus a positive local model advantage cannot be explained simply by using more same-genus training species than the external model; it still may be original trait-provider/detection/source opportunity bias.

Every original focal island must have ≥5 original Direct trait carrier species for a particular endpoint. For the strongest **matched-congener direct-source** evaluation, ≥5 of those carriers must have exact name-disjoint donors in BOTH source worlds after full original target-island name removal. The percentage of all original trait carriers surviving these strict steps is ALWAYS disclosed; no pooled "all 4,505 are externally confirmed" claim.

### Three different proper/aggregate score contrasts

All tested real plant categories and focal islands are IDENTICAL between compared models.

- **Local minus external genus:** `logloss(external genus) − logloss(local genus)` per plant in an original focal genus with documented donors in both regions. Positive means local-exclusive species genus states generalize better.
- **External genus minus genus-blind reference:** `logloss(external source category prior) − logloss(external matched genus)` on EXACT the same supported original target species. Positive indicates plant genus identity still contains useful SOURCE-WIDE categorical information when even the original source donor species are region-external AND plant-name-disjoint.
- **Focal island category-frequency squared-loss comparison:** same observed directly documented island category-frequency vector, two balanced genus source predictions; pair frequency D² gain is a different estimand from per-species proper log/Brier. The report preserves positive and negative signs and an alternate original Brier-score source check.

Source physically indexed 10° group-equal and independently grouped **first-original-GIFT-native potential mainland source** 999-resample conditional intervals are provided separately, not phylogenetic or statistically independent causal replication. Original full island study samples share genera, floras and Direct trait labels, and the sixteen nested source/status/region comparisons are NOT independent evolutionary events.

## Primary original all-recorded region results, physical 10° block-equal log loss

**A. Regional exclusive source advantage over matched across-region source** (positive if local source predicts better; unit nat/plant in matched original target genera):

| True Direct endpoint | Northern midlatitude (n original eligible matched floras) / local−external 10° block equal, 95% | Tropical (n) / local−external 10° block equal, 95% |
|---|---|---|
| Primary flower colour | **1,597**, **+0.2612 [+0.2275,+0.2985]** | **712**, −0.0168 [−0.0365,+0.0032] |
| Floral architecture (form) | **1,436**, **+0.3537 [+0.3129,+0.3937]** | **483**, **−0.1195 [−0.1510,−0.0888]** |
| Floral symmetry | 1,595, **−0.0035 [−0.0054,−0.0017]** | 756, **+0.0168 [+0.0101,+0.0257]** |
| Flower-size class | 1,327, +0.0732 [+0.0534,+0.0924] | 558, +0.1322 [+0.1057,+0.1593] |
| Ordered floral tube-depth class | 1,290, +0.0513 [+0.0359,+0.0662] | 441, −0.0069 [−0.0221,+0.0077] |
| Inflorescence display class | 1,631, +0.0068 [−0.0058,+0.0194] | 744, +0.0456 [+0.0273,+0.0646] |
| Original self-incompatibility categories | 1,463, +0.0427 [+0.0338,+0.0516] | 406, −0.0290 [−0.0448,−0.0149] |
| Autonomous selfing capacity | **369**, +0.0505 [+0.0368,+0.0646] | **82**, **−0.0536 [−0.0825,−0.0200]** |
| Mating system | 875, +0.0759 [+0.0574,+0.0935] | 38, −0.2825 [−0.3176,−0.2444] |
| Cleistogamy | 139, +0.0000 [0.0000,0.0000] | **no congener-supported target** |
| True same-plant colour×form | 1,102, +0.6803 [+0.5961,+0.7656] | 218, −0.0717 [−0.1868,+0.0332] |
| True same-plant colour×symmetry | 1,336, +0.3411 [+0.2958,+0.3884] | 467, −0.0708 [−0.0959,−0.0469] |
| True same-plant form×symmetry | 1,150, +0.4263 [+0.3673,+0.4901] | 306, +0.0451 [+0.0241,+0.0672] |

**Do not select positive rows post hoc as a new confirmed hypothesis family.** These are many trait/region/source/status contrasts, paired conditional group-resampling ranges, and the outcome is a source-annotated species-level category not a measured plant population change. Large raw category log-loss differences across axes with different category vocabulary size K are NOT directly comparable in biological effect magnitude. In addition, cross-region exclusive donor support is nonrandom and may preferentially include genera with stable or easy-to-describe floral forms.

**B. Does EXTERNAL-only genus information improve prediction above a shared genus-blind outside-region prior?** (positive = transferable genus category information on previously unseen accepted plant names; same target islands and matching):

| True Direct endpoint | Northern midlatitude (block equal proper log-loss improvement) | Tropical (block equal proper log-loss improvement) |
|---|---:|---:|
| Flower colour | **+0.4184 [+0.3912,+0.4449]** | **+0.3191 [+0.2997,+0.3402]** |
| Floral architecture | **+0.4379 [+0.4017,+0.4748]** | **+0.7858 [+0.7496,+0.8218]** |
| Floral symmetry | **+0.3462 [+0.3389,+0.3536]** | **+0.3153 [+0.3048,+0.3258]** |
| True colour×form | +0.2320 [+0.1533,+0.3167] | +0.7171 [+0.6316,+0.8037] |
| True form×symmetry | +0.7692 [+0.7227,+0.8145] | +1.1740 [+1.1250,+1.2221] |
| Self-incompatibility categories | **−0.0220 [−0.0331,−0.0116]** | +0.0108 [+0.0034,+0.0189] |
| Autonomous selfing | **−0.0587 [−0.0801,−0.0381]** | +0.0037 [−0.0227,+0.0269] |

The **separate original first-GIFT-source-group analysis** shows the same qualitative external genus-prediction advantage for *floral symmetry* (north +0.3394 [+.3326,+.3463], 96 original first-source groups; tropical +0.3143 [+.3034,+.3251], 78 groups) and *floral architecture* (north +0.4524 [+.4231,+.4811], tropical +0.7972 [+.7610,+.8332]). **Northern original autonomous selfing external genus vs blind remains NEGATIVE** with source-group mean −0.0498 [−0.0710,−0.0305]. Tropical autonomous selfing first-source group mean +0.0003 [−0.0267,+0.0271], not identifiable as positive under its documentation support. These group intervals are conditional resamplings of an observational source reference, NOT phylogenetically independent genus evolutionary histories.

### Direct trait-provider opportunity: do not confuse estimated value with representativeness

Strict original matched congener support is unequal:

| Endpoint | Northern-midlatitude distinct focal floras with ≥5 carriers with BOTH exclusive source groups / proportion of ALL original Direct target carrier records matched | Tropical |
|---|---|---|
| Floral symmetry | **1,595 / 46.6%** | **756 / 42.3%** |
| Floral form | 1,436 / 34.0% | 483 / 32.0% |
| True colour×form | 1,102 / 30.0% | 218 / 28.0% |
| Strict original self-incompatibility | 1,463 / 31.3% | 406 / 22.9% |
| Strict original autonomous selfing | **369 / 7.6%** | **82 / 8.8%** |

For symmetry among the matched original target Direct carriers, an **observed MULTI-CATEGORY reference genus** occurs in only **5.17% local /0.93% external** carrier source opportunities in northern midlat, and **4.72% local /7.49% external** in tropical. Northern-high-latitude has only **34** exclusively local Direct documented symmetry donor species and **zero** multi-state matched-genus variation, even though it has 190 donor-matched floras with ≥5 target congeners and external genus vs blind original log-loss improvement +0.2297 under physical-block equalization. **High predictive accuracy in coarse floral symmetry thus may stem largely from genuinely or textually conserved singleton-category genus identities, not source-confirmed repeated ecology or regional island evolution.**

Autonomous-selfing original Direct carrier support is especially narrow; although tropical all-recorded n5 direct island scope is **610**, just **82** have ≥5 potential balanced genus-supported carrier records after strict region/source species-name withholding. Strict source-native n5 islands have still lower overall opportunity. Therefore the lack of transferable autonomous selfing effect is NOT evidence of universal island breeding-system plasticity, selfing absence, or no pollinator-driven selection.

### Southern extratropical and northern highlatitude do not disappear

The ORIGINAL complete 4,505 population accounting includes two further geographic regions and their exact 13×2 trait/status supports, with no posthoc omission:
- Original southern extratropical floral symmetry has **160** eligible original all-recorded matched-source floras; external genus vs blind 10° block-equal logloss +0.2924 [+0.2686,+0.3160], and local−external +0.0017 [−0.0001,+0.0033]. Floral form has 119 and external predictive +0.6458 [+0.5268,+0.7523], local−external −0.0843 [−0.1992,+0.0406]. The genuine Direct same-plant form×symmetry source has only 82 matched floras and leaves substantial uncategorized source opportunity.
- Original northern-high-latitude floral symmetry has 190 matched original all-recorded floras and external genus-vs-blind +0.2297 [+0.2068,+0.2537], but its 34 exclusively local available Direct symmetry species have no matched-genus source variation. Northern-high-latitude floral form 129 matched floras, external benefit +0.4119 [+0.3292,+0.5039], whereas local−external is **−0.1332 [−0.2637,−0.0112]**. Strict original high-lat autonomous trait has **no** source-congener support across both regions, and true colour×form only 17 matched floras, BELOW the gate. The single unresolved source geography flora is always counted but is not assumed to be an independent region donor/test.
- The original strict-native nested subsamples are fully included in the machine-readable anonymous 100 reported source/status/trait/region rows, but many lack ≥30 matched original floras or ≥10 independent physical blocks; these report **NOT IDENTIFIABLE**, not zero biology. The separate original 425 geology-typed oceanic/continental source flora analysis is a different nonrandom source population and is NOT the 424 northern-high-latitude GSHHG category.

## Scientific interpretation and publication limit

**Positive descriptive source-backed generalization:** even when original observed accepted source plant species occur **only in other nonfocal GSHHG regions** and every original target flora accepted plant name is excluded from training, original botanical GENUS identity provides strong category prediction of the flower *symmetry* and *architecture* of previously unseen accepted plant names relative to a genus-blind external source profile. It is plausible to discuss coarse botanical morphological conservatism across regions. Whether this is biological conserved phylogeny, a morphology-category definition, regional source/niche ecology or similar botanical description provider is NOT identified without independent provider and phylogeny.

**Transferability is trait-dependent:** original STRICT reproductive categories, particularly autonomous selfing, do not reproduce the strong cross-region genus-specific proper prediction improvement in the same source study, but available matched founder plant evidence is much sparser. The north-midlat floral forms differ in source regional-specific shape composition (local beats external), while tropical forms transfer better from outside-exclusive congeneric plants. Those opposite signs mean geographic source conditions and original provider composition are meaningful rival drivers; they cannot be interpreted as adaptive morphological change after island colonization.

**Do NOT claim** that geographical source species' original Direct trait labels were independently measured on the mainland, that cross-region geographical contrast demonstrates gene flow, independent phylogenetic drift, adaptive change, pollinator-syndrome loss, pollinator transition, autonomous selfing selection or actual pollen delivery. These are DOCUMENTED species-wide morphological categories assembled from overlapping original trait literature and taxonomic names. Original 4505 island biodiversity, per-trait Direct support and ~425 independently geology-typed source floras remain separately labeled.

**Next actual identifying data:** original trait provider-aware and external-source independent phylogenetic split, genus multistate annotation-quality audit, more independent mainland botanical morphological studies and island plant–pollinator interaction/fitness data. A geographically exclusive source list is NOT the same as phylogenetically independent botanical trait measurement. Keep fd PR #16 and island PR #326 **DRAFT**.

[Protocol fixed before original external-region source results](../../docs/INTERREGIONAL_GENUS_TRANSFER_FROZEN_DESIGN.md) · [Original source-protected runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/fd_interregional_genus_name_exclusive_transfer.py) · [Pure synthetic source-name exclusion and proper scoring tests](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/tests/test_fd_interregional_genus_name_exclusive_transfer.py) · [Original CI](https://github.com/zuizui0223/island/actions/runs/43700248327). Full source species, genus and location IDs remain protected in original GitHub CI, only anonymous regional/trait aggregates published.
