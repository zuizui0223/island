# 2026-11-16 — Original 4,505 island floral-source supply × independently qualified GBIF guild support

**Status: COMPLETED real source execution; observational composition and support audit, NOT an identified ecological guild effect.** Original protected [island Actions #43834081663](https://github.com/zuizui0223/island/actions/runs/43834081663) SUCCESS (including 11 synthetic gate tests and source-hash-frozen botanical/GIFT/GBIF joins). Original aggregate-only artifact [#14209568425](https://github.com/zuizui0223/island/actions/runs/43834081663/artifacts/14209568425), JSON SHA256 `c96ca00a4ddd69d885a3a5129d543fc13a3580717b27a38e7c1b858c7f721747`. No original taxon, island, geography, source-region, occurrence or observer identifiers were published. Original runner [island script](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_source_guild_support.py) and [fd strict gate](../../src/fd/source_guild_support.py).

## Exact design and nonidentification

Global original GIFT flora **4,505 island checklists**, including unresolved geographic/status cells in the denominator, **748** original native GIFT mainland source entities, **five** original directly documented morphological sectors (zygomorphic; tubular/bell; open-radial; deep tube; actual SAME accepted species zygomorphic×tubular/bell), **five** independently qualified GBIF guild channels, **112,625 original private island×guild×sector cells**, **125 region×guild×sector potential comparisons** (including one geography-unresolved stratum). Separate reproductive categories were NOT treated as flower morphology. Each sector uses its own >=5 directly coded island species and >=10 independently directly coded *mainland candidate-only* species from five nearest GIFT regions. All focal island names are excluded from each candidate-only pool. These are plausible native mainland floristic source checklists, not actual seed arrivals or historical ancestors.

**Critical fixed data semantic correction:** Original qualified GBIF `disrupted` = source-available guild that was adequately *not detected* under GBIF documentation effort, NOT ecological guild extinction, historical decrease, or a bird/bee pollinator shift. `retained` = documented guild presence. `unresolved` never means true absence; `structurally_absent` is a separate source-state. Previous code had incorrectly expected a string `adequately_nondetected` and was repaired before this source run.

Exploratory *support gate*: >=10 retained and >=10 source-available adequately nondetected original islands; >=10 10-degree physical blocks EACH state; >=3 first-GIFT-source signatures containing BOTH; nondegenerate overlapping source-sector candidate frequencies. A stricter requirement >=5 blocks containing both states was transparently added **after the first support counts were read but before any status-effect comparison**. It is not preregistration. For descriptive matched comparisons, the actual same island plant trait denominator, source choice, guild qualification, and region are held fixed.

## Main discovered support: three slots, one underlying guild contrast

**3 / 125** source×guild×sector slots pass the strict source/block/n5 source10 gate; **122 / 125** do NOT. These are three correlated morphological outcomes for **the SAME tropical flower-visiting bird channel**, not three independent ecological confirmations.

| tropical morphological sector | original GBIF retained | adequately nondetected | retained blocks | nondetected blocks | physically mixed blocks | first GIFT-source signatures with both |
|---|---:|---:|---:|---:|---:|---:|
| zygomorphic symmetry | 593 | 13 | 79 | 11 | 10 | 8 |
| tubular or bell | 504 | 12 | 76 | 10 | 9 | 8 |
| open-radial form | 504 | 12 | 76 | 10 | 9 | 8 |

Of 125 original strata, gate failures overlap: 122 insufficient state-specific blocks; 120 insufficient *mixed* blocks; 119 insufficient both state counts; 100 insufficient mixed source signatures; 81 no positive source-sector frequency overlap. A failure is not proof that ecological effects are zero.

## Measured same-island, source-referenced, mixed-block contrasts do not identify a guild effect

For each of the three admitted tropical visitor-bird slots, define actual island documented sector share `y` and independently documented GIFT-source candidate-only share `s`. Contrast `(y-s)_retained - (y-s)_adequately_nondetected` using **only blocks that contain both guild states**; 1,999 mixed-cluster bootstrap draws, whole-block resampling, not random guild assignment. Positive values mean source-relative morphological share higher among documented bird-present islands.

| sector | mixed physical blocks | physical-block mean, percentage points | conditional 95% bootstrap, pp | first-source-group mean, pp and CI | joint first-source×block mean, pp and CI |
|---|---:|---:|---|---|---|
| zygomorphic | 10 | +2.06 | [−7.88, +9.31] | +1.77 [−3.89, +8.41] | +4.20 [−7.29, +13.90] |
| tubular/bell | 9 | −3.74 | [−13.81, +6.14] | −3.32 [−12.07, +3.97] | −3.73 [−14.70, +7.03] |
| open-radial | 9 | +1.46 | [−5.38, +8.45] | −0.58 [−6.28, +5.93] | +2.71 [−5.52, +10.10] |

**All 9 conditional, mutually nonindependent source-relative intervals cross zero**. The archive also separately reports the actual island-share and the source-share component for every grouping, exact diagnostic block sign flips, and leave-one-mixed-group ranges. These post-selection block sign diagnostics are not randomized treatment p-values; no source×guild *interaction coefficient* was fit and no claim of absence of ecological effect is licensed. Three correlated selected endpoints among 125 screened slots preclude picking nominal unadjusted p-values as discovery. Nondetected islands are too sparse for a reliable multi-guild or historical mechanism interpretation.

## Descriptive SOURCE-ONLY morphological supply versus island composition (all recorded provenance)

For context, category share difference = `island - GIFT candidate source-only`, original 10-degree **block-equal** recorded-species proportions. These are a new direct category-frequency comparison, not previous marginal/joint Rao Q, not a source-listing-breadth/genus-preserving causal filtering null, and not independent plant-population phenotypes:

| region | zygomorphic pp [95%] | tube/bell pp [95%] | true deepest tube pp [95%] | open radial pp [95%] |
|---|---|---|---|---|
| Northern midlatitude | −5.00 [−6.12, −3.85] | +0.15 [−1.32, +1.71] | +1.88 [+0.45, +3.46] | −0.35 [−2.04, +1.37] |
| Tropical | −8.44 [−10.14, −6.69] | +2.89 [+0.90, +4.73] | +9.96 [+8.11, +11.95] | +0.15 [−1.50, +1.84] |
| Southern extratropical | −10.43 [−13.55, −7.29] | −12.59 [−17.60, −7.68] | +1.79 [−2.13, +5.71] | +8.74 [+4.35, +13.54] |
| Northern high latitude | −15.22 [−18.56, −11.63] | +1.41 [−1.31, +4.07] | −3.81 [−6.15, −1.12] | +6.43 [+2.41, +10.39] |

Cohorts differ by individual trait: tropical zygomorphic **954 islands / 101 blocks**, tropical tubular/bell **779 / 96**, tropical deep tube **733 / 92**; southern form **157 / 37**, southern symmetry **179 / 37**; high-latitude symmetry **260 / 54**. High-latitude observed source-relative morphological composition is NOT invariant. Southern extratropical is not ecologically uniform and includes geographically dissimilar mainland sources. These descriptive intervals do NOT correct for multiple-testing, genus/source listing-breadth, incomplete Direct label propensity or native-provenance selection.

## Evidence and claim ceiling

What is supported: source-only candidate floral-category composition and actual documented island flora differ in **regionally distinct DIRECTIONS**; the full original 4,505 flora framework permits a narrowly supported observational bird-availability comparison without reducing the main paper to historic 928 complete cases; that narrowly supported contrast shows NO robust source-referenced association under physical-block, donor-source or joint blocking.

What is NOT supported: a general reduction of morphological space on isolated islands; a causal source abundance × guild extinction mechanism; any historical change of local effective visitors, lineage adaptation or pollen-delivery/seed-set mechanism. The comparative donor pool uses observed species literature-wide state labels; native/introduced provenance is mostly unresolved; unrecorded mainland candidate taxa are not measured dispersal failures. Existing source genus/finite hypergeometric null and strong differential trait-documentation bias remain more serious rivals than bird guild status.

Next stringent analysis should compare the identical valid island cohort under **geo5 vs climate10** donor definitions; separately test strict-native/endemics and Direct trait-documentation filters; seek **actual same-species island/mainland flower-visitor observations** rather than broaden a noisy GBIF guild proxy. No unobserved pollinator or floral traits should be filled to force a guild interaction.

Related full source-null context: [direct source finite-sampling genus null](../20261115_global_island_functional_beta/DIRECT_SOURCE_FINITE_SAMPLING_NULL_GENUS_FILTER_4505.md) and [independently qualified historic guild/pollination limitations](../../docs/SOURCE_GUILD_INTERACTION_IDENTIFICATION_CONTRACT.md). Both original fd and island PRs remain DRAFT.
