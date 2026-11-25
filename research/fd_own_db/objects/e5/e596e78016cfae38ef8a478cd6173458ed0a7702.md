# 13 Nov 2026 | Experimental GloPL pollen limitation × island floral FD

**Original SHA-verified 2018 GloPL and Direct recorded floral traits — completed original real-data execution, no confirmed mechanism.**

Original source execution: [island Actions #43357788907](https://github.com/zuizui0223/island/actions/runs/43357788907) **SUCCESS**, anonymous [artifact #14007176415](https://github.com/zuizui0223/island/actions/runs/43357788907/artifacts/14007176415). Direct flower source and original GSHHG matched original source hashes. Public GloPL `Output/GloPL.csv`, original commit `abbe193770981b2c61713e75a0f01e004bdded8e`, SHA256 `264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de`.

## Scientific question and two non-interchangeable levels

*Are island floras with lower source-relative floral functional diversity exposed to stronger experimentally measured pollen limitation?* This is a cross-sectional **island-level ecological association**, **not** the effect of historical pollinator decline, effective visitation identity, within-species floral evolution, or experimental causal mediation.

The prerequisite [outcome-blind GloPL flora co-location](https://github.com/zuizui0223/island/actions/runs/43356660466) previously found 31 of the original 928 eligible Direct floral islands with at least one exact-geolocated GloPL experiment; **zero exact same-island same-species plant names** were common to both the original morphological and experimental response sources. Therefore **NO individual-plant pollination-morphology effect is identifiable**, and the analysis cannot claim floral traits reduced their own pollen limitation.

The new runner tested an **independent island-level null** with source-type labels, original island area/distance and broad region controls. Source-available Direct floral species >=30 on each *island and its five geographically nearest GIFT candidate mainland source regions* were uniformly rarefied to exactly 30 species in each of 30 deterministic draws. Morphology: colour, floral form and symmetry from original species-direct high/medium single-state records, without genus imputation. Outcomes: floral marginal Rao Q, **actual recorded triple-profile joint Rao Q**, open-radial flower prevalence, each as island minus its separate candidate source flora. Geography and trait outcomes were formed and the **25-island / 15-spatial-block / 20-experiment-study / 10 northern-mid / 10 tropical** support gate passed *before reading any GloPL effect sizes*.

GloPL response `PL_Effect_Size` is a published log hand-pollen-treatment vs natural response, not a count of missing bees. Effects within the same original **publication × island** are averaged, and each island's mean across available publications is the independent ecological response unit; one island receives one regression row, regardless of how many published measurements exist. The model controls region, original log mainland distance and log island area, with whole 10-degree geographic-block resampling and leave-one-block sensitivity. Original source and GloPL source record-level identities are not exported.

## Real-data results

Frozen trait/source & GloPL data yielded:
- **30** matched source-Q, finite GloPL-response islands, across **26** physical spatial blocks; **94** distinct experiment publication keys and **102** publication × island cells.
- Of 30 islands: **13 northern midlatitude, 12 tropical, 5 southern extratropical**. The earlier 31-island *outcome-blind geolocated coverage* differs by one island from the final source-Q/finite-outcome model; never silently identify 31 as 30.
- Per SD difference in original recorded *island-minus-source* floral functional metric, adjusted slope of the island's mean experimental pollen-limitation response:

| Source-relative island floral FD exposure | Adjusted original per-SD slope | 1,999 whole-spatial-block conditional 95% interval | Share of leave-one-block-out estimates with original sign |
| --- | ---: | ---: | ---: |
| **Marginal categorical Rao Q** (primary) | **−0.00436** | **[−0.50081,+0.28724]** | 65.4% |
| **Joint colour×form×symmetry Rao Q** (primary) | **+0.06373** | **[−0.24912,+0.31854]** | 96.2% |
| **Open-radial share** (secondary) | **−0.14428** | **[−0.41494,+0.17966]** | 100.0% |

**All three block intervals include zero**. The first primary response is near-zero point estimate; the second joint-Q positive point estimate is imprecise and cannot establish an inverse functional-decline gradient. The open-flower negative point estimate also lacks uncertainty support. Leave-one-block sign agreement alone is **not statistical validation**.

The pre-existing GloPL full global experiment finding of higher pollen limitation with mainland distance is a different all-site model with 1,248 sites, not demonstrated here as a mediator of any fd floral functional island-versus-source difference.

### Stronger experimental-definition sensitivity: hand-pollen supplementation only

Restrict the original experiment rows to the GloPL `PL_Effect_Size_Type2='Sup'` source coding. The eligible cohort declines to **19 islands / 18 spatial blocks**, **59 original publication keys / 66 publication-island cells**, with **10 northern midlatitude, 5 tropical and 4 southern**. It fails the frozen minimum 25 islands and 10 tropical contexts, so **no pollen-only rescue regression is reported**. This reinforces the sensitivity of the 30-island general-effect model to experimental endpoint definitions.

### Unresolved inferential issues

The original GloPL measurement definitions and supplement protocols vary across publications; the simple exploratory FD model averages finite reported responses instead of the parent GloPL model's original measurement-type fixed effects and publication-total weighting. Its whole-block bootstrap is *conditional*, not historical geographic or botanical source-clade independence. GIFT five nearest potential mainland pools are NOT measured actual ancestors or known counterfactuals. The sampled 30 islands are a highly selected intersection of floral trait completeness and pollen experiment availability; 30/928 is not a global random island sample.

The absence of a strong ecological correlation is not evidence that **pollinator loss has no effect**: real loss, change in flower-specific visits, evolutionary trait response and time-matched effective pollen delivery were not measured. Conversely, existing direct island flora morphology/SC/selfing source-relative traits, 11-island actual EuPPollNet recorded visit networks and GloPL experimental effects cannot be assembled into a confirmed causal ecosystem-service chain without common plants/sites/times and study-method control.

## Decision and next falsification

**Observed positive fd findings retained:** source-relative all-recorded floral marginal/joint diversity reduction and open flower prevalence in well-sampled northern midlatitude islands; on identical plant species, distinct source-relative SC vs selfing mating-system representation; documented link vulnerability to hypothetical all-bee disappearance in original European networks. These remain separately supported **different ecological estimands**, not three proven causal steps.

**GloPL result classification:** `ecological_pollen_limitation_to_flora_FD_association_not_supported` under the source-pinned data and given sample/measurement constraints. Do not fit a bespoke mediation path to force the result. The next genuinely ecological test requires matched native plant floras and **same plant population** measured visitation, floral traits, and pollen addition; or multiple genuinely independent island-wide standardized service survey populations.

## Auditable source

Reusable [fd analytical module](../../src/fd/glopl_floral_fd_association.py) and [synthetic validation](../../tests/test_glopl_floral_fd_association.py), original rights-safe [island source runner](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/scripts/run_fd_glopl_flora_association.py) + [workflow](https://github.com/zuizui0223/island/blob/analysis/fd-functional-guild-impact-20261113/.github/workflows/run-fd-glopl-flora-association.yml) + [source Actions #43357788907](https://github.com/zuizui0223/island/actions/runs/43357788907). The public original-source output is a single anonymous aggregate JSON. All individual experimental species, site/study IDs, island IDs, GIFT mainland sources and georeferences remain within the protected input process. This exploratory extension remains Draft and is not an established biological effect of actual guild loss.
