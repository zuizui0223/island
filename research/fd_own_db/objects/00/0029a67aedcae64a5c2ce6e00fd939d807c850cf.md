# Family-conditioned rare floral category support on 928 islands — 13 November 2026

**Status:** original-source test completed; exploratory categorical morphology, not demonstrated within-species evolution or pollinator-mediated selection. Source [island Actions #43406742054](https://github.com/zuizui0223/island/actions/runs/43406742054) **SUCCESS**, anonymous [artifact #14026919625](https://github.com/zuizui0223/island/actions/runs/43406742054/artifacts/14026919625). Scientific [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and protected source [island PR #326](https://github.com/zuizui0223/island/pull/326) remain Draft.

## Source and eligibility

The original 928-island flower-colour × floral-form × symmetry data comprise **2,234 original high/medium Direct species** and **100 distinct ACTUALLY recorded three-trait categories**. An audit of 297 original progressive source CSV/TSV files found **no direct species→family table**. An independently frozen Catalogue of Life Extended Release source, Island Plant Database 2.0.0-alpha2 taxonomy, original run 34668274725, original taxonomy crosswalk SHA256 **079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3**, provides accepted-name EXACT matching candidates; candidate family was required to agree with the original reported GBIF family and genus.

**2,190/2,234** Direct source species have strict double-consistent candidate family mappings, covering **119 families and 882 genera**. This source explicitly marks all automatic taxonomy matches **unreviewed/not promoted**. We do not promote candidates as verified taxonomy or reconstruct a phylogenetic tree.

The 44 unverified species are unevenly distributed: **one** of 40 rare-profile source species, **seven** of 460 middle-profile carriers, and **36** of 1,734 common-profile carriers. **One of 29 globally rare recorded profiles has no family-qualifying source carrier**, whereas none of 25 common-profile categories is entirely missing from mapped species. The family-filtered island sample and its fixed 100-profile measurement reference therefore cannot be assumed identical to the original all-2,234 species cohort.

## Original design and what it does not identify

Fixed reference categories: 29 rare categories with 1–2 source Direct carrier species each, 46 middle, 25 common with ≥26 carriers. Metric is **coverage of ACTUALLY observed original categories within one original categorical colour/form/symmetry mismatch**, NOT PCoA 17D Gaussian KDE or pollen function.

Each original island was compared to original five-nearest GIFT candidate mainland native source floras, after excluding **every original island accepted species** from the source. The **family-identity n5** experiment selects five EXACT same family IDs on island and source, one independent accepted source species/family per side, for 40 rarefactions; n8 sensitivity. A stricter source control requires **different genus identities inside each chosen same family**, even though family identities are held fixed. A later **post-hoc** difference-of-estimands samples five families per side with equal **family count** but independently selected family IDs on the SAME original island. All sample sizes and family source IDs are outcome-independent. Reference k-null holds EXACT observed number k of distinct source-profile types and uses 399 uniform or original species-frequency-proportional k-category draws; 1,999 equal physical block and separately first-nearest-GIFT-signature bootstraps.

Original all-recorded sample: tropical **106 islands/42 blocks/36 GIFT donor signatures** for n5 family-identity match; northern midlatitude **743/53/66**. Each result below is a fraction of original globally recorded category-reference support, **island minus potential mainland source**; negative is a category-neighborhood shortfall.

## Tropical 106 islands: is any independent family-conditional rare category deficit present?

Primary response: **rare-minus-common reference category neighborhood coverage, corrected for k**, block-equal mean and 95% conditional interval:

| Original source control | Islands | Block mean [95%] | Source-signature mean [95%] |
| --- | ---: | --- | --- |
| **Same five family IDs** | **106** | **−0.01188 [−0.03015,+0.00591]** | −0.01193 [−0.03152,+0.00621] |
| Same five families, but DIFFERENT genus IDs in every matched family | **106** | **−0.01065 [−0.02986,+0.00851]** | −0.00817 [−0.02978,+0.01343] |
| Same eight family IDs | 102 | **−0.00507 [−0.02305,+0.01244]** | −0.00488 [−0.02522,+0.01592] |
| Same eight families, different genus ID within family | 99 | −0.00107 [−0.02100,+0.01963] | −0.00474 [−0.02601,+0.01731] |
| **Equal FIVE family units, family IDs free to differ** | 106 | **+0.00119 [−0.01794,+0.02065]** | +0.00271 [−0.01794,+0.02397] |
| Equal EIGHT family units, IDs free | 106 | −0.00124 [−0.02422,+0.02100] | +0.00101 [−0.02427,+0.02519] |

**No tropical conditional 95% interval excludes zero.** On exactly the **same 106 physical islands**, matched-five-family-ID minus family-count-only gap is **−0.01307 [−0.03279,+0.00643]** (block), **−0.01463 [−0.03588,+0.00648]** (source). The eight-family paired difference **−0.01126 [−0.02937,+0.00601]** also includes zero. Thus **neither a residual within-family trait deficit NOR the isolated contribution of family identity itself is established**. Non-rejection alone is not evidence of equivalence or proof of family-driven filtering.

Tropical rare-only coverage in the same-five-family n5 experiment is **+0.00754 [−0.00390,+0.01846]**, whereas common-only coverage is **+0.01975 [+0.00702,+0.03252]**; the negative rare−common gap can therefore reflect stronger common-profile coverage and a reference ceiling rather than absolute missing rare traits.

### Family-qualified 30-species Shapley accounting

For each n30 island and independently sampled source, divide the reference coverage into a shared-family component and a component from families recorded on **only one side**; overlap is split exactly 50/50 in an additive two-player Shapley allocation, NOT a causal lineage model.

| Region | Eligible family-qualified islands/blocks/sources | Total rare support difference | Shared-family part | One-sided-family part |
| --- | --- | --- | --- | --- |
| **Tropical** | **105/41/35** | **−0.01413 [−0.03612,+0.00893]** | −0.01904 [−0.04925,+0.01067] | +0.00491 [−0.02968,+0.04094] |
| **Northern midlatitude** | **729/50/63** | **−0.02593 [−0.04742,−0.00546]** | **+0.10394 [+0.07641,+0.13060]** | **−0.12987 [−0.15780,−0.10077]** |

The tropical family-qualified n30 result is not comparable one-to-one with the earlier all-species genus Shapley result on all 106 islands: mapping 44 species (including a source-rare profile) out of the family-conditioned estimand and requiring family-qualified source floras changes its support. In northern midlatitude the one-sided-family contribution remains strongly negative, but source lists have greater family richness than the island flora, so this is source-repertoire accounting, not historical family disappearance.

## Northern midlatitude: limited conditional residual

On the **743** broadly recorded northern-midlatitude islands, same-five-family n5 k-adjusted gap is **−0.01708 [−0.03495,−0.00149]** block and **−0.01482 [−0.02968,−0.00184]** source-equal. Same-eight-family comparison is **−0.01852 [−0.03814,−0.00162]** block and −0.01672 [−0.03148,−0.00368] source-equal.

However, requiring different **genera within those same families** makes the k-adjusted gap uncertain: n5 **−0.00983 [−0.02927,+0.00702]**; n8 **−0.00733 [−0.02114,+0.00690]** block; source-equal n5/n8 intervals also contain zero. Post-hoc SAME-ISLAND paired same-five-family-ID minus merely same-five-family-count gap is **−0.01280 [−0.02397,−0.00130]** under geographic block equality but **−0.00584 [−0.01587,+0.00395]** under first-GIFT-source equality, so the effect is not independently source-level stable. This is not a globally uniform within-family morphological shift.

Original source-native-only tropical has **17 eligible island floras** for n5/n8 family comparisons; all corresponding k-adjusted rare−common intervals span zero. Recorded native-status incompleteness and independent ancestral source uncertainty remain major constraints.

## Scientific interpretation and stop rules

The all-recorded **rare-tail category support difference** between island and potential mainland candidate floras is real under the earlier fixed profile-count and genus-richness reference models, but the new family-qualified and identity-controlled original source analyses show it is **not a confirmed family/phylogeny-independent functional contraction**. A signal requiring changing genus names, family names or source lists can reflect species composition, source region flora inventory, study/provider taxon coverage, taxonomic recognition and geographic history—not necessarily island evolutionary phenotype changes. Neither family name nor genus is a branch-length phylogeny. Actual pollen deposition, fecundity, visitor histories and pollinator community shifts were NOT measured.

**Do not merge Draft PRs as causal island syndrome evidence.** The next defensible experiment needs phylogenetically/authority-validated flowering-plant taxonomic coverage with genuinely comparable island/mainland strictly-native floras and independently sampled groups, and repeated actual flower visitors or reproductive endpoints for a pollination mechanism. The present paper can report original source-controlled categorical floral **floristic composition** results, with source and taxonomy caveats.

## Reproducibility

- Scientific code: [family identity and rarity estimator](../../src/fd/family_conditioned_rare_volume.py), [COL XR strict family source map](../../src/fd/colxr_family_opportunity.py), and [source-independent synthetic safeguards](../../tests/test_family_conditioned_rare_volume.py).
- Original protected run: [island original runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_family_conditioned_rare_volume.py), [protected workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-family-conditioned-rare-volume.yml), completed [#43406742054](https://github.com/zuizui0223/island/actions/runs/43406742054), anonymous [aggregate artifact #14026919625](https://github.com/zuizui0223/island/actions/runs/43406742054/artifacts/14026919625).
- The original aggregate JSON SHA256 is **02d75a8ad3be98b7f586f262971380712869ef2142a8f1ae28b5b8d2553a4f22**. Rights-restricted accepted species/island names, family/genus names and GIFT region entities remain within the protected original-source job. No user-visible raw dataset exported.
