# Data sources, analysis provenance and remaining publication gates

## Verified external source citations

The following references were rechecked on 24 November 2026 against publisher or official database pages. They identify the research sources, while the pinned internal source artifacts identify the exact analyzed copies. A citation does not authorize redistribution.

| Source | Published reference | Exact analytical copy and role |
|---|---|---|
| Plant floras | [Weigelt et al. 2020](https://doi.org/10.1111/jbi.13623); [Denelle et al. 2023](https://doi.org/10.1111/2041-210X.14213) | Project-specific GIFT reconstruction: 4,505 botanical island floras. The public GIFT region count is not the analytical count. Original release, retrieval date, original primary checklist references and permissions need reconciliation. |
| Physical geography | [Wessel & Smith 1996](https://doi.org/10.1029/96JB00104); [GSHHG 2.3.7 official release](https://doi.org/10.5281/zenodo.7007502) | 8,264 reference entities, reconstructed hash crosswalk. Frozen source archive SHA256 8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41. Official release is a reference; do not assert this exact file version until original archived manifest is reconciled. |
| Bees | [Marshall et al. 2026, Nature](https://doi.org/10.1038/s41586-026-11325-6); [Data](https://doi.org/10.5281/zenodo.23222636) | Published 18 Nov 2026. Published native richness 4,128 species/305 islands. Data.zip MD5 372c22e2c24e203600a29d7f874a907d; Data/model_table.csv SHA256 a7ddae4cdb41453313d20811371cf64bd38b42b788fc89ae96c284d8d667ddd7. Use FinalRichNat, original polygons, AreaModel/DistModel and other source distance variables. Raw name union is not published accepted species richness. |
| Bird checklist | [Matthews et al. 2023, Ecology Letters](https://doi.org/10.1111/ele.14203); [DARs pinned source](https://github.com/txm676/DARs/tree/8b381ff26e3d6f4da17730dda0c5ab7dbad12eed) | Contemporary true-island matrices in Data/Island_datasets/True_island_datasets. Eleven selected matrices (200 island columns; 132 tropical candidates and 68 Ryukyu sensitivity columns). All footer richness checks reproduced. The 38 tropical bird–bee matches are provisional name/area candidates. Preserve source-level attribution and restrictions for each constituent checklist. |
| Bird diet | [Wilman et al. 2014](https://doi.org/10.1890/13-1917.1) | Figshare file 5631081 SHA256 97216eb1797da077169ebb1ebea275db293b09fc62f8bb8911f9beb98c50d321. Diet.Nect fraction; ≥50, partial, 0. Not measured local pollination. |
| Bird morphology | [Tobias et al. 2022](https://doi.org/10.1111/ele.13898) | Secondary AVONET BirdTree CSV mirror, [commit 757d1659](https://github.com/stephchia/bird-morphology-tda/tree/757d16598fa00ce0e5696a22a54807ad2d299e01), data/raw/AVONET_birdtree.csv. Cite original source; mirror is the analyzed copy. |
| Butterflies | [Yau et al. 2025](https://doi.org/10.1038/s41597-025-05333-w); [Figshare dataset](https://doi.org/10.6084/m9.figshare.25037645) | Revised occurrence file 66997676 (project archive calls it revised v1.1); SHA256 508125561dea6723ae90faf678353f5f27d6182010a89e5f24f03f04b0364c79, MD5 b99c4f94e3781da5cae25035703c8b45, 747,825 rows. Earlier paper has 730,190 rows. File ID/hash rather than an unverified DOI version suffix identifies the actual input. Record data are positive-only. |
| Butterfly traits | [Shirey et al. 2022](https://doi.org/10.1038/s41597-022-01473-5); [LepTraits pinned source](https://github.com/RiesLabGU/LepTraits/tree/f73317d28b38f80923ad8e02cc06d35f5c16ed1f) | consensus/consensus.csv SHA256 6ec35b8a31e96c971aeaa229a48aae9f107c40c33695f0d470aa4382ca6d635b; 12,448 input rows. Species match and duplicate/ambiguity audit precede trait analysis; documented host families are not actual island host use. |

## Protected original inputs (do not release rows)

Restored from [island run 34232450884](https://github.com/zuizui0223/island/actions/runs/34232450884) artifact chapter1-progressive-analysis-34232450884. Relative files and checksums:

| Archived input | SHA256 |
|---|---|
| fixed/canonical/input/chapter1_status_flora.csv.gz | 6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a |
| fixed/source/current/evidence/direct_only_trait_ledger.csv.gz | 9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745 |
| fixed/source/result/island_source_assignments.csv.gz | 8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44 |
| fixed/source/source/gift_native_mainland_flora.csv.gz | 3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef (logged on original run, not originally pre-pinned) |
| corrected_geography_covariates.csv, island commit d4a5dafab5f0737ea6c1c1c0263fc20decf6035c | 08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7 |
| public animal/geometry bridge JSON | 467f111ca14be773bc7dbe6cd76880c443d5f8d6a81c0dad55ca6b8799302776 |
| public bird covariates JSON | 872956d06661bc0fcc8e555beed07c910077423c59bb0d6da58d66570144892c |
| public bird incidence JSON | 6ae6aecf0a2850f91147cea78f0f3e5b6c348f4e1104fa78bf12a43372d7feab |

The original island PR #346 source scripts and anonymous output reports preserve frozen data access without exposing individual botanical rows. Do not use historical public Direct name-unit aggregates as a substitute for original same-species form/depth combinations.

## Main figure and analysis replay

| Result | Verified report / execution | Exact analysis code |
|---|---|---|
| Global recurrence and geographical reassembly | [Main evidence ledger](../../results/20261118_positive_evidence/README.md); [original tropical covariance](https://github.com/zuizui0223/island/actions/runs/44235884651) | Original archived source analyses; global methods and resampling definitions still require source-code reconciliation. |
| Source join and geographical contrasts | [45358282275](https://github.com/zuizui0223/fd/actions/runs/45358282275); [45356938723](https://github.com/zuizui0223/fd/actions/runs/45356938723) | scripts/tropical_pollinator_public_flora_bridge.py; scripts/tropical_three_guild_candidate_audit.py |
| Discovery diet coefficients and family/genus null | [45364977926](https://github.com/zuizui0223/fd/actions/runs/45364977926); [45372672995](https://github.com/zuizui0223/fd/actions/runs/45372672995) | scripts/tropical_specialist_partial_flora.py; scripts/tropical_public_flora_guild_association.py; scripts/tropical_bird_diet_label_null.py |
| Protected original source-relative and form/depth | [45365112213](https://github.com/zuizui0223/island/actions/runs/45365112213); [45372589127](https://github.com/zuizui0223/island/actions/runs/45372589127) | island scripts/fd_original_source_relative_pollinator.py; scripts/fd_original_tubular_depth_diet_split.py |
| Original dual-null | [45376455820](https://github.com/zuizui0223/island/actions/runs/45376455820) | island scripts/fd_original_dual_plant_bird_null.py |
| Exact taxonomy, origin and bounds | [45384576074](https://github.com/zuizui0223/island/actions/runs/45384576074), [45385821139](https://github.com/zuizui0223/island/actions/runs/45385821139), [45386194167](https://github.com/zuizui0223/island/actions/runs/45386194167) | island scripts/fd_original_tubular_taxonomy_decomposition.py; scripts/fd_original_flower_origin_decomposition.py; scripts/fd_original_native_origin_bounds.py |
| Figure data + vector rendering | [Figure captions](FIGURE_CAPTIONS.md), [successful replay 45395209126](https://github.com/zuizui0223/fd/actions/runs/45395209126) | fd scripts/render_integrated_manuscript_figures.mjs; figure_data.json |

## Submission blockers

- GIFT requires the exact GIFT version, original database reference, and all relevant primary source references (official policy: https://gift.uni-goettingen.de/about). Check which restricted checklist contributors require permission; do not presume all raw records may be deposited.
- Confirm actual original GSHHG file version, archived regional regime definitions and original global geographic-block bootstrap/weighting. Exact values are recorded, but methods must be traced to the original implementations.
- Complete bibliographic credit for constituent bird checklists and source trait evidence; verify rights for the secondary AVONET mirror and the revised butterfly occurrence archive (including source-specific GBIF terms). Preserve file hashes and original published references.
- Supplementary figures S1–S6 still need production and review; main SVG replay passed but visual journal layout and accessibility remain to be finalized.
- All potential-partner models remain exploratory, with provisional bird geography, sparse independent groups, nonrandom floral documentation and unresolved plant provenance. No evidence of a causal guild replacement or robust native floral–diet association.
- This manuscript pass changes text and provenance only. No new statistical test, new source reconstruction or independent replication was performed.

## Original botanical pair-method reconciliation (25 November 2026 JST)
[Run 43857181963](https://github.com/zuizui0223/island/actions/runs/43857181963) SUCCESS (job 133264343391), anonymous artifact 14221701164, JSON SHA256 478e89e0dbdde12c967c1611c70b4b6fc18f500eb1fc5b8832e53da0f7cbaf54. [Archived verified report](https://github.com/zuizui0223/fd/blob/analysis/discrete-functional-hypervolume-20261113/results/20261116_true_joint_pollination_syndrome_source/README.md) specifies both-source joint-carrier support (island ≥5, both potential sources ≥10), documentation-balanced threshold (≤10 percentage-point difference), and physical-block versus first-source-group aggregation. [Genus covariance report](https://github.com/zuizui0223/fd/blob/analysis/discrete-functional-hypervolume-20261113/results/20261118_original_genus_flower_pair_covariance/README.md) specifies 999 geographic block resamples. These global-method details are now included in the manuscript. Complete source acquisition and primary references remain open. Do not confuse the 96 documented-native botanical symmetry/depth floras with the 20 animal-matched native-supported form/depth floras.
