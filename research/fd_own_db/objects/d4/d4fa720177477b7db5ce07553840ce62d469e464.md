# 14 November 2026 — 3 published dated backbones and a stricter within-family tree-distance audit

**Source-backed observational macroecology; no direct causal evidence of flower evolution or pollinator loss.** Original source workflows: [3 independently named original backbones, 3/3 successful #43448315502](https://github.com/zuizui0223/island/actions/runs/43448315502) and stricter within-family cutoff diagnostic [#43454436980](https://github.com/zuizui0223/island/actions/runs/43454436980), both 3/3 successful. Only anonymous regional aggregates are uploaded; original island, taxon, family, original source entity and precise coordinate records remain private.

## Original frozen source and inferred targets

The input remains the original 928 island floras, 2,234 single-state Direct species, 100 observed colour×form×symmetry profiles. Tree source is Jin & Qian V.PhyloMaker2, original source commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9 and the original published TPL (74,529 tips), LCVP (73,420), World Plants (72,570) backbones, each with finite positive branch lengths. Of 2,234 directly categorized plants, respectively 1,187/1,162/1,160 have strict published tree-tip name, genus and source-family agreement. These original published tree tips are not necessarily independently sequenced by this study; genus/family-inserted V.PhyloMaker hypotheses are excluded.

Original five-nearest GIFT candidate mainland native flora sources and exact same focal-island accepted species-name exclusion were fixed. The floral outcome is fraction of ORIGINAL rare-profile (29 types with 1–2 source species carriers) minus common-profile (25 types with ≥26 source species carriers) one-category-mismatch neighborhood coverage, conditional on exact recorded profile count k. This is categorical floristic support, not pollination service or genuine within-plant flower evolution.

## Fifteen published tree-tip plants vs 15 different nearest eligible donor plants

| Original backbone | Tropical islands / spatial blocks / first GIFT source signatures | Island minus optimally nearest GIFT donor k-adjusted rare−common support [95% block] | Same-island phylo-nearest minus random donor difference [95% block] |
| --- | ---: | --- | --- |
| TPL | 103 / 42 / 36 | **−0.04965 [−0.07491, −0.02595]** | **+0.00058 [−0.01130,+0.01350]** |
| LCVP | 103 / 42 / 36 | **−0.04463 [−0.07071, −0.01926]** | **−0.00675 [−0.02145,+0.00780]** |
| WP | 103 / 42 / 36 | **−0.05384 [−0.08010, −0.02913]** | **−0.00313 [−0.01850,+0.01400]** |

All three nearest-source matched gaps are negative; all three same-island **phylogenetic source selection changes** have 95% intervals crossing zero. The three plant trees share most of their underlying megatree provenance and are not three independent island evolutionary radiations. North-midlatitude ~739–741 tested islands have source gap intervals crossing zero; the northern-high and southern first-nearest GIFT sources do not meet the original predeclared support gate.

**Critical phylogenetic closeness qualification:** tropical 15-tip 'optimized nearest phylogenetic donor' still matches only **8.8–9.2%** in the same genus and about **54.5–55.3%** in the same independently matching CoL XR candidate family. Real original dated source pairwise distances average ~116 in matched sets versus ~228 under a random mainland source draw; 25.7–26.1% of matched pairs have original branch-path distance ≤50 and about 53.0–53.4% distance ≤100. No category below is an inferred actual ancestral lineage. Thus a broadly negative source-flora contrast is NOT proof of an independent effect within near relatives.

## Post-hoc stricter source-family plus dated-branch length ceiling

On original frozen tree tips and *within the same CoL candidate family*, obtain a one-to-one original observed source plant pairing, before inspecting phenotypes, whose source distance is below a fixed path threshold. Independently repeat 20 times with 5 original distinct plant pairs. A separate same-original-island 5-tip control samples source plant pairs matched on family but without imposing the tree-distance bound. Both source/tip branches are exact original published tips, not taxonomic grafts. These 50 and 100 thresholds are **post-hoc**, a descriptive stress-test rather than a preregistered selection experiment.

| Tropical same-family constraint | TPL eligible islands; gap [95% block] | LCVP eligible islands; gap [95% block] | WP eligible islands; gap [95% block] |
| --- | --- | --- | --- |
| Tree path ≤100 | 106; **−0.06598 [−0.08910,−0.04372]** | 105; **−0.06285 [−0.08710,−0.03916]** | 105; **−0.06084 [−0.08435,−0.03647]** |
| **Tree path ≤50** | 90; +0.00224 [−0.01818,+0.02380] | 90; +0.00120 [−0.01907,+0.02217] | 88; −0.00188 [−0.02588,+0.02088] |
| Same genus; tree path ≤100 | 40; +0.02013 [−0.01159,+0.05037] | 37; +0.01413 [−0.02466,+0.05000] | 39; +0.02091 [−0.00960,+0.05259] |

Under the source-signature equal bootstrap, the ≤100 same-family gap remains below zero in all three; ≤50 and within-genus intervals overlap zero in all three. At cutoff ≤100, TPL and LCVP matched-minus-random-same-family donor gap has negative conditional 95% intervals but WP crosses zero, showing source-tree sensitivity.

**Most important limitation:** the ≤50 source-qualified subset contains fewer and nonrandom original island floras (88–90 compared with 105–106 at ≤100), and its n5 species sample is more phylogenetically narrow; therefore apparent loss of signal at ≤50 is NOT by itself an independent causal impact of evolutionary separation. To distinguish this source-cohort attrition from tree-distance restriction, a new **paired exact-same-island 50-versus-100** analysis is coded and flagged explicitly as post-hoc, without claiming identical species on either side of the two thresholds. Its result must be read separately and cannot prove phylogenetic adaptation.

## Interpretation ceiling

These original-source results support a documented source-relative **floristic rare floral categorical repertoire** pattern. The pattern can be altered materially by source **species phylogenetic closeness, family, taxonomic provenance and selected island eligibility**; absence of a 95% interval-separated contrast among strict close relatives does not establish biological equivalence. Source tree branches are a published evolutionary hypothesis; one accepted species has a fixed recorded phenotype across all sites. No population-level flower phenotype evolution, true historical island colonization, visitor guild decline or pollen deposition/seed set was measured.

## Code and rights-safe audit

- Original dated tree distance algorithm: [true_backbone_phylo_distance.py](../../src/fd/true_backbone_phylo_distance.py).
- Actual plant-tree/GIFT donor analysis: [true_backbone_phylo_source_match.py](../../src/fd/true_backbone_phylo_source_match.py).
- Same-family and tree-distance-bounded plant-pair analysis: [true_tree_close_family_source.py](../../src/fd/true_tree_close_family_source.py).
- Protected original island executable: [run_fd_lcvp_phylogenetic_source_match.py](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_lcvp_phylogenetic_source_match.py).
- Source original tree runs: [#43448315502](https://github.com/zuizui0223/island/actions/runs/43448315502), [#43454436980](https://github.com/zuizui0223/island/actions/runs/43454436980). All source data rights stay in upstream authorized CI; no original accepted species or exact source region identifiers or copyrighted tree files redistributed.

Keep fd PR #16 and island PR #326 **Draft**; no adaptive island syndrome causal claim.

## Additional post-hoc exact-SAME-ISLAND threshold contrast (original 3 trees, all successful)

The separate [original workflow #43455223738](https://github.com/zuizui0223/island/actions/runs/43455223738) ran the 0–50 and 0–100 within-family n5 donor pair tests and computed a directly paired contrast on the EXACT SAME original physical island IDs and GIFT source signatures, avoiding the different-island denominator. This does NOT fix the individual focal plant species sampled in each threshold comparison; those may still be different.

| Tropical same-island comparison | TPL 90 islands / 38 blocks / 33 sources | LCVP 90 / 38 / 33 | WP 88 / 38 / 31 |
| --- | --- | --- | --- |
| **≤50 minus ≤100 k-conditioned rare−common gap**, block mean [95%] | **+0.05892 [+0.03297,+0.08728]** | **+0.06045 [+0.03590,+0.08654]** | **+0.06166 [+0.03579,+0.08939]** |
| Same island **family-only random source** comparison, block | +0.03918 [+0.01115,+0.06665] | +0.04068 [+0.01820,+0.06235] | +0.04906 [+0.02584,+0.07190] |
| Threshold-bounded minus family-only response, block | +0.01974 [+0.00469,+0.03398] | +0.01976 [+0.00575,+0.03318] | +0.01260 [−0.00078,+0.02607] |
| Threshold-bounded minus family-only response, source-signature | +0.00601 [−0.01059,+0.02273] | +0.01666 [+0.00076,+0.03101] | +0.01123 [−0.00373,+0.02530] |

The same-physical-island contrast does show greater recorded rare-category support in the <50 condition than in <100 on selected original tropical islands. However, even a source-matched FAMILY-only comparator shows a difference and **2/3 source-signature conditional intervals for the threshold-specific incremental effect cross zero**. Selection of WHICH focal species and GIFT donor species can still confound a branch-length effect.

An additional posthoc source comparator requiring the EXACT SAME FIVE original island plant tips under non-overlapping donor path-distance windows (≤50 versus >50 and ≤100, same original accepted CoL candidate family) is executed and treated as a separate support gate; until completion no effect is inferred. Source family labels are unreviewed CoL candidate matches, the published tree tips are not necessarily newly sequenced, and no pollination service was measured.
