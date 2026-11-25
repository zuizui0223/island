# 14 November 2026 — Identical island plant taxa under disjoint near/mid dated-tree mainland source windows

**ORIGINAL protected fd/island macro results, all 3 original dated phylogenetic source workflows completed successfully.** Source [island Actions #43456704114](https://github.com/zuizui0223/island/actions/runs/43456704114) original study input matched across TPL/LCVP/WP, each source original tree and original CoL XR crosswalk SHA verified; anonymous original [TPL artifact #14055115216](https://github.com/zuizui0223/island/actions/runs/43456704114/artifacts/14055115216), [LCVP #14054311779](https://github.com/zuizui0223/island/actions/runs/43456704114/artifacts/14054311779), [World Plants #14055125284](https://github.com/zuizui0223/island/actions/runs/43456704114/artifacts/14055125284). The original source accepted plant taxa, GIFT native mainland source IDs, family names, exact original geographic coordinates and tree-tip identities remain PRIVATE.

## Research question and original source observation constraints

An earlier real-source analysis using 15 tree-tip plants on each side showed tropical island-minus-GIFT source rare-category functional coverage gaps across all three trees, but only about **9%** of those nearest published-tree source pairs were same genus and about **55%** same family. A more stringent same-family cutoff trial yielded a strong negative original tropical gap at <=100 dated-tree path units but no directional gap at <=50. Matching identical physical island IDs for these two cutoff trials gave differences of +0.0589 TPL, +0.0604 LCVP and +0.0617 WP in k-controlled rare-minus-common reference coverage, but those comparisons still changed WHICH local and source plant species were selected.

This new result keeps the **EXACT SAME FIVE focal island accepted plant species** fixed within each draw, all of which are true tips in the published original backbone, while independently choosing two disjoint source windows: a different species of the SAME original CoL XR candidate family with TRUE tree-path distance <=50, and a different donor species of that SAME family with tree-path distance (50,100]. Both independent 5-donor bands use one-to-one actual source taxon matching. Only origin and family/source name and original tree branch paths determine eligibility; flower phenotype states do NOT determine source matching. All focal accepted plant species are excluded from the five-nearest GIFT candidate mainland source floras. No originally absent taxon is interpreted as a true extinction. This posthoc original source support test is not a controlled evolutionary experiment, and neither tree family taxonomy nor GIFT geographic proximity establishes actual common ancestry of an island colonist.

Responses retain original 2,234 source Direct plant species, 100 actually recorded colour/form/symmetry triple categorical reference profiles. Of them 29 rare profiles (1–2 global original documented plant carriers each), 25 common (>=26 carriers each). Outcomes are one-trait-mismatch reference neighborhood coverage and rare-minus-common gap corrected by a fixed original exact occupied-profile-count k reference null. Original source 20 independent repeated five-source taxon draws where feasible; 999 original physical 10-degree spatial-block and nearest-GIFT-source-signature-equal conditional bootstraps. Source-support gate >=20 islands, >=10 geographic blocks and >=10 candidate source signatures.

## Tropical, all-recorded original floras: source distances 0–50 versus 50–100

| True dated tree | Eligible distinct tropical islands / 10° blocks / GIFT donor groups | Difference in rare-profile coverage: medium-source minus near-source [95% BLOCK] | Same exact focal species: k-profile-number-adjusted difference [95% BLOCK] |
| --- | ---: | --- | --- |
| TPL | **66 / 29 / 23** | **+0.01888 [+0.00638,+0.03095]** | +0.02305 [−0.00041,+0.04486] |
| LCVP | **61 / 27 / 21** | **+0.02120 [+0.00660,+0.03586]** | +0.01918 [−0.00901,+0.04363] |
| WP | **62 / 28 / 23** | **+0.01965 [+0.00522,+0.03398]** | +0.01474 [−0.01093,+0.03876] |

**Sign convention:** each difference = (island minus NEAR source) minus (island minus MEDIUM source) = MEDIUM source coverage minus NEAR source coverage for the EXACT SAME focal island plant taxa. A positive raw rare-profile difference indicates the farther within-family source donor group covers more rare categorical neighborhoods than the near donor group. All three rare-only contrasts have whole-block and separately equal-first-GIFT-source conditional 95% intervals wholly positive. In the TPL group, source-equal +0.01729 [+0.00105,+0.03564]; LCVP +0.02280 [+0.00080,+0.04694]; WP +0.02369 [+0.00447,+0.04493].

**Critical falsification:** after conditioning on the observed number k of different exact floral triple profiles in each donor population, **ALL THREE** whole-block and first-source-signature intervals for the same-exact-focal-taxon near-minus-middle category gap include zero. This indicates the raw rare-reference contrast is not independently established beyond differences in documented floral category occupancy number. Consequently **NO causal or phylogenetically independent island flower-space contraction** follows.

Mean original dated path lengths among the paired source taxa were roughly 34 units near and 79–81 units medium in tropical. Each run used ORIGINAL published tree positive lengths, no genus/family-inferred tips and no v.PhyloMaker hypothetical graft. The three tree backbones are related phylogeny hypotheses, not three independent molecular field experiments or separate island colonizations.

## Northern midlatitude: contrast not robust to weighting or tree

The same-five-focal-taxa comparison covers **698 TPL, 686 LCVP, 683 WP** original source-qualified midlatitude islands. Raw rare-reference near-minus-medium differences are +0.00096 TPL, +0.00164 LCVP and −0.01233 WP; the last has a block interval excluding zero, unlike the first two. The exact-k rare-minus-common gap near-minus-middle difference has block intervals crossing zero for TPL and LCVP and negative interval for WP; equal-GIFT-source intervals are negative in all three. These different normalizations/phylogenies do NOT establish a universal globally transferable near-tree-tip effect on floral phenotype space.

## What the result supports and what remains unresolved

**Supported empirical category-space description:** within selected tropical island floras, original published true-tree-tip candidate mainland plants from the same original family but two distinct phylogenetic distance bands differ in how many globally rare morphological category neighborhoods they document, even after fixing the five focal island species per draw. This is a source-composition dependence of documented original floral phenotype types, not observed loss of their real effective pollinators.

**Not supported independently:** distance-based change in true morphological FUNCTIONAL diversity after exact observed-profile richness k adjustment, within-species floral adaptation, lineage-specific historical island extinction/colonization, causal selection, pollinator-guild decline or replacement, pollen limitation, seed set, UV/nectar function, or true verified geographically ancestral mainland sources. The original CoL XR family labels are strict exact **unreviewed classification candidates**, not a gene-sequence phylogenetic tip validation; the published dated tree labels may not all represent individually sequenced taxon placements.

## Original source code

- [Dated published tree LCA source algorithm](../../src/fd/true_backbone_phylo_distance.py).
- [True source-tip island vs donor morphology source estimator](../../src/fd/true_backbone_phylo_source_match.py).
- [Bounded same-family and SAME-focal-plant independent donor-distance-band tests](../../src/fd/true_tree_close_family_source.py).
- [Synthetic identical-focal-plant source separation and cutoff safeguards](../../tests/test_true_tree_close_family_source.py).
- [Protected original three-dated-tree workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-lcvp-phylo-source-match.yml).

All source result evidence is anonymous region-level aggregate. Keep [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and [island PR #326](https://github.com/zuizui0223/island/pull/326) **Draft**, without converting posthoc conditional category-occupancy differences into causal pollinator-shift biology.
