# 2026-11-16 — Fixed-island, two-source EXACT four-term floral genus accounting across 4,505 GIFT floras

**Original protected-data numerical result SUCCESS:** [island Actions #43852272780](https://github.com/zuizui0223/island/actions/runs/43852272780) verified synthetic finite signed identity, source-locked accepted plant Direct phenotypes and original GIFT botanical floras, independently rebuilt GSHHG 10° geographic blocks. Original anonymized [artifact #14218297756](https://github.com/zuizui0223/island/actions/runs/43852272780/artifacts/14218297756), extracted JSON SHA256 `f6d049c7f96db17dd1a7efa43aa16ab1f90b162536678581844878730c2ea70d`. No private island, species, genus name, GIFT-source, geometry or GBIF identifiers exported. [Runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_source_genus_exact_modes.py), [exact signed aggregate library](../../src/fd/source_genus_exact_partition_modes.py), [independent pure genus math](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/fd_direct_source_genus_accounting.py).

**Independent implementation equivalence:** Original optimized source-genus-count replay [island Actions #43853678261](https://github.com/zuizui0223/island/actions/runs/43853678261) SUCCESS with 5 fd synthetic tests plus 12 island original exact-genus/equivalence tests. Original anonymous [artifact #14221290295](https://github.com/zuizui0223/island/actions/runs/43853678261/artifacts/14221290295), SHA256 `2d2f44fcfce7813b709c33a2355a72d73a9ce8d0c479613cb29dfb411d66fb09`. Across all 100 anonymous report panels and nested CI summaries, the optimized cached-donor calculation reproduces the previously released exact decomposition: maximum numeric absolute difference **2.22e−16** (pure floating-point rounding), with 959 scalar slots differing only at rounding level. This is a deterministic code-equivalence audit on the SAME source archive, not an independent biological replication.

## Why this is a new, stronger check

The prior [actual **different accepted mainland species** congener-only reference](../20261116_source_congener_reference/README.md) is useful but changes the focal plant cohort: only 29–50% of tropical original Direct island plants often have valid source-only congeners. A small congener-source deviation therefore does NOT numerically decompose the original whole-flora difference. Here we **hold each original island, plant status, Direct morphological category, and its accepted focal species names fixed**, and partition the *same original observed species-fraction difference* in FOUR algebraically exact terms, separately for actual `geo_k5` and `geo50_climate10` source-candidate pools. All 4,505 original accepted botanical island floras remain the study population; 4 original nested status scopes ×5 measured morphological sectors ×4 regional regimes+1 unresolved = 100 reported panels / 90,100 original private cells. Not the old 928 triple-complete / n30 subset.

For a true recorded island Direct category and a distinct-name mainland candidate source, define genus weights \(w_{Ig},w_{Mg}\), documented within-genus frequencies \(p_{Ig},p_{Mg}\), and genus sets \(G_I,G_M\). The exact identity is:

\[
D=\sum_{g\in G_I\cap G_M}
\left[\underbrace{(w_{Ig}-w_{Mg})(p_{Ig}+p_{Mg})/2}_{\text{shared genus weight}}+
\underbrace{(w_{Ig}+w_{Mg})(p_{Ig}-p_{Mg})/2}_{\text{shared genus distinct species states}}\right]+
\sum_{g\in G_I\setminus G_M}w_{Ig}p_{Ig}
-\sum_{g\in G_M\setminus G_I}w_{Mg}p_{Mg}.
\]

Every botanical source species name recorded in the entire original focal island checklist is removed from BOTH mainland candidate-only pools. A zero shared-genus component is not evidence for biological evolutionary equality; source/plant detection and taxon histories differ. The four SIGNED terms can strongly cancel, so their shares are **not causal variance fractions**. Mainland 'candidate' means GIFT native mainland checklist, NOT actual seed arrival, ancestor, locally sampled phenotype or independent selection trial.

## Full original tropical all-recorded morphology on SAME source/plant scope

10-degree physical-block-equal ORIGINAL category percentage points. The four terms EXACTLY sum to island−candidate source recorded species-category prevalence on each fixed original island before aggregation. The two source modes independently select different botanical candidate accepted species.

| Morphology | Original both-source island floras / blocks | Mode | Shared genus weight | Within-shared-genus distinct species states | Island-only genus | Candidate-mainland-only genus | Exact observed difference |
|---|---:|---|---:|---:|---:|---:|---:|
| **Zygomorphic** | 945 /99 | geo5 | +10.38 | **+0.15** | +14.44 | −33.35 | **−8.37** |
| Zygomorphic | 945 /99 | climate10 | +13.14 | **+0.11** | +12.18 | −34.00 | **−8.58** |
| **Tube/bell** | 773 /94 | geo5 | +10.86 | **−0.03** | +34.47 | −42.10 | **+3.20** |
| Tube/bell | 773 /94 | climate10 | +14.09 | **+0.01** | +31.39 | −41.77 | **+3.73** |
| **Deepest tube** | 728 /90 | geo5 | +9.44 | **+0.91** | +17.49 | −17.93 | **+9.91** |
| Deepest tube | 728 /90 | climate10 | +12.10 | **+0.95** | +15.27 | −18.32 | **+10.00** |

A strong and coherent **observational taxonomy-linked floristic reassembly** is evident: the large all-plant bilateral deficit and tube/deep-tube positive shares arise alongside large but cancelling presence/weight contributions from botanical genera. The actual genus-shared distinct-accepted-species state component is far smaller than the full differences. This does NOT statistically identify a lineage filter or quantify its causal contribution. In particular, the very negative candidate-only genus contribution reflects the GIFT candidate botanical checklist and recording, NOT observed propagule arrivals or demonstrated establishment failures.

## Robustness and fragility of the within-shared-genus component

Conditional group-resampling 95% intervals (all comparisons observational; no family-wide multiplicity correction):

| Tropical all-recorded, within-shared-genus distinct species component | geo5 pp [physical-block CI] | climate10 pp [physical-block CI] | GIFT source/plant DIRECT coverage <=10 pp difference cohort |
|---|---|---|---|
| Zygomorphic, n945 | +0.15 [+0.01,+0.33] | +0.11 [−0.01,+0.24] | **n106/28 blocks:** +0.04 [−0.11,+0.17] / +0.01 [−0.16,+0.21] |
| Tube/bell, n773 | −0.03 [−0.57,+0.46] | +0.01 [−0.54,+0.55] | **n304/53:** +0.25 [−0.48,+1.13] / +0.33 [−0.39,+1.14] |
| Deepest tube, n728 | **+0.91 [+0.49,+1.32]** | **+0.95 [+0.45,+1.46]** | **n338/67:** +0.35 [−0.22,+0.88] / +0.35 [−0.18,+0.85] |

Thus a nominal positive within-genus deepest-tube component is evident on all-recorded floras, but its apparent directional support **does not survive the original <=10pp Direct documentation gap screen**. Do not call it provenance-invariant natural selection or universal within-lineage evolutionary change. The congeners-only *selected cohort* statistic of +3.24 geo5 /+2.89 climate10 pp estimates a different weighting; it is NOT in conflict with the exact full-cohort +0.91/+0.95 pp signed term.

### Narrow positive candidate in STRICT native nonendemic plants (posthoc and selected)

On 107 tropical original strict known-native nonendemic island floras with real candidate sources, the whole-status cohort within-genus term is **geo5 +1.20 [+0.12,+2.29] vs climate10 +0.78 [−0.29,+1.75]**, physical blocks; first GIFT-source grouping intervals cross zero for both. In an **observed Direct-documentation-gap <=10 percentage point** subset, only **36** original island floras remain across **24 physical blocks and 18 first GIFT source groups**. There the same exact signed within-genus term is **geo5 +1.79 [+0.30,+3.36] and climate10 +1.45 [+0.33,+2.49]** under physical blocks; grouped by first GIFT source, geo5 +1.28 [+0.13,+2.52], climate10 +1.15 [+0.32,+2.00].

This is an *interesting but posthoc*, small, nested and geographically/source-shared **distinct-species state** signal. It was found among numerous correlated region×scope×morphology×source×coverage checks; no family-wise correction or independent provider validation has been applied, and congeneric accepted species are not known ancestral/sister lineages. Do not turn its conditional intervals into confirmed island floral adaptation or a measured pollinator replacement mechanism. It is a candidate for an independent future validation on actual local morphology/source relatives.

## Contrasts outside tropics reject a universal functional contraction narrative

- Northern-midlatitude original all-recorded **zygomorphic**: n1750, geo5 raw −4.99 vs exact within-shared-genus **−0.33 [−0.48,−0.19]** pp; climate10 raw −7.12 vs within **−0.39 [−0.58,−0.22]**. A much larger source membership/presence component accompanies the signed trait difference.
- Northern high-latitude original all-recorded **zygomorphic**: n258, geo5 raw **−15.17** vs within **−0.11 [−0.35,+0.17]**; climate10 raw −15.68 vs within −0.20 [−0.46,+0.05]. High latitude cannot be called 'island source-invariant,' but neither can documented negative bilateral category shares be equated with independently repeated within-genus loss.
- Southern extratropical original all-recorded **open-radial**: n157, geo5 raw **+8.74** vs within **+0.30 [−0.33,+0.94]**; climate10 raw **+10.76** vs within +0.16 [−0.50,+0.83]. Original Direct annotation-gap <=10pp leaves only 31 islands; southern mainland ecosystems are ecologically heterogeneous.

Original observational phenotypes are one species-wide documented category, not location-specific island flower images, within-plant state changes or fitness. Geological origin labels have independent but sparse original coverage; cannot conflate continental isolation with lost bee/bird channel.

## Decision / strongest defensible manuscript claim

**Main supported observation**: The global island flower-assemblage differences relative to explicitly defined GIFT native mainland candidate floras are **directionally regional and morphology-specific**, and their exact signed source representation is mainly located in botanical genus membership and relative taxonomic weights, not a universal unidirectional within-shared-genus morphological shift. This is consistent with nonuniform source-conditioned floristic assembly and falsifies a simplistic homogeneous floral 'syndrome loss' narrative.

**Not supported**: a numerical causal 'genus filtering explains X%' conclusion; actual arrival vs establishment separation; ubiquitous within-genus island evolution; effective pollinator-guild transitions; a confirmed species-level deep-tube selection response. The earlier original source-breadth fixed-genus finite null and strong flower Direct documentation imbalance remain independent interpretation limits.

**Next methodological priority:** externally validate the narrow tropical strict-native nonendemic deep-tube signal under a new source/trait-provider, phylogenetically independent lineage groups, complete native geographic provenance and measured tubular depth where possible. Maintain [fd PR #16](https://github.com/zuizui0223/fd/pull/16) DRAFT. Preserve original full 4,505-island universe and present separate marginal/joint Rao Q results, rather than replacing source-category reassembly with historical 17D KDE volume.
