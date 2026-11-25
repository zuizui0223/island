# 13 November 2026 — Does the original island rare flower-profile gap reflect genus richness or genus identity?

**Source-back original-data full 928-island macro audit, completed.** Original source CI [island Actions #43387213850](https://github.com/zuizui0223/island/actions/runs/43387213850) **SUCCESS**; anonymous source aggregate [artifact #14019535508](https://github.com/zuizui0223/island/actions/runs/43387213850/artifacts/14019535508). Scientific source [fd PR #16](https://github.com/zuizui0223/fd/pull/16) and protected original runner [island PR #326](https://github.com/zuizui0223/island/pull/326) remain Draft. Original Direct taxa, source GIFT entities, native-status rows, island IDs and exact geographic coordinates were not released. **No actual pollinator decline or fitness change is measured.**

## Primary falsifiable alternative after the genus Shapley accounting

The immediately preceding [original full-928 genus Shapley analysis](GENUS_SHAPLEY_MACRO_928.md) found that the difference in rare-category floral support between island and donor floras is attributed disproportionately to genus labels **observed on only one side**, not the genera observed in both floras. But this might be a *trivial richness artifact* because the original candidate GIFT mainland lists frequently contain many more distinct plant genera than the island.

We therefore measure a new distinct support estimand that fixes both **species count AND the number of observed distinct genera** in every island-versus-donor comparison, while deliberately leaving **genus IDENTITY unrestricted**:
- Fixed source input: original 2,234 Direct high/medium unambiguous species, 100 *actually observed* three-trait colour×form×symmetry profiles. Global rare bins = 29 profiles recorded for 1–2 species, common = 25 profiles recorded in ≥26 species; 46 intermediate.
- Original fixed observational island sample: 928 island floras with >=30 directly characterized accepted species; source is the five-nearest geographically assigned GIFT candidate native mainland region flora, not phylogenetic ancestor.
- **Strict source exclusion:** remove all exact accepted focal island plant species names from every candidate mainland GIFT donor list before sampling.
- For each island and source independently, randomly select **15 distinct botanical genera** from the directly characterized species roster, and **one accepted Direct plant species per chosen genus**. Thus both floral samples have EXACTLY **15 different genera and 15 different species**, but the genus identities usually differ. Repeat 40 source- and island-independent draws. A separate fixed **10 distinct genus** group is a predeclared sensitivity.
- Outcome: reference support fraction for original *rare profile neighborhoods minus common profile neighborhoods*, measured on the globally frozen 100-profile reference using Hamming distance at most **one original trait mismatch**. Report raw gap, and gap after 399 exact-k profile-count conditional reference subsets, both uniform and original species-frequency weighted. Do not call the latter the globally species-frequency-weighted **volume**; it is a *reference null for k profile centers*, not a changed measure of reference volume.
- Conditional block- and first-nearest-source-signature-equal bootstrap (1,999 draws) for every eligible region. Required ≥10 islands, ≥5 blocks and ≥5 first-nearest-source signatures. This is observational taxon/source composition conditional uncertainty, not historical phylogenetic replication.
- The genus identities are **not** matched in this experiment. It is distinct from the already executed **identical-genus 5/8-unit control**, which asks about replacement of plant species **within the same observed genus names**.

## Original source-backed matched-genus-richness results

**n=15 distinct botanical genera and exactly 15 species on each side**. Main response is *source-relative rare-minus-common category-neighborhood coverage, adjusted for exact observed profile richness k* (unit = fraction of source reference profiles); figures are block-equal mean, conditional 95% block interval and a separate source-signature-equal interval.

| Scope and original region | Actual eligible islands / physical blocks / first source signatures | Rare-minus-common gap adjusted for exact k, block | First-source-signature equal |
| --- | --- | --- | --- |
| **All-recorded tropical** | **106 / 42 / 36** | **−0.03669 [−0.05308,−0.01896]** | **−0.03972 [−0.05605,−0.01981]** |
| **All-recorded northern midlatitude** | **742 / 52 / 65** | **−0.02694 [−0.04257,−0.01153]** | **−0.02227 [−0.03617,−0.00902]** |
| All-recorded southern extratropical | 26 / 13 / 11 | −0.00879 [−0.03587,+0.01868] | −0.00903 [−0.03550,+0.01608] |
| Source-native northern midlatitude | 67 / 16 / 23 | +0.01276 [−0.01511,+0.03674] | +0.01608 [−0.00848,+0.03961] |
| Source-native tropical | 17 / 10 / 9 | −0.00865 [−0.04995,+0.02809] | −0.01418 [−0.05275,+0.01949] |

The northern high latitude cohort has 53 eligible islands but **only three first-nearest GIFT source signatures**, so it is suppressed from geographic/source-independent directional inference even though block values exist.

The tropical **UNADJUSTED** rare coverage difference under this equal-15-distinct-genus cohort is −0.03581 [−0.05245,−0.01873], but ordinary common-reference coverage difference +0.00202 [−0.00320,+0.00695] is near ceiling. The k-adjusted rare-minus-common gap **−0.03669** remains directional under both source vs block weighting. A more restrictive uniform-unique-profile k-reference but with global species-frequency-proportional original reference draws has tropical **−0.03193 [−0.04651,−0.01707]**, separate source-equal **−0.03581 [−0.05141,−0.01872]**.

### Independent sensitivity: 10 instead of 15 distinct genera

- **Tropical 106 / 42 blocks / 36 source signatures:** k-adjusted **−0.03160 [−0.04770,−0.01538]** block and −0.02932 [−0.04545,−0.01253] equal-source.
- **Northern midlatitude 742 / 52 / 65:** k-adjusted **−0.01913 [−0.03423,−0.00516]** block and −0.01756 [−0.02899,−0.00727] source-equal.
- Strict-native tropical n10 retains 17 eligible islands and still has a confidence interval containing zero.

Thus a full **106-island** tropical source-relative deficit in the rare end of the original Direct category-space persists even when every island and GIFT candidate source contains the **same number of plant genera and same number of characterized plant species**. Mere **genus richness imbalance is insufficient** to account for this specific observed reference-measure contrast.

## Triangulation: different genus identities versus different species within identical genera

These are different source-selection/population/metric estimands, so they should not be compared as an exact causal mediation or one-to-one attenuation decomposition:

| Source-qualified analysis | Tropical eligible islands | Recorded rare/profile category finding |
| --- | ---: | --- |
| Original 30-species category source contrast, all taxa | 106 | Negative rare-minus-common k-controlled reference gap |
| 30-species **exact Shapley genus-group accounting** | 106 | Rare support: **total −0.02970**, shared-genus contribution **+0.00945**, one-sided-genus contribution **−0.03914**; only the last contribution has conditional block & source intervals fully below zero |
| **15 different botanical genera on island AND source** | **106** | **k-controlled rare-minus-common −0.03669 [−0.05308,−0.01896]**; preserving genus richness does NOT eliminate documented gap |
| **Exactly same 5 botanical genus identities, one different accepted species each** | **81** | **+0.00261 [−0.01579,+0.01982]**, no directional negative support |
| Same **8 genus identities** | 50 | +0.00914 [−0.01475,+0.03532], no directional negative support |
| Source-native original Direct tropical flora | 17 | n15 and n10 equal-genus-richness intervals include zero |

**Interpretation:** it is the **identity/composition of recorded botanical genera** — not just how many genera are in a sample, and not a demonstrated phenotypic contraction among source-comparable member species of the same genus — that is most consistent with the observed *all-recorded tropical* rare category reference-neighborhood shortfall.

This statement remains **association and descriptive accounting**, NOT demonstrated that named genera died out, failed to colonize, evolved smaller flowers or lost their pollinators. The identical-genus 5/8 plant comparisons select genera shared by two checklists, which is a nonrandom subset and has fewer species than the all-genera n15 analysis.

## Key limitations before claiming universal island syndrome

- The original Direct floral traits are **fixed at one accepted species-wide single category each** across the mainland and all islands; within-species evolution and plasticity cannot be inferred by design.
- The original source pool is geographically close but not empirically established as a true phylogenetic ancestral/propagule supply pool. Different genus identities may represent geography/biogeography, detection, genus-level data availability, provenance status, climate or floral/ecological filtering.
- Even though each n15 comparison has **equal number of species and distinct genera**, local taxonomy **family/phylogenetic structure** and source-provider study/trait classification opportunity remain unmatched. Confidence intervals from source and physical block weights are not independent evolutionary replicates.
- Strict source-native Direct coverage is extremely sparse in the tropics (17 eligible islands) and lacks directional evidence. A general island-adaptation story cannot be based on broad recorded provenance alone.
- **Common reference** may be nearly 100% supported under one-mismatch neighborhoods (a mathematical ceiling). The primary indicator is the rare-reference coverage itself, while rare-minus-common gap and k-referenced residual are inferential controls.
- Original source 29 rare profiles are represented by **only 1–2 Direct recorded species**, sometimes potentially an observation/trait-label rarity, not a confirmed ecological specialization or low-pollinator-redundancy function.
- No source-pinned original record on actual flower visiting, pollen transfer, fruit/seed set or historic pollinator decline is joined to these genus-normalized 928 floras. Do not infer causal role of pollinator shifts from this contrast.

## Reproducibility

[Source code](../../src/fd/genus_richness_rare_volume.py), [independent synthetic tests](../../tests/test_genus_richness_rare_volume.py), [original-rights-safe source runner](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_genus_richness_rare_volume.py), [original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-genus-richness-rare-volume.yml) and [successful original run #43387213850](https://github.com/zuizui0223/island/actions/runs/43387213850), [source-anonymous aggregate artifact #14019535508](https://github.com/zuizui0223/island/actions/runs/43387213850/artifacts/14019535508). Source original Direct accepted taxa, GIFT mainland entities, island identifiers, source genus names and coordinates are processed only in protected source job; the public artifact contains regional data-support statistics and conditional group summaries only. No experimental loss results generated.

**Decision:** promote as **source-relative rarity-sensitive functional morphological *assemblage composition***, and prioritize a botanical family/phylogenetic and observation-status source-availability-conditioned null before any claim of evolutionary phenotype filtering.
