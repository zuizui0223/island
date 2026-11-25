# 15 November 2026 — Regional floral genus-reference identification contract

## Why this audit is necessary

The protected original 4,505-island GIFT botanical frame includes 2,324 floras with ≥5 high/medium **species-direct** colour+form records on the **same original accepted plant species**. The prior all-recorded, across-block northern-midlatitude **global** within-genus profile-swap residual was positive on 138,117 species-disjoint island pairs, but almost all those pairs are >2,000 km apart; independent strict-native same-island comparisons greatly attenuated the result. Neither result establishes independent island floral adaptation or a pollinator syndrome.

The geographical audit's pooled genus-reference donor set contains Direct-colour+form species recorded across all studied regions. Thus an observed-versus-global-genus reference contrast can reflect **regional genus-level source opportunity** rather than independent local plant phenotypic selection. The next study conditions the null *reference* on the focal broad botanical region **without changing the focal pairs**.

## One dataset, two reference expectations

For exactly the same focal island pair `(i,j)`, define `p_i(c)` as its fraction of accepted species with an original directly documented flower-colour×form combination `c`, and

`D²_ij = sum_c [p_i(c) - p_j(c)]²`.

Keep the original island-by-accepted-species incidence matrix and the original entire recorded two-trait profile attached to each accepted plant species. Compare:

- **Global genus null:** analytically expected `D²` when those complete colour+form profiles are swapped within their botanical genus among all globally Direct-measured accepted species. Original genus name parsing and fixed incidence masks match `fd_island_functional_beta.py`.
- **Region-local genus null:** analytically expected `D²` when swapping complete colour+form profiles within genus among only trait-measured accepted species represented in eligible island botanical floras of the **same broad region** as the focal pair. A species's observed trait is never changed in the real data; reference shuffling is counterfactual. Regional donor selection includes the focal plants, so this is explicitly an **in-sample descriptive diagnostic, not prospective out-of-island or mainland prediction**.

The measured `D²` and focal pair graph must be numerically identical between the two calculations. The exact paired comparison is `(D²-E_region)-(D²-E_global) = E_global-E_region`. A shift in this reference-only quantity does not identify the fraction of ecological geographical filtering explained.

## Fixed focal support and no dyad pseudoreplication

- Keep the original 4,505 island botanical population, split all-recorded vs strict original source-native plant-name scopes. **Never** require a third phenotype, a KDE space, a 30-complete-species list, or a known real pollinator.
- Each focal island has ≥5 directly known colour+form species in its own scope; other individual traits have independent support and remain outside this endpoint.
- Eligible focal pairs share **zero accepted species in the FULL original botanical species lists**, even when the measured endpoint uses a native-only subset; are from the SAME pre-existing broad region; and are in DIFFERENT independently derived 10° spatial blocks.
- Within n5, separate <500, 500–2,000 and ≥2,000 km distances, plus an all-distance row. n≥10 Direct carriers and ≥10% paired-trait documentation coverage are independently reported sensitivity gates, not matched populations.
- Suppress output for strata with <100 pairs, <30 unique focal islands or <10 spatial blocks. Geographic isolation and trait coverage exclusions are determined without reading floral outcomes.
- Report pair-equal, island-equal and block-equal contrasts plus leave-one-block extreme values. Reused dyads, numbers of numerical seeds, source bootstrap resamples and overlapping All/Native subsets are **never independent biological replicates**.
- Retain anonymous source metadata: regional number of available genus groups, number of genera with actual multiple observed whole colour+form profiles, and fraction of eligible distinct species falling in such variable genera.

## Pre-registered interpretation gates

**An in-sample regional-genus zero does not demonstrate ecological equivalence.** If almost every local genus is recorded with one species or one floral profile, the regional within-genus null becomes deterministic and reproduces the observed category frequencies mechanically. This is a loss of null exchangeability, not evidence against adaptation.

- A stable positive or negative `region-relative residual` across at least ten physical blocks is a **descriptive region-conditioned morphological assemblage difference**. It is not geographic source lineage-independent evolution.
- If the global sign reverses or attenuates when the donor pool is region restricted, report the change and the changing within-genus donor support. Do not convert this into causal percentages or phylogenetic p-values.
- If all <2,000 km region strata fail the prespecified support gates, no **local** region-conditioned ecological effect can be inferred from abundant distant pairs.
- If strict-native strata shrink or select a different island-pair sample, do not interpret All-minus-Native contrasts without an identical-dyad design.
- No local pollen effective visitor network, past guild decline, true plant population changes, phylogenetic branches, or native ancestral distributions are measured by this source; leave pollination and adaptive selection claims unmade.

## Reproduction and data protection

Original trait and botanical membership data are restored inside rights-protected `island` GitHub Actions and SHA256 checked. The independent source GSHHG geography is rebuilt from its audited commit. The original runner is `island/scripts/fd_island_beta_regional_genus.py`, with synthetic tests `tests/test_fd_island_beta_regional_genus.py`, and isolated workflow `.github/workflows/run-fd-island-beta-regional-genus.yml`. Only grouped and anonymous region/scope/support/distance scores leave the source runner; never publish original accepted species, island names, coördinates, mainland source entity IDs, or private study keys.
