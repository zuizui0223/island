# Directional deficits in island floral functional category space

**2026-11-15 pre-outcome protocol (exploratory, not preregistered).**

## Population and estimands
The primary population stays all 4,505 source-original island botanical floras. Each single trait uses its own high/medium species-direct labels, and each joint trait uses only real same-species documentation for exactly those traits. Do not require the historic three-trait/30-species 928-island intersection. Original species-level traits are constant across islands; this concerns floristic composition, not measured within-species adaptation.

We distinguish: (1) direction of the share of directly recorded flora in a predeclared morphological sector; (2) directional deficit of that sector's recorded categories conditional on same available regional state frequencies and trait documentation; and (3) asymmetry of near versus far island profile-set nestedness. Cross-sectional nestedness does not establish historical local extinction.

## Biological sector definitions, fixed before outcome inspection

- Open architecture: floral_form exactly open_radial.
- Bilateral architecture: floral_symmetry exactly zygomorphic, NOT the broader term symmetrical and NOT inferred bee pollination. Radial symmetry is actinomorphic.
- Tube-or-bell: floral_form in the fixed set tubular, bell_campanulate. Bell-only is descriptive sensitivity.
- Highest depth bin: derive only from a separately verified original ordinal tube_depth_class dictionary. Never equate tubular form and measured tube depth; if ordinal classes are not verified, abstain.
- A putative bilateral×tubular sector must use the same accepted species with BOTH direct form and symmetry. Bilateral×deep×bell may not be inferred by multiplying marginal category frequencies on unrelated plants.
- Autonomous selfing capacity and self-incompatibility are different reproductive characters; selfing must not be added to a floral morphology hypervolume, but is a separate directional response if source labels are interpretable.
- Do not label any morphological sector an effective bee, moth or bird pollination guild.

## Geographic and support plan, fixed in advance
Use independently joined, SHA-pinned original GSHHG analysis_regime and distance_to_continent_km (NOT distance between islands), area_km2, climate PCs, island latitude, geology and 10-degree physical blocks. Original unresolved geography remains in the 4,505 population denominator but is not assigned false isolation distances. Primary region groups are northern_midlatitude, tropical, southern_extratropical; absolute latitude >=60 degrees is an OVERLAPPING nested high-latitude sensitivity, never an independent fourth replication.

Directional user hypotheses before inspecting isolation gradients: northern-midlatitude bilateral symmetry, tubular/bell and high depth share decrease with distance, open radial does not decrease; tropical bilateral symmetry decreases, tubular/bell is retained; high latitude has weak slope (testing equivalence requires a fixed practical-effect margin, not simply p>0.05). These are predictions, NOT empirical findings.

Require >=5 documented Direct species per focal island for each endpoint, and >=30 qualified islands from >=10 physically distinct original 10-degree blocks within region/response to claim an estimable slope. Other source observation gates n>=10 and >=10% Direct documentation, actual source-native status and original typed geology are sensitivity cohorts only. Use exactly one island as a unit, with geographic-block-equal bootstrap uncertainty and leave-one-block-out direction checks; never treat millions of reused dyads or bootstrap draws as independent evolution. Adjust a predeclared main five categorical slopes (north 3, tropical 2) with Holm multiplicity; if upper-depth qualifies, include it as a sixth test rather than hunt for a favourable category. Each raw effect is explicitly observational.

Use log1p(distance_to_continent_km) predictor. Correct possible confounding descriptively for log(area), direct-recording fraction, direct carrier counts, and optionally climate PCs. If covariate gates remove many islands, report attrition and sensitivity rather than claim causal adjustment. Distance must be finite, nonnegative, and must not be substituted with nearest GIFT candidate mainland region distance.

## Mathematical caution: directional nestedness cannot use equal k on both islands

For any two finite category sets, **if far and near floras are forced to have the same number k of distinct category profiles**, the statement “all profiles of the more isolated flora are a subset of the nearer flora” is equivalent to **the two sets being IDENTICAL**. It cannot be a directional shrinkage/subsetting phenomenon. Therefore the fixed-k null for where absences are concentrated is a different analysis from directional island-pair subset loss.

For a future near/far island comparison, retain the actual possibly different `k_far` and `k_near`, condition a reference on BOTH richness counts, regional category frequencies and genus support, and report the excess probability of true `profiles_far ⊆ profiles_near` only when `k_far < k_near`. This is an analytical frequency/source-composition control, not a historical island colonization-loss direction. A scalar Baselga category-set nestedness-resultant beta component is symmetric under swapping island names; its sign by itself cannot be identified as “more isolated island lost the specialist profile.” Geographical proximity and original donor/source-history matching must be checked, and original dyad reuse clustered by islands/archipelagos rather than treating millions of pairs as independent replications.

## Missing category space is not total volume
For region-specific possible reference profile categories, compute actual coverage of predeclared sectors and compare to a fixed-k reference with k the island's actually occupied distinct profiles. A uniform fixed-k reference preserves k but NOT original category frequency or genus; it is a negative-control visualization only. Add a region/category-prevalence and carrier-count null, and separately genus-aware availability when identifiable. Available source trait category support comes from observed taxa, not a directly observed ancestral or independently field-measured mainland assemblage. Absent source documentation is not biological absence; no made-up cross-products, Gaussian 17D KDE niche volume or hybrid pollination syndrome.

Stop if trait states cannot be mapped unambiguously, original native status/provenance cannot be verified, or source/space support fails; publish the coverage/identification limitation and retain the global 4,505-island primary descriptive results. Keep original taxon/island/source IDs within protected island CI and PRs Draft.
