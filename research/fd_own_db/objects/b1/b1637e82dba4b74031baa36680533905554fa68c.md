# 2026-11-15 — Original GIFT source-range LISTING opportunity matched direct native floristic-sector contrast

This is a **posthoc sensitivity designed before inspecting new original breadth-matched outcomes**, not a causal establishment propensity model. Global source population still all 4,505 botanical island floras, independently reported 4 northern/tropical/southern/highlatitude regimes and one geographically unresolved island; use the original single Direct morphological/reproductive trait sectors individually, not the 928 three-trait complete-case.

## Risk to the initial direct source comparison
Original direct comparison [2026-11-15 protected workflow #43702039443](https://github.com/zuizui0223/island/actions/runs/43702039443) compares strictly native and nonendemic plant species documented on the island with candidate native GIFT mainland source plant species NOT documented anywhere in the island's FULL accepted flora. This preserves genus-level assembly differences but may be strongly confounded by different plant range sizes, GIFT checklist representation, botanical record effort and the sampling of plants that have Source Direct flower traits.

The archive's number of **distinct mainland GIFT source entities reporting each species** is a usable, directly observed documentation/source-availability **proxy**, NOT an independent biological geographic range area or actual seed arrival. It is still useful to compare source plants at similar data-driven mainland recording opportunity.

## Fixed category matching
For each island and strict native_nonendemic species also listed in one of its own 5 candidate GIFT native mainland regions, and every other candidate native mainland species NOT recorded in the island's entire full flora, restrict to accepted species with an unambiguous original high/medium Direct category for the SAME trait. Precondition >=5 recorded strict native source-listed species and >=10 unrecorded source candidate species on that axis.

Assign accepted species to 5 **predeclared** bins by total number of distinct GIFT native mainland entities reporting the species across all 748 original source regions:
- exactly 1
- 2–3
- 4–10
- 11–30
- 31+.

For each bin b with ≥1 directly labelled species in BOTH the recorded strict native list and candidate-unrecorded source list, let w_b = min(n_native_b,n_source_b). Compute `Δ_b` as the actual mean original binary trait state among the native names in b minus the actual mean among source candidate-only names in b. The breadth-matched direct difference is `Σw_b Δ_b / Σw_b`. **No genus conditioning** is applied: genus composition is the direct candidate sorting signal being studied. For each SAME island, also compute the unadjusted native-minus-source trait frequency contrast on the identical eligible flora and their difference. Report every original island's matched bin count and overlap total; require `Σ w_b≥5` for inclusion and ≥30 geographically supported target islands across ≥10 original 10° blocks to report geographically blocked conditional 999-resample intervals. Also compare first-GIFT-mainland-source-group weighting and retain negative/null results.

**This does not equalize genus composition, flowering phenology, dispersal opportunity, propagule pressure, native plant introduction/establishment effort, source-specific studies or phylogenetic histories.** GIFT mainland source entity counts may be biological geographic breadth but also reflect source region size and compiled publication effort. This observational control cannot identify a colonization event, isolate seed transport from establishment, or prove preferential survival and reproduction.

As a separate direct comparison on identical original island floras with ≥5 strict native_nonendemic Direct carrier species and ≥5 original introduced Direct carrier species, compare introduced minus native_nonendemic category fractions; record denominators and block/first-source-group uncertainty. Introduced taxonomy must be original `origin_status=introduced`, not inferred from native-nonendemic or GIFT/WCVP presence, and source-introduced may include horticulture/planting/naturalization stages, so it is a documentary community contrast, not measured establishment. Do not treat repeated island/species/site documents as independent biological evolution.

Strict native endemic status remains independently crosschecked with original floristic and origin labels; the original direct endemic same-island same-genus test only becomes inferentially informative if **≥10 distinct actually state-variable genus IDs**, ≥30 original eligible source islands and ≥10 original physical blocks; original taxonomic 1-value invariant-genus zeros are unidentifiable, not evidence for no evolution.

All raw source island, plant/taxon, original status and GIFT source region labels stay inside protected `island` Actions. Only anonymous numerical summaries and source fingerprint validation exported to `fd` draft PR. A next causal level requires independent native-range/phylogeny/pollen visitation/seed-set observation.
