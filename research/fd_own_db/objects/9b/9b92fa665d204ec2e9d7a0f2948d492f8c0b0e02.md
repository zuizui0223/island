# Genus- and regional-availability-conditioned floral repertoire test

**Research home:** `zuizui0223/fd` independent of `island` Chapter 1.  
**Status:** exploratory follow-up, methods written after prior geological-origin comparisons. Not a prospectively preregistered test or causal inference.

## Why this gate matters

The existing [37 matched oceanic–continental island pairs](../results/20261108_matched_flower_repertoire/README.md) differ in **exact recorded** flower-colour × form repertoire after matching island area, climate, distance and recording coverage (original All/Direct BY q≈0.00042). The original trait-distance-aware kernel test does **not** reject. This remains compatible with **taxonomic sorting**: oceanic and continental islands recruit different genera and species, which already have different floral traits.

The next question is **conditional**: does the repertoire difference remain unusual if we fix the observed island × taxon-name incidence and the distribution of flower profiles among taxon names within genera, plus geographic availability?

## Source checkpoint and strict reconstruction

Use frozen, access-controlled historical GitHub Actions source archives (they are accessible to the linked GitHub account but **NOT redistributable**):

- `zuizui0223/island` run 34191508045, artifact 10042377593: reviewed Direct trait ledger and the 106,295-name axis population.
- `zuizui0223/island` run 34232450884, artifact 10058653212: frozen All/Direct trait snapshots, complete island floristic-status table.
- `island` PR #321 source commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`, exact ontology and selection contract.

Raw frozen gzip member SHA-256 expected:
```
direct_species_trait_ledger.csv.gz  7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577
species_axis_coverage.csv.gz       a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a
snapshot_all.csv.gz                2d95136fc7df1dcce92aab738ce0364cdb1434746e1f7b794c529900192cbf4d
snapshot_direct.csv.gz             d7d6a0372ce7f3ebe958047a222ae0504f01564e364f93bab360273235c6661b
chapter1_status_flora.csv.gz       6b5006322d15dd87aa3aad3b44d050d4c6420485521be314ae511f39e05b859a
```

The file hashes above were verified against the frozen `island` analysis manifest before this exploration. The runner must fail closed if any actual source member has a different digest.

Validate exact PR #321 reconstructed statistics before testing: 2,847 recovered reviewed Direct cells, 106,124 analysis-eligible taxon-name units, 16,739 All and 8,455 Direct complete colour/form names, 496 All / 467 Direct distinct profiles. Compare each observed island's complete-name count and distinct-profile count with `fd` G0 output (all 1,367 north-mid+tropical eligible-support islands in All, 1,292 in Direct). Stop on any discrepancy. One multistate value is a recorded set, not an invented colour×form cross-product.

## Fixed cohort and null hypotheses

Reconstruct the exact existing 37 continental–oceanic island pairs via `fd.repertoire.match_origin_islands`, using the original fixed geology crosswalk (not a new flower-trait match). Twenty-three pairs are northern-midlatitude; fourteen tropical. Both All and Direct share pairs; the scopes overlap.

The statistic matches the earlier `fd` repertoire metric exactly:
```
For island i: p_i(a) = 1/n_present_profiles_i if profile a is recorded, else 0.
For pair j: delta_j = p_oceanic(j) - p_continental(j).
For sampled region S: mean_delta = average_j(delta_j).
L2(S) = mean_delta^T I mean_delta
K(S) = mean_delta^T [1 − colour/form state-set Jaccard dissimilarity] mean_delta
```
The paired **global** 37-island-pair statistic is the primary endpoint, with northern-midlatitude and tropical calculations as diagnostics. Unequal genus abundance does not generate multiple profile points by itself.

**Conditioning:**

- **Fixed island × taxon identity matrix:** the same recorded name occurs on the same islands across all randomizations, preserving name geographic range, island genus composition and entire completeness/missingness pattern.
- **Fixed reviewed All/Direct profile population:** shuffle *whole existing profiles* of complete taxon names, never shuffle colour and form separately.
- **Primary null: genus × four-regime geographic incidence signature** (one signed four-regime availability pattern per name, determined from observed floras): shuffle existing profiles only among eligible names sharing a genus and that same set of observed broad climatic/geographic regions. One random profile assignment per name is reused at all islands and among the 37 pairs. Thus the regional taxon profile multiset and geographic availability are fixed by construction. The signature is **observed source-region opportunity**, **not** a reconstructed mainland colonization pool.
- **Sensitivity null: genus-only global permutation** preserving one globally coherent taxon profile and island incidence. It does not preserve within-genus regional trait distribution.
- Optional native-only assessment: restrict observations to explicitly native entries, require ≥20 complete native name units per endpoint, keep the original matches and report selection attrition. Do not treat unresolved status as introduced or absent.

Use deterministic independent seeds, at least 1,999 conditional draws, and finite tail `p=(1+#(T_null≥T_obs−tolerance))/(B+1)`. Adjust BH and BY over four **global primary tests** (two metrics × All/Direct) for the primary regional-signature null. Per-regime and alternative-null tests are **secondary sensitivities**, with their own declared family. The original matched-pair label-swap test and the new genus-permutation test are different reference distributions; results must not be interchanged.

## Supported interpretation

If observed global L2 returns within the genus+availability reference, that supports **consistency with sorting among existing genus-linked floral repertoires**; it **does not prove equivalence**, nor demonstrate how source assembly operates. If global difference remains unusual, there is **residual taxon-profile sorting beyond genus and the coarse regional availability**; it still may be driven by native/introduced status, phylogeny within genera, imperfect recorded traits, and geological match uncertainty, not necessarily selection.

Null features **not conditioned**: fully resolved phylogeny, true colonist source pools and arrival histories, trait detection process, flower-specific functional pollinator measurements, and within-species changes. The analyses are post-hoc and use nested All/Direct evidence; avoid declaring independent external replication or a causal geological origin effect.

## Distribution constraints

Original plant taxon×trait rows, source metadata, and island×taxon presence data remain in private, transient workspace. Publish only generic code, input digests, aggregate support counts, observed/expected statistics, Monte Carlo tail estimates, and non-identifying tests. Do **not** commit locally reconstructed species-to-profile mappings, selected island pairs containing identifiers, private source checkpoint ZIPs or the geological source-ID crosswalk. Public `fd` will retain aggregation-only results and enough source hashes to enable authorized future reruns.

## Replay from locally authorized private checkpoints

These commands are **not a public Actions workflow**. Authenticate to the original `island` Actions artifacts on the user's own computer. Reconstruct `species_profiles.csv.gz` through the original `island` PR321 `scripts/run_floral_space_exploration.py` with its frozen hashes and ontology, or supply the same private scoped mappings already generated by that runner. Never add these inputs to the `fd` git index.

```bash
python -m fd.genus_source_null \
  --private-all /secure/island/floral_space_20261106/all/species_profiles.csv.gz \
  --private-direct /secure/island/floral_space_20261106/direct/species_profiles.csv.gz \
  --flora /secure/island/chapter1_status_flora.csv.gz \
  --crosswalk /secure/fd/island_origin_crosswalk.csv \
  --g0 work/fd_g0 \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --draws 1999 --output results/run_private_genus_conditioned
python -m fd.genus_source_null \
  --private-all /secure/island/floral_space_20261106/all/species_profiles.csv.gz \
  --private-direct /secure/island/floral_space_20261106/direct/species_profiles.csv.gz \
  --flora /secure/island/chapter1_status_flora.csv.gz \
  --crosswalk /secure/fd/island_origin_crosswalk.csv \
  --g0 work/fd_g0 \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --draws 999 --native-only --min-native-names 20 \
  --output results/run_private_genus_native
```

**Important:** All paths containing species profile maps, the island-flora source-status table or crosswalk details must remain outside the public repository. The runner emits only two aggregate-only outputs in the requested directory: `genus_conditioned_summary.csv` and `manifest.json`.

## Independent GIFT native-mainland source-repertoire sensitivity

The earlier `island` source-pool study archived GIFT 3.2 native-mainland floras, and outcome-blind island→mainland source-region assignments under `geo_k5`, `geo_k10`, `geo_k20`, `geo50_climate10`. All 74 matched island endpoints have all four modes. This provides a finer **potential source-flora profile baseline** than broad four-regime signatures, but is still not reconstructed ancestral dispersal.

For each source mode and matched island: take the set-union of all resolved colour/form profiles of exact-match taxon names native to its assigned mainland entities. Normalize that unique-profile vector to sum 1. Subtract it from the observed island unique-profile vector, then evaluate the paired origin residual with sign swaps, preserving the same frozen geographic pairs. Do not confuse this *mainland spectrum subtraction* with either genus-permutation null. The source spectrum can be richer than the island spectrum; there is no colonist sampling model that equalizes them.

```bash
python -m fd.mainland_source_residual \
  --private-all /secure/island/floral_space_20261106/all/species_profiles.csv.gz \
  --private-direct /secure/island/floral_space_20261106/direct/species_profiles.csv.gz \
  --flora /secure/island/chapter1_status_flora.csv.gz \
  --source-assignments /secure/island/island_source_assignments.csv.gz \
  --mainland-native /secure/island/gift_native_mainland_flora.csv.gz \
  --crosswalk /secure/fd/island_origin_crosswalk.csv \
  --g0 work/fd_g0 \
  --reference work/island_public/floral_space_20261106/all/profile_coordinates.csv \
  --permutations 9999 --output results/run_private_gift_source
```

The GIFT calculation writes only aggregate comparisons, source-profile overlap summaries and hashes; no restricted source entities, species names, geocoordinates or taxon-to-profile data are published. The `geo_k5` global residual 2 scopes × 2 endpoints are the stated four-test exploratory family; other modes/regions are unadjusted sensitivity checks. Neither a successful source subtraction nor a nonsignificant genus conditional null proves ecological selection, actual ancestral sources, or no functional differentiation.

