# fd v0.1 — scientific and computational analysis contract

Date: 2026-11-08 JST. Source: [island PR #321](https://github.com/zuizui0223/island/pull/321), frozen aggregate inputs from source commit `d4a5dafab5f0737ea6c1c1c0263fc20decf6035c`.

## 1. Question and estimands

The primary question is whether island floral assemblages show **local compression and between-island differentiation of recorded colour/form phenotype space**. These are logically distinct properties and may coexist.

- **H1 local restriction:** Within-island taxon-weighted functional dispersion is lower than a source-/lineage-appropriate conditional expectation. Target quantity is Rao's Q, not a Gaussian KDE volume.
- **H2 spatial differentiation:** Between-island matching distances exceed an incidence-appropriate, species/source-aware expectation. Primary distance target is symmetric nearest-profile dissimilarity; exact-profile Jaccard is a comparator. Turnover/nestedness of exact labels cannot substitute for distance-aware differentiation.
- **H3 spatial drivers:** Isolation and island area explain systematic patterns in H1/H2 *after* source pool, genus/lineage, region, geography, detection/trait coverage, floristic status and spatial dependence are considered. Test isolation × area on continuous scales; do not dichotomize arbitrarily.
- **H4 functional interpretation (separate gate):** Observed phenotype-space differences can reflect pollination-access and signal architecture. They are not measured pollination services. Attribution to channel loss, selection, adaptation, evolution or functional replacement requires independent data.

The old claim 'islands necessarily converge to one simple floral syndrome' is a rival hypothesis, not a prior. Observed higher turnover by itself does not establish environmental selection.

## 2. Observation unit and profile

- Unit: retained **taxon-name unit** on each island, not established reconciled biological species. Some names include known hybrids; taxon range size and observational uncertainty can affect results.
- Existing island datasets include native, introduced and unresolved origins. Any native-only inference needs a different, audited input.
- Profile: `(flower_primary_color state set, floral_form state set)`. In the multistate scope the entire state set remains one taxon profile; the Cartesian product is NEVER treated as observed combinations.
- Evidence scopes: `all`, `direct`, `all_single_state`, `direct_single_state`. The four scopes share the same 496-profile All reference; Direct is a source-evidence class, not solely field sampling.
- Eligible islands: >=50 observed complete name units, consistent with the upstream descriptive comparison. Missing evidence means no observation, **not** a biologically absent profile.
- Pair IDs: upstream `pairwise_observed_and_null.csv.gz`, preserving the pair selection. These are repeated, nonindependent pairs. No independent-pair sample-size claim.
- Source-pool/ancestral state is not identified by a broad geographic region alone.

## 3. Exact dissimilarity; no ill-supported hypervolume

For profile `x` and `y`, let `C_x,F_x` be nonempty recorded state sets.

```text
d(x,y) = 0.5 * [1 - |C_x intersect C_y| / |C_x union C_y|
              + 1 - |F_x intersect F_y| / |F_x union F_y|]
```

Thus single-state colour or form mismatch contributes 0.5; mismatch in both contributes 1.0. A multistate record is one point defined by its sets. There is no trait-scale calibration: equal weights are a stated choice, to be sensitivity-tested.

With `p_{i,a}=n_{i,a}/N_i` as the proportion of taxon-name carriers of observed profile `a`:

```text
RaoQ(i) = sum_a sum_b p(i,a) p(i,b) d(a,b)
Q_unbiased(i) = N_i / (N_i-1) * RaoQ(i), if N_i>1
Q_equal(i) = mean_{a,b in present profiles(i)} d(a,b)
```

The finite-name formula is a conventional without-replacement correction to the name-unit average; it does not repair trait missingness, effort or source composition.

For islands A/B, score directions from each recorded profile to the nearest profile in the other island:

```text
NN(A->B) = sum_a p(A,a) min_b d(a,b)
NN_sym(A,B) = [NN(A->B) + NN(B->A)]/2
cross(A,B) = sum_a sum_b p(A,a) p(B,b) d(a,b)
```

`NN_sym` is a bounded 0..1 distance between *supports* (modulated by carriers), not a Wasserstein distance, an abundance-compositional metric or a taxon-excluded comparison. Identical profile supports give zero regardless of differing frequencies. `cross` is affected by within-island dispersion; it is secondary. Exact Jaccard ignores how dissimilar two nonmatching profiles are.

PCoA is visualization only. Its 2 axes cover 15.35% of positive inertia and 17 axes are needed for 80%; no ecological hypervolume is inferred from the 2-D projection. Gaussian KDE volumes in these truncated axes are not substitutes for complete categorical trait breadth.

## 4. Sequential inference gates

**G0 (implemented):** Reconstruct D on the 496 state-set profiles from public `profile_coordinates.csv`, ignoring coordinates. Check reference identity, name counts, profile counts, island metadata and matched-pair identity. Compute Q and NN on All/Direct and single-state cohorts, record checksums and missingness.

**G1 (next):** Fixed row/column **profile incidence** null (Curveball) on same matched islands and region/scope; preserve per-island unique profile counts and regional profile prevalence. This can test *equal-profile* Q and NN support distances. It does **not** preserve original profile carrier counts, lineage, source or colonization, and cannot be used as the null for carrier-weighted RaoQ without an extra count-preserving design. Explicitly report MCMC mixing and approximate tails and control multiplicity across declared families.

**G2:** Source/genus-conditioned null using *authorized private inputs*, maintaining one regional taxon-profile reassignment across islands, all taxon membership and genus/state frequencies; assess both local Q and island-pair NN under an explicitly fixed cohort. This cannot be reconstructed from the currently public aggregate. The upstream shared-name-exclusive PR321 inference belongs to a different estimand and is not reproduced by G0/G1.

**G3:** Fit isolation, source access, area, climate, regional/regime, status, lineage and coverage models to the G2 effects; account for island repeat degree in pair comparisons (e.g., island-block bootstrap or multi-membership effects) and spatial dependence. Freeze analysis support and use source-matched distances if available. Do not reinterpret the current generic continent-distance zeros as real water gaps without auditing geography.

**G4:** Interpret differences in recorded colour/form architecture alongside independent pollinator guild, functional access, tube depth, symmetry, flower size, reproduction and observed interaction data. Separate signal, access, reward and selfing; do not infer pollinator identity from floral templates or claim causal adaptation.

## 5. Falsification and reporting requirements

A defensible positive conclusion requires that directions persist across support-matched All/Direct and a biologically interpretable primary scope, against an appropriately conditioned null and with effect size/uncertainty. A positive fixed-margins result can be real without a positive genus-conditioned result; the two test different questions. 'Significant in one region, not in another' is not a test of regional difference. Any between-region contrast must itself be tested.

Required audit fields: frozen cohort, number of island units and *distinct* pairs, shared-taxa policy, null constraints, tail and dependence correction, trait coverage, status, genus/source pool, missingness and source hashes. Overlapping evidence scopes are sensitivity checks, not 4 independent discoveries.

Existing [PR321](https://github.com/zuizui0223/island/pull/321) positive results remain evidence from that exploration and are not automatically confirmed by these new function-distance statistics. G0 is descriptive, so no p-values or 'supported' labels are produced.

## 6. Data redistribution boundary

`island/docs/CHAPTER1_DATABASE_RIGHTS_AUDIT.md` is controlling. Only explicitly public **aggregated island×profile counts** and references can enter `fd`. Species-to-profile rows, provider source ledgers, exact identity-level audits and restricted raw source extracts must not be pushed. The importer pins one immutable public commit and logs exact local SHA256. Full reconstruction from original evidence requires the restricted private checkpoint, not a public repo copy.
