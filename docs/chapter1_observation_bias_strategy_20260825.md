# Chapter 1 observation-bias robustness strategy — 2026-08-25

## Why this matters

The canonical Chapter 1 question is **when and where isolation-associated floral/reproductive filtering is detectable**. That question can be distorted if heavily surveyed islands, trait-rich floras, or heavily studied archipelagos contribute much more information than poorly surveyed systems.

The current pipeline already removes one major source of effort bias: raw GBIF record multiplicity is collapsed to unique `island_id × accepted_species` incidence before floral composition is built. Repeated records of the same species therefore do not increase the species numerator or denominator.

However, the grouped-binomial likelihood currently uses the number of trait-resolved species (`trials`) as information weight. Across the current all-native direct-evidence table, this varies strongly by island. In the frozen 2026-08-25 inputs, median / maximum trial counts are approximately:

- flower colour: 102 / 2193 species;
- floral form: 33 / 861 species;
- self-incompatibility: 40 / 1104 species.

Therefore a small number of very information-rich islands can have substantially more likelihood weight than a typical island even though record multiplicity itself has already been removed.

Hawaii should **not** be manually excluded merely because it is intensively studied. The analysis should instead be designed so that Hawaii, Taiwan, Japan, large Mediterranean islands, Indonesia, or any other well-sampled system cannot dominate the global inference by observation effort alone.

## Four distinct observation processes

1. **Occurrence multiplicity** — repeated GBIF records of the same species on an island.
   - already neutralized in the primary response by unique species incidence.

2. **Flora completeness** — poorly surveyed islands may be missing true resident species.
   - requires observation-effort diagnostics and coverage sensitivity.

3. **Trait-resolution completeness** — only some observed species have direct trait evidence.
   - requires island × trait direct-evidence fraction audits and weighting sensitivity.

4. **Geographic / archipelago study intensity** — entire island groups can be systematically better surveyed than others.
   - requires within-archipelago / spatial-block contrasts and leave-one-block-out influence tests.

These processes must not be collapsed into one generic `n_records` covariate.

## Bias-resistant analysis ladder

### A. Species-incidence response — canonical and already implemented

Each species has at most one contribution to an island trait denominator. Raw record abundance is not treated as ecological abundance.

### B. Observation-quality gate — required before manuscript freeze

Reconnect the existing `observation_diagnostics.py` layer to the canonical when/where workflow. At minimum retain and report:

- `n_records`;
- `n_species`;
- `n_datasets`;
- preserved-specimen fraction;
- cultivated-record fraction;
- coordinate uncertainty;
- latest observation year;
- accepted-angiosperm species and record coverage where reviewed scope is available.

The existing config thresholds remain an inclusion / diagnostic sensitivity, not a way to declare biological absence.

### C. Information-weight robustness — required sensitivity

The canonical grouped-binomial model uses trait-resolved species count as statistical information. Repeat the response-vector omnibus under predeclared alternative effective-trial rules while preserving each island's observed trait proportion:

1. canonical species-count weighting;
2. capped effective trials at 100;
3. capped effective trials at 50;
4. capped effective trials at 20;
5. equal-island weighting (`effective_trials = 1`).

For a row with observed share `p = successes / trials`, use

`effective_successes = p × effective_trials`.

These are robustness pseudo-likelihoods, not alternative biological abundance models. Their purpose is to test whether the when/where result requires highly trait-resolved islands to carry disproportionately large likelihood weight.

### D. Trait-resolution coverage robustness — required sensitivity

For each `island × stratum × trait`, calculate

`direct_trait_fraction = n_direct_trait_species / n_observed_stratum_species`.

Report its distribution by region and test whether it correlates with isolation. Repeat the when/where omnibus after either:

- restricting to a predeclared minimum direct-trait fraction; or
- adding a trait-specific coverage term in a declared observation model.

Do not infer that high `direct_trait_fraction` means the island flora itself is complete; this addresses trait-resolution bias only.

### E. Within-block / archipelago robustness — strongest protection against geographic observation bias

Cluster-robust standard errors alone do **not** remove a block-level confounder; they adjust uncertainty after fitting. Therefore add a sensitivity in which isolation effects are identified from **within spatial block / archipelago contrasts** rather than between highly different geographic regions.

Preferred implementation:

- include response-specific spatial-block fixed intercepts where estimable; or
- estimate block-level isolation responses and combine them at the regional level;
- require at least two (preferably three) islands with isolation variation per block.

Any fixed block-level observation process then cancels from the isolation contrast.

### F. Leave-one-block-out influence gate

Repeat the headline northern-midlatitude and tropical WHERE tests and the northern-vs-tropical BETWEEN-WHERE test while excluding one spatial block at a time.

Report:

- fraction of leave-one-block-out runs preserving each conclusion;
- the most influential block;
- whether deletion of a Hawaii-containing block changes inference;
- whether any single block changes support classification.

The manuscript-level result should be downgraded if it depends on one high-information archipelago/block.

## Coverage-based rarefaction: useful but not a magic correction

Coverage-based standardization is preferable to equal sample-size rarefaction when comparing biodiversity samples of unequal completeness. But GBIF occurrence data are opportunistic and spatially aggregated; ordinary Good–Turing coverage assumptions can be violated. Therefore record-based sample coverage should be used as an observation diagnostic / sensitivity layer, not as proof that missing species are missing at random.

## Manuscript claim gate

The Chapter 1 when/where conclusion is considered observation-robust only if:

1. northern-midlatitude and tropical WHERE remain supported under capped/equal-island weighting;
2. northern-vs-tropical BETWEEN-WHERE remains supported under the same weighting sensitivities;
3. the result persists in native non-endemic flora;
4. observation-effort / direct-trait-coverage sensitivity does not reverse the result;
5. leave-one-block-out analysis shows that no single archipelago/spatial block is necessary for the conclusion.

If one of these fails, report the dependence explicitly and downgrade the corresponding where/when claim rather than adding more mechanism-specific data to rescue it.

## Interpretation of Hawaii

Hawaii is a useful stress test, not a special-case exclusion.

- Its many occurrence records should not increase species-level response weight directly: record multiplicity is already collapsed.
- Its trait-resolved species count can still affect grouped-binomial information weight.
- Its shared survey history with the Hawaiian archipelago can create block-level observation structure.

Therefore the correct test is **equal/capped island information + leave-Hawaii-block-out + within-block inference**, not `remove Hawaii because it is unusual`.
