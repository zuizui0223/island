# Flora status, endemism, lineage, and model-support plan

## Decision

The acquisition target is no longer global completion of `species x axis` cells.
The analysis is reorganized around the final estimands and the biological stage at
which isolation can change flora composition.

## Primary flora strata

Every `island x species` row must carry source-backed floristic status:

- `native_nonendemic`
- `endemic`
- `introduced`
- `unresolved`

`all_native` is computed from origin status and therefore also includes native
species whose endemic status is unresolved. Status is never inferred from the
number of GBIF island occurrences. Conflicts among sources fail closed to
`unresolved`.

The first public status source is GIFT because its checklist system explicitly
retains native, naturalized and endemic status when the underlying flora provides
it. Islands or taxa not covered by a status source remain unresolved.

## Analyses unlocked by the status layer

1. `Isolation -> endemic_share_of_native` is a separate response, not a nuisance
   covariate.
2. Every focal trait is summarized in parallel for:
   - all native species;
   - native non-endemics;
   - endemics;
   - introduced species as a secondary negative-control stratum.
3. A shift present among native non-endemics is compatible with colonization /
   establishment filtering. A shift concentrated among endemics is compatible
   with endemic lineage turnover or in-situ differentiation. Neither pattern alone
   proves the process; regional flora, environment and lineage remain alternatives.
4. Introduced species are not mixed into the primary native composition. They can
   be used secondarily to test whether human transport weakens the isolation
   gradient.

## Lineage sensitivity

Two operations have distinct meanings and must not be merged:

- **genus-fixed null:** preserve each island's observed genus composition and ask
  whether the observed trait composition departs from expectation under within-genus
  trait reassignment. This diagnoses lineage-composition confounding.
- **genus probability inference:** propagate missing trait states from empirical
  within-genus distributions as a measurement sensitivity. It does not create
  confirmatory direct evidence.

The repository already contains probabilistic genus-inference machinery. A later
step will connect it to the status strata rather than duplicating it here.

## Model-specific support replaces the 50% three-axis gate

For each final outcome (`flower_colour`, `floral_form`, `self_compatibility`) and
flora stratum, report:

- number of species in the stratum;
- number and fraction with direct species-level evidence;
- number of islands with at least 30 direct species;
- isolation 5/25/50/75/95% quantiles and range;
- spatial-block count where available;
- region-specific support.

The operational continuation rule is:

- `>=50` directly supported islands: continue to modelling before more acquisition;
- `30-49`: targeted acquisition only;
- `<30`: do not start a broad acquisition campaign automatically. First test whether
  recoverable missing species can materially expand isolation support.

These island counts are support diagnostics, not biological thresholds.

## Targeted acquisition objective

For outcomes below the support target, rank unresolved direct evidence by:

1. whether one new species record immediately moves an island over the declared
   per-island direct-species minimum;
2. inverse gap to that minimum;
3. number of affected islands;
4. the isolation range of those affected islands;
5. evidence recoverability in the acquisition ledger.

This replaces the previous objective of maximizing the number of ready islands by
collecting widespread common species irrespective of where those islands lie on the
isolation gradient.

## Interpretation boundary

The status split improves process diagnosis but does not by itself prove
colonization filtering or evolutionary change. Final models still require climate,
island area, spatial dependence, region/floristic context, and genus/lineage
sensitivity. The key improvement is that those alternatives can now be separated
instead of being folded into a single whole-flora trait slope.
