# Own GBIF floral joint diversity: exact local margin partition

Source: [successful protected run 45566580429](https://github.com/zuizui0223/island/actions/runs/45566580429), original commit `17db797d5cb5d0b820455d3ae214dbe316a3d872`. [All 96 panels recovered from the job log](verified_panels.json); artifact metadata digest recorded, not a ZIP byte verification or rerun.

## Actual inputs (not a GIFT island flora)
- `fixed/canonical/input/chapter1_status_flora.csv.gz`: own GBIF incidence roster with attached floristic status.
- `fixed/source/current/evidence/direct_only_trait_ledger.csv.gz`: own integrated species-Direct high/medium single-state evidence.
- Corrected geography, no mainland donor lists. Native subset depends on appended status. Historical "GIFT4505" labels in source are preserved for replay only. Ownership: island.

## Results
96 planned panels, 43 meet geographic/sample/fit support gates. These are NOT 43 significant slopes. Descriptive mean across recorded island units with >=5 carriers, before geographic fit:

| Trait pair | Islands | Observed joint Q | Local fixed-margin expectation | Genus × grade component | Within-stratum component |
|---|---:|---:|---:|---:|---:|
| flower_primary_color × floral_form | 2323 | 0.900499 | 0.918255 | -0.018907 | 0.001151 |
| floral_form × floral_symmetry | 2610 | 0.837824 | 0.876407 | -0.038555 | -0.000028 |
| floral_form × tube_depth_class | 2274 | 0.932197 | 0.939559 | -0.007230 | -0.000133 |
| floral_symmetry × tube_depth_class | 2368 | 0.771239 | 0.793105 | -0.021858 | -0.000007 |

The average departure from local fixed margins is negative for all four pairs in this particular broad cohort. Most numerical departure is reproduced by conditioning on genus and evidence grade. It is an exact combinatorial partition, not a causal explained fraction; morphology axes can share ontology. The conditional within-stratum remainder near zero is not proof that species-level evolution is absent.

Only 1/43 joint-Q isolation slopes have a pointwise bootstrap interval excluding zero. The within-stratum term has 9/43 such slopes; the repeated form × symmetry results are small, nested, and not multiplicity-controlled. No positive broad universal joint-Q decline or pollinator mechanism claim. All negative and insufficient-support panels retained.

Keep outside the frozen GEB H1–H4 paper. Future tests should use a coherent species-level resampling design and examine whether ontology-linked traits inflate apparent association. Official GIFT-only work belongs to fd, and must not cite this result as independent external replication.
