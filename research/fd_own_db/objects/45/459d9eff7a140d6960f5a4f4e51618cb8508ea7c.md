# Exact plant taxonomy and source-only donor decomposition — 2026-11-24
[Protected run 45384576074](https://github.com/zuizui0223/island/actions/runs/45384576074) SUCCESS; anonymous artifact 14913864601 SHA256 3254b0adee1d4ad7bdc9c1138698e8c6bb751cf2d1728bbe2dcdbbba2c5302bb. Original 4505 botanical /1558 tropical /8264 physical island parents; 748 mainland regions. All 24 panels estimable; analytic conditional expectation (block means) exactly decomposes observed diet contrast coefficient; all 199-draw Monte Carlo mean replay checks pass. Same original distinct-species target and source pool, not causal mediation or within-species evolution.

## Same 30 islands, same joint form+depth scored taxa, all-recorded min5
| Candidate source / endpoint | Observed island-minus-source coefficient | Genus x grade expectation | Within-block residual |
|---|---:|---:|---:|
| geo5 tubular | .019803 | .113496 | -.093692 |
| geo5 deep | .310323 | .114743 | .195580 |
| geo5 tubular and deep | .169000 | .204551 | -.035551 |
| climate10 tubular | -.004555 | .069979 | -.074534 |
| climate10 deep | .338377 | .150697 | .187680 |
| climate10 tubular and deep | .156366 | .186584 | -.030217 |

Form records on the same 30 islands but different scored taxa: tubular geo5 observed .255493, expected .341604, residual -.086111; climate10 .283316/.351079/-.067762. The positive broad tubular association is compatible with plant composition and depends strongly on trait-record roster. The positive deep residual is not independent of selected scope: all min30 joint 21 islands residual +.053655/+ .056833 but observed - .349249/-.202876; native min5 joint 20 islands residual -.102319/-.148903, observed -.246593/-.228290; high-only joint min5 28 islands residual +.295297/+ .308073. Different island/trait/grade rosters and controls make these non-equivalent estimands. Previous dual-null 0/160 Holm support remains unchanged.

## Strict independent mainland donor audit
All 48 panel×genus(/grade) donor checks insufficient coverage; at least 2 distinct mainland donors in the same block, excluding **every recorded plant name from every matched target island**, require >=5 eligible names and >=50% coverage per island. Largest qualifying cohort only 4 islands (all-recorded min5 climate10 form-only genus); joint 30-island geo5 median coverage .226844 genus/.200000 genus-grade, 0 eligible; climate10 .307255/.255435, one eligible island. No independent mainland coefficient identified. Do not relax donor support to select favourable result.

## Next exact origin decomposition
Protected island commit ef73602c61e09693bf058ec62d03c4e63a200b92, tests native-status sign reversal using 16 panels (depth only and jointly scored form-depth; high-medium/high-only; full-scored threshold 5/30; two source definitions). Fixed original source roster and genus-grade expectation across the full and common native cohort. Exact slope change from all-recorded to literal native = island selection + change in scored-count control + roster change, the last exactly separated into introduced vs unresolved-origin contributions. Source labels native/introduced/unresolved retained as annotations, no inference of alien plants from missing native label. No result claimed until successful execution.
