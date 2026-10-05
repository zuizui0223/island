# Far-tail genus attribution diagnostic — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

Parent signal: `yellow_orange__butterfly_deep_tube_given_colour` above the support-derived
783.017 km hinge.

## What was added

The diagnostic reconstructs the species-level numerator/denominator behind the frozen
island-level coupling counts, ranks genera in the >783 km tail, removes one genus at a
time and refits the same hinge model, and repeats the hinge fit on floristic-status
strata.

Successful workflow run: **37302624291**.

## Headline result

### Northern mid-latitude

Direct all-observed baseline:

- 1,514 complete model islands;
- 62 >783 km islands in 17 spatial blocks;
- post-hinge slope = **+0.4589**, P = **0.0212**.

Top direct far-tail deep-tube contributors:

| genus | deep occurrences | deep islands | share of deep occurrences |
|---|---:|---:|---:|
| Oenothera | 23 | 13 | 10.6% |
| Aloe | 17 | 11 | 7.8% |
| Lonicera | 12 | 12 | 5.5% |
| Oxalis | 12 | 12 | 5.5% |
| Tecomaria | 9 | 9 | 4.1% |
| Rhododendron | 8 | 8 | 3.7% |
| Jacobaea | 8 | 8 | 3.7% |
| Calystegia | 7 | 7 | 3.2% |
| Hemerocallis | 7 | 7 | 3.2% |
| Alcea | 6 | 6 | 2.8% |

The first ten genera account for about **50%** of direct far-tail deep occurrences.

The classic candidate genera proposed as a possible explanation are not the source of
the signal: **Isoplexis, Canarina and Clermontia are absent; Lotus occurs but contributes
zero deep-tube successes.**

Leave-one-genus-out shows concentration but not a single-genus effect:

- remove Oenothera: slope **+0.2979**, P = 0.115;
- remove Aloe: slope **+0.3507**, P = 0.076;
- remove Rhododendron: slope **+0.3532**, P = 0.060.

More importantly, the direct far-tail numerator is dominated by unresolved or introduced
floristic status:

- unresolved: **150/218 = 68.8%** of deep occurrences;
- introduced: **49/218 = 22.5%**;
- confirmed native: **19/218 = 8.7%**.

When the same hinge is reconstructed on confirmed native flora, the northern positive
tail does not persist:

- all-native Direct: 221 complete islands, **23 tail islands / 8 blocks**,
  slope **-0.6673**, P = 0.556;
- native-nonendemic Direct: 221 complete islands, **23 tail islands / 8 blocks**,
  slope **-0.6923**, P = 0.494.

Given only eight native far-tail blocks, these are not precise negative-effect estimates.
The important diagnostic fact is that the positive all-observed northern tail is not
recoverable in the confirmed-native stratum.

### Tropical

Direct all-observed baseline:

- 855 complete model islands;
- 349 >783 km islands in 53 spatial blocks;
- post-hinge slope = **+0.2015**, P = **0.00094**.

Unlike the northern result, the positive tail persists after restricting to native flora:

- all-native Direct: 120 complete islands, **67 tail islands / 26 blocks**,
  slope **+0.8074**, P = 3.53e-6;
- native-nonendemic Direct: 119 complete islands, **66 tail islands / 26 blocks**,
  slope **+0.9844**, P = 1.22e-19.

Confirmed-native direct far-tail deep occurrences are led by:

| genus | deep occurrences |
|---|---:|
| Cordia | 49 |
| Thespesia | 42 |
| Ipomoea | 30 |
| Physalis | 20 |
| Decalobanthus | 15 |
| Tecoma | 9 |
| Crescentia | 8 |
| Hibiscus | 6 |
| Abelmoschus | 5 |
| Brunfelsia | 5 |

These genera span multiple plausible pollination/reproductive strategies, so the raw
`butterfly-associated` label must not be interpreted as realized butterfly pollination.

## Interpretation

This diagnostic changes the earlier post-hoc reading.

1. The northern far-tail positive response cannot currently be used as evidence that the
   tropical signal is a generic >783 km isolation phenomenon. It is highly dependent on
   the all-observed/status-unresolved layer and disappears in the confirmed-native flora.
2. The tropical far-tail response remains positive after restricting to native and
   native-nonendemic strata, but its biological interpretation is unresolved. Several
   leading genera are characteristic pantropical/coastal lineages, so dispersal filtering
   (including coastal/ocean-dispersed floras) is a live alternative explanation that must
   be tested before attributing the pattern to pollination.
3. The tropical response is not attributable to the proposed classic Macaronesian/Hawaiian
   bird-pollination genera. Any pollinator interpretation must be made genus/species by
   genus/species and remain secondary to the raw floral-architecture result.
4. Native-stratum P values should not be treated as strong inferential evidence: only 26
   far-tail spatial blocks contribute to the tropical native fits, and the current
   cluster-robust normal approximation can be anticonservative with few clusters. The
   native samples are also support-imbalanced around the 783 km knot, so pre-hinge slope
   estimation is weak.
5. This entire hinge family remains post-hoc and outside frozen Chapter 1 inference.
