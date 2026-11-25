# Geographic neighborhood reference (N1c) substantially reduces apparent floral packing

**9 November 2026 — public-source result, observational only.** [Validated source-pinned full workflow #42584462875](https://github.com/zuizui0223/fd/actions/runs/42584462875) produced the complete `results/local_geography_ses/` aggregate artifact. The published method and original-distance reference are [N1c protocol](../../docs/LOCAL_GEOGRAPHIC_NULL.md) and [implementation](../../src/fd/local_geography_ses.py).

## Question

After recording frequency of flower-colour/form category profiles is accounted for, is a focal island's within-flora profile-nearest-neighbor distance still unusually low **relative to its nearest observed island floras** rather than the entire broad region?

N1c builds a leave-focal-out pool from the **10 closest other eligible islands in the same broad analysis regime**, regardless of geological class. It samples fixed observed-profile richness, giving each profile a weight equal to the number of those 10 donor islands with records. N1b uses the same rule but every other regional island. The original distance D is mean set-Jaccard dissimilarity for colour and form. The endpoint is average distance from each unique focal profile to its nearest OTHER focal profile (exclude the diagonal). A negative SES means nearer recorded categories than the null expects. N1c is descriptive only; it does **not** condition on genus, independently verified source ancestry, true plant absences or detection.

All values below are from one validated CI artifact, not stitched across distinct Monte Carlo runs.

## Main results: uniform 30-distinct-profile comparison

| All scope, 30-profile observed thinning | Continental, 130 islands | Oceanic, 130 islands |
|---|---:|---:|
| Median observed nearest-profile D | 0.290527 | 0.275905 |
| **N1b: regional weighted reference, median D** | 0.286164 | 0.290890 |
| **N1b: median SES** | **+0.165571** | **−0.554704** |
| N1b: fraction SES < 0 | 42.3% | 84.6% |
| **N1c: nearest-10 weighted reference, median D** | 0.294432 | 0.282709 |
| **N1c: median SES** | **−0.050883** | **−0.231046** |
| N1c: fraction SES < 0 | 52.3% | 70.8% |
| N1c: median observed-minus-null D | −0.001216 | −0.005534 |

For oceanic islands, the median negative distance departure shrinks from −0.013566 to −0.005534 under the nearest-10 reference (a descriptive ~59% reduction in the magnitude of the reported median departure). The two SES medians differ because their null reference distribution shifts; they are NOT estimates of a causal fraction explained by geography. Continental islands' median SES is approximately zero under the geographically local comparator.

The nested Direct scope gives a consistent *descriptive* pattern for equal30 nearest D: N1b SES +0.1634 continental (123 islands) versus −0.5432 oceanic (119), and N1c SES −0.0483 continental versus −0.2583 oceanic, with 52.8% versus 68.9% of N1c scores negative.

Original within-island **mean distinct-pairwise D** tells a similar but weaker story: All/equal30 N1c SES continental −0.0075, oceanic −0.2205. These endpoints are dependent measurements of the same categorical repertoire, not independent replication.

## Same-island paired audit: reference shift is not population replacement

[Source-pinned exact same-island comparison CI #42584985224](https://github.com/zuizui0223/fd/actions/runs/42584985224) exported `local_geography_matched_reference_comparison.csv`, evaluated on the **intersection of targets with estimable N1b and N1c**. Every included island retains its original observed distance and 30-profile thinning; only the null reference changes. No new geological p-value is generated.

| Equal30 nearest-profile D, paired-median statistic | All continental (n=130) | All oceanic (n=130) | Direct oceanic (n=119) |
|---|---:|---:|---:|
| **Within-island N1c SES − N1b SES**, median | **−0.165075** | **+0.285840** | **+0.305889** |
| Within-island N1c null mean − N1b null mean (original D), median | **+0.004094** | **−0.007306** | **−0.007886** |
| Negative N1b SES fraction | 42.31% | 84.62% | 84.87% |
| Negative N1c SES fraction | 52.31% | 70.77% | 68.91% |

For oceanic All, nearby recorded floras make the *expected* nearest-profile original D smaller by a median 0.00731 relative to broad-regional prevalence; hence the observed 0.276 nearest-profile D is less exceptional. In contrast, nearby mainland-derived island floras yield a slightly *larger* expected D than their region-wide source reference. This supports **different empirical source-neighborhood compositions** as a descriptive driver of the null-score change. It does not identify underlying taxa, dispersal or the causal geological origin effect.

## Geographic and support limitations

- The median distance to the **10th nearest eligible donor** is **386.68 km** on continental islands versus **1149.21 km** on oceanic islands. Thus nearest-10 does **NOT** hold physical spatial range constant between origins. The oceanic reference is often geographically much broader.
- All/equal30 has usable nearest-10 support for **130/130 continental and 130/130 oceanic** islands. Full observed-n support is unavailable for many continental islands (**111/130** potentially estimable, versus **128/130** oceanic); do not interpret their full-n comparison as the same population. The paired 30-profile design is the least sample-selected descriptive comparison.
- Direct/equal30 has **123 continental and 119 oceanic** scoreable islands; one oceanic island lacks sufficient scored 30-profile support.
- Nearest island floras may share **identical plant species** or strong genus history with the focal island. Neither an observational N1b nor a nearby-island N1c reference controls repeated species, phylogenetic supply or the data-reporting process.
- N1c does **not** generate new origin-group significance p values. The independent 10-pair geological matching and old 5-spatial-component audits still do not establish a multiplicity-supported geological effect.

**Scientific statement supported now:** recorded floral phenotypes cluster relative to a broad regional profile-frequency sample in oceanic islands, but this excess is **substantially smaller under geographically nearest recorded floral source mixtures**. This is consistent with, but does not prove, spatially structured lineage assembly and sampling opportunities as a major reason for the prior negative SES. Do NOT say 'geography explains 59% of biological selection', 'oceanic adaptation proven', or 'functional hypervolume contracted'.

## Next falsification gate: taxon-name-held-out N2

The [private-only source/genus N2](../../docs/N2_SOURCE_DISTANCE.md) uses observed genus counts and candidate GIFT native-mainland floral donors, with all evaluated taxon names held out across donor entities. Existing N2 code additionally checks heldout-name predictive accuracy, donor coverage, expected-versus-observed distinct richness and per-genus carrier-level generative simulations, since a fixed-n weighted profile sample is not an exact conditional trait permutation. **There are no empirical N2 results yet; restricted checkpoint data are needed.**

## Reproduction and data rights

```bash
python -m fd.local_geography_ses \
  --input work/island_public --g0 work/fd_g0 \
  --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
  --output results/local_geography_ses \
  --neighbors 10 --draws 999 --thinnings 99
```

The CI artifact contains non-identifying grouped `local_geography_ses_summary.csv`, `local_geography_support.csv`, and method/source-hash manifest; no focal island identities, pair membership, neighbor coordinates, or restricted species-to-trait and GIFT source mappings are published.
