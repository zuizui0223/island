# Original GIFT 4,505 botanical island floras × DoPI v2: independent Curveball versus edge-swap positive visitor-network null — 2026-11-18

**Status: original source-protected numerical workflow SUCCESS.** This is a posthoc algorithm-sensitivity comparison on **identical source study × published site/grid × year × month networks**, using only positively documented European flower–visitor records for accepted plant species listed somewhere on the original global 4,505 GIFT island floras. It is NOT a direct original island visitation dataset, actual effective pollination, historical guild loss or evidence of in situ flower adaptation. The **global parent 4,505 floras remain unchanged**.

## Reproducibility and exact source evidence

- Protected [original island Actions #44221323430](https://github.com/zuizui0223/island/actions/runs/44221323430) **SUCCESS**, island workflow commit \`4db093f7b29dc5c0c097c62109f7c110fe59c5b8\` pinning tested fd source \`c1a41c42762458b47e1aef46023ff9ade5cda397\`. Independent [fd Curveball and edge-swap synthetic source tests #44221112205](https://github.com/zuizui0223/fd/actions/runs/44221112205) **SUCCESS**.
- Original SHA-frozen anonymous GitHub artifact **#14391635619**, file \`anonymous_original_dopi_gift_guild_motifs.json\`, verified 2,677,216 bytes, **SHA256 \`d6274a471b786a3121c112065e378e332191309aee91951af9411e192c827195\`**. Original botanical flora, independent species-Direct flowering morphologies, and published DoPI v2 archive hash/size were checked inside original source CI. Original site IDs, publication IDs, exact accepted plant and insect names remain within protected \`island\` Actions only.
- [Pure independent Curveball / matched-source algorithm](../../src/fd/dopi_curveball_paired_degree_null.py), [Curveball row/column-degree preservation and matched-source synthetic tests](../../tests/test_dopi_curveball_paired_degree_null.py), [original secure GIFT/DoPI driver](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_dopi_gift_guild_motifs.py), [pinned protected original workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-original-dopi-guild-motifs.yml).
- Methodological references: [Strona et al. 2014, Nature Communications 5:4114](https://doi.org/10.1038/ncomms5114) original Curveball trades, and [Carstens 2015, Physical Review E 91:042812](https://doi.org/10.1103/PhysRevE.91.042812) on stationary distribution and retaining self-loop/failed trade opportunities. This implementation randomly samples two plant rows and redistributes their unique visitor neighbors while holding common visitors unchanged; failed proposals remain Markov transition opportunities.

## Unchanged biological and source estimand

The response is the difference in the number of *positively documented INSECT GUILD categories* between accepted plants with versus without each GIFT Direct species-wide flower morphology (open radial, tubular/bell, zygomorphic, deep tube), in **one original source publication × published geographic site/grid × observation year/month**. Each plant must document >=2 true identified visiting insect species; every botanical morphology comparison needs >=2 true Direct motif-positive plants and >=2 known motif-negative plants. The source guild categories and exact accepted species names were held fixed. DoPI unreported guilds are NOT ecological absences; source visitation evidence grade1/2/3 is not necessarily effective pollen deposition. These are not necessarily field observations ON any of the global original 4,505 islands.

Two actual binary-incidence randomization algorithms independently preserve **each plant's real insect species degree** and **each real insect species' original number of documented plant hosts** exactly:

1. Previous **checkerboard edge-swap**: rewire two original plant–insect species edges if cross-edges are absent.
2. Independent **Curveball**: select two distinct original plant rows and trade unique visitor species while retaining shared visitors and original row sizes.

Each runs **four independently seeded chains with eighty null snapshots**, with actual accepted trade/swap and chain-mean-disagreement screening. A source network/flower profile is included only if BOTH algorithms pass all existing source-prespecified computational eligibility and convergence-proxy gates. Hence **the same original sites and the same original source papers are compared within each flower morphology**, not two selection-biased algorithm-specific cohorts. An immutable, genuinely unshufflable positive source network is NOT_IDENTIFIABLE, not a no-effect result. There is no claim that finite MCMC burn-in proves perfect uniform mixing.

For each supported original source site, compute observed flower-positive minus negative guild count, the expected difference under each reference and the respective residuals, then the Curveball minus swap residual **on that exact original source network**. Average within original source publication and resample whole distinct publication IDs 999 times; reported 95% intervals are **descriptive exploratory**. The four flower outcomes were selected after an older 27-morphology screen, so no preregistered or original 27-family multiplicity control.

## Same-source, actual original numerical comparison

All units: difference in NUMBER of positively documented insect visitor guild categories; original 95% equal-source-paper bootstrap intervals. No guild probability, no local flower visitation rate, no selection coefficient.

| True Direct flower morphology | Original matched SOURCE papers / source site×month networks | Observed motif-positive minus negative | Edge-swap residual [95% interval] | **Curveball residual [95% interval]** | Curveball minus edge-swap residual [95% interval] |
|---|---:|---:|---|---|---|
| **Open radial** | **59 / 339** | **−0.17572** | **−0.18457 [−0.28640,−0.08211]** | **−0.18377 [−0.28323,−0.08170]** | **+0.00080 [−0.00579,+0.00692]** |
| Tubular/bell | 43 / 220 | +0.08068 | +0.12633 [+0.02167,+0.22328] | +0.12961 [+0.02650,+0.22538] | +0.00328 [−0.00278,+0.00909] |
| Zygomorphic | 69 / 491 | +0.04438 | +0.05784 [−0.02693,+0.15449] | +0.05435 [−0.03078,+0.15085] | −0.00349 [−0.00825,+0.00202] |
| Deep tube | 15 / 22 | −0.19343 | +0.04358 [−0.12293,+0.19138] | +0.05675 [−0.10675,+0.20353] | +0.01317 [−0.00265,+0.03292] |

**Interpretation:** the open-radial negative source record association survives this stringent **algorithm change**. On the very same 59 original European source papers and 339 exact published site×month networks, source residuals differ by **only +0.00080 guild categories** between algorithms, their paired 95% study-bootstrap interval crosses zero narrowly (−0.00579 to +0.00692), and **42/59 Curveball publication-level residuals are negative** versus **43/59 edge-swap**. Only **one** original source paper changes residual sign between methods. For the Curveball result, leaving out each original publication in turn yields **[−0.20644,−0.16997]**, always negative.

**Other three motifs:** zygomorphic and deep tube remain inconclusive. Tubular/bell is positive under both **marginal, posthoc** bootstrap intervals on 43 papers, but this was previously NOT supported under a **selected-four-motif simultaneous** bootstrap and it does not have a lineage-independent demonstrated effective pollinator mechanism. Do NOT relabel it as confirmatory or island adaptive specialization based on two algorithms giving similar marginal results; they use the same data.

### Original source computational eligibility limits

- Original source public-study positive visitor index comprises **16,757** study × published site × observation year/month networks with relevant validated real insect binomials, 193 literature source keys.
- For the paired algorithms **16,097** site networks failed preflight before comparable morphology/degree-randomization; a further **66** candidate site×motif combinations failed at least one of the computational chain/support gates. The latter count is a motif-level count and may include the same network for multiple profiles; neither is a biological non-visitation/absence count.
- Source comparisons were performed ONLY on those sites qualified for BOTH methods, with 59/43/69/15 original papers for the four distinct morphology panels. A source paper can occur in more than one panel, so those counts are NOT independent repeated evolution events.
- Curveball cross-chain spread and accepted trades provide diagnostics; do not infer perfect uniform sampling or valid group-exchangeability of plant lineages from them. The source data remain selected literature observational records, with unknown original plant community flower abundance, observation effort and visitor detection probabilities.

## Consequences for pollination syndrome and island biogeography claims

**New supported analytical statement:** the **negative association** between open radial flower morphology and positively documented European insect visitor-guild repertoire is robust to which of two degree-preserving bipartite network null generators is used, when comparing identical study/site/month observational networks.

**Still NOT identified:** true source-local visitation preference at equal plant abundances; effective pollen movement, stigma contact, fruit/seed production; any GIFT island's current pollinator-guild repertoire or loss; phylogenetic or source-genus-independent physiological flower specialization; or same-species island-versus-mainland adaptive phenotype change. Most strict within-genus floral comparison is unsupported in DoPI; source original GIFT accepted plants can belong to very different lineages.

**Paper decision:** retain 4,505-island functional beta, source-relative true same-species morphological pairs and exact signed genus accounting as botanical main results. Use this as a bounded supplementary validation and *counterexample to automatic open-flower generalism*. A change of network-null algorithm **is not independent ecological replication**. The null here is a computational conditional reference only, and the two algorithms derive from the SAME original observational data. Keep fd PR #16 draft; do not upgrade the main title to pollinator syndrome shift or island adaptive evolution.
