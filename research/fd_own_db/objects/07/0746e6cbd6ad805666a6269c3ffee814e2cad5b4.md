# Are island flowers experiencing real pollinator shifts, or only a change in documented guild opportunities?

**fd PR #14 · 12 Nov 2026 · original-data evidence audit, not evolutionary shift certification.**

## Two different evidence problems, now resolved separately

**A. Global source-documentation fragility.** A SHA256-verified 2026-06-04 GloBI v0.9 raw original run [island Actions #43034251795](https://github.com/zuizui0223/island/actions/runs/43034251795) succeeded. It retested the **SAME 344 botanical genera / 77 unambiguous families** before and after excluding the dominant `sourceArchiveURI`-first group (95.9584% of citation-backed genus×visitor-channel reference units). Both predictions remove each tested family from model training and refit the GloBI-reference-effort-only baseline. All results concern **global DOCUMENTED references**, not pollinator occurrence or effective flower service.

| Documented GloBI reference guild | Family-equal morphology MSE gain before | Family-equal gain after | Paired after-minus-before gain (4,999 family-bootstrap CI) |
| --- | ---: | ---: | ---: |
| Bombus | −0.000629 | −0.000974 | −0.000345 [−0.002078,+0.001148] |
| Non-Bombus bees | **+0.004348** | +0.001087 | **−0.003262 [−0.006630,−0.000129]** |
| Lepidoptera (butterflies and moths combined) | +0.001609 | −0.000069 | **−0.001679 [−0.003375,−0.000039]** |
| Flower-visiting birds | **+0.002243** | **+0.002621** | +0.000378 [−0.000508,+0.001350] |
| Diptera | +0.001274 | +0.000771 | −0.000503 [−0.002169,+0.001177] |
| **Five-channel mean** | **+0.001769** | +0.000687 | **−0.001082 [−0.002177,−0.000107]** |

The source omission **attenuates the pooled flower-open-form incremental predictive gain** on matched taxa/lineages. Bird documentation remains positive *after omission*, with post-omission family-equal 95% interval **[+0.000879,+0.004586]**. The paired bird-specific change is NOT distinguishable from zero in its conditional interval. Non-Bombus bee and Lepidoptera documentation sensitivity are negative, but multiple unregistered guild outcomes and dependency of the five shares prohibit independent evolutionary significance claims.

The dominant `sourceArchiveURI` label might contain many studies; **archive-label exclusion is NOT yet confirmed leave-one-independent-provider-study-out**. The 77 families can share reference providers and phylogenetic deep ancestry. A GloBI row classified as `pollination` also need not measure pollen deposition or seed set.

**B. Real island pollinator channel opportunity.** Independently join the **928 islands** supporting >=30 fully measured species-direct floral colour/form/symmetry profiles to the original five-channel GBIF exact-island observation/source qualification (artifact 10314984302). This dataset's statuses are `retained` = documented GBIF animal occurrence, `disrupted` = source region available but island adequately nondetected, and `unresolved` = no justified inference. Neither positive animal occurrence nor adequate nondetection establishes species-specific plant–pollinator interaction or historical replacement.

- **22/928** morphology-complete islands have a qualified Bombus nondetection. At those islands, **22** also have non-Bombus bees, **22** Lepidoptera, **22** Diptera and **21** flower-visiting bird occurrence records.
- **The SAME 21 Bombus-undetected / bird-recorded islands also have all three other alternative guilds recorded.** Thus the data cannot distinguish *bird takeover* from broadly multi-channel opportunities.
- The 21 candidate bird/Bombus co-recording condition spans **14 geographic spatial blocks**: **15 tropical islands in 8 blocks**, **5 northern midlatitude islands in 5 blocks**, **one southern extratropical island in 1 block**. No northern-high-latitude island satisfies that particular condition.
- Broad regional comparison fails overlap/positivity: tropical Bombus `disrupted=16` but `retained=2`; northern midlatitude `disrupted=5` vs `retained=480`. The strict gate requiring >=10 islands and >=10 blocks of EACH state within a region–guild yields **zero qualified cells** for the five-guild dataset. The original N1 channel isolation heterogeneity gate was unsupported.
- **Zero certified plant-specific bird-for-Bombus evolutionary pollinator shifts can be extracted from the available *table schema***, because there are **no plant-species × visiting-animal links**, mainland/island paired plant visitation observations, effective pollen-transfer measurements or ancestral reconstruction records in the five animal channel occurrence status columns. This is **lack of the required data**, NOT a discovery that actual islands have no pollinator shifts.

A completely separate island-level original morphology result remains valid: marginal trait Rao Q and colour×form/symmetry coordination are distinct; **928** islands have 30+ Direct trait-complete species, and **909/928** show positive colour–architecture coupling relative to island trait-marginal label reshuffles. This morphological excess is highly taxonomically conditioned and does NOT itself measure realized pollinator guild matching.

## How to test competing mechanisms

Three distinct non-nested hypotheses must be recognized:

1. **Historical lineage pollinator shift:** the SAME plant taxon (or phylogenetically supported sister lineage) uses a sufficiently documented mainland effective pollinator guild A, but on an island its **effective pollen delivery** is dominated by guild B, with A absent or unimportant under matched sampling. This is an actual ecological/evolutionary shift; it needs plant-linked and time/space-indexed visitor and pollen-transfer records.
2. **Community-level compensatory generalization:** when one potential visitor guild is poorly represented, *multiple* others visit/pollinate the same island plant taxa. A phenotype may remain morphologically conservative while realized interactions change. The 21 islands with four alternative animal-guild occurrence records are only **candidates to survey** for this hypothesis, not proof of same-plant sharing.
3. **Flora assembly/lineage sorting:** plant families and genera with preexisting open flowers colonize/are recorded in certain island environments. Genus-conditioned morphology and GloBI global source documentation can covary WITHOUT any within-lineage trait evolution or ecological pollinator switch.

**Predeclared measurement gate for a positive island shift claim:** retain plant species/genus ID and indigenous/origin status; identify actual mainland or plausible source plant visitor interaction partners and an independently sampled island visitor network with flower visits quantified and effort-standardized; require taxon A and B evidence and appropriate nondetection sampling; compare **source and island with the SAME lineage**; use effective pollen deposition or seed success to interpret actual pollination. Include bird, flies and bees, split **moths from butterflies** if primary taxonomic records support it, and allow potential reptile/lizard pollination without imputing missing records.

## Manuscript interpretation

**Supported:** Floral morphological disparity can coexist with nonrandom trait coordination. Coarse open-radial flower frequency contains a small source-dependent predictive association with GloBI global documentation shares; bird references are positive after deleting the largest archive grouping, while bee/Lepidoptera gains attenuate. Several island floras have Bombus adequately nondetected while multiple other animal groups are recorded.

**Unproven:** birds actually replace bees on those islands, floral morphology evolved in response to replacement, generalist functional compensation maintains seed set, moth-vs-butterfly switching, or island-level true pollinator syndrome alignment. Interpret shifts as **testable channel-gated mechanism candidates**, not confirmed ecological results.

Relevant independent published biological anchors: [Canary Islands bee-to-bird colour transition, functional explanation (2024)](https://doi.org/10.1111/1365-2435.14493); [oceanic-island pollinator syndrome mismatch, Wang et al. (2020)](https://doi.org/10.1038/s41598-020-70954-7); [island plant mutualism filter, 2024](https://doi.org/10.1038/s41586-024-07110-y); [Ogasawara co-invasive networks, 2026](https://pubmed.ncbi.nlm.nih.gov/42032836/). They are **independent context**, not additional fd source-specific empirical confirmations.

## Rights and replication

- `src/fd/pollinator_shift_opportunity.py` has fail-closed five-channel input support and outputs only anonymous pair-of-channel candidate counts and the explicit **zero certified species-level shift available in input** guard.
- `tests/test_pollinator_shift_opportunity.py` synthetic missing/unresolved/duplicate and false-bird-replacement prevention tests.
- `paired_provider_omission_77families_fiveguilds.csv`: six original source-pinned, aggregate-only paired family gain contrasts from completed island workflow #43034251795 and artifact #13864221183.
- `bombus_nondetection_multiguild_candidate_support.csv`: four non-identifying island-level candidate guild co-recording conditions.
- Original precise island IDs, species names, original GloBI provider/reference identifiers and spatial-block membership are NOT released. Both GitHub PRs stay Draft pending a true plant×visitor×pollen-transfer bridge.
