# First matched animal-guild/public floral state-set association — 2026-11-24
[Actions 45363899725](https://github.com/zuizui0223/fd/actions/runs/45363899725) SUCCESS; artifact 14906555345 SHA256 f327f2b892757c93da9b0440e8df34988099286dd80f3b7dfbde645f1cf5395c. Fast pair-specific workflow reads only verified source bridge/bird/butterfly artifacts and pinned public floral carriers, no repeated original bee archive download. The broader bee-only 130-island fit has NOT been completed in this run; separate full archive retrieval 45363436156 timed out. Data: PR321 historical colour/form taxon-name state sets, not original 4505 radial x deep phenotypes. All 40 exploratory panels and all eight outcomes are archived, no adjusted significance claimed. Within-archipelago group-equal conditional OLS, log1p animal species counts, log2 area, mainland distance and floral recorded names; bee count held constant. Butterfly positives require >=10 original records, no inferred zero species. Bird match provisional.

## Bird endpoints
| Scope / minimum floral names | Candidate islands | Informative groups | Specialist >=50% tubular membership coefficient | Extended-access shape coefficient | Extended-access group resampling 2.5–97.5% |
|---|---:|---:|---:|---:|---|
| Direct / 5 | 32 | 6 | +.064516 | +.106748 | [-.002705,.692284] |
| Direct / 50 | 24 | 5 | +.045312 | +.060550 | [-.025115,.079488] |
| Direct single-state / 5 | 29 | 5 | +.037713 | +.102051 | [.046375,.148535] |
| Direct single-state / 50 | 18 | 4 | +.008833 | +.066053 | [.034168,.085938] |

Extended shape is recorded membership in tubular/salverform/funnel_trumpet/spurred; not directly measured depth, accessibility or syndromes. Tubular coefficient corresponds to state-set membership, not exclusive single tubular flowers except in single-state scope. Both direct multistate specialist point estimates stay positive on group omission. However, the group numbers are sparse and percentile bootstrap ranges are not calibrated CIs. Broad nectar >0 extended shape coefficients (same four nested panels) +.001728,+.009305,+.016801,-.003270, after including nonnectar bird count -.047033,-.013428,-.010022,+.003190. Red/pink x extended membership lacks consistent positive support and can be negative. Specialist coefficient has NOT been tested in a direct paired model with nonnectar/partial feeders; do not interpret as evidence of guild-specific compensation or morphological specialization. Next explicit specialist-vs-partial feeder contrast on same island; original protected joint floral endpoints are more relevant.

## Butterfly Q and source sensitivity
| Direct min5, source | Islands | Raw Rao Q coefficient | Record-count-adjusted |
|---|---:|---:|---:|
| all | 54 | +.013974 | +.005055 |
| uncertainty <=1km | 40 | +.002370 | -.000708 |
| no GBIF ID with reference | 37 | +.011206 | +.022416 |
All/54 raw resampling [.000209,.029677]; non-GBIF-ID/ref 37 raw [.001520,.035419], adjusted [.006207,.066535]. Those leave-one-group coefficients positive. BUT Direct n50 adjusted coefficients for all/precision/ref n40/34/26 are -.002604,+.002432,-.008865; latter range [-.022471,-.008572] and leave-one group negative. Single-state n50 adjusted all/precision/ref n23/23/19 coefficients +.020771,+.026842,-.040736, latter resampling [-.042827,-.018138]. Therefore no robust Q improvement, no evidence of functional compensation. Butterfly shape effects also unstable by mask and threshold. Preserve all outputs rather than selected positive panels; records not independent true butterfly richness estimates.

## Original protected endpoint
Public-only animal geometry/record bridge independently exported in [45363392892](https://github.com/zuizui0223/fd/actions/runs/45363392892), JSON SHA256 467f111ca14be773bc7dbe6cd76880c443d5f8d6a81c0dad55ca6b8799302776. Protected [island draft PR #346](https://github.com/zuizui0223/island/pull/346) created to test same-species actual radial x deep, bilateral x deep and joint category Q using original 4,505 parent, original native/high-quality and genus x quality nulls; [run 45363746484](https://github.com/zuizui0223/island/actions/runs/45363746484) started. No original-data result claimed until execution succeeds. Only anonymous summaries may leave protected runner.
