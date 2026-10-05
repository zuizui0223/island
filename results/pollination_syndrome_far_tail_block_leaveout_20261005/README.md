# Pollination-syndrome far-tail spatial-block leaveout — 2026-10-05

Status: **post-hoc diagnostic; not part of the frozen Chapter 1 submission inference.**

> **Superseded biological interpretation (2026-10-05; PR #313).**  
> The coefficients below are retained as a historical post-hoc diagnostic trail. Later
> species/status attribution showed that the northern mid-latitude positive far tail
> disappears in confirmed native flora. In tropical native flora, the positive tail
> disappears after removing literature-qualified ocean-dispersed coastal species
> (native-nonendemic **+0.9844 → -0.2330**; **-0.2294** after recomputing
> `selfing_core`), whereas tube-depth coding flags alone do not remove it. These
> outputs therefore should **not** be interpreted as a pollination-syndrome mechanism
> or as a shared very-remote-island biological response across regions. Current
> follow-up: `results/tropical_far_tail_dispersal_20261005/README.md`.


This diagnostic asks whether the newly identified remote-island positive
`yellow/orange × butterfly-associated deep-tube` response is driven by one
archipelago/spatial block.

## Design

The parent model is unchanged:

- corrected mainland distance;
- frozen raw `architecture | colour` island counts;
- selfing-core, island area and climate PC1–4 adjustment;
- fixed hinge at **783.017 km**;
- spatial-block cluster-robust uncertainty.

For northern mid-latitude and tropical floras separately, and for both All and
Direct evidence, every spatial block containing at least one >783 km observation
was omitted in turn and the hinge model was refit.

## Result

The post-hinge slope remained **positive in every single leave-one-block fit**:

- northern mid-latitude All: **17/17**;
- northern mid-latitude Direct: **17/17**;
- tropical All: **54/54**;
- tropical Direct: **53/53**.

### Northern mid-latitudes

All-evidence post-hinge slope across block deletions ranged from **+0.265 to
+0.600**. Direct ranged from **+0.348 to +0.550**.

Nominal P < 0.05 was retained in 13/17 All and 14/17 Direct deletions. The
important robustness result is sign and magnitude stability, not significance
counting.

### Tropics

All-evidence post-hinge slopes ranged from **+0.122 to +0.271** and remained
positive in all 54 deletions. The All layer remains imprecise, as in the parent
fit.

Direct-evidence slopes ranged from **+0.186 to +0.266** and remained nominally
supported in **53/53** leaveouts. The hinge model remained AIC-preferred in every
Direct tropical deletion.

## Interpretation

The remote-island positive deep-tube pattern is **not a one-archipelago artefact**
in either northern mid-latitudes or the tropics.

Combined with the direct north-mid versus tropical far-tail contrast, the most
defensible reading is now:

> **The yellow/orange deep-tube response is a distributed very-remote-island
> pattern shared across at least northern mid-latitude and tropical floras,
> rather than a tropical-specific response driven by one island group.**

This does not mean butterflies are the realized visitors. The response remains a
phenotype–syndrome concordance pattern.

## Submission boundary

Do not add these 2026-10-05 coefficients or leaveout counts to the current
Chapter 1 manuscript or SI. Their role is to prevent over-interpretation of the
full-range tropical association as a unique tropical pollinator mechanism.

## Reproducibility

- workflow run: **37295699203**
- artifact: **11337988422**
- script:
  `scripts/audit_pollination_syndrome_far_tail_block_leaveout.py`
- machine-readable summary:
  `results/pollination_syndrome_far_tail_block_leaveout_20261005/far_tail_block_leaveout_summary.csv`
- full leave-one-block rows are retained in the workflow artifact and regenerated
  by the workflow.
