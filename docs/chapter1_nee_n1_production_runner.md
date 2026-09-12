# Chapter 1 NEE one-shot N1 production execution

## Purpose

This workflow is the final executable surface of N1. It runs only after the separate
channel-qualification workflow has completed and produced the frozen projected channel
states plus qualification receipt.

The trigger file contains one explicit upstream run ID:

```text
qualification_run_id=<run id>
```

Without that trigger, the workflow is a no-op.

## Inputs

The workflow requires exactly one artifact named
`chapter1-nee-channel-qualification-<run id>` and reads only:

- `nee_projected_channel_state.csv`;
- `nee_channel_qualification_receipt.csv`.

It independently downloads the canonical Chapter 1 geography/covariate artifact:

- artifact ID `8270544465`;
- ZIP SHA-256 `695f35b97bae07e81b05deab537dd73fa687b2d99e9efdb1cb3babd2fa12dfb6`;
- `results/purpose_shortest_island_data.csv`.

No pollinator occurrence collection, source-state construction, floral trait data, plant
assemblage outcome, or N2 evidence is read here.

## Execution

The merged `chapter1_nee_n1_run` wrapper is called exactly once. If fewer than three
channels are confirmatory, it writes `N1_not_evaluable_for_NEE_gate` and does not fit a
rescue model. Otherwise it runs the frozen common and heterogeneous retention models,
the global channel-by-isolation Wald test, leave-one-spatial-block-out robustness and
leave-one-raw-source-region-out robustness.

## Stop/pass boundary

A failed or non-evaluable N1 must write:

`stop_before_N2_and_keep_frozen_Chapter1`

A passing N1 may write:

`N2_may_open`

but this workflow **does not launch N2 automatically**. The N1 artifact must be inspected
and locked before any N2 execution begins.

## Output

The workflow uploads `chapter1-nee-n1-<run id>` containing the complete frozen N1 output
surface: model support, coefficients, channel slopes, global test, model comparison,
both deletion robustness tables, run metadata and the gate receipt.
