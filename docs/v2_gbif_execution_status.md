# v2 GBIF acquisition status

The frozen island universe (GSHHG-derived, 8,265 islands >= 5 km2) is packed
into **103 regional GBIF request blocks**, tracked in
`config/gbif_full_acquisition_v2.json`. GBIF credentials are configured as the
`gbif-downloads` environment secrets `GBIF_USERNAME` / `GBIF_PASSWORD`.

## Self-driving campaign

Two scheduled workflows drain the 103 blocks unattended:

- **`submit_gbif_frozen_full_acquisition`**  -  every 15 min. Regenerates the
  frozen block manifest, selects pending blocks not yet recorded with a real
  download key, and submits up to `batch_size` (default 3). The submitter is
  capacity-aware: it stops the moment GBIF reports its 3-simultaneous-download
  limit, records the overflow as transient `deferred_capacity`, and those
  blocks are re-selected on the next tick.
- **`poll_gbif_full_acquisition`**  -  every 15 min (offset). Advances each active
  download's status (`submitted` -> `running` -> `succeeded` / `failed`) and
  records the DOI.

Both share the `gbif-full-acquisition` concurrency group (so they never run
concurrently against the shared ledger) and rebase onto `main` before pushing.

**Throughput is set by GBIF, not the cron interval:** only 3 downloads can be
active at once, so the campaign advances ~3 blocks per GBIF completion cycle.
With ~95 blocks left this is on the order of a day of wall-clock, depending on
per-download processing time and GitHub's schedule punctuality.

To pause: disable the `submit_gbif_frozen_full_acquisition` workflow in the
Actions tab (in-flight downloads keep completing and can still be polled).

## Remaining manual / not-yet-automated step

`collect`  -  once blocks reach `succeeded`, download the SIMPLE_CSV archives and
assign returned coordinates back to the **original exact island polygons** (not
the buffered query catchments). This is the next module to build; it is not part
of the submit/poll loop.
