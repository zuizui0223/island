# v2 GBIF acquisition status

The code, configuration, test coverage, and manual GitHub Action are present on `v2-channel-architecture`.

The first global run is blocked only by unavailable external inputs in the repository runtime:

1. source island geometry (`GlbIslands.gdb` / `BigIslands` or a documented replacement);
2. GBIF download credentials configured as `GBIF_USERNAME` and `GBIF_PASSWORD`.

Once those are available, run the GitHub Action in order: `prepare`, `requests`, repeated bounded `submit`, repeated `poll`, then `collect`.
