# Europe PMC exact article-level metadata audit

2026-11-26 JST. The historical source family has 2,024 unique, parseable article URLs: MED 1,302, AGR 537, PPR 128, CBA 25, ETH 21, PAT 10, PMC 1. These are historical identifiers, not necessarily 2,024 independent studies. They include literature databases, preprints, theses and patents.

The frozen rights audit marks 18,826 cells with this family as a blocker (15,295 reproductive); only 883 cells have it as the sole blocker (781 reproductive). Of these affected cells 8,431 have Validated-Low upstream lineage recovery. The exact lineage crosswalk uses all stored direct tokens plus the recovered upstream tokens; every affected cell has at least one exact Europe PMC identifier, none has an unmapped Europe PMC token. The 883 are the maximum if every required Europe PMC lineage is cleared; this is not a permission forecast.

## Protocol

Use [Europe PMC official REST service](https://dev.europepmc.org/RestfulWebService), exact SRC and EXT_ID (PMCID for PMC) batches, strict source/ID checks and allowlisted license metadata. The [copyright policy](https://dev.europepmc.org/Copyright) distinguishes accessibility from reuse permission and requires source-specific terms. Do not bulk scrape site pages. The former NCBI OA web service was retired in August 2026; it is not used ([NCBI announcement](https://pmc.ncbi.nlm.nih.gov/tools/oa-service/)). No full text or excerpts retained.

A current exact article metadata record with OA, DOI, PMCID and explicit CC BY, CC0 or CC BY-SA label in MED/PMC is a candidate for item-level license/edition verification, not an authorized trait row. Noncommercial, no-derivatives, missing, conflicting, unavailable and ambiguous version/source records remain separate. License versions and third-party source restrictions still need verification. Evaluate cell-level ceilings only when every frozen direct/derived Europe PMC token is in the candidate set and no other source family blocks it. Never weaken the frozen rights gate. The runner publishes only aggregate counts; optional private metadata output can be generated in authorized storage. No publication, editor email, source policy change or new scientific fit.

Local pinned-archive roster audit completed. The current environment's direct API request could not connect; the workflow performs the official metadata requests on its runner. Numerical license results await successful execution.
