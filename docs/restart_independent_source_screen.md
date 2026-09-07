# Independent-source acquisition screen, 2026-09-07

Baseline and formal increment remain Run 34093932899 and Run 34096408983:
net +35 species-axis. No rows were admitted by this screen. Reviewer: Codex,
automated source-grounded screening, not independent human biological review.

## Castanopsis and Meyer source inventory

Existing direct SI for C. carlesii, C. eyrei and C. fargesii all points to
the same Meyer Dryad compilation (10.5061/dryad.cc2fqz6hr). The repository's
`Output_Data_MS.csv` has species, family and mating-system category only;
it cannot establish three independent original studies.

The official Dryad API was read successfully:

- https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.cc2fqz6hr
- https://datadryad.org/api/v2/versions/422585
- https://datadryad.org/api/v2/versions/422585/files

This is version 4, published 2026-02-02, NOT a newly discovered dataset release.
Its input workbook is `Input_Data_MS.xlsx`, file ID 4581599, advertised size
684829 bytes and SHA256
`6df359d9f1b52d9548c1b695974aaa981837d35896a4b7a88a9b30ca808c982e`.
Whether it contains original-study citations is not yet established.

The API download requires a bearer token. The ordinary public download
https://datadryad.org/downloads/file_stream/4581599 returned 4309 bytes of
HTML validation content, SHA256
`5eb0cb334e1fda5921bc0d4dfa6104ca2c81ce0551055311206ba78a06ba357d`.
It was renamed locally as an HTML response, not retained as a supposed XLSX.
The normal in-app browser remained at the automatic validation page; no
authentication or protection was bypassed. Workbook retrieval and analysis
remain incomplete. Do not retry the same failed route without a changed state.

Targeted searches for C. carlesii SI and C. eyrei breeding system mostly
returned studies of other species in forests containing Castanopsis. Those
co-occurrence mentions are not evidence for Castanopsis. An additional
Spermacoce query excluding the already reviewed three species did not yield
a usable independent species-level source. Existing duplicate study remains
excluded from new acquisition counts.

## Torenia: distinguish seed dispersal from self-fertilization

The current Torenia direct reproductive records all share
DOI 10.1104/pp.106.083832. A different study was located:
Wu, Qin and Zhao (2008), Plant Growth Regulation 55:137-148,
https://link.springer.com/article/10.1007/s10725-008-9268-5 . The publisher
provides an abstract about pollen-tube hormone responses, but the relevant
self-compatible pollination statement was not available in the publisher
body. Search indexing is not an exact-source receipt; no claim was promoted.

The Australian regulator's original monograph was inspected:
https://www.ogtr.gov.au/sites/default/files/files/2021-07/the_biology_of_torenia.pdf
The Biology of Torenia spp., July 2008, sections 4.1-4.3 (printed pp.11-12).
Its statement “T. fournieri may self seed each year” concerns reseeding and
cites Gilman and Howe (1999); it does not establish autonomous selfing or
self-compatibility. Artificial pollination and seed production likewise do
not identify self versus cross pollen. Genus/family pollination statements,
hybrid sterility and named cultivar traits were not transferred to wild species.
Plant Pono and this monograph must not automatically be counted as independent
experimental sources. The monograph provides leads, not a new Torenia Low rule.

Seven discovery queries were used in this screen (Castanopsis 2, Spermacoce 1,
Torenia 2, exact Wu article title 1, OGTR monograph 1). Search billing is not
exposed; cost is unknown, not zero. No accepted claims, new cells, or Low
invalidations were produced. The source-deficit queue remains active; this
specific Dryad download route is paused pending access, while other sources
and the already validated-rule candidates remain available for further work.
