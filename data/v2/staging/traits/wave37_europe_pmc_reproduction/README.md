# Wave37 Europe PMC reproduction acquisition

Wave37 freezes 500 high-information `genus x reproductive trait` tasks derived
from the formal Wave36 coverage. The priority order is reproduction first, with
support-two tasks followed by high-impact support-one tasks. Search uses the
official, credential-free Europe PMC API; snippets are discovery only and all
candidate evidence is extracted from fetched full-text XML.

Formal lineage anchors:

- integrated all-evidence baseline: Run `32932103226`, artifact ID `9593698966`
- lossless coverage baseline: Run `33137984367`, artifact ID `9672774830`
- fixed denominator: 106,295 species and 318,885 species x axis cells

Tracked inputs:

- `wave37_frozen_top500_queue.csv.gz`: exact task queue, SHA-256
  `db9e7e3149f0a2e5110e9521ed9765870dee5c608c59306aaf91c42cae4b17ad`
- `wave37_candidate_review_decisions.csv`: one frozen decision for every one of
  the 44 candidates in the reviewed pilot acquisition
- `albizia_saman_two_backbone_snapshot.json`: WFO and GBIF exact-synonym
  consensus used to map *Albizia saman* to *Samanea saman*

The checkpoint rejects unknown or missing candidate IDs, verifies the two-
backbone synonym decision, verifies all manually corrected quotations against
the fetched XML, and promotes only 11 reviewed direct species x trait rows.
Validated Low is not created here. The shared trait-specific all-evidence audit
rebuilds rules at `genus x axis x trait_name`, including both leave-one-species-
out and leave-one-source-lineage-out validation. Family inference and global
fallback are excluded.

Run `.github/workflows/build-wave37-europe-pmc-reproduction.yml` to reacquire
the frozen queue, build the candidate audit, rebuild the shared rules, apply the
lossless formal overlay, and publish all raw searches, fetched XML, review
artifacts, and final coverage as one immutable GitHub artifact.
