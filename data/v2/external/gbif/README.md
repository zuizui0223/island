# GBIF v2 retention policy

`request_manifest.csv` is the authoritative state ledger for island-level GBIF downloads.

Each row retains the exact WKT hash, payload, request status, GBIF download key, DOI, and archive link. A failure is recorded as `submit_failed` or with a polling error; it is not converted to an empty flora.

The first raw species table is only a candidate list and must undergo taxonomic and establishment-status review before trait collection.
