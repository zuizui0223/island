# Wave42 external bulk reproductive support

Wave42 reuses already acquired reproductive databases before issuing another
broad Web campaign.  Species outside the fixed 106,295-species island universe
are admitted only as species-direct training support for trait-specific
`genus x axis x trait_name` rules.  They never enter confirmatory direct
coverage.

## Frozen identity gate

- 1,671 source names were absent from the Wave41 query set and the fixed target
  universe.
- Exact species names were reconciled against the World Flora Online Plant List
  release 2026-06 (archive SHA-256
  `75f1ad1f371978c9e46f3044152c07ed276fe57be9fb9a15b3621b19cf231987`).
- WFO accepted species and family then had to agree exactly with GBIF species
  matches.  Fuzzy, ambiguous, family-conflicting, backbone-disagreeing, and
  fixed-target mappings fail closed.
- 1,453 GBIF responses were reused from immutable earlier checkpoints and only
  64 new GBIF requests were needed.  Query cost was USD 0.

## Evidence policy

The packet contains 1,390 evidence rows for 1,299 external species and 1,318
external species-trait cells: Ferrer 2024 (650), Meyer 2026 (669),
Razanajatovo 2016 (70), and Goodwillie 2005 (1).  Self-incompatibility,
autonomous selfing capacity, and cleistogamy remain separate traits.

All four secondary databases use one conservative shared source lineage,
`wave42-secondary-reproductive-databases:shared-redistribution-guard`.
This prevents overlapping original studies redistributed by multiple providers
from being counted as independent evidence.  It intentionally sacrifices some
power rather than inflating genus-rule support.

The formal workflow validates this packet, rebuilds all-evidence trait-specific
genus rules together with Wave41 support, and overlays only newly eligible Low
cells onto the formal Wave41 baseline.  Family inference and global fallback are
absent, and all previously filled cells must be retained without downgrade.
