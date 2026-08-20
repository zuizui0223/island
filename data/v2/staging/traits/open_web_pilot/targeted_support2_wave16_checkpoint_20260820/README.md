# Targeted support-two Wave 16 checkpoint

Wave 16 adds 27 individually reviewed species-direct trait rows to the shared
all-evidence input. It does not implement genus inference. Validated Low is
rebuilt only by the common trait-specific `genus x axis x trait_name` path.

The checkpoint contains:

- 25 rows pinned to three completed acquisition artifacts, including audited
  PlantNET NSW Flora morphology and two NC State Extension rows;
- `Rhamnus alaternus` mating-system evidence from the University of Porto
  Natural History and Science Museum species page; and
- `Viscum coloratum` mating-system evidence from the peer-reviewed species
  review archived by PubMed Central.

Every row stores an exact excerpt, stable URL, retrieval time, source lineage,
content fingerprint, and an individual accept decision. The committed source
snapshot records the originating run and artifact IDs, so rebuilding does not
repeat Web acquisition.

Fail-closed exclusions include seven descriptions of a solitary **spikelet**
(not a solitary flower display), coarse NC State `<1 inch` size categories,
and any cross-trait substitution. `pollen_vector_mode` and `reward_type` remain
independent traits and are absent from the strict three-axis inputs here.

Regenerate the deterministic files with:

```powershell
$env:PYTHONPATH = "src"
py -3 -m island_v2.targeted_support2_wave16_checkpoint
```

Formal direct and Validated Low gains are measured only by the review/promote
workflow artifact; support-two queue ceilings are not reported as coverage.
