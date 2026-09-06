# Wave39 FloraWeb/BIOLFLOR synonym recovery

Wave39 recovers FloraWeb/BIOLFLOR records whose historical scientific name
did not exactly match the fixed 106,295-species universe.

The frozen WFO 2026-06 Darwin Core release is used only to prefilter possible
routes.  Promotion then requires independent agreement between the live WFO
exact matcher and the GBIF strict species matcher, species rank, Plantae,
family agreement, and membership in the fixed universe.  Fuzzy, ambiguous,
single-backbone, family-conflicting, and backbone-disagreeing matches remain
in the audit and are not promoted.

The directory stores the complete prefiler, all 446 WFO and 446 GBIF response
payloads, the strict mapping audit, extraction and provenance audits, direct
three-axis evidence, and independent-trait evidence.  The independent ledger
keeps pollen vector, reward, dichogamy, and sex system outside the strict
three-axis coverage.

The formal secondary overlay is based on:

- Wave33 Run `32932103226` / artifact ID `9593698966`;
- Wave37 Run `33143109604` / artifact ID `9674790380`;
- Wave38 Run `33149235490` / artifact ID `9677027199`.

New direct evidence and newly eligible `genus x trait_name` rules are added
losslessly.  Counterexamples that would invalidate already-published Low rules
are retained as a shadow invalidation audit; they are not silently deleted
from the append-only formal secondary ledger.  The scientifically rebased
diagnostic remains separately reproducible for future replacement or direct
repair acquisition.
