# Restart: affected genus-trait frontier

2026-09-07. This is a non-promoting diagnostic, not an acquisition result.

Formal baseline: Run 34093932899 / artifact 10007872352. Latest integrated
Run 34096408983 / artifact 10008766395: net +35 species-axis (46 new direct
species-trait), 222,410 filled and 96,475 unresolved. The +1,000 target remains
active; 965 additional net cells are still required. Historical Low is not
revalidated by this diagnostic.

The 46 added source records affect 37 genus-trait pairs. The common
`all_evidence_trait_audit.py` implementation at commit
ff8249c6de682eaa70dc5a953c5ed420d85e1f98 was used, without a genus-axis join
or a new local inference algorithm. Expected baseline source-lineage keys
for these groups all have raw receipts (zero missing keys). Receipt presence
alone does not prove independent upstream studies or resolve conflicts.

| Setting | Eligible rules | Newly eligible | Newly ineligible | Unresolved axis candidates |
| --- | ---: | ---: | ---: | ---: |
| Current min-species 3 | 5 | 5 | 0 | 21 |
| Mild relaxation, min-species 3 | 5 | 4 | 0 | 21 |

Current eligible genus-trait pairs: Aegilops/mating_system;
Anacyclus/inflorescence_display; Anemonastrum/inflorescence_display;
Briza/inflorescence_display; Capsella/floral_symmetry. These are NOT adopted
Low rules. Added records retain genus-training=false; the diagnostic deliberately
models their possible use, without authorizing promotion. Counts are restricted
to affected groups and are not a full-genus sensitivity report. Newly ineligible
rules=0 is not a claim of zero invalidations throughout the historical Low ledger.

## Decision

Relaxing thresholds produces no additional candidates in this batch. Do not
relax them merely to chase the +1,000 target. Audit upstream citations for the
21 candidates, but move the main acquisition effort toward unresolved
reproductive genus-trait groups and original multi-species tables. Check source
DOIs against the existing source inventory before fetching the same studies.
The support-two queue is a search aid, not a promised number of unlocked cells.

## Reproduction

Run `scripts/diagnose_restart_genus_frontier.py` with `--common-root` pointing
to the pinned common repository, `--recovered-raw` to the recovered raw source
inventory, `--baseline` to the recovered baseline artifact, `--integrated` to
the latest integrated artifact, and a new `--output` directory. Outputs include
before/after rules, unpromoted frontier, missing-lineage receipt keys, and a
summary with all input SHA256 values and formal source run/artifact IDs.
The raw reconstruction is locally available but is not yet packaged as a
standalone formal artifact; this diagnostic therefore remains local and must
not be described as a formally reproduced new Low artifact.

## Further source screen

Three discovery queries were used for Malpighia breeding systems, Portulaca
comparative reproduction, and the twelve-species Malpighiaceae study. Search
billing is unavailable, not assumed zero. No snippet was admitted as evidence.

- Wigandia urens, DOI 10.1007/s11258-026-01671-w: original publisher page was
  accessible. Current ledger already has autonomous selfing and mixed mating
  from Europe PMC preprint PPR1257894. Treat publication/preprint as a potential
  lineage alias, not new cells or independent support.
- Portulaca amilis, PMC9157154: reproductive axis is unresolved, but direct PMC
  retrieval returned a browser challenge. Search content is a lead only; no
  new evidence was promoted.
- PMC11805948: original retrieval returned a browser challenge. Repository
  Wave42 already lists supplementary dataset `mcae056_suppl_supplementary_data.xlsx`
  and DOI 10.1093/aob/mcae056. Resolve article identity and prior source coverage
  before downloading or claiming a new bulk dataset.

No additional direct or inferred cells were promoted by this diagnostic/screen.
