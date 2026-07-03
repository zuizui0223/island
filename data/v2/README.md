# v2 data layout

This directory starts intentionally empty of biological results. v2 is rebuilt from traceable inputs rather than importing v1 derived tables.

```text
templates/
  taxa_template.csv                 accepted taxon batch for web-search work
  accepted_evidence_template.csv    one accepted human-reviewed evidence record per row

staging/                            ignored: raw LLM candidate exports and review queues
curated/                            versioned, human-adjudicated trait evidence and decisions
external/                           versioned input snapshots and acquisition manifests
```

## Minimal workflow

1. Create a species batch from `taxa_template.csv`, using a fixed taxonomy backbone.
2. Run the manual `LLM trait web search` GitHub Action.
3. Download the artifact: it contains raw API responses, cited source URLs, a prompt hash, ontology hash, and candidate CSV.
4. Review candidates; retain accepted and rejected decisions in the curated evidence table.
5. Use `island-v2-hierarchical` to create probability-preserving genus/family candidates only from accepted direct evidence.
6. Analyse direct-conservative, direct-broad-web, expanded-genus, and expanded-family tracks separately.

No action run can silently edit `curated/` or merge values into analysis inputs.
