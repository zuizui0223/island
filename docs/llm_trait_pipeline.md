# LLM trait acquisition pipeline (v2)

## Why v2 does not require papers only

Global plant trait coverage is too uneven for a literature-only approach. Species descriptions in institutional and curated web sources often provide the only obtainable colour, gross floral-form, or phenology information for a large part of the island flora.

v2 therefore uses **broad retrieval with provenance**, rather than restricting the database to journal sources or pretending that missing traits are biologically absent.

## Four evidence layers

```text
Layer 1  species direct, formal sources
         primary studies, floras, monographs, taxonomic treatments

Layer 2  species direct, curated web sources
         herbarium, botanic-garden, government, university, museum,
         curated biodiversity portal, or trait database

Layer 3  species direct, specialist web sources
         identifiable taxonomic scope but weaker editorial guarantees

Layer 4  declared hierarchical inference
         genus or family distribution built from reviewed Layers 1–3 records
```

All layers are useful. They are never silently pooled as if equally certain.

## Why genus inference remains necessary

A narrow direct-only database would exclude too much of the global flora and can itself create taxonomic and geographic bias. Genus inference is retained under two conditions:

1. it names the supporting congeneric taxa and their accepted evidence records; and
2. it remains a probability distribution rather than an unmarked hard trait code.

For example, a focal species without direct colour evidence may receive:

```text
scope: genus_inference
support: 8 accepted congeners
P(white): 0.625
P(red_pink): 0.250
P(blue_purple): 0.125
review status: pending
```

The analysis samples or integrates over that distribution in expanded-coverage sensitivity analyses.

## The controlled web-search action

`src/island_v2/trait_extraction.py` is a manual batch runner. It uses the Responses API web-search tool, requests structured JSON, and writes:

- raw API responses;
- URLs returned by the web-search response;
- source-cited trait candidate rows;
- prompt and ontology hashes; and
- a run manifest.

It **does not** update curated evidence or analysis data. Review is a separate explicit step.

## Review and analysis tracks

The reviewer stores accept/reject decisions separately. After review, `src/island_v2/hierarchical.py` can derive genus/family distribution candidates from accepted direct evidence.

The resulting analyses are always reported as parallel tracks:

```text
1. direct-conservative
2. direct-broad-web
3. expanded-genus-distribution
4. expanded-family-distribution
```

The v2 conclusion is strongest when the sign and regional pattern persist across these tracks. A result that appears only after broad hierarchical inference remains informative, but is reported as coverage-sensitive rather than definitive.
