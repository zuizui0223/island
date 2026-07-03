# Hierarchical trait coverage in v2

## Purpose

The v2 database is designed for broad global coverage. Species-direct papers are insufficient for many taxa, so curated web descriptions and transparent genus/family inference are part of the design.

The key rule is: **inference increases coverage but does not become an unmarked species observation.**

## Evidence tracks

| Track | Permitted records | Role |
|---|---|---|
| Direct-conservative | species direct/indirect, source reliability A or B | main conservative sensitivity analysis |
| Direct-broad-web | direct species evidence including credible specialist webpages (C) | broad direct-evidence analysis |
| Expanded-genus | direct-broad-web plus genus probability distributions | coverage-oriented sensitivity analysis |
| Expanded-family | expands only when genus support is absent and family thresholds are met | widest, most uncertain sensitivity analysis |

No single track is declared "the truth." The manuscript reports whether a conclusion is stable across them.

## Inference is a distribution, not a forced code

For a missing species trait, v2 does not write a single guessed category whenever there is genus support. It calculates the category distribution among accepted supporting congeners.

Example:

```text
Accepted direct flower colour records in genus G:
white = 3 species
red_pink = 1 species

Candidate for focal missing species:
source scope = genus_inference
P(white) = 0.75
P(red_pink) = 0.25
status = pending human review
```

Island-level models can then be repeated across draws from those distributions. A pattern supported only by a hard modal-category fill but not by distribution-aware draws is not treated as robust.

## Source hierarchy

1. Primary studies, floras, monographs, and taxonomic treatments.
2. Curated databases and institutional web sources: herbaria, botanic gardens, government and university flora portals, museums, and maintained biodiversity databases.
3. Specialist web sources with identifiable taxonomic scope.
4. Weakly curated web sources, retained only as low-confidence leads and requiring review.

A species-level institutional web description is direct evidence about that species. It is not discarded just because it is not a journal article.

## Guardrails by trait domain

- Colour and gross floral form may use genus distributions at relatively modest support thresholds.
- Pollination guild and reproductive traits may use genus/family distributions only at higher support thresholds and remain sensitivity-track values.
- Self-compatibility does not imply autonomous selfing; the latter cannot be inferred from the former.
- Missing occurrence evidence is never used as evidence of Bombus absence.

The thresholds and accepted source grades are versioned in `config/hierarchical_inference.yml`.
