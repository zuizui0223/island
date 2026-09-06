# Chapter 1 canonical trait-coverage contract

## One canonical coverage number

Chapter 1 uses one scientific trait-coverage denominator:

- analysis-applicable angiosperms: **106,295 species**;
- strict analysis axes: **3** (`flower_colour`, `floral_structural_complexity`, `reproductive_assurance`);
- fixed species × axis denominator: **318,885 cells**.

The current working strict checkpoint after private TRY adjudication and Wave55 Batch 5 is:

- resolved strict species × axis cells: **222,759 / 318,885**;
- strict coverage: **69.8556%**;
- unresolved strict cells: **96,126**.

This is the only percentage that should be described as **Chapter 1 trait coverage** during the continuing acquisition campaign.

The private TRY request is not redistributed through this repository. The public repository retains the adapter, provenance/conflict contract, and reproducible public evidence layers; the working checkpoint above therefore records the current private-plus-public analysis state rather than claiming that the private TRY rows are publicly reproducible.

## What is not Chapter 1 coverage

The `all-master` acquisition inventory uses a different population and a different purpose:

- 115,328 master names, including 9,033 taxa outside the strict angiosperm analysis universe;
- many individual trait fields rather than the three strict analysis axes;
- candidate/source reach upstream of strict adjudication;
- optional machine-candidate and inference inventory channels that do not automatically become strict accepted evidence.

Percentages such as `any-trait broad`, five-field completion, or other 115,328-species all-master rates are therefore **secondary acquisition-inventory metrics**. They must never be compared numerically with 69.8556% or presented as the current Chapter 1 coverage.

Use the all-master inventory only to answer operational questions such as:

- which sources still contain potentially useful records;
- which species/traits have no configured source candidate;
- which provider contributes unique candidate reach;
- which records need taxonomic or ontology reconciliation.

## Promotion rule

A newly acquired record changes the canonical Chapter 1 percentage only after it passes the strict evidence contract and resolves a previously unresolved `accepted_species × analysis_axis` cell.

Machine-only candidates, global/family fallback fills, unresolved direct conflicts, out-of-scope taxa, and candidate presence without strict adjudication do not increase canonical coverage.

## Reporting rule

Every acquisition checkpoint should report, in this order:

1. strict resolved cells / 318,885 and percentage;
2. strict resolved counts by the three axes;
3. net newly resolved strict cells and any losses/downgrades;
4. remaining unresolved cells by axis and recoverability class;
5. secondary source-inventory metrics only in a clearly labelled appendix/diagnostic section.

This separation prevents broad source-discovery reach from being mistaken for analysis-ready evidence.