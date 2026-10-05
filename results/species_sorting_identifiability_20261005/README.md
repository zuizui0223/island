# Chapter 1 source-matched assembly decomposition — 2026-10-05

Status: post-baseline mechanistic diagnostic; does not replace the final traitwise H1.

## Question

Can the source-evaluable part of the current seven-trait H1 be separated into
pre-existing source-pool composition, genus-level assemblage structure, species sorting
within represented genera, and a residual within-species evolutionary term?

The current data identify the first three levels. They do **not** contain
population-specific trait states, so a within-species term is missing rather than zero.

## Exact decomposition

For each island × outcome × source-mode row:

```
raw H1 mean
= source-species expectation
+ genus/species-structure enrichment
+ within-genus species-sorting enrichment
```

The source-species expectation preserves source-species prevalence across the assigned
mainland source entities. The genus/species expectation additionally preserves the genus
of each observed source-evaluable species slot. Thus the second increment is the
contribution of which genera structure the island assemblage, whereas the final increment
asks whether species within those represented genera are additionally sorted by H1 state.

The row identity closes to machine precision (maximum absolute error
1.11e-16 in the validated refined run).

## Result

### Northern midlatitudes

The raw source-evaluable H1 response is positive in all four frozen source modes.

All-analysis:
- raw H1 mean isolation slope: +0.01862 to +0.01992; all four source-mode vectors FDR-supported;
- genus/species-structure enrichment: **+0.01627 to +0.02214**, positive and FDR-supported in all four modes (maximum q = 1.49e-5);
- within-genus species sorting: **-0.00105 to +0.00161**, not four-mode FDR-supported (maximum q = 0.136).

Direct-only:
- raw H1 mean: +0.02247 to +0.02315;
- genus/species-structure enrichment: **+0.01842 to +0.02333**, positive and FDR-supported in all four modes (maximum q = 5.54e-6);
- within-genus species sorting: **-0.00095 to +0.00232**, not supported (maximum q = 0.130).

The northern-midlatitude source-matched signal is therefore localized primarily to
**which genera compose the island assemblage**, not to additional selection among
species within already represented genera.

Atomic effects are heterogeneous. In the all-analysis scope, plain colour,
selfing mating system and shallow/open tubes make the clearest positive
genus-structure contributions; self-compatibility is also positive in all four source
modes. No atomic result is promoted on its own from this post-baseline decomposition.

### Tropics

The source pool and realized assemblage move in opposite directions.

All-analysis:
- source-species expectation: **+0.01853 to +0.01949**, positive and FDR-supported in all four modes;
- genus/species-structure enrichment: **-0.02018 to -0.01474**, negative and FDR-supported in all four modes (maximum q = 0.00722);
- within-genus species sorting: **+0.000004 to +0.000404**, not supported (maximum q = 0.703);
- total species sorting relative to the source expectation: -0.02000 to -0.01455.

Direct-only:
- genus/species-structure enrichment: **-0.01836 to -0.01466**, negative and FDR-supported in all four modes (maximum q = 0.0498);
- within-genus species sorting: -0.00101 to +0.00051, not supported.

Thus increasingly isolated tropical source pools tend to contain more H1-oriented trait
states, but realized genus composition shifts in the opposite direction and nearly
cancels that expectation. The antagonistic response is again localized above the
within-genus species level.

### Unsupported contexts

Northern high latitudes and southern extratropics do not meet the frozen
minimum-50-islands-per-outcome source-evaluable gate. They remain not testable; the
threshold is not relaxed.

## What this resolves

The previous statement that contemporary island-flora data cannot distinguish species
sorting from within-lineage change was too coarse for the source-evaluable subset.

The present result supports a narrower and stronger statement:

> Where source support is adequate, the isolation-associated assemblage response is
> concentrated in genus-level composition, while additional sorting among species
> within represented genera is not supported.

This is positive evidence for an **assemblage-filtering / compositional route**. It
does not identify whether the responsible demographic process is arrival, establishment,
persistence or extinction.

## What remains unresolved

The active H1 state is one value per accepted species. Both all-analysis and Direct-only
state tables contain no population/locality dimension and no duplicate
accepted-species × trait state. Consequently, within-species trait change cannot be
estimated from the current H1 input.

The audit identified 10,751 trait-resolved species that occur on at least two
source-backed native-nonendemic islands and at least one mainland source entity. These
are candidates for a true island-versus-source population design, not evidence that
within-lineage change occurred.

## Population-level bridge

A recoverable route exists for mating system. Whitehead et al. (2018) compiled
population-level multilocus outcrossing rates (tm), whereas the Chapter 1 ledger retained
only species-level mating states. The hash-pinned Whitehead source was audited against the 10,751 repeated-lineage
candidates. It contains **743 population rows across 105 species with numeric species
means**; **36 Chapter 1 candidate species overlap, contributing 255 population rows with
numeric tm**. This is a real recoverable within-lineage data bridge.

However, the Whitehead population table contains no explicit latitude, longitude,
locality, site name or island/source label. Therefore the island-versus-source contrast
is **not yet directly estimable** from that table alone. The missing object is geography,
not the mating-system measurement.

The required next gate is strict:
1. recover study-specific locality metadata for the 36 overlapping species;
2. retain population-specific tm;
3. require at least one island and one mainland/source population for the same species;
4. preferably require both to come from the same underlying study;
5. fit the island/source or isolation contrast within species before pooling lineages;
6. make no claim of heritable evolution unless genetic/common-garden evidence separates
   evolutionary divergence from environmental plasticity.

## Provenance

Final validated implementation:
- workflow: Chapter 1 species sorting identifiability
- run: **37254767442**
- commit: **a5a76ed0052c71eee1f719af3601e4af7b387efd**
- artifact: **11322366778**

The final implementation retains all trait-resolved observed species when reporting the
source-overlap fraction while keeping the modeled decomposition explicitly
source-evaluable. A sparse-matrix optimization temporarily broke source membership
matching because CSR product column indices were not guaranteed sorted; the final code
sorts those indices before binary search and includes a regression guard. The validated
results above match the pre-optimization refined estimates.

Claim boundary:
- assemblage composition localized: yes, on testable source-matched contexts;
- colonization versus establishment/persistence/extinction: no;
- within-species phenotypic divergence: not yet;
- within-lineage evolution: not yet.
