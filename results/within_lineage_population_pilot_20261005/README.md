# Within-lineage population pilot — 2026-10-05

Status: feasibility / boundary-condition analysis. This does not replace the Chapter 1 assemblage inference and is not promoted as a global within-lineage evolutionary result.

## Question

Does the first recoverable population-level evidence support a simple rule that island populations become more selfing within species?

The pilot deliberately requires population-specific mating measurements rather than the species-level trait states used by Chapter 1.

## Arabidopsis lyrata

Source: Willi & Määttänen (2010), *Journal of Evolutionary Biology*, DOI 10.1111/j.1420-9101.2010.02073.x.

The source reports 18 North American populations with coordinates and multilocus outcrossing rate, `tm`. Three sites are explicitly named freshwater islands:

- Beaver Island: tm = 0.878
- Isle Royale: tm = 0.134
- Apostle Islands: tm = 0.903

The remaining 15 populations are treated only as the non-island comparison set for this feasibility test.

Observed means:

- island mean tm = **0.63833**
- non-island mean tm = **0.74407**
- island − non-island tm = **−0.10573**
- equivalently, island selfing fraction is higher by **0.10573** on average.

An exact 3-of-18 island-label permutation test gives:

- two-sided P = **0.67034**
- one-sided lower-tm-on-islands P = **0.28554**
- exact labelings = **816**

Thus the three island populations do not provide evidence for a uniform shift toward lower outcrossing. The island values themselves are strongly heterogeneous: one is predominantly selfing while two are predominantly outcrossing.

This comparison is not a clean oceanic-isolation experiment. The island sites are Great Lakes islands, geography and phylogeographic history are structured, and the test is retained as a feasibility / negative-control analysis rather than a Chapter 1 effect estimate.

## Phormium tenax

Source: Howell & Jesson (2013), *New Zealand Journal of Botany*, DOI 10.1080/0028825X.2013.772904.

The source directly compares two latitude-matched island-mainland pairs:

1. Little Barrier Island vs Tāwharanui Regional Park
2. Tiritiri Matangi Island vs Shakespear Regional Park

The paper reports **no significant difference in seedling multilocus outcrossing rates in either island-mainland comparison**. All four populations were highly outcrossing.

The repository retains that source-reported inference rather than digitizing Figure 4 into pseudo-precision.

## Interpretation

The first population-level bridge does **not** support a simple universal rule that island populations become more selfing within species.

That negative result is compatible with the Chapter 1 source-matched assembly result:

- assemblage filtering can be strong at genus composition;
- additional within-genus species sorting can be weak;
- within-species population responses can still be heterogeneous rather than uniformly directional.

This does not show that within-lineage evolution is absent. Two species are far too few, and neither analysis separates heritable evolutionary divergence from plasticity or other local ecological effects.

## Next gate

Expand only with studies that satisfy all of the following:

1. the same accepted species has multiple population-specific reproductive measurements;
2. at least one population is explicitly an island and at least one is mainland/source;
3. site identity or coordinates are recoverable from the original study;
4. the mating-system metric is comparable within the study;
5. population contrasts are not reconstructed from species means;
6. evolutionary claims require independent genetic/common-garden evidence beyond phenotypic divergence.

Whitehead et al. (2018) provides 36 Chapter 1-overlapping species and 255 population `tm` rows, but its compiled table lacks site geography. Original-study locality recovery is therefore the next scalable route.

## Provenance

Validated workflow:
- workflow: Within-lineage population pilot
- run: **37258014953**
- commit: **83bd8df336f13fc85dd7f712a42fd6a0dcbabed4**
- artifact: **11323152646**

Curated inputs:
- `data/v2/curation/within_lineage_population_pilot_arabidopsis_lyrata_20261005.csv`
- `data/v2/curation/within_lineage_population_pilot_phormium_pairs_20261005.csv`

Implementation:
- `src/island_v2/within_lineage_population_pilot.py`
- `scripts/analyze_within_lineage_population_pilot.py`
