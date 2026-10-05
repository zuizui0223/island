# Eichhornia paniculata: founder filtering versus post-colonization evolution test

Status: candidate design frozen before any new genotype reanalysis.

## Why this system is highest priority

*Eichhornia paniculata* is unusually close to the data structure required to separate
within-species founder filtering from post-colonization evolution.

Existing work provides:

- a mainland source region in north-eastern Brazil and Caribbean island populations in
  Jamaica and Cuba;
- population-level morph structure, selfing-variant frequency and mating-system data;
- glasshouse-grown material demonstrating persistent floral/mating-system differences;
- multilocus population-genetic evidence linking Caribbean populations to Brazil;
- a coalescent estimate of Caribbean colonization on the order of 125,000 years ago;
- genetic mapping showing that independently derived selfing variants can be controlled
  by different mating-system modifier loci;
- modern reference-genome and sequencing resources.

The existing literature also makes an important founder-filter prediction: Jamaica and
Cuba are consistent with descent from a shared long-distance dispersal event, and a
self-pollinating M-morph founder from Brazil is a leading historical scenario. Thus this
system should not be treated as prior proof of post-colonization evolution.

## Competing models

### F — founder-genotype filtering

The Caribbean was founded by one or more Brazilian genotypes already carrying the
selfing-capacity allele/haplotype.

Predictions:

- the Caribbean selfing-modifier haplotype is present, or ancestrally represented, in
  Brazil;
- the modifier lineage predates Caribbean colonization;
- Caribbean genomes are nested within / derived from a source ancestry carrying the same
  functional state;
- no island-specific derived causal mutation is required to explain selfing capacity.

### E — post-colonization evolution

The founder lineage lacked the derived Caribbean selfing modifier and the functional
change arose after island establishment.

Predictions:

- the functional Caribbean allele/haplotype is derived and monophyletic within the
  Caribbean lineage;
- it is absent from dense ancestry-matched Brazilian source sampling, beyond a
  prespecified detection bound;
- the inferred age of the causal modifier is younger than the Caribbean colonization
  split;
- ancestry-matched mainland/source haplotypes retain the outcrossing architecture.

### M — mixed route

A selfing-capable founder passed the colonization filter, followed by additional derived
Caribbean changes in selfing syndrome or floral architecture.

This is biologically plausible because historical crossing studies indicate that the
genetic basis of stamen modification differs among regional selfing variants.

## Primary temporal estimand

Let

- (T_C) = Caribbean colonization time;
- (T_M) = origin time of the focal Caribbean selfing-modifier haplotype.

The primary discriminator is the posterior ordering of these times.

```
Founder filtering favoured:
P(T_M > T_C | data) >= 0.95
or the functional allele is observed in the source population.

Post-colonization evolution favoured:
P(T_M < T_C | data) >= 0.95
AND the derived functional allele is Caribbean-monophyletic
AND dense ancestry-matched source sampling excludes that allele.

Otherwise:
unresolved / mixed.
```

The 0.95 rule is a design target, not a retrospectively chosen significance threshold.

## Required analyses

1. Reconstruct the Brazil–Caribbean ancestry graph using genome-wide neutral markers.
2. Identify the functional/selfing-modifier haplotype(s) associated with Caribbean
   semi-homostyly.
3. Genotype those loci densely across Brazilian, Jamaican and Cuban populations.
4. Estimate modifier-haplotype genealogy and mutation/origin time jointly with the
   population split.
5. Compare (T_M) with (T_C) under explicit uncertainty.
6. Repeat the test separately for primary selfing capacity and downstream floral-syndrome
   modifiers; do not force both into one event.

## Public genomic resources already identified

- population-genetic history: Ness, Wright & Barrett 2010,
  DOI 10.1534/genetics.109.110130;
- genetic architecture of tristyly/selfing modifiers: Arunkumar et al. 2017,
  DOI 10.1111/mec.13946;
- older sequence resources: BioProject PRJNA266681;
- mapped/genome resources associated with PRJNA310302 / PRJNA310303;
- chromosome-scale *E. paniculata* reference and current genomic resources:
  BioProject PRJNA863567.

These resources establish feasibility; they do not by themselves provide the final
founder-versus-evolution classification.

## Stop rules

- Do not label a Caribbean-private SNP as causal without functional or mapping evidence.
- Do not treat absence from small Brazilian samples as ancestral absence.
- Do not compare allele age and colonization time without propagating uncertainty.
- If causal modifiers remain unidentified, report only ancestry-resolved founder
  compatibility, not level-3 evolution.
- If the functional allele is found in source populations, level 2 founder filtering is
  supported even if subsequent island-specific floral modifiers also evolved.

## Relation to Chapter 1

This test addresses the remaining mechanistic gap after the global assemblage
decomposition:

```
global H1
  -> genus-level assemblage filtering         [identified in testable contexts]
  -> founder-genotype filtering within species [not globally quantified]
  -> post-colonization evolutionary change     [not globally quantified]
```

The *Eichhornia* analysis is therefore a mechanistic validation system, not a replacement
for the global Chapter 1 estimand.
