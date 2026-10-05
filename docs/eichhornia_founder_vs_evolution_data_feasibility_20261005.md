# Eichhornia founder-versus-evolution data feasibility audit — 2026-10-05

## Decision

The decisive founder-filtering versus post-colonization test is **not yet executable from
one existing public dataset**.

The required ingredients exist, but they are split across different sampling designs:

1. broad population ancestry / colonization history;
2. causal or mapped mating-system modifier loci;
3. dense source-population allele frequencies at those causal loci.

The first two exist. The third is the missing bridge.

## Existing data layers

### Population-history layer

Ness, Wright & Barrett (2010), DOI 10.1534/genetics.109.110130:

- 225 individuals from 25 populations;
- north-eastern Brazil, Jamaica and Cuba;
- 10 EST-derived nuclear loci;
- population morph structure and selfing-variant frequency;
- strong regional genetic structure;
- Caribbean populations linked to a Brazilian source;
- coalescent colonization estimate on the order of 125,000 years.

This layer is suitable for population history, but the ten loci were designed as broadly
distributed nuclear markers rather than the later-mapped causal selfing architecture.

### Functional / mapping layer

Arunkumar et al. (2017), DOI 10.1111/mec.13946:

- independently derived L- and M-morph semi-homostylous selfing variants;
- 462 backcross progeny;
- 1,450 GBS markers;
- QTL mapping of style length and anther height;
- selfing variants governed by different modifier loci.

This layer identifies genetic architecture but is not a dense Brazil–Caribbean
population-frequency survey of the mapped modifier haplotypes.

### Transcriptome / comparative-genome layer

Public sequence resources include:

- BioProject PRJNA266681 / SRP049636: transcriptome/genomic comparisons including
  Brazilian outcrossing, Jamaican selfing and Nicaraguan selfing material;
- PRJNA310302 / PRJNA310303: mapping/reference resources associated with the earlier
  genetic-architecture work;
- PRJNA863567: modern chromosome-scale *E. paniculata* reference and related genomic /
  transcriptomic data.

These resources make candidate-locus identification and read remapping feasible. They do
not provide dense source-population sampling at the functional modifier loci.

## Why the missing bridge matters

A Caribbean-private present-day allele is not enough to infer post-colonization
evolution. It may simply have existed at low frequency in the Brazilian source and been
sampled by the founder event.

Conversely, finding a functionally equivalent selfing allele in Brazil would directly
support founder-genotype filtering even if subsequent Caribbean-specific modifiers also
evolved.

Therefore the key missing object is:

> **source-population allele-frequency and haplotype data at the causal / mapped selfing
> modifier regions, sampled across the ancestry-matched Brazilian source range.**

## Minimal new data product

The smallest decisive dataset is not another global trait table. It is a targeted panel:

- ancestry-matched Brazilian source populations from the 2010 sampling frame;
- Jamaican and Cuban populations from the same ancestry cluster;
- targeted capture / amplicon / low-pass sequencing covering mapped selfing-modifier
  intervals and flanking neutral haplotypes;
- enough individuals per source population to bound an unobserved source allele
  frequency.

The analysis should retain two distinct outputs:

1. **presence / frequency test:** is the Caribbean functional allele already present in
   source populations?
2. **genealogy / age test:** does the functional haplotype coalesce before or after the
   Caribbean colonization split?

## Detection bound

If a functional allele is not observed among (n) independent source chromosomes, the
absence claim should be expressed as a frequency bound rather than literal absence.

For zero observations under simple random sampling,

```
upper 95% source-frequency bound = 1 - 0.05^(1/n)
```

Examples:

- n = 20 chromosomes -> upper bound about 0.139;
- n = 50 -> about 0.058;
- n = 100 -> about 0.030;
- n = 200 -> about 0.015.

Thus a convincing “Caribbean-specific” claim requires substantially deeper source
sampling than a handful of Brazilian genotypes.

## Final identification rule

### Founder filtering supported

Any of the following is sufficient to block a pure post-colonization claim:

- the same functional allele / haplotype occurs in the ancestry-matched Brazilian source;
- the modifier genealogy predates the Caribbean split;
- Caribbean functional haplotypes are nested inside source variation in a way consistent
  with founder sampling.

### Post-colonization evolution supported

Require all of:

1. functionally mapped Caribbean-derived allele or haplotype;
2. dense ancestry-matched source sampling with a prespecified low upper frequency bound;
3. Caribbean monophyly / derived state at the functional locus;
4. posterior support that modifier origin postdates the Caribbean split;
5. robustness to gene flow and incomplete lineage sorting.

### Mixed route

A source selfing allele passes the colonization filter, but additional island-specific
modifier changes postdate establishment.

## Current blocker

No currently identified public dataset contains both:

- broad Brazil–Caribbean population sampling, and
- the causal / mapped modifier loci at sufficient population depth.

Therefore the next empirical step is **targeted population genotyping at the modifier
regions**, not another reanalysis of species-level traits.

This is the natural stop point for the current public-data-only route.
