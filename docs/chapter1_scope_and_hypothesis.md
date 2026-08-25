# Chapter 1 scope: biogeographically contingent island floral filtering

## Decision

Chapter 1 is the **pattern and boundary** chapter. Its direct hypothesis is:

> **Biogeographically contingent island floral filtering hypothesis.** Island
> isolation does not produce one universal floral and reproductive syndrome.
> Instead, isolation-associated floristic and trait composition can vary among
> biogeographic and floristic contexts after area, climate, observation support,
> floristic status, and lineage composition are represented.

`Channel-gated island assembly` is retained one level above Chapter 1 as the
dissertation-wide conceptual explanation to be evaluated across chapters. It is
not the direct Chapter 1 estimand because Chapter 1 does not observe pollinator
functional deficit or effective replacement.

The direct question is therefore:

```text
Where and under what biogeographic and floristic contexts does isolation leave
a floral or reproductive compositional signal beyond status and lineage turnover?
```

The safer term is **filtering**, not **syndrome**. A syndrome should be named only
if multiple prespecified trait axes retain a coherent, coverage-robust direction.

## Direct hypotheses

### H1 — Universal syndrome test

Test whether pooled isolation associations form one consistent plain,
generalized, self-compatible, or otherwise simplified trait vector. A pooled
coefficient is not sufficient if region, floristic status, or lineage explains
the pattern.

### H2 — Biogeographic and floristic contingency

Test whether isolation associations differ among `Northern`, `Northern high
latitude`, `Tropical`, and `Southern` regions and among `all_native`,
`native_nonendemic`, `endemic`, and secondary `introduced` strata.

Endemic share is a separate response:

```text
Isolation -> endemic share of resolved native flora
```

It is not inserted as a nuisance covariate into the same response that it may
help generate.

### H3 — Category decomposition after lineage control

For each supported trait category, preserve the observed island genus
composition and analyse:

```text
observed direct category share - same-genus null expectation
    ~ distance + area + climate
```

Category results require the declared island-support boundary, spatial-block
uncertainty, false-discovery-rate control, and direct-coverage sensitivity.

## Evidence hierarchy

1. **Floristic status:** does isolation predict endemic or other status turnover?
2. **Broad residual composition:** do plain colour, generalized form, or SC remain
   associated with isolation after status and genus composition are represented?
3. **Fine residual categories:** do individual colour, form, architectural, or
   reproductive categories retain supported associations?
4. **Multivariate coherence:** only then ask whether the retained vector warrants
   the word `syndrome`.
5. **Mechanistic concordance:** only in Discussion, compare the retained vector
   with literature-defined pollinator-functional expectations.

This order prevents a raw whole-flora association from being promoted over a
status- or lineage-based explanation.

## Bombus and alternative pollinators

Bombus occurrence, absence, applicability, and hypervolume are not direct
Chapter 1 explanatory variables. The existing hypervolume is retained as an
exploratory climatic-opportunity diagnostic and as a recorded falsification: it
did not recover the predicted Northern pathway after status and genus
composition were represented.

If a coherent Northern trait vector is eventually recovered, it may be described
as **concordant with** documented loss of large/long-tongued pollinator function.
If a tropical or Southern vector is recovered, it may be compared symmetrically
with bird-, butterfly-, hawkmoth-, fly-, or other-bee expectations. Neither a
region nor a single floral trait identifies the pollinator channel.

No additional Bombus occurrence campaign is required to rescue Chapter 1.
Direct occurrence, visitation, trait matching, pollination effectiveness, and
functional replacement belong to the mechanism chapter.

## Current empirical boundary

The status/lineage workflow verified in PR #133 currently shows:

- a robust Northern isolation-associated increase in regional endemicity;
- approximately null broad plain/generalized/SC residual slopes among native
  non-endemics after genus composition is fixed;
- one small, comparatively coverage-robust Northern all-native decrease in
  red/pink residual share;
- a tropical non-endemic red/pink decrease that collapses at higher direct
  coverage; and
- no floral-form category with at least 50 islands surviving FDR.

The current-main post-pilot port in PR #135 adds that tropical endemic
generalized form is unsupported at the frozen 30-island acquisition boundary
(`29` complete-case islands, beta `-0.0432`, p `0.122`). Northern endemic form
and SC retain nominal coefficients but only on 20 and 12 islands, respectively,
and remain below the pilot support boundary.

Thus Chapter 1 currently establishes a strong **boundary result**: isolation is
associated most clearly with floristic/endemic-lineage turnover, while a broad
residual floral syndrome has not been recovered. This is a falsification of the
simple syndrome claim, not a failed analysis.

## Chapter boundaries

```text
Chapter 1 — Pattern and boundary
  isolation x biogeographic/floristic context -> residual composition

Chapter 2 — Mechanism
  pollinator identity vs functional diversity vs effective service
  vs replacement vs reproductive assurance

Chapter 3 — Phenotypic realization
  within-lineage divergence of measured floral components

Dissertation synthesis
  integrate the three levels to evaluate channel-gated island assembly
```

Chapter 1 cannot identify historical Bombus loss, functional replacement,
causal floral evolution, or assembly versus in-situ evolution. Those limits are
part of the result and must remain explicit.
