# Thesis positioning — Chapter 1

## Role in the dissertation

This repository is the **Chapter 1 / macroecological pattern-and-boundary** component of the dissertation.

The shared dissertation-level question is:

> **How does geographic isolation alter plant reproduction and floral phenotype, and why do different phenotype components respond differently across islands and lineages?**

`island` addresses the first step:

> **Does island isolation generate one coherent floral/reproductive syndrome, or does it reorganize particular phenotype components independently? Where are those component responses detectable, and do regional isolation effects genuinely differ?**

The companion repository [`zuizui0223/izu-core`](https://github.com/zuizui0223/izu-core) is Chapter 2 and addresses the mechanistic follow-up: which interaction states make responses branch, propagate or buffer. [`zuizui0223/shimahotarubukuro`](https://github.com/zuizui0223/shimahotarubukuro) is Chapter 3 and supplies a directly measured within-lineage floral phenotype for *Campanula microdonta* across the Izu island series.

## Chapter 1 headline hypothesis: component-specific floral reorganization

The current Chapter 1 contrast is:

### H1 — coherent island-syndrome hypothesis

> **Isolation produces a coordinated floral and reproductive syndrome whose broad components move together and remain after floristic-status and lineage safeguards.**

A classical version predicts coordinated shifts in broad summaries such as reproductive assurance / self-compatibility, floral generalization and reduced or simplified floral signalling.

### H2 — component-specific floral reorganization hypothesis

> **Isolation alters particular floral components without necessarily producing one coherent syndrome; broad summaries can cancel, weaken or disappear when multistate categories and inherited lineage composition are preserved.**

Current frozen evidence favors **H2 over H1**. See [`docs/chapter1_frozen_result_20260825.md`](docs/chapter1_frozen_result_20260825.md).

## Secondary boundary test — biogeographic contingency

Regional context remains scientifically important, but it is no longer the headline hypothesis.

The formal test is:

```text
M0 = baseline + biogeographic-context main effects
M1 = M0 + common isolation slope
M2 = M1 + isolation × biogeographic-context interactions
```

The key inferential rule is:

> **A significant isolation slope in one region and a nonsignificant slope in another is not evidence that the regional slopes differ.**

Regional heterogeneity is claimed only when the joint cluster-robust interaction test survives the declared multiplicity correction and has adequate context support.

The current frozen run finds only one FDR-supported joint contingency category (`endemic × yellow_orange`), based on 32 northern-midlatitude and 30 tropical islands. Because both counts are pilot-level, Chapter 1 currently has **no confirmatory biogeographic-contingency result**.

## What the current data do support

The current evidence supports isolation-associated **atomic trait reorganization**.

In the northern-midlatitude context, repeated FDR-supported signals in `all_native` and `native_nonendemic` strata include:

- `bell_campanulate` — positive;
- `composite_head` — positive;
- `funnel_trumpet` — positive;
- `spurred` — positive;
- `other_described` — negative;
- `papilionaceous` — negative;
- `salverform` — negative;
- `red_pink` — negative.

Sixteen FDR-supported within-context slopes meet the count component of confirmatory support, and all sixteen occur in northern-midlatitude islands. This is a strong **within-context concentration of signal**, but it is not by itself proof of a statistically different isolation effect from the tropics or Southern Hemisphere.

## Lineage guardrail and rejection of a simple classic syndrome

M3 uses a genus-composition-preserving null to test broad outcomes after inherited lineage structure is constrained.

In confirmatory-support northern-midlatitude and tropical subsets, the broad outcomes

- `generalized_form`;
- `plain_colour`;
- `self_compatibility`

do not retain FDR-supported isolation slopes in the main `all_native` or `native_nonendemic` strata.

Therefore Chapter 1 should not be summarized as a demonstrated `generalized + plain + SC` island syndrome.

The stronger conclusion is:

```text
isolation
→ trait-component-specific floral reorganization
≠ one coherent floral/reproductive syndrome
```

## Pollination-syndrome concordance is a discussion layer

Chapter 1 does not require a global Bombus-deficit variable for its primary inference.

The northern-midlatitude component vector may be compared with Bombus / large- or long-tongued-bee literature for **partial concordance or mismatch**. Because the broad classic syndrome is not recovered after M3, Chapter 1 must not reinterpret the northern signal as proof of a Bombus-loss syndrome.

Likewise, tropical or southern bird-, butterfly-, moth-, hawkmoth- or alternative-bee interpretations are permitted only if the frozen trait vector actually matches those literature-defined expectations. The current formal interaction result does not justify inventing a contrasting tropical/southern syndrome simply because northern slopes are stronger.

The allowed logic is:

```text
observed component-specific trait direction
<-> literature-defined pollination-syndrome expectation
= discussion-level concordance or mismatch
```

not:

```text
trait pattern -> inferred pollinator guild -> causal mechanism
```

## Status of the Bombus data products

Existing Bombus applicability, environmental-compatibility and occurrence-diagnostic products are retained for provenance, exploratory diagnostics and sensitivity analyses. They are not required to define the Chapter 1 primary model.

In particular:

- climatic compatibility is not realized occurrence or pollination service;
- opportunistic non-detection is not historical loss;
- a Bombus record is not equivalent to effective pollination;
- Bombus absence is not assumed to cause SC, floral generalization or colour change; and
- weak or null Bombus proxies are not rescued by additional acquisition simply to preserve a preferred mechanism.

## What Chapter 1 identifies

The primary estimand is **comparative assemblage reorganization**, not within-lineage evolutionary causation.

The analysis asks which floral/reproductive categories vary with isolation after accounting for the strongest available alternatives, including:

- island area and climate;
- biogeographic-context mean composition;
- floristic status / endemicity;
- lineage composition through a genus-preserving safeguard; and
- declared support thresholds and evidence resolution.

The intended inference is:

```text
baseline biogeography / lineage
→ common isolation-associated component effects
→ formal test of regional slope heterogeneity
→ category-preserving interpretation
→ candidate ecological concordance or mismatch
```

## Handoff to Chapter 2

Chapter 1 should end with an unresolved mechanistic question rather than answering it with proxy data:

> **Why does isolation reorganize some floral components but not produce one coherent syndrome, and which ecological interaction states determine which components respond?**

Candidate explanations include pollinator identity, long-tongued pollinator function, pollinator functional diversity, trait matching, effective pollination service, functional replacement, reproductive assurance, network context and non-pollination geography/history.

Chapter 2 (`izu-core`) is the appropriate place to distinguish those alternatives because it applies a stricter mechanistic standard and explicitly separates visitor identity, functional exposure, effective service, reproductive dependency and downstream response.

This makes the chapter transition stronger:

```text
Chapter 1: WHICH phenotype components reorganize under isolation?
            and is regional heterogeneity actually supported?
                         |
                         v
Chapter 2: WHY do particular response components branch / propagate / buffer?
                         |
                         v
Chapter 3: WHAT multidimensional phenotype actually diverges within one lineage?
```

## Relationship to Chapters 2 and 3

| Chapter 1 — `island` | Chapter 2 — `izu-core` | Chapter 3 — `shimahotarubukuro` |
| --- | --- | --- |
| global / multi-island comparative scale | mechanistic response architecture | one focal lineage across five Izu islands |
| asks **which components reorganize and where signals are detectable** | asks **how and why responses branch** | asks **what phenotype actually diverges** |
| tests coherent syndrome vs component-specific response | distinguishes candidate interaction mechanisms and response states | directly measures multidimensional floral divergence |
| regional heterogeneity requires formal joint interaction support | pollinator function is treated mechanistically | causal pollinator contrasts remain outside the morphometric pipeline |
| pollinator syndromes are discussion-level concordance/mismatch | separates visitor identity, service and dependency | does not equate Pst with Qst or selection |
| output is frozen atomic trait-response structure | output is mechanism / branching / propagation / buffering evidence | output is site-corrected trait differentiation and Pst / pairwise Pst |

## Claim boundary

Chapter 1 should not claim that:

- all isolated islands share one floral syndrome;
- a strong regional slope proves that the slope differs statistically from another region;
- current data demonstrate confirmatory biogeographic contingency;
- northern component signals constitute a demonstrated Bombus-loss syndrome;
- bird- or butterfly-like floral traits prove functional replacement;
- a pollination syndrome is a direct observation of pollinator use; or
- cross-sectional assemblage composition alone demonstrates within-lineage floral evolution.

The Chapter 1 contribution is narrower and stronger: **to show that isolation-associated floral change is component-specific rather than one coherent syndrome, to identify where those component signals are strongest, and to reserve claims of regional heterogeneity and ecological mechanism for tests that actually support them.**
