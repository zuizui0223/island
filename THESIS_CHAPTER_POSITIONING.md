# Thesis positioning — Chapter 1

## Role in the dissertation

This repository is the **Chapter 1 / macroecological boundary-condition** component of the dissertation.

The shared dissertation-level question is:

> **How does geographic isolation alter plant reproduction through changes in ecological interactions, and why do those changes produce different floral outcomes across islands and lineages?**

`island` addresses the first half of that question:

> **Where and under what biogeographic conditions does island isolation generate floral and reproductive filtering?**

The companion repository [`zuizui0223/izu-core`](https://github.com/zuizui0223/izu-core) is Chapter 2 and addresses the mechanistic follow-up: how a changed pollination environment is translated into reproductive and floral responses within lineages and field systems. [`zuizui0223/shimahotarubukuro`](https://github.com/zuizui0223/shimahotarubukuro) is Chapter 3 and supplies a high-resolution within-lineage empirical phenotype for *Campanula microdonta* across the Izu island series.

## Chapter 1 hypothesis: channel-gated island assembly

The working hypothesis is the **Channel-gated island assembly hypothesis**.

> Island isolation should not generate a universal floral syndrome. Instead, reproductive and floral trait filtering should emerge when isolation disrupts a pollination channel that is regionally available and biologically relevant to the source-pool context. Structural absence of a pollinator guild is therefore not equivalent to its loss. Comparable geographic isolation may produce different floral assemblages depending on whether the relevant channel is retained, disrupted, or not regionally applicable.

A broader conceptual implication is that pollinator functional groups can differ in their susceptibility to oceanic barriers. Consequently, a classical floral island syndrome should be strongest where isolation removes a regionally available pollination channel, and should weaken, disappear, or take a different phenotypic direction where that channel is structurally absent and other pollination functions can persist.

## Operational test in this repository

The comparative analysis deliberately separates three Bombus states:

1. **retained** — Bombus is regionally applicable and supported on the island;
2. **expected-but-deficient** — Bombus is native to the source-region context and the island is environmentally compatible, but island evidence indicates a reduced or missing Bombus channel under adequate observation;
3. **structurally absent** — Bombus is not native to the relevant source-region context, so absence is not interpreted as channel loss.

The primary Chapter 1 prediction is not `Bombus absent -> island syndrome`. It is:

```text
regionally applicable Bombus channel
+ island isolation
+ expected-but-deficient Bombus state
-> increased representation of classical island-syndrome traits
```

The focal response domains are reproductive assurance / self-compatibility and floral generalization, with plain or inconspicuous colour retained as a secondary phenotype rather than a universal endpoint.

## Falsification and counter-pattern

The structurally absent regime is not treated as a failed Bombus projection. It is an explicit falsification condition.

If Bombus absence alone caused the classical syndrome, structurally absent regions should converge on the same reproductive and floral pattern. The alternative prediction is that these regions need not show the same generalization. Existing category-preserving floral traits are therefore used to test whether conspicuous colours and tubular, deep, or otherwise specialized floral forms remain represented.

Such a counter-pattern may be discussed as **consistent with** pollination syndromes associated with birds, Lepidoptera, hawkmoths, other bees, or other alternative functional groups. This repository does **not** infer that those pollinators actually caused the pattern unless independent interaction/effectiveness evidence exists.

## What Chapter 1 identifies

The primary estimand is **comparative assemblage filtering**, not within-lineage evolutionary causation.

The analysis asks whether floral and reproductive composition differs conditionally across island regimes after accounting for the strongest available alternatives, including:

- island isolation, area, and climate;
- floristic status / endemicity;
- source-pool and lineage composition where available;
- establishment / reachability structure; and
- observation and evidence-coverage processes.

The intended inference is therefore:

```text
biogeographic filtering first
-> residual interaction-associated filtering second
```

A Bombus-associated residual pattern is evidence compatible with a channel-gated assembly mechanism. It is not, by itself, evidence that historical Bombus loss caused floral evolution within an extant lineage.

## Relationship to Chapters 2 and 3

The three chapters answer different levels of the same problem.

| Chapter 1 — `island` | Chapter 2 — `izu-core` | Chapter 3 — `shimahotarubukuro` |
| --- | --- | --- |
| global / multi-island comparative scale | Izu / lineage / plant / network scale | one focal lineage across five Izu islands |
| asks **where and when** filtering appears | asks **how and why** biological responses are generated | asks **what phenotype actually diverges** |
| identifies boundary conditions for assemblage filtering | identifies reproductive and interaction mechanisms | directly measures multidimensional floral divergence |
| Bombus is a regional channel indicator / natural experiment | Bombus identity alone is insufficient; effectiveness and dependency must be identified | causal Bombus contrasts are deliberately excluded |
| output is trait-composition / assemblage association | output is reproductive, pollen-transfer, dependency, and response-mode evidence | output is site-corrected trait differentiation and Pst / pairwise Pst |
| does not claim within-lineage evolution from cross-sectional composition | separates local reproduction `F(z)` from establishment / reachability `E(z)` | measures phenotype but does not equate Pst with Qst or selection |

The handoff is:

```text
Chapter 1: Geography -> channel regime -> assemblage filtering
                         |
                         v
Chapter 2: channel regime -> plant-specific functional/reproductive response
                         |
                         v
Chapter 3: focal lineage -> realized multidimensional floral divergence
```

Chapter 1 establishes whether a pollination-channel boundary is a plausible condition for island floral filtering. Chapter 2 tests why the same broad functional change can yield continuous erosion, threshold reproductive assurance, interaction rewiring, no shared morphology response, or alternative-guild specialization. Chapter 3 then anchors those broader ideas in a directly measured focal phenotype, showing which *C. microdonta* floral axes actually differ among island populations without treating the pattern itself as causal proof.

## Claim boundary

This repository should not claim that:

- all isolated islands share one floral syndrome;
- simple Bombus absence is equivalent to Bombus loss;
- a floral syndrome identifies the pollinator that caused it;
- bird- or butterfly-like floral traits prove functional replacement; or
- cross-sectional assemblage composition alone demonstrates within-lineage floral evolution.

The Chapter 1 contribution is narrower and stronger: **to identify the biogeographic conditions under which a classical island floral syndrome appears, fails to appear, or takes a different phenotypic direction.**
