# Chapter 1 frozen result — 2026-08-25

## Canonical run

This checkpoint freezes the current Chapter 1 scientific interpretation after the formal nested interaction test was added.

- workflow run: `32833362756`
- artifact: `chapter1-context-main-32833362756`
- artifact id: `9557730792`
- artifact digest: `sha256:e571f7117736e101e4f0d79f635183aadd75cc525d522d2e5a2715346f7dc9f9`
- direct-trait input run: `32702160934`
- floristic-status input run: `32559322028`
- isolation/context input run: `29228212586`
- M3 genus-fixed null draws: `1000`

The fitted Chapter 1 design contains no Bombus or other pollinator predictor.

## Primary conclusion

> **The current data do not support a confirmatory claim that the effect of island isolation differs among biogeographic contexts. They do support trait-component-specific floral reorganization associated with isolation, with the strongest within-context signals concentrated in northern mid-latitude islands.**

This distinction is mandatory. A significant isolation slope in one region and a nonsignificant slope in another is not, by itself, evidence that the two slopes differ.

## Formal contingency test

The canonical nested ladder is:

```text
M0 = baseline covariates + biogeographic-context main effects
M1 = M0 + one common isolation slope
M2 = M1 + isolation × biogeographic-context interactions
```

For each atomic trait category, M2 uses a cluster-robust joint Wald test of all `isolation × context` interaction coefficients. Joint p-values are BH-corrected across atomic states within each `stratum × trait` family.

Frozen result:

- fitted atomic categories: **43**
- context-specific slope rows: **95**
- FDR-supported within-context isolation slopes: **17**
- FDR-supported joint contingency categories: **1**

The sole FDR-supported joint contingency category is:

- stratum: `endemic`
- trait: `flower_primary_color`
- state: `yellow_orange`
- joint Wald p = **0.026617**
- joint FDR q = **0.026617**
- northern-midlatitude slope = **-0.7224** log-odds per SD, n = **32** islands
- tropical slope = **-0.1617**, n = **30** islands

Both context counts are only **pilot-level (30–49 islands)**. Therefore this result does **not** establish confirmatory biogeographic contingency.

## Trait-component isolation signals

Of the 17 FDR-supported within-context slopes, **16 meet the count component of confirmatory support and all 16 occur in the northern-midlatitude context**. The remaining slope is the pilot endemic `yellow_orange` result above.

Repeated northern-midlatitude signals in both `all_native` and `native_nonendemic` strata include:

- `bell_campanulate` — positive
- `composite_head` — positive
- `funnel_trumpet` — positive
- `spurred` — positive
- `other_described` — negative
- `papilionaceous` — negative
- `salverform` — negative
- `red_pink` — negative

These are **component-specific changes**, not a coherent classical island-syndrome vector.

The context-aware M1 common-isolation model also recovers FDR-supported atomic effects in `all_native` and/or `native_nonendemic` strata, including positive bell/funnel/composite/spurred states and negative `red_pink` / `other_described` states. Thus the data contain isolation-associated floral reorganization even when confirmatory slope heterogeneity among regions is not established.

## M3 lineage guardrail

M3 compares broad outcomes against a genus-composition-preserving null. In confirmatory-support northern-midlatitude and tropical subsets, the broad outcomes

- `generalized_form`
- `plain_colour`
- `self_compatibility`

do **not** retain FDR-supported isolation slopes for `all_native` or `native_nonendemic` strata.

Low-support exceptions exist in some northern-high-latitude or endemic subsets, but their island counts are below the declared confirmatory threshold and they cannot define the manuscript conclusion.

Therefore the current data do not support a manuscript-level claim of a coherent classic `generalized + plain + SC` island syndrome after lineage composition is constrained.

## Chapter 1 hypothesis status

The previously stated **Biogeographically contingent floral island syndrome** hypothesis is too strong for the frozen evidence and is demoted from the Chapter 1 headline hypothesis.

The new headline contrast is:

### H1 — coherent island-syndrome hypothesis

Isolation produces a coordinated floral/reproductive syndrome whose broad components move together and survive status/lineage safeguards.

### H2 — component-specific floral reorganization hypothesis

Isolation is associated with changes in particular floral components, while broad syndrome summaries may cancel, weaken, or disappear after category preservation and lineage control.

Current evidence favors **H2 over H1**.

A regional-contingency test remains important as a boundary-condition analysis, but the current joint-Wald result does not justify making confirmed biogeographic heterogeneity the headline finding.

## Pollination-syndrome interpretation boundary

Pollination syndromes are considered only after the above result is frozen.

The northern-midlatitude component vector may be compared with Bombus / large- or long-tongued-bee literature for **partial concordance or mismatch**. However, because the broad `generalized + plain + SC` signature is not recovered after M3, Chapter 1 must not describe the result as a demonstrated Bombus-loss syndrome.

Likewise, tropical/southern bird-, butterfly-, moth- or hawkmoth interpretations must not be forced: the current frozen result does not provide a confirmatory contrasting regional syndrome.

The permitted statement is:

> observed component-specific trait response ↔ literature-defined pollination-syndrome expectation = discussion-level concordance or mismatch

not:

> trait pattern → inferred pollinator guild → causal mechanism

## Chapter 1 → Chapter 2 handoff

Chapter 1 now ends with a sharper mechanistic question:

> **Why does isolation reorganize some floral components but not produce one coherent reproductive/floral syndrome, and which ecological interaction states determine which components respond?**

That question is handed to Chapter 2 (`izu-core`), where pollinator identity, functional diversity, trait matching, effective service, reproductive assurance, replacement and network context can be distinguished mechanistically.
