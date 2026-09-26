> **HISTORICAL / SUPERSEDED — pre-corrected Chapter 1 surface.** Retained for provenance/replay only. The current submission is selected by `config/chapter1_submission_current.json`; use `submission/chapter1_current/MANUSCRIPT.md` and `docs/PAPER_PIPELINE.md` for current results.

# Chapter 1 manuscript v8 reframing — 2026-09-13

## Status

This note is the canonical writing target after the merged hierarchical-depth and response-geometry audits (PR #212). It does not replace the frozen H1–H5 analysis contract. It changes the manuscript product: what the paper leads with, what is secondary, and what is explicitly reported as a failed or non-promoted explanation.

The working question is no longer **when and where**. No island age, colonisation time, divergence time or other temporal exposure is in the Chapter 1 analysis. The manuscript question is:

> **Where, whether, and at what assembly level does an isolation-associated floral island response emerge?**

A fourth, explicitly secondary question is:

> **Does the response depart from a smooth/monotonic geographic gradient in an identifiable way?**

The answer to the fourth question is currently negative under the frozen geometry gate.

---

## 1. Working title

### Preferred EL-facing title

> **Biogeographic context reshapes the direction and hierarchical expression of island floral assembly**

This title abstracts one level above “the island floral syndrome is not one syndrome” without claiming a universal taxonomic-depth law. “Hierarchical expression” covers source/family/genus/beyond-genus representation while allowing the new matched depth result to remain bounded.

### Safer GEB fallback

> **The island floral syndrome is not one syndrome: biogeographic branching, source filtering and lineage assembly across global islands**

### Do not use

- “when and where” in title/abstract/Figure 1;
- “assembly versus evolution” as a solved dichotomy;
- “pollinator-loss threshold” or “non-monotonic island syndrome” as a Chapter 1 headline;
- “endemic versus non-endemic separates in-situ evolution from filtering”.

---

## 2. One-sentence manuscript claim

> **Isolation is associated with different floral and reproductive response directions in different biogeographic contexts, and the hierarchical level at which those responses are represented is itself context dependent under bounded conditions; neither a universal area mechanism, a universal nonlinear threshold, nor a pollinator-channel mechanism is promoted by the current global data.**

The shortest public-facing version is:

> **The ecological response to island isolation changes not only in direction, but also in where it is expressed within assemblage hierarchy.**

The second sentence is a generalisation/interpretive headline, not a claim that every context has a different taxonomic depth.

---

## 3. Canonical evidence hierarchy for v8

### Primary result A — direction is context dependent

H1/H2 remain first.

- one universal floral/reproductive syndrome is not recovered;
- Palearctic source separation is associated with increased accessibility/generalisation and increased reproductive assurance;
- Tropical direct-only assemblages show increased reproductive assurance while accessibility/generalisation decreases;
- therefore the two major components are not an obligatory serial syndrome.

This remains the first empirical result because it is the preregistered/frozen primary branching architecture.

### Primary result B — the strongest Palearctic response is genus-level assembly

Wave52 is canonical, not the older Wave36 V2 checkpoint.

Wave52 Palearctic primary vector:

- observed: 4/4 source modes;
- after family expectation: 4/4;
- after source-matched genus expectation: 0/4;
- same classification in all-analysis/direct-only and all-native/native-nonendemic.

Canonical wording:

> **The broad Palearctic response persists beyond measured family composition but is compatible with source-matched genus-level lineage assembly.**

Do not retain the obsolete Wave36 sentence that the broad Palearctic response survives genus adjustment.

### Primary result C — hierarchical expression varies among contexts, but only as a bounded secondary re-expression

PR #212 added a matched context × taxonomic-stage audit on the frozen final dataset using the same four-component plant-architecture response before versus after genus adjustment.

Native non-endemic assemblages show direct context × stage differences under all four source modes in both evidence scopes:

- all-analysis q range: 0.00462–0.02166;
- direct-only q range: 0.00972–0.04921.

The all-native result is weaker:

- all-analysis: 0/4 source modes pass the frozen family;
- direct-only: 3/4 pass, with the remaining source mode q≈0.0598.

Publication classification is therefore `bounded_hierarchical_depth_signal`.

This result supports the broader conceptual statement that the hierarchical expression of an isolation response can depend on biogeographic context, but it must be described as a **matched secondary re-expression of a frozen dataset**, not as a newly preregistered primary test.

Do not infer in-situ evolution from a beyond-genus residual. Such a residual can combine within-genus species sorting, unmeasured source composition, other ecological filters and genuine within-lineage change.

### Primary negative D — area does not earn mechanism promotion

H4 remains a useful competing explanation rather than a side robustness result.

- 0/16 cells pass the frozen heteroskedastic-null promotion gate;
- 10/16 show apparent small-island amplification across the fitted sensitivity modes;
- therefore area is retained as a measurement-sensitive modifier, not promoted as a capacity/founder/pollinator-persistence mechanism.

Canonical wording:

> **Apparent small-island amplification was not distinguishable from measurement-sensitive heterogeneity strongly enough to promote an area/capacity mechanism.**

Do not write “area had no effect”.

### Primary negative E — global nonlinear geometry was testable but not promoted

The nonlinear analysis was staged specifically to prevent post-hoc breakpoint hunting.

1. Candidate geometries were frozen before observed nonlinear inspection: G0 flat, G1 cline, G2 step, G3 hinge, G4 reversal.
2. V1 showed naive AICc selection was unsafe under the observed spatial-block structure, so observed geometry stayed closed.
3. V2 used independent calibration and validation simulation seeds while preserving the real island exposure distribution, trial counts, baseline covariates, spatial blocks and response-specific missingness.
4. V2 cross-scope qualification: step 11/12 cells, hinge 0/12, reversal 4/12.
5. Observed geometry was then opened once under the frozen V2 gate.
6. Final observed result: **12/12 cross-scope cells = `monotonic_or_unresolved`; identified step = 0; identified reversal = 0.**

Two scope-specific nonlinear signals crossed their calibrated gate but failed replication across evidence scopes and were not promoted.

Canonical wording:

> **A prospective geometry audit found no cross-evidence support for a common step or reversal in the broad floral, structural or reproductive contrasts. Isolation therefore reorganises response direction and hierarchical expression without evidence for one global nonlinear transition.**

This belongs in the main Results as a short bounded negative, not as a new headline figure unless reviewers request it.

### Primary negative F — independent pollinator-channel heterogeneity failed prospectively

The N1 chain stays closed.

- joint isolation × channel Wald test: W=1.6187, df=3, p=0.65516;
- leave-one-spatial-block-out and leave-one-source-region-out checks remained non-significant;
- the preregistered action was stop before N2;
- no reduced-channel rescue, threshold retuning or N2 opening is allowed.

The descriptive channel slopes may be reported only as descriptive compatibility, not mechanism evidence.

Canonical wording:

> **An independently frozen pollinator-channel test was executable but did not support channel-specific isolation responses; the preregistered chain therefore stopped before the mechanistic N2 stage.**

This is a transparency asset and a claim-ceiling result, not a reason to hide the pollinator hypothesis.

---

## 4. Abstract v8 draft

Island floras are often expected to converge toward a coherent floral “island syndrome”, combining reproductive assurance with reduced floral specialisation or attraction. Yet isolation filters regional source pools before it can act on assemblage-level phenotype, and the same oceanic barrier need not have the same biological meaning across biogeographic contexts. We asked whether isolation produces one floral/reproductive response, where departures from that response occur, and at what hierarchical level those departures are represented.

We analysed a fixed universe of 8,265 islands and 106,295 accepted angiosperm species under a progressive contract that froze hypotheses, support gates, source-pool safeguards, model order and claim ceilings while trait evidence improved. The final snapshot resolved 222,688 of 318,885 species-by-axis cells across flower colour, floral structural complexity and reproductive assurance. Mainland distance was treated as a composite gradient of source separation, connectivity and accessibility rather than as a mechanistically pure treatment.

A universal floral/reproductive syndrome was not recovered. In the Palearctic, increasing separation was associated with greater accessibility/generalisation and reproductive assurance, whereas Tropical direct-only assemblages combined increasing reproductive assurance with decreasing accessibility/generalisation. The strongest Palearctic response persisted after family adjustment but disappeared after source-matched genus adjustment across all four source definitions in both evidence scopes and both all-native and native-nonendemic strata, placing the primary response at the level of genus assembly rather than a robust beyond-genus residual.

A matched secondary audit further showed that the change in plant-architecture response across taxonomic stages differed between northern-midlatitude and Tropical contexts across all four source definitions in native non-endemics in both evidence scopes, although this result weakened when island endemics were included. By contrast, continuous island area did not pass the prespecified mechanism-promotion gate, and a separately calibrated response-geometry audit promoted no global step or reversal across colour, structure or reproductive contrasts. An independently frozen pollinator-channel test also failed its prospective heterogeneity gate and stopped before downstream mechanism testing.

These results indicate that isolation does not simply strengthen one island floral syndrome. Instead, response direction and hierarchical expression depend on biogeographic context, with lineage assembly forming a major part of the global pattern. The global data delimit rather than identify pollination mechanism; direct functional tests of regime transitions belong at the within-system scale.

---

## 5. Introduction: final paragraph replacement

Replace the current final prediction paragraph with the following logic:

> Our central prediction is not that distant islands should be uniformly selfing, generalised or inconspicuous. We instead ask three nested questions. First, **whether** source separation produces one coherent floral/reproductive response or different response directions among biogeographic contexts. Second, **where** those contrasting responses are expressed geographically. Third, **at what assembly level** a supported response remains after source-matched family and genus expectations are removed. A universal island syndrome predicts common response direction; a context-dependent assembly model predicts decoupled response components and context-specific hierarchical expression. We additionally tested, under a separately calibrated identifiability gate, whether broad response components exhibit reproducible nonlinear transitions rather than smooth or unresolved geographic change. Mechanistic attribution to pollinator loss, replacement or effective service was reserved for independent pollinator data rather than inferred from plant phenotype.

Do not use “when” in this paragraph.

---

## 6. Methods order for v8

Keep the existing data/trait/support sections, then order the inferential sections as follows:

1. **H1: universal-syndrome rival** — unchanged.
2. **H2: biogeographic branching** — unchanged; still the primary pollinator-name-free vector.
3. **H3: source and lineage assembly** — explicitly declare Wave52 as the canonical final depth result.
4. **Matched context × taxonomic-stage audit** — new short subsection; label as secondary re-expression, not outcome-blind preregistration.
5. **H4: area moderation** — keep as a competing modifier/mechanism test.
6. **Response-geometry identifiability audit** — new subsection. Report V1 failure, V2 calibration/validation, then one-shot observed opening. Do not list observed breakpoint estimates for non-promoted cells in the main text.
7. **Pollination-associated architecture concordance** — secondary.
8. **Prospective independent pollinator-channel gate (N1)** — replace the old future-looking H5 section with the result that N1 was executed and failed; N2 remained closed.
9. Robustness/missingness and stopping rule.

Historical code/module names containing `when_where` can remain for reproducibility. The manuscript prose should not use “when” as a biological axis.

---

## 7. Results order for v8

### Result 1 — one universal syndrome is rejected

Lead with the direct multivariate H1 result.

### Result 2 — response direction branches by biogeographic context

Palearctic versus Tropical is the clearest ecological contrast. Keep effect sizes and intervals/FDR, not just support counts.

### Result 3 — the Palearctic primary response resolves to genus-level assembly

This is the decisive H3 result. Use effect attenuation/interval plots where possible; keep `4/4 → 4/4 → 0/4` only as a compact robustness summary, not as the main evidential language.

### Result 4 — hierarchical expression itself is context dependent, but bounded

Present the new matched interaction audit here. The biological message is not “Tropical = evolution”. It is:

> **Genus adjustment changes the same plant-architecture response differently across biogeographic contexts in widespread native assemblages.**

Then state immediately that all-native support is weaker, preventing a universal taxonomic-depth claim.

### Result 5 — default alternatives do not earn mechanism promotion

Put H4 area, nonlinear geometry, and N1 in one coherent falsification section:

- area amplification does not pass the heteroskedastic-null gate;
- global nonlinear geometry promotes 0/12 steps/reversals;
- independent pollinator-channel heterogeneity p=0.655 and N2 remains closed.

This creates a strong narrative: the paper does not reach its conclusion by ignoring standard alternatives; it tests and bounds them.

### Result 6 — secondary pollination-associated architecture

Retain V4 common factor and residual results as ecological interpretation, not as pollinator identity.

---

## 8. Figure architecture

### Figure 1 — conceptual hypothesis tree

Title inside figure:

> **Where, whether, and at what assembly level does an island floral response emerge?**

Left-to-right flow:

`source separation / connectivity`
→ `regional plant source availability`
→ `island assemblage`
→ `floral + reproductive response`

Overlay four questions rather than four mechanistic arrows:

- **Whether?** universal vector vs context-dependent branching.
- **Where?** analysis regime / realm.
- **At what level?** observed → family-adjusted → genus-adjusted → residual.
- **How?** mechanism deferred to independent pollination evidence / Chapter 2.

Add a small side box:

`response geometry: flat / cline / step / reversal`
→ `global audit: no promoted nonlinear transition`

Do not draw a solid causal arrow `distance → Bombus loss → floral trait`.

### Figure 2 — response direction

Effect-size forest/arrow plot for Palearctic and Tropical primary axes. Show all-analysis and direct-only side by side. This is the H1/H2 figure.

### Figure 3 — hierarchical attenuation

Plot effect estimates with intervals at observed / after-family / after-genus stages rather than only support votes. Facet by biogeographic context. Native non-endemic should be visible because the matched stage × context signal is strongest there.

### Figure 4 — source/lineage assembly details

Genus entry/loading or equivalent H3 decomposition, whichever best shows what “genus-level assembly” means biologically.

### Figure 5 — falsification / claim-ceiling panel

Compact three-part panel:

- area: heteroskedastic-null gate 0/16;
- response geometry: identified step/reversal 0/12;
- independent pollinator N1: p=0.655, stopped before N2.

This can be Supplement if the main journal has a strict figure limit, but the logic should remain in the main text.

Secondary syndrome-template common-factor material can move to a later main figure or Supplement depending space.

---

## 9. Discussion order

### Paragraph 1 — general result

Do not open with Bombus. Open with hierarchy:

> Isolation-associated ecological responses are not defined only by effect direction. In the global island flora, both response direction and the hierarchical level at which response is represented depend on biogeographic context.

Immediately bound this by saying the direct matched depth result is strongest in native non-endemics and weaker in all-native assemblages.

### Paragraph 2 — island syndrome as assembly outcome

Explain Palearctic genus-level assembly. Emphasise that lineage composition is not a nuisance covariate but part of the biological response.

### Paragraph 3 — why Tropical matters

The Tropical assurance↑ / accessibility↓ combination is the cleanest contradiction of a single opportunistic-generalisation pathway. Discuss possible retained/redirected specialised floral architecture without assigning visitor identity.

### Paragraph 4 — what the negative alternatives teach us

H4, geometry and N1 all fail to earn promotion. This is not three “nulls” to hide. It narrows what can be claimed:

- no universal area/capacity pathway;
- no global step/reversal geometry;
- no independently demonstrated channel-specific isolation response.

### Paragraph 5 — double-filter hypothesis remains a framework, not a demonstrated mechanism

Isolation can in principle filter both plant source pools and interaction partners. Chapter 1 directly demonstrates the plant/source/lineage side. The pollinator side remains compatible but unverified.

### Paragraph 6 — Chapter 2 handoff

Use the old `izu-core` insight explicitly:

> Global data did not promote a shared nonlinear transition. At the within-system scale, however, different biological response channels can still follow different geometries. The Izu design therefore tests whether continuous floral/outcrossing change and threshold reproductive assurance co-localise with directly measured pollination-regime transitions.

This is the cleanest Chapter 1 → Chapter 2 bridge.

---

## 10. What changes in the EL argument

The EL argument should not be:

> “island syndrome = genus assembly.”

That is strong but field-specific.

The broader argument is:

> **The same broad environmental gradient can generate different response directions and different hierarchical expressions depending on biogeographic context.**

The global island system is useful because the same exposure, source controls, trait ontology and hierarchical decomposition can be applied across thousands of islands.

The current data support this only with boundaries:

- primary direction heterogeneity is strong;
- Palearctic genus-level assembly is strong;
- matched context × stage heterogeneity is strong in native non-endemics but weaker in all-native assemblages;
- endemic strata are not promoted as an evolutionary time axis;
- global nonlinear geometry is negative;
- pollinator mechanism is not identified.

Therefore **EL first remains a rational high-risk/high-reward submission**, with GEB as the natural fallback. Do not inflate the abstract to claim a universal hierarchical-depth law.

---

## 11. Required v7 → v8 textual corrections

1. Replace all biological uses of “when and where” with “where, whether, and at what assembly level” or equivalent prose. Historical module names remain unchanged.
2. Use Wave52, not Wave36, as the canonical H3 depth result.
3. Update V4 shared-factor values to Wave52 values: **86.85% all-analysis; 86.44% direct-only**.
4. Add the bounded matched context × stage audit after H3.
5. Add the calibrated nonlinear-geometry negative result: **12/12 monotonic_or_unresolved**.
6. Replace future-looking H5 language with the actual prospective N1 result and explicit stop before N2.
7. Keep area as a failed mechanism-promotion gate, not “area absent”.
8. Keep syndrome templates as plant-architecture concordance; never infer realised pollinator guild.
9. Do not use endemicity as a temporal proxy or as proof of in-situ evolution.
10. Prefer effect sizes and intervals in figures/text; retain mode counts only as robustness shorthand.

---

## 12. Chapter 1 / Chapter 2 division after v8

### Chapter 1 — `island`

**Where / whether / at what assembly level?**

- global response direction;
- biogeographic contingency;
- source and lineage assembly;
- bounded hierarchical-depth differences;
- explicit failures of area, nonlinear geometry and independent channel heterogeneity to earn mechanism promotion.

### Chapter 2 — `izu-core`

**How / why within a biological system?**

- direct pollinator regime;
- trait matching and effective pollen transfer;
- reproductive outcome;
- response geometry across biological channels;
- cline versus threshold under direct functional evidence.

The two chapters therefore share the response-geometry idea without requiring the same empirical geometry at both scales: the global chapter tests whether a common nonlinear form generalises and finds that it does not; the focal-system chapter tests how specific functional channels can nevertheless generate different response shapes.
