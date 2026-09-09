# Chapter 1 submission hypothesis framework — H1–H5

## Purpose

This is the paper-facing hypothesis framework for the current Chapter 1 manuscript. It preserves the frozen PR142 analysis contract while clarifying the causal hierarchy implied by the double-geographic-filter interpretation.

The paper should present H1–H5 as a sequence of increasingly specific ecological questions:

1. **H1 — Is there one universal floral island syndrome?**
2. **H2 — If not, where do source-distance responses branch?**
3. **H3 — At what plant-assemblage depth is the supported response generated?**
4. **H4 — Does island area modify the strength of that filtering?**
5. **H5 — What pollination-channel mechanism could generate the remaining or lineage-level filter?**

The analysis itself remains the frozen PR142 pipeline. The H5a/H5b distinction below is a paper-level causal interpretation and prospective extension; it does not retroactively alter the fitted PR142 models.

---

## Central geographic exposure

The common exposure is increasing geographic separation, implemented formally as `log1p_distance_to_continent_km` and interpreted as a composite gradient of separation, connectivity and source accessibility.

The ecological limiting picture is:

```text
regional source system
        |
      d -> 0
        |
near-source islands
        |
intermediate separation
        |
remote islands
        |
stronger opportunity for differential filtering
```

Distance is not treated as a direct causal force on floral phenotype.

---

## H1 — Universal floral island syndrome

### Question

Does increasing source separation push island floras along one coherent floral/reproductive response vector across biogeographic contexts?

### Prediction under H1

A universal syndrome predicts similar multivariate responses across regions, for example increased reproductive assurance together with increased accessibility/generalization or loss of specialized floral architecture.

### Required evidence

- within-context multivariate response support;
- direct between-context vector comparison;
- no inference of heterogeneity from `significant here / nonsignificant there` alone.

### Current result

**Not supported.** The response to distance differs among contexts.

### Biological meaning

The island floral syndrome is not one deterministic global rule.

---

## H2 — Biogeographic branching

### Question

If H1 fails, do different biogeographic contexts follow different source-distance trajectories?

### Primary response

The primary response is deliberately pollinator-name-free:

- `accessibility_generalization`;
- `reproductive_assurance`.

### Current result

**Supported at the pattern level.**

Palearctic source separation is associated with increasing accessibility/generalization and reproductive assurance.

Tropical source separation can instead combine increasing reproductive assurance with decreasing accessibility/generalization, i.e. maintenance or strengthening of specialized/attractive architecture.

### Biological meaning

The classical island syndrome decomposes into at least two partially separable components:

```text
source separation
      |\
      | \
      |  -> reproductive-assurance component
      |
      -> floral-architecture component
```

The direction and coupling of these two branches depend on context.

### Claim ceiling

The regional pattern is supported, but measured-climate-independent causal effects of realm identity are not established.

---

## Pollination-associated concordance — interpretation bridge, not H6

After the H2 plant-side response is established, fixed `large_bee_like`, `butterfly_like`, and `bird_like` floral templates are evaluated.

Their role is:

```text
H2 plant response
      -> floral-architecture concordance
      -> candidate functional-channel hypothesis
      -> independent H5 test
```

They are not pollinator classifiers.

Current pattern:

- Palearctic: large-bee-like and butterfly-like architecture decline with isolation;
- Tropical: all three sampled templates increase;
- approximately 87% of the template variation is carried by one shared plant-architecture factor.

Therefore template concordance constrains interpretation but cannot identify pollinator retention, loss, mobility, visitation, effectiveness or replacement.

---

## H3 — Source and lineage assembly

### Question

At what level of the plant assemblage is the H2 distance response represented?

### Formal decomposition

```text
observed response
      -> after family composition
      -> after source-matched genus composition
      -> beyond-genus residual
```

### Current Palearctic result

`observed 4/4 -> after family 4/4 -> after genus 0/4`.

### Biological meaning

The broad Palearctic floral/reproductive response is compatible with differential representation of source-available genera.

H3 therefore identifies **where in the plant assemblage hierarchy the filter is expressed**.

### What H3 does not identify

Genus-level assembly does not prove plant propagule dispersal alone caused the filter.

Differential genus representation may arise from:

- plant dispersal or propagule supply;
- establishment and persistence;
- habitat filtering;
- demographic history;
- interaction dependence;
- pollination-channel retention or loss;
- combinations of these.

Thus disappearance after genus adjustment is not evidence against a pollination mechanism.

---

## H4 — Area / capacity moderation

### Question

Does island area modify the strength of the source-distance filter?

### Formal form

`response ~ distance + area + distance × area`

with continuous area and no post-hoc small/large island threshold.

### Current result

All 16 frozen primary V3 classifications remain `retain_area_as_measurement_sensitive_modifier_only`; zero heteroskedastic-null promotion gates pass.

### Biological meaning

Area can modify the apparent magnitude or precision of filtering, but the current data do not identify whether this reflects founder filtering, habitat capacity, demographic persistence or pollinator persistence.

H4 is therefore a boundary condition, not a closed mechanism.

---

## H5 — Pollination-channel mechanism

### Core question

Can the pollinator side of the geographic barrier explain why source-distance filtering produces the observed plant assemblage response?

Pollinator mobility belongs upstream of channel state:

`E_poll(g,d) = f(source availability, distance, flight/dispersal ability, ocean crossing, establishment, habitat, realized community)`.

The same geographic distance can therefore generate different channel states for different functional pollinator groups.

### H5a — channel-dependent lineage filtering

H5a asks whether pollinator-side geography helps generate H3 itself.

```text
source pollination channel
        -> retention / disruption with distance
        -> effective service available to colonists
        -> differential plant establishment / persistence
        -> genus entry / representation
        -> H3 assemblage pattern
```

The future estimand is lineage entry/loading conditional on independently measured channel state and lineage functional dependency.

A positive H5a result would show that a pollination-channel filter is one mechanism producing genus-level assembly.

### H5b — beyond-genus residual channel mechanism

H5b corresponds to the original strict PR142 residual gate.

```text
source / genus composition fixed
        -> retained / disrupted pollination channel
        -> visitation
        -> single-visit effectiveness
        -> effective service
        -> residual floral/reproductive response
```

A positive H5b result would identify pollinator-side effects that remain after plant lineage composition is accounted for.

### Current status

**Not established in Chapter 1.**

The plant database does not contain the complete independent source-channel -> retention/disruption -> visitation -> single-visit effectiveness -> effective-service chain required for causal promotion.

### Important logical relation

`H5b null` does not falsify `H5a`.

A pollination channel may affect which genera establish or persist, in which case its effect is absorbed by H3 genus composition and no beyond-genus residual is expected.

---

## One-line interpretation of each hypothesis

| Hypothesis | Paper-facing question | Current answer |
|---|---|---|
| H1 | Is there one universal syndrome? | No |
| H2 | Where do responses branch? | Biogeographically/contextually |
| H3 | Where in the plant assemblage is the main filter expressed? | Genus-level assembly in the Palearctic broad response |
| H4 | Under what island capacity conditions is filtering stronger? | Area modifies, but mechanism remains unidentified |
| H5 | Why are those lineages / residual responses filtered? | Pollination-channel mechanism remains prospective |

---

## Final synthesis

The current Chapter 1 result is best stated as:

> Increasing source separation does not impose one universal floral island syndrome. It produces context-dependent assemblage trajectories in which reproductive assurance and pollination-associated floral architecture can decouple. The strongest Palearctic response is expressed primarily through source-matched genus assembly. This identifies the plant-side level of filtering but not its ultimate cause. Pollinator-channel mobility and retention remain candidate upstream mechanisms that could generate lineage assembly itself, while independent effective-service data are required to test any beyond-genus pollinator mechanism.
