# Island floral syndrome — v2

This repository supports **Chapter 1 / macroecological when-and-where analysis** of island floral and reproductive filtering.

The current Chapter 1 question is:

> **When and where is isolation-associated floral/reproductive filtering detectable, and where does the multivariate response to isolation differ?**

The primary result is now a **response-vector** result, not a count of significant atomic traits and not a Bombus-deficit model.

> **Isolation-associated filtering is confirmatorily detectable in both northern mid-latitude and tropical island floras, persists in native non-endemic assemblages, and the northern-versus-tropical multivariate isolation-response vectors differ significantly. Current data are insufficient for equivalent confirmatory conclusions in northern high-latitude or southern-extratropical floras.**

See:

- [`docs/chapter1_when_where_frozen_result_20260825.md`](docs/chapter1_when_where_frozen_result_20260825.md) — canonical when/where result;
- [`THESIS_CHAPTER_POSITIONING.md`](THESIS_CHAPTER_POSITIONING.md) — dissertation role and claim boundary;
- [`docs/chapter1_pollination_syndrome_concordance_20260825.md`](docs/chapter1_pollination_syndrome_concordance_20260825.md) — post-freeze Discussion-only pollination-syndrome audit.

## Canonical result

Formal workflow run: `32837335384`.

### WHERE

Confirmatory response-vector tests detect isolation-associated filtering in:

- **northern mid-latitude** islands;
- **tropical** islands.

The same result is recovered in `all_native` and `native_nonendemic` flora.

### BETWEEN-WHERE

Using 17 atomic responses with confirmatory support in both regions, northern-mid-latitude and tropical isolation-response vectors differ jointly:

- all native: joint Wald χ² = **69.78**, df = 17, q = **2.35e-08**;
- native non-endemic: joint Wald χ² = **61.02**, df = 17, q = **7.13e-07**.

Thus the current evidence supports **biogeographically structured isolation filtering between northern mid-latitude and tropical island floras**.

### WHEN

The filtering signal persists within **native non-endemic** flora in both northern mid-latitude and tropical contexts. It therefore is not confined to endemic-lineage turnover.

Endemic-only vector tests do not currently meet the confirmatory support requirement, so endemicity is not claimed to be either necessary or irrelevant.

Northern high-latitude and southern-extratropical floras are currently **unresolved confirmatorily**, not demonstrated null regions. At the 30-island pilot threshold, southern-extratropical vectors are nonzero but are supported by only a few atomic responses.

## Canonical analysis surface

The primary implementation is:

- `src/island_v2/chapter1_when_where_omnibus.py`
- `config/chapter1_when_where_omnibus.yml`
- `.github/workflows/run-chapter1-context-main.yml`

The inference sequence is:

```text
frozen island universe
→ direct trait evidence
→ floristic-status strata
→ area + climate + isolation + frozen biogeographic context
→ WHERE: within-context response-vector joint test
→ BETWEEN-WHERE: pairwise response-vector difference test
→ WHEN: persistence across floristic-status strata
→ lineage / genus-composition guardrail
→ atomic trait decomposition
→ pollination-syndrome interpretation only after freeze
```

Atomic category models remain important, but only to explain **what creates the supported regional vectors** after when/where is established.

## Evidence tiers

Complete fill is not complete evidence. Trait resolution remains visible.

- **Confirmatory:** direct source-backed species evidence; vector responses require >=50 contributing islands per retained atomic response.
- **Pilot:** >=30 and <50 islands per retained response.
- **Secondary robustness:** taxonomic inference at genus/family level.
- **Sensitivity only:** global fallback or other declared inference layers.

## Bombus and pollination-syndrome boundary

Bombus, bird, butterfly, moth, hawkmoth and other pollinator labels do **not** enter the primary when/where model.

Existing Bombus applicability, environmental-niche and occurrence-diagnostic products remain provenance-preserving diagnostic/sensitivity products. Climatic compatibility is not realized occurrence; occurrence is not visitation or service; opportunistic non-detection is not historical loss.

After the when/where result is frozen, the northern-versus-tropical vector difference may be compared with pollination-syndrome literature for **concordance or mismatch only**. Trait vectors do not identify a causal pollinator guild.

## Lineage guardrail

The genus-composition-preserving M3 layer remains mandatory before broad syndrome interpretation. It does not recover a coherent confirmatory `generalized + plain + SC` syndrome in the main northern/tropical strata.

Therefore Chapter 1 can conclude that regional isolation-response vectors differ without naming either vector a classical pollination syndrome.

## Thesis handoff

```text
Chapter 1 — island
WHEN / WHERE is isolation filtering detectable?
WHERE do multivariate responses differ?
        ↓
Chapter 2 — izu-core
WHY do those contexts generate different response architectures?
        ↓
Chapter 3 — shimahotarubukuro
WHAT phenotype axes actually diverge within one focal lineage?
```

## Repository layout

- `src/island_v2/` — reusable v2 data and analysis utilities
- `analysis/v2/` — statistical analysis scripts and execution notes
- `config/` — frozen contracts, ontology, and artifact locks
- `data/v2/` — external/staging/curated/template data layers
- `docs/` — scientific design, data policy, methods, and reproducibility notes
- `.github/workflows/` — active validation and canonical analysis workflows
- `legacy/v1/` — frozen v1 provenance only

## Reproducibility rule before submission

GitHub Actions artifacts are temporary and are not a permanent supplement. The manuscript release must archive all critical inputs and outputs durably with checksums, identify one canonical workflow per main analysis, and report attrition from the frozen 8,265-island universe to every fitted model.
