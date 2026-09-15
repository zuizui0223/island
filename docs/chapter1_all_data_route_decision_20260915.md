# Chapter 1 all-data route — integrated decision, 2026-09-15

Status: **exploratory route resolved enough for a claim hierarchy; do not replace frozen v11 native-assembly manuscript yet.**

## 1. What the expanded route was designed to test

The expansion separates three previously entangled restrictions:

1. trait-evidence quality;
2. probability-model form;
3. floristic-status identity.

The broad data layer uses every trait-evaluable observed island flora rather than shrinking the main pattern analysis to source-backed native-status records.

Evidence hierarchy:

- primary: `all_analysis_eligible` = High + Medium + trait-specific Validated Low;
- sensitivity: `direct_only` = High + Medium.

Primary probability model:

- beta-binomial logit regression on island-level trait counts;
- six atomic outcomes;
- no arbitrary minimum species-per-island cutoff;
- distance + island area + climate PC1--PC4;
- spatial-block cluster-robust covariance;
- direct North--Tropical joint Wald test of distance x context interactions.

## 2. Broad observed-flora pattern is real as an observational result

Primary beta-binomial run `34961775336`, artifact `10394245237`:

- all-analysis: **3,626 islands / 187 blocks**, North--Tropical p = **0.0007633**;
- direct High/Medium: **3,569 islands / 187 blocks**, p = **0.004502**.

Matched grouped-binomial model-form sensitivity run `34963020944`, artifact `10393863668`:

- all-analysis: **p = 0.001385**;
- direct: **p = 0.004501**.

Therefore the broad multivariate difference is neither a Validated-Low artefact nor a beta-binomial-likelihood artefact.

## 3. The old two-axis score fails for a construction reason, not simply because six dimensions were reduced to two

Historical species-level composite route on the expanded all-observed flora:

- all-analysis: p = **0.73998**;
- direct: p = **0.06142**.

But a linear projection of the same six beta-binomial interaction estimates onto two predeclared families retains the difference.

Response-compression audit run `34964112334`, artifact `10394930700`:

### all-analysis

- full six-axis vector: **p = 0.0007633**;
- two family-level coefficient contrasts: **p = 0.0016467**;
- discarded within-family contrasts: **p = 0.5085**;
- Euclidean signal retained by the two-family subspace: **78.22%**.

### direct-only

- full six-axis vector: **p = 0.0045020**;
- two family-level coefficient contrasts: **p = 0.0120921**;
- discarded within-family contrasts: **p = 0.1438**;
- retained Euclidean signal: **64.82%**.

The minimum-informative-trait rule is not the explanation either. Syndrome support-intersection audit run `34964326934`, artifact `10394368256`:

- historical `minimum_informative_traits = 2`: p = **0.73998** all-analysis / **0.06142** direct;
- relaxed `minimum_informative_traits = 1`: p = **0.50197** all-analysis / **0.85119** direct.

A separate trait-first family refit also remains weak rather than reproducing the coefficient-space projection. Run `34964919955`, artifact `10394986981`:

- family averages using available components: p = **0.50765** all-analysis / **0.25947** direct;
- requiring complete three-component support: p = **0.11072** all-analysis / **0.06374** direct.

Decision: **the primary estimand should remain the six-atomic multivariate slope vector.** Accessibility/generalisation and reproductive assurance can be shown as coefficient-level linear summaries of that fitted vector, but the historical species-level weighted concordance scores should not define the primary expanded response.

## 4. Floristic status is the biological claim boundary

On source-backed native records only, the six-atomic North--Tropical difference is not supported:

- all native: p = **0.2491** all-analysis / **0.3018** direct;
- native non-endemic: p = **0.05097** all-analysis / **0.2079** direct.

Holding the island frame fixed at 375 status-supported islands does not rescue this as a power explanation:

- all observed species on those islands: **p = 0.0005607** all-analysis / **0.0003959** direct;
- native species only on the same islands: **p = 0.2491** / **0.3018**.

The outcome-free status-resolution audit run `34963956873`, artifact `10394441857`, shows that origin-status resolution is geographically nonrandom within contexts, but the North--Tropical difference in the distance slope of status resolution is absent (**p = 0.849**). Thus the discrepancy is not explained by a simple context difference in whether status was recorded.

## 5. WCVP greatly expands regional-native-compatible coverage, but it does not restore the North--Tropical difference

WCVP regional-native compatibility run `34968461262`, artifact `10396382401`.

WCVP coverage:

- fixed island-flora species target: **115,328**;
- exact accepted WCVP species: **101,101**;
- species with a non-doubtful extant native TDWG-L3 range: **100,887**.

For unresolved source-status rows, an observation was admitted as *regional-native-compatible* only when the focal island had an unambiguous accepted TDWG-L3 mapping and that region was included in the species' WCVP native range. This is a compatibility sensitivity, not exact island-native proof.

Coverage gain:

- unresolved rows upgraded to regional-native-compatible: **362,789**;
- islands receiving at least one upgrade: **2,281**;
- resulting regional-native-compatible rows: **513,320**;
- resulting islands: **2,372**.

North--Tropical six-atomic test after this expansion:

- all-analysis: **1,848 analysis islands / 178 blocks**, p = **0.07954**;
- direct: **1,827 islands / 178 blocks**, p = **0.17472**.

So the broad all-observed difference still does **not** transfer to the much larger regional-native-compatible flora.

## 6. Partition audit locates where the broad observed-flora contrast resides

The same WCVP run partitioned source-status unresolved records without relabelling them as introduced.

### all-analysis evidence

- source-backed native: p = **0.24915**;
- WCVP-compatible unresolved: p = **0.06562**;
- WCVP-incompatible unresolved: p = **0.002010**;
- WCVP-unclassifiable unresolved: p = **0.01368**;
- source-introduced + WCVP-incompatible unresolved: p = **0.001102**;
- regional-native-compatible combined: p = **0.07954**.

### direct High/Medium

- source-backed native: p = **0.30184**;
- WCVP-compatible unresolved: p = **0.35251**;
- WCVP-incompatible unresolved: p = **0.003065**;
- WCVP-unclassifiable unresolved: p = **0.002687**;
- source-introduced + WCVP-incompatible unresolved: p = **0.009484**;
- regional-native-compatible combined: p = **0.17472**.

Interpretation: the broad all-observed North--Tropical contrast is concentrated in records that are **not currently defensible as regional-native-compatible**. This does not prove that every incompatible or unclassifiable record is introduced; it does mean that the broad global contrast cannot currently be promoted as native island assembly.

## 7. Claim hierarchy going forward

### Level A — broad observational pattern

Supported:

> Across thousands of contemporary observed island floras, isolation-associated floral and reproductive composition differs among biogeographic contexts at the six-atomic-trait level.

This is the correct use of the all-data analysis.

### Level B — evidence and model robustness

Supported:

- High/Medium-only direct evidence retains the broad result;
- grouped-binomial and beta-binomial models agree;
- no arbitrary 30-species-per-island threshold is required.

### Level C — native island-assembly interpretation

**Not promoted globally.** The North--Tropical six-atomic contrast weakens when restricted either to source-backed native taxa or to the much larger WCVP regional-native-compatible set.

The frozen v11 Palearctic genus-structured native result therefore remains the current defended native-assembly core. The all-data route complements it; it does not supersede it.

## 8. Paper/poster consequence

For Q1, use a layered presentation rather than calling all 3,626 islands a native island-syndrome test:

1. **global observed-flora pattern:** thousands of islands, six-atomic probability vector;
2. **direct-evidence/model sensitivities:** result survives;
3. **native-status gate:** global context contrast attenuates;
4. **defended assembly result:** Palearctic native response is genus-structured in frozen v11.

This creates a stronger and more accurate message:

> **Global island floras show context-dependent floral/reproductive composition, but the part that can currently be defended as native island assembly is narrower and lineage-structured.**

The next data-improvement target is no longer trait coverage. It is island-specific floristic origin/status resolution for the records that are WCVP-unclassifiable or regionally incompatible.
