# Chapter 1 all-data probability result — 2026-09-15

Status: **exploratory expansion completed; do not replace the frozen v11 submission result yet**.

## Execution provenance

Primary probability run:

- branch: `ch1-all-data-primary`;
- workflow: `Run Chapter 1 all-data probability analysis`;
- successful run: **34961775336**;
- artifact: `chapter1-all-data-probability-34961775336`;
- artifact ID: **10394245237**;
- digest: `sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b`;
- frozen upstream Chapter 1 progressive artifact: run `34232450884`.

All fitted beta-binomial optimizers converged.

## Analysis hierarchy

### Primary probability model

- beta-binomial logit regression on island-level trait counts;
- six predeclared atomic outcomes;
- island-specific denominator = number of species with a resolved value for that outcome;
- **no arbitrary minimum species-per-island cutoff**;
- mean model: distance + island area + climate PC1--PC4;
- spatial-block cluster-robust sandwich covariance;
- North--Tropical difference = direct joint Wald test of distance x context coefficients.

### Evidence

1. `all_analysis_eligible` = High + Medium + trait-specific Validated Low, primary;
2. `direct_only` = High + Medium, evidence-quality sensitivity.

### Flora

1. `all_observed` = broad contemporary observed-flora composition;
2. `all_native` = status-resolved sensitivity;
3. `native_nonendemic` = stricter status-resolved sensitivity.

## 1. Expanded all-observed result

### All-analysis evidence

Within-context six-outcome vector support:

- northern mid-latitude: **2,164 islands**, 86 blocks, p = **1.06e-6**;
- northern high latitude: **409 islands**, 63 blocks, p = **1.10e-6**;
- tropical: **1,462 islands**, 110 blocks, p = **1.56e-5**;
- southern extratropical: **299 islands**, 41 blocks, p = **2.09e-13**.

Formal North--Tropical interaction-vector test:

- **3,626 unique islands**;
- 187 spatial blocks;
- chi-square = **23.1011**, df = 6;
- **p = 0.0007633**.

Thus the broad observed-flora North--Tropical response difference is supported under the beta-binomial model.

### High/Medium direct sensitivity

The formal North--Tropical test uses **3,569 islands**, 187 blocks and remains supported:

- chi-square = **18.8072**, df = 6;
- **p = 0.004502**.

The broad observed-flora result is therefore not created solely by Validated-Low trait evidence.

## 2. Model-form sensitivity

A matched grouped-binomial analysis retained the same six atomic outcomes, island sets, covariates and spatial-block robust inference.

Provenance:

- run: **34963020944**;
- artifact ID: **10393863668**;
- digest: `sha256:dbd5b20b09659404470c5e706449fe90bb870878aed3fc3716037b7fedcd6507`.

North--Tropical six-outcome test:

- all-analysis: **p = 0.001385**;
- direct-only: **p = 0.004501**.

The broad result is therefore **not a beta-binomial artefact**.

## 3. Floristic-status restriction changes the biological interpretation

The all-observed result does not transfer cleanly to status-resolved native flora.

All-analysis evidence:

- `all_observed`: **p = 0.000763**;
- `all_native`: **p = 0.2491**;
- `native_nonendemic`: **p = 0.05097**.

Direct-only:

- `all_observed`: **p = 0.004502**;
- `all_native`: **p = 0.3018**;
- `native_nonendemic`: **p = 0.2079**.

This is not merely loss of islands. On the same 375 status-supported islands:

- all observed species: **p = 0.0005607** all-analysis / **0.0003959** direct;
- native species only: **p = 0.2491** all-analysis / **0.3018** direct.

Unresolved-status species alone remain informative on that same frame:

- all-analysis: **p = 0.01914**;
- direct-only: **p = 0.03044**.

Introduced-only data do not satisfy the 50-island-per-outcome six-vector support rule. These diagnostics do not imply `unresolved = introduced`; they show that origin-status uncertainty materially changes the inference.

## 4. Status-resolution observation process

A separate outcome-free audit tested whether origin-status *resolution itself* is geographically structured.

Provenance:

- workflow run: **34963956873**;
- artifact ID: **10394441857**;
- digest: `sha256:09d0bd8c7d0d69de8e716ac389e54d544405f837e4642f405ff11277c2921586`.

Raw resolved fractions of observed island-species rows:

- northern mid-latitude: **0.1403**;
- northern high latitude: **0.0112**;
- tropical: **0.2206**;
- southern extratropical: **0.1882**.

Status resolution increases with distance within each fitted regime, but the direct North--Tropical difference in that distance slope is unsupported:

- 3,748 islands / 190 blocks;
- interaction estimate `+0.0996`;
- SE `0.5231`;
- **p = 0.8490**.

Thus the native/all-observed discrepancy cannot be reduced to a simple North--Tropical difference in the distance gradient of status availability. The remaining issue is which unresolved species are actually native, not merely which islands have more status information.

## 5. Response-definition audit: dimensionality is not the main problem

The historical species-level two-axis concordance route (`generalized_accessible + selfing_core`) is weak on the expanded all-observed flora:

- all-analysis: p = **0.73998**, q = **0.88797**;
- direct-only: p = **0.06142**, q = **0.13094**.

That initially suggested that reducing six traits to two dimensions destroyed the signal. A matched linear-contrast audit of the **same six beta-binomial interaction estimates and their full covariance** rejects that explanation.

Provenance:

- response-compression run: **34964112334**;
- artifact ID: **10394930700**;
- digest: `sha256:960802ffc47a2287849693fc37788466c6778beae6388823453bf39eb51782fc`.

### All-analysis evidence

- full six-axis vector: **p = 0.0007633**;
- two equal-weight family means: **p = 0.0016467**;
- four within-family contrasts discarded by the two-family projection: **p = 0.5085**;
- Euclidean signal retained by the two-family subspace: **78.22%**.

### Direct-only

- full six-axis vector: **p = 0.0045020**;
- two family means: **p = 0.0120921**;
- discarded within-family contrasts: **p = 0.1438**;
- retained Euclidean signal: **64.82%**.

Therefore **two-dimensionality itself is not the explanation**. A simple two-family projection of the atomic slopes remains supported. The discrepancy lies in the historical species-level concordance construction and/or its support weighting.

The historical composite route differs because it:

- first creates a multivariate concordance score for each species;
- requires at least **two informative traits per species**;
- ignores missing traits in the within-species denominator;
- averages the resulting species scores within each island;
- then fits a continuous island-score model.

The atomic probability route instead uses every species resolved for each individual trait outcome and carries the outcome-specific species denominator into the likelihood.

A dedicated support-intersection audit is therefore the next diagnostic; the specific test is `minimum_informative_traits = 2` versus 1 with all other syndrome definitions held fixed.

## Decision

The expanded route now establishes that:

1. a broad contemporary observed-flora North--Tropical response difference is estimable across thousands of islands;
2. it survives High/Medium-only trait evidence;
3. it survives matched grouped-binomial model-form sensitivity;
4. the difference can still be represented by two family-level linear contrasts, so the old composite failure is **not** simply a six-versus-two dimensionality problem;
5. floristic-status identity remains the main biological claim boundary.

The current broad claim is therefore:

> **Across thousands of observed island floras, isolation-associated floral and reproductive composition differs among biogeographic contexts. This pattern is robust to trait-evidence quality and probability-model form, but its interpretation as native island assembly is limited by unresolved floristic status.**

Do not promote the all-observed result to native colonisation filtering, in-situ evolution or historical pollinator-driven adaptation. Frozen v11 remains canonical until the status-identity problem and the species-level syndrome-support discrepancy are resolved.
