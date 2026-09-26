> **HISTORICAL / SUPERSEDED — pre-corrected Chapter 1 surface.** Retained for provenance/replay only. The current submission is selected by `config/chapter1_submission_current.json`; use `submission/chapter1_current/MANUSCRIPT.md` and `docs/PAPER_PIPELINE.md` for current results.

# Chapter 1 submission freeze — P1-defended — 2026-09-15

Status: current submission-facing surface after P1 assembly-depth defense. Historical v8 and 2026-09-14 freezes remain provenance records.

## Canonical manuscript

`docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md`

Title:

> **A floral island syndrome can emerge from hierarchical lineage assembly: biogeographic contingency across 8,265 islands**

The manuscript begins from the island question—why island flowers are expected to become generalized or self-reliant—and reaches the broader trait-syndrome inference problem only after the island evidence is resolved.

## Canonical publication-facing claim

> **Geographic isolation does not impose one floral island syndrome. The strongest Palearctic floral/reproductive response is strongly structured by real genus composition: true genus boundaries attenuate the response more than arbitrary within-family partitions of identical grouping complexity, while the exact incremental family-to-genus attenuation remains spatially imprecise.**

General conceptual statement:

> **An apparent macroecological trait syndrome can be carried by non-random lineage assembly rather than one repeated organismal response; taxonomic localization and its uncertainty should be established before mechanism is assigned.**

This statement arises from the floral-island problem and is not used as a top-down premise.

## P1 final defense

Integrated lock:

- `config/chapter1_p1_final_decision_result_lock.json`

### P1a — same-support safeguard

- run `34935183075`;
- artifact `10382749051`;
- digest `sha256:4c4689d5d8e0eaca7e3ac66e1e866a48942b2f7220e962f7b1bb7e7ac1d43af2`;
- verdict: observed/family/genus stages use the same focal observations and information weights; stage-specific sample loss is not a viable attenuation explanation.

### P1c — matched-complexity genus null

Frozen permutations:

- source run `34936193944`;
- 40/40 shards successful;
- 2,000/2,000 valid permutations;
- no permutations regenerated after aggregation wiring failure.

Final aggregate:

- run `34941827774`;
- artifact `10385820775`;
- digest `sha256:861d18fe87b9190619f925a2446be5fd4d460b818825578930883257c0a6ed16`;
- true-genus statistic `0.7206615`;
- matched pseudo-genus null median `0.2243375`;
- 57/2,000 null values ≥ observed;
- one-sided randomization `p=0.0289855`;
- verdict: `true_genus_exceeds_matched_complexity_null`.

### P1d — paired spatial uncertainty

- run `34939113182`;
- artifact `10383754545`;
- digest `sha256:822c69e9a1be4eb2c340391c1fb84dfba91b45b8a765549bb8ad1ec6379ef480`;
- 2,000 paired spatial-block draws;
- total genus attenuation remains large;
- 0/8 direct-only primary profiles have a 95% interval for the *additional* family→genus attenuation entirely above zero.

Interpretation: genus-specific taxonomic structure is supported; a precisely estimated family→genus increment is not.

## Locked main figures

### Figure 1

`config/chapter1_v8_figure1_result_lock.json`

Hierarchy-of-assembly conceptual map. Its role is unchanged; v9 legend language no longer presents the family→genus increment as precisely estimated.

### Figure 2

`config/chapter1_v8_figure2_result_lock.json`

Primary biogeographic branching and component-decoupling figure.

### Figure 3 — replaced by P1-defended version

`config/chapter1_v9_figure3_p1_defense_result_lock.json`

- run `34942779753`;
- artifact `10386555245`;
- digest `sha256:93957ed4e12193fd2630aa3c274a99f004af4c8e40a6f63985e46a77c4cdc239`;
- visual review passed;
- includes both favourable matched-null evidence and adverse paired-uncertainty evidence.

The old v8 Figure 3 lock is historical and no longer defines the submission Figure 3.

### Figure 4

`config/chapter1_v8_figure4_result_lock.json`

Observation, response-geometry and H5 claim-boundary figure.

## Existing robustness retained

- V5 trait-resolution MNAR: Palearctic primary vector survives the finite frozen grid; unrestricted missingness remains outside the claim.
- V6 species-list detection: Palearctic accessibility survives `99/100` baseline-supported surfaces and all `80/80` scenarios where list completeness declines with isolation; tropical accessibility is substantially more sensitive.
- observed nonlinear geometry: `0/12` promoted shapes under calibrated cross-scope gates.
- N1: channel-specific isolation interaction unsupported (`p=0.65516`); N2 remains unopened.
- source-side GloBI breadth: `0/4` promoted cells.
- H5c: independent biotic-vs-wind specificity not supported in the sole qualified cell (`p=0.41221`).
- H5d: `0/8` distributed-threshold design cells pass identifiability; observed lineage-threshold distributions remain closed.

## Claim ceiling

The paper may state:

1. one universal floral/reproductive island syndrome is not recovered;
2. response composition and direction differ among biogeographic contexts;
3. the positive Palearctic branch is substantially more observation-robust than the tropical accessibility component;
4. the Palearctic response is strongly structured by true genus composition, beyond arbitrary matched fine grouping;
5. the broad Palearctic response is not a robust beyond-genus response repeated uniformly across the flora;
6. several simple global area/geometry/pollination explanations fail explicit promotion gates.

The paper must state:

- the exact additional family→genus attenuation is spatially imprecise.

The paper must not state:

- that a precise family→genus taxonomic breakpoint has been estimated;
- that genus attenuation proves dispersal-only assembly;
- that within-lineage evolution is absent;
- that pollinator loss caused the Palearctic response;
- that H5c proves pollinators are irrelevant;
- that tropical accessibility is defended as strongly as the Palearctic branch;
- that endemicity is a time axis;
- that a smooth global cline rules out local thresholds.

## Chapter handoff

Chapter 1 answers:

> **Where, in which components and contexts, and at what lineage-assembly level is the floral island response represented?**

Chapter 2 / `izu-core` answers:

> **How and why does a resolved interaction → service → reproduction → phenotype chain change within one island system, and is its response geometry clinal or threshold-like?**

The global failure to identify a distributed-threshold generator is not a weakness to rescue; it defines why the local chapter is necessary.
