# Chapter 1 submission freeze — 2026-09-15

## Canonical narrative surface

The current Chapter 1 manuscript remains:

- `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`

The mandatory figure-reference / legend synchronization layer is:

- `docs/chapter1_v8_submission_figure_sync_20260915.md`

Until the next full-text copyedit commit, manuscript interpretation and figure numbering must follow these two files together. Earlier v7 and 9/14 figure sketches remain historical development surfaces.

## Canonical main figures

### Figure 1 — hierarchy-of-assembly inference map

- result lock: `config/chapter1_v8_figure1_result_lock.json`
- run: `34862804485`
- artifact: `10355018602`
- digest: `sha256:551e5f9738ca0193083232c27ba1abbdccc92e718c9553f175d68c1cde4c5e66`
- status: rendered and visually checked

### Figure 2 — primary biogeographic branching

- result lock: `config/chapter1_v8_figure2_result_lock.json`
- run: `34810404509`
- artifact: `10334412842`
- digest: `sha256:c177e28ac75a8bccf00e9f0e1a467f2958907ba3c3eae87334a62e5924780675`
- status: rendered and visually checked

### Figure 3 — source-matched assembly depth

- result lock: `config/chapter1_v8_figure3_submission_result_lock.json`
- run: `34816668470`
- artifact: `10336875587`
- digest: `sha256:853fd1c590d86341d6302c1da8e998af93b0b5fa3b35c464dd3762378385bd4d`
- status: corrected submission layout rendered and visually checked

### Figure 4 — falsification and claim boundaries

- result lock: `config/chapter1_v8_figure4_result_lock.json`
- run: `34817349228`
- artifact: `10337380746`
- digest: `sha256:f41a90089ddda1d285186ca441d416d71cf1c276e742a7da917dbb572f39fbb0`
- status: rendered and visually checked

## Frozen paper-level evidence sequence

1. **Context branching:** one universal floral/reproductive response is not recovered.
2. **Component decoupling:** reproductive assurance and floral accessibility can move differently across biogeographic contexts.
3. **Assembly depth:** the strongest Palearctic response is retained after family adjustment but not after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`); genus attenuation is `78.8–85.9%` of the observed vector and conditional family→genus attenuation is `70.6–79.1%`.
4. **Observation defense:** the Palearctic accessibility branch survives `99/100` baseline-supported V6 bias surfaces; the tropical component is more sensitivity-bounded (`35/75` survive), while the North–Tropical multivariate contrast survives `70/75`.
5. **No promoted common nonlinear geometry:** `0/12` broad atomic cells pass the frozen all-analysis + direct-only nonlinear promotion rule.
6. **No identified global pollinator-specific mechanism:** N1, sampled source breadth and H5c do not promote a causal channel interpretation; H5c gives `distance × biotic = +0.06495`, `p=0.41221` in the sole prospectively qualified cell.
7. **Distributed thresholds not identifiable globally:** H5d qualifies `0/8` cells, so observed genus-level threshold distributions remain closed.

## Publication-facing conceptual claim

> **A floral island syndrome can be a hierarchically assembled phenotypic syndrome: a visible community-level trait pattern whose direction and taxonomic depth depend on biogeographic context.**

Generalized inference:

> **Trait syndromes observed across environmental gradients need not be repeated organismal adaptations; identifying the assembly depth of a syndrome is a prerequisite for mechanistic interpretation.**

## Claim ceiling

The submission must not state that:

- pollinator loss caused the Palearctic response;
- genus attenuation proves dispersal alone or absence of within-lineage evolution;
- H5c `p=0.412` proves pollinators do not matter;
- tropical accessibility is as observation-robust as the Palearctic branch;
- a smooth global response proves local systems lack thresholds;
- distributed thresholds are supported by Chapter 1;
- endemicity is a direct time/evolution axis.

## Chapter 1 / Chapter 2 handoff

Chapter 1 resolves:

- response direction;
- component decoupling;
- biogeographic context;
- taxonomic assembly depth;
- observation and mechanism claim boundaries.

Chapter 2 (`izu-core`) resolves the within-system chain:

`interaction state -> effective service -> reproductive outcome -> phenotype`

and can prospectively test cline, step, shared-breakpoint and channel-specific response geometry. The global H5d non-identifiability result is the empirical reason this local mechanistic resolution is needed.
