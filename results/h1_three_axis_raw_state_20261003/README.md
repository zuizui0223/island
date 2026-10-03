# H1 three-axis raw-state reanalysis — 2026-10-03

## Purpose

This analysis returns H1 to the original Chapter 1 measurement design: **106,295 species × three biological axes** rather than treating seven binary contrasts as seven complete traits.

The frozen database contains 318,885 possible species×axis cells and **222,688 resolved cells (69.83%)**:

- flower colour: **82,556**;
- floral structural complexity: **91,635**;
- reproductive assurance: **48,497**.

Every resolved cell contributes through whatever ontology-valid component states are observed. Missing component traits are not set to zero, complete cases are not required, and raw multistate memberships are retained.

Formal inference fits trait-specific state prevalences with their own denominators and tests the full raw-state isolation-response vector within each biological axis and geographic region.

## Frozen inputs and validation

- workflow run: **37117607718**
- artifact: **11272053673**
- artifact digest: `sha256:b40721ca0f20b3c1af3e92782c3466d50bb687e9dea17dbe8b9613e9aad574ed`
- trait database: **222,688 / 318,885 resolved species×axis cells**
- minimum formal raw-state support: **30 unique species**
- formal raw states retained across evidence scopes: **97**

### Ontology audit

Direct use of the species×axis cells exposed an old low-confidence bookkeeping defect: 900 all-analysis cells from the historical validated-low layer contained within-axis trait-label permutations. All were uniquely recoverable using only the set of trait slots already declared in the same cell plus the frozen trait ontology; no species identity, geography or outcome was used.

After this audit:

- all-analysis structural cells: **91,635 / 91,635 ontology-valid**, including 860 uniquely repaired cells;
- all-analysis colour cells: **82,556 / 82,556 ontology-valid**;
- all-analysis reproductive cells: **48,497 / 48,497 ontology-valid**, including 40 uniquely repaired cells;
- Direct-only cells required **no repairs**.

Thus all **222,688 resolved species×axis cells** remain ontology-valid. After the 30-species raw-state support gate, **222,687 / 222,688 all-analysis cells** still contribute at least one state to the formal H1 tests (99.9996%). Direct-only retains 102,924 / 102,932 resolved cells (99.9922%).

## Primary all-observed result

### All-analysis

| Axis | Northern mid-latitude | Northern high latitude | Tropical | Southern extratropical |
| --- | ---: | ---: | ---: | ---: |
| Reproductive assurance | q=1.09e-20 | q=7.81e-18 | q=2.16e-7 | q=2.15e-6 |
| Floral structural complexity | q=3.21e-95 | q=1.26e-82 | q=2.02e-34 | q<1e-300* |
| Flower colour | q=4.57e-4 | **q=0.155** | q=5.45e-8 | q=9.57e-7 |

### Direct-only

| Axis | Northern mid-latitude | Northern high latitude | Tropical | Southern extratropical |
| --- | ---: | ---: | ---: | ---: |
| Reproductive assurance | q=1.81e-16 | q=4.65e-53* | q=5.43e-7 | q=4.37e-15 |
| Floral structural complexity | q=8.77e-55 | q=6.46e-176 | q=1.17e-27 | q<1e-300 |
| Flower colour | q=2.94e-9 | **q=0.189** | q=2.06e-6 | q=1.63e-12 |

*One raw-state optimizer flag persists in each marked primary axis test. Removing the failed state entirely yields a fully converged axis test with the same decision: southern all-analysis structure remains supported (drop-failed q < 1e-300) and northern-high Direct reproductive assurance remains supported (drop-failed q = 2.41e-28).

**Primary conclusion:** reproductive assurance and floral structure are recurrent isolation-associated response domains in all four regions and both evidence scopes. Colour composition is geographically contingent, with no supported broad northern-high-latitude axis response.

## Raw-state interpretation

The axis tests are not direction-free black boxes.

Across the broad flora, reproductive assurance repeatedly shifts toward assurance/selfing states: examples include increasing self-compatibility or autonomous/selfing states and decreasing self-incompatibility or predominantly outcrossing states. The exact component carrying the response varies among regions.

Structural composition also reorganizes in every region, but the phenotype is not universal. Northern-high and tropical floras show increases in open/radial and actinomorphic states with decreases in several restricted forms, whereas southern extratropical islands combine increasing open/radial structure with increasing deep tube and decreasing shallow tube. The recurrent result is therefore **structural reorganization**, not universal simplification.

Colour is the least recurrent domain. It changes in three of four broad regions but is unsupported in northern high latitudes under both evidence scopes.

## Floristic-origin / introduced-species audit

### Strict source-backed native records

Only northern mid-latitudes and the tropics have sufficient source-backed native support.

In both regions, **reproductive assurance and structural complexity remain strongly associated with isolation**. In the tropics, all-analysis q values are:

- reproduction: **5.46e-4**;
- structure: **1.29e-213**;
- colour: **0.0508**.

Direct-only tropical native results are strongly supported for all three axes.

This changes the interpretation of the earlier seven-indicator strict-native tropical q=0.206: that result showed that the selected directional seven-indicator vector was not supported; it did **not** show absence of native reproductive or structural reorganization when the full raw measurement axes are retained.

### Regional-native-compatible sensitivity

Reproduction and structure remain supported in **4/4 regions** in both evidence scopes. Colour is 3/4 in all-analysis and 4/4 in Direct-only.

On identical Level-3-area-complete support, reproduction and structure remain 4/4 supported before adjustment and after adding log TDWG-Level-3 area. The geographic resolution of WCVP therefore does not explain these two recurrent domains.

### Native-compatible versus incompatible/introduced response vectors

The complementary regionally incompatible/introduced flora also changes with isolation, but generally **not through the same raw-state vector**.

Formal status-by-isolation interaction tests show:
- structure differs between the two provenance partitions in **4/4 regions** in both evidence scopes;
- reproductive assurance differs in **4/4 regions** in both evidence scopes;
- colour differs in 3/4 all-analysis regions and 4/4 Direct-only regions.

Thus the broad observed-flora signal is not a single response reproduced identically by native-compatible and non-native-compatible plants; it is a mixture of provenance-specific assembly trajectories.

### Strict known-origin comparison

Source-backed native versus source-backed introduced records are jointly testable only in the tropics.

All-analysis:
- colour vector difference: q=**0.441**, unsupported;
- structure: q=**1.45e-57**, strongly different;
- reproductive assurance: q=**2.36e-11**, strongly different.

Direct-only:
- colour: q=**0.0231**;
- structure: q=**6.98e-28**;
- reproductive assurance: q=**1.53e-10**.

The reproductive-assurance response is especially distinct: cosine similarity between native and introduced vectors is **-0.758** in all-analysis and **-0.853** in Direct-only. Known introduced tropical flora therefore do not mimic the native reproductive-assurance isolation response; their response is approximately opposed in multivariate direction.

## Claim ceiling

The three-axis result supports contemporary assemblage reorganization, not a historical evolutionary mechanism. Source-backed native data are insufficient for northern-high and southern-extratropical four-region replication, and regional-native compatibility is not exact focal-island nativeness.

The defensible H1 is:

> **Geographic isolation is associated with recurrent reorganization of reproductive assurance and floral structure across all four regions, while colour composition is more geographically contingent. These responses persist in regional-native-compatible floras and are not reproduced as the same response vector by introduced/incompatible floras.**

This does not distinguish colonization filtering, persistence/extinction, species sorting or within-lineage evolution.
