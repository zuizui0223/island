# Trait-specific source-genus prediction: decision memo (2026-11-17)

**Status:** Reanalysis of already verified original protected-source aggregates, not a new original-source run. Primary source: [source-only real joint congener prediction](../results/20261116_source_only_congener_real_joint_prediction/README.md). Keep PR #16 DRAFT.

## A stronger, falsifiable biological contrast than generic island floral simplification

The original same-accepted-species morphology pair test finds that knowing the **observed island's botanical genus composition** can improve source-only prediction for some floral combinations and worsen it for others. This is an outcome-specific **predictability heterogeneity** finding; it does not establish a causal genus assembly effect or an effective pollinator guild.

All-recorded tropical source-only same-target comparison, geographically equal 10-degree-block Brier gain (positive = genus-weighted mainland donors predict better than source matched-genus pooled donors):

| True same-species flower pair | Qualified islands | Genus weighted vs source pooled Brier gain, geo5 [block 95%] | Extra actual paired source state vs within-genus product of separate trait marginals |
|---|---:|---|---|
| Actinomorphic × deepest tube | 331 | **+0.11872 [+0.10286,+0.13290]** | +0.000054 [−0.000016,+0.000156] |
| Zygomorphic × deepest tube | 331 | **+0.02712 [+0.01804,+0.03639]** | +0.000005 [−0.000193,+0.000135] |
| Tubular/bell × deepest tube | 260 | **−0.03975 [−0.05418,−0.02659]** | +0.00120 [−0.00194,+0.00300] |

These are **different endpoint-specific eligible cohorts**, not a paired significance test of difference between traits. The morphology-specific sign reversal is the descriptive observation. Never label this a statistically established pairwise trait interaction without common-island paired resampling and a predeclared contrast.

Independent source-donor stringency (at least two accepted donor species per genus in both mainland modes) preserves the descriptive contrast: radial/deep +0.12987 (225 islands), tube/bell/deep −0.01553 (152 islands). In northern midlatitude, radial/deep +0.04583 (1227 islands), but no same-cohort tropical-vs-midlatitude contrast has been computed.

### Identification ceiling

The source-only model *uses actual observed island genus membership*. It is a conditional morphology prediction, NOT a model that predicts which genus arrives on an island. The target accepted species names are removed from the mainland donor pool, but whole archipelago and literature provider were not necessarily held out. Source GIFT candidate lists are not measured propagule arrivals. Trait recording is strongly differential: only **6/331** tropical radial/deep cases have island-vs-source documentation gap <=10 percentage points; only **1** case also has variable source genus support for both traits. Thus even the large +0.11872 predictive gain is not confirmed in a balanced independent cohort. Pollinator guild mechanism is not identified.

### Required next real-data evaluation

1. Build a **common-target/common-island** radial/deep vs tubular/deep evaluation where all true same-species Direct traits exist on each target species, both source modes and both genus donor gates pass. Preserve source- and genus-coverage counts; report N original botanical islands, independent 10-degree blocks, first source entities, archipelagos and providers. Do not restrict the global study population to these complete-profile islands: this is only a specific mechanistic falsification subset.
2. On exactly these matched islands calculate per-island paired `gain_radial_deep - gain_tubular_deep`, block-equal mean and group-resampled interval. If fewer than 30 islands/10 blocks or too few independent providers, status NOT_IDENTIFIABLE.
3. Require actual held-out provenance records: target accepted names excluded; target archipelago and phenotype provider excluded from training; original source region provenance not leaked. Merely setting booleans to True is not proof. Recompute donor models inside each fold; do not relabel already trained in-sample source predictions as cross-validated.
4. Recompute for strict native, independently documented flower phenotypes, and no overlapping literature sources. If support fails, retain descriptive conditional source-genus predictability and do not claim cross-region adaptation or guild shifts.
5. Only after morphology robustness, test whether independently documented bird, butterfly, Sphingidae, Bombus and other bee effort-qualified availability improves heldout predictions above climate, geographical access, source species availability and genus composition. Missing visitor records remain missing, never guild absence.

## Main-paper consequence

Retain global 4,505-island floristic assembly, functional beta and regional real same-species morphological pairs as main results. Present **trait-specific source-genus predictability** as an exploratory secondary result with the explicit negative tube/bell control, not a one-direction claim that genus always explains flower architecture. No 17D KDE volume, no unverified guild names, and no in-situ island evolutionary inference.
