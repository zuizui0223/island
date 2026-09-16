# Chapter 1 H3A broad observed-flora taxonomic depth — result

Status: **completed; the broad H2 contrast is not erased by source-free family/genus residualization.**

## Provenance

- branch: `ch1-all-data-primary`
- architecture contract: `chapter1_hypothesis_architecture_v3`
- H3A contract: `chapter1_h3_observed_taxonomic_depth_v1`
- workflow: `Run Chapter 1 H3 observed taxonomic depth`
- successful run: **35043562731**
- artifact: `chapter1-h3-observed-taxonomic-depth-35043562731`
- artifact ID: **10426327347**
- digest: `sha256:05c23e93a7954b76399efae1df494f0fc9a78e58a09e1166e1a971b0a391562d`
- frozen Chapter 1 input run: `34232450884`

The H3A design was created to fill the missing cell: all-observed H2 response x taxonomic-depth decomposition. GIFT source assignments are deliberately absent because all-observed flora contains introduced and unresolved-status records.

## Design

For each of the six atomic outcomes, scored species were joined to fixed family/genus taxonomy. Family and genus expectations were estimated as leave-one-species-out means from the scored species pool. A species was retained only if both its family and genus contained at least one other scored species for that outcome. The same retained island-species observations were then used for all three stages:

1. `observed_score`;
2. `after_family_residual = state - LOO family mean`;
3. `after_genus_residual = state - LOO genus mean`.

No per-island species cutoff was imposed. The residual-stage model used equal island weight, climate PC1-4 and island area controls, and spatial-block cluster-robust covariance. A paired 300-draw spatial-block bootstrap quantified attenuation uncertainty.

Because equal-island linear decomposition changes the model form relative to H2, a separate gate re-fitted the **same beta-binomial H2 model** on the exact H3 common species support before any attenuation interpretation.

## H2 common-support gate

The original North-Tropical six-atomic H2 contrast remains strongly supported after restricting to species eligible for both family and genus LOO decomposition:

- all-analysis evidence: **3,560 islands / 187 blocks, p = 8.39e-05**;
- direct High/Medium: **3,487 islands / 187 blocks, p = 1.21e-04**.

Therefore loss of taxonomically isolated species does not remove the broad H2 response. The H3A decomposition is addressing the same broad response, not a different support-defined signal.

## Source-free taxonomic residualization

### All-analysis evidence

Equal-island linear decomposition:

- observed common-support vector: `p = 0.05452`, norm `0.05387`;
- after family residualization: **`p = 0.01043`**, norm `0.05281`;
- after genus residualization: **`p = 0.01150`**, norm `0.02993`.

Point attenuation:

- family: `1.96%`;
- genus: `44.44%`.

Paired spatial-block bootstrap genus attenuation:

- median `44.17%`;
- 95% interval `[-2.36%, 68.57%]`.

The attenuation magnitude is therefore imprecise, but the post-genus North-Tropical residual vector remains supported.

### Direct High/Medium evidence

Equal-island linear decomposition:

- observed common-support vector: `p = 0.23247`, norm `0.05257`;
- after family residualization: **`p = 0.02497`**, norm `0.05021`;
- after genus residualization: **`p = 0.004593`**, norm `0.04802`.

Point attenuation:

- family: `4.49%`;
- genus: `8.65%`.

Paired spatial-block bootstrap genus attenuation:

- median `14.02%`;
- 95% interval `[-94.12%, 58.94%]`.

Again, the exact attenuation fraction is not identified, but genus residualization clearly does not erase the broad context difference.

## Which components remain after genus residualization?

The equal-island post-genus residual is not a simple scaled copy of the original vector. The clearest retained components are:

- `generalized_form`: positive Tropical-minus-North residual difference (`p = 0.00363` all-analysis; `0.00327` direct);
- `self_compatibility`: positive residual difference (`p = 0.0212` all-analysis; `0.00768` direct).

Other components are weaker or change contribution after residualization. H3A therefore indicates reorganization below the raw genus-composition level rather than one genus-fixed copy of the original six-trait vector.

## Comparison with defended native H3B

This result is structurally different from frozen v11 native assembly evidence. In the defended Palearctic native analysis, the response passed the observed and family-adjusted gates but disappeared after source-matched genus adjustment (`4/4 -> 4/4 -> 0/4`), with approximately `78.8-85.9%` attenuation across broader frozen profiles; the matched-complexity pseudo-genus test also showed true genus structure stronger than arbitrary grouping (`p = 0.02899`).

H3A instead shows that the broad all-observed North-Tropical H2 contrast remains detectable after source-free genus residualization.

Therefore the two results must not be written as one causal ladder:

> **Broad contemporary-flora context dependence is not simply the same genus-assembly signal seen in the defended native Palearctic subset.**

The similarity is only that taxonomy matters in both layers; the depth and inferential meaning differ.

## Final structural decision

The Chapter 1 architecture remains a nested two-panel design:

- **Panel A — global contemporary flora:** H1 universal response -> H2 biogeographic branching -> H3A source-free taxonomic representation depth; broad H2 persists below genus.
- **Panel B — defended native assembly:** narrower native Palearctic response -> H3B source-matched family/genus decomposition; strong genus structuring.

H4 is restored to **area moderation** and remains measurement-sensitive. H5 remains the independent-mechanism cross-examination.

The broad and native panels can be compared, but H3B must not be used to claim that the global all-observed H2 response is explained by native genus assembly.
