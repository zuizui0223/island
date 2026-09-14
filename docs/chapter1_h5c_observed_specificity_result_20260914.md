# Chapter 1 H5c independent pollination-mode specificity — observed result

## Decision

The only prequalified observed H5c cell does **not** support the predicted biotic-versus-wind specificity of the Palearctic accessibility response.

This result is deliberately narrow. It is stronger than floral-syndrome concordance because the pollination-mode labels were independently documented, but it remains a negative-control comparison rather than a direct measurement of visitor loss or effective pollination service.

## Prospective sequence

Observed outcomes were kept closed while GIFT `pollen_vector_mode` coverage and the realized Chapter 1 island/mode design were audited. The qualification run recovered 5,771 unambiguous independent pollination-mode species (4,538 biotic; 1,233 abiotic wind). Of eight scope × stratum × context cells, only one passed both the frozen support gate and the power/type-I gate:

- direct-only evidence;
- native non-endemics;
- Palearctic;
- 132 islands with biotic support and 60 with wind support;
- 22 and 18 represented spatial blocks, respectively;
- type-I error `0.088`;
- recovery at the frozen 0.5-SD interaction `1.00`.

All other observed cells remained closed. Tropical cells failed wind-support gates; all-analysis Palearctic cells and direct-only all-native Palearctic exceeded the frozen 0.10 simulation type-I ceiling.

Qualification receipt:

- run `34803493307`;
- artifact `10332810251`;
- digest `sha256:074a52323ba835201043a3a5980cc2b2d0ff79878e83bc9be19511ae83a1c0e8`.

## Observed result

Canonical observed run:

- run `34803837463`;
- artifact `10332871066`;
- digest `sha256:fb8927bcf2b6e1c762b3ecbe59c42360faac679cf37b1ad0ed31386631a0f405`.

The frozen test was the cluster-robust `distance × biotic` interaction in the direct-only native-nonendemic Palearctic subset, with wind-pollinated plants as the reference mode and the standard area/climate covariates.

- interaction estimate: `+0.06495`;
- cluster-robust SE: `0.07921`;
- 95% CI: `[-0.09030, 0.22020]`;
- p = `0.41221`;
- classification: `no_pollination_mode_specificity_support`.

Mode-specific distance slopes were:

- abiotic wind: `-0.05393`, 95% CI `[-0.20672, 0.09886]`;
- biotic: `+0.01102`, 95% CI `[-0.01042, 0.03246]`.

The fitted design contained 133 unique islands and 23 spatial blocks. Species-support sums were 8,773 for biotic and 131 for wind rows, exactly matching the pre-outcome qualification support identity.

## Scientific consequence

This test does not provide the favorable specificity evidence that would be expected if the Palearctic accessibility gradient were clearly stronger among independently documented animal-pollinated plants than among wind-pollinated plants.

Accordingly, Chapter 1 should **not** promote a pollinator-specific mechanism from the global data. This conclusion is concordant with, but logically distinct from, the failed N1 channel-heterogeneity gate and the null GloBI source-breadth extension.

The result does not show that pollinators are irrelevant. The wind comparison is sparse relative to biotic support, only one design cell passed the prospective qualification gate, and `pollen_vector_mode` is a coarse dependency label rather than effective-service data. Its proper role is therefore an adverse but bounded piece of evidence in the mechanism case.

## Court-style evidentiary interpretation

Evidence supporting a robust plant-side case still includes:

1. direct biogeographic response-vector heterogeneity;
2. Palearctic persistence across evidence scopes and native non-endemics;
3. V5 trait-missingness and V6 species-detection defenses;
4. strong family-to-genus attenuation;
5. rejection of a common global breakpoint.

Evidence **against** a strong Chapter 1 pollinator-causation claim now includes:

1. N1 prospective channel heterogeneity not promoted;
2. GloBI source-side sampled channel breadth not promoted;
3. the independently documented biotic-versus-wind specificity test not supported;
4. H5d distributed-threshold generator not identifiable against heterogeneous smooth clines.

The defensible synthesis is therefore that the floral/reproductive island response is real, context-dependent and strongly assemblage-structured, while the pollinator-side cause remains unidentified at global scale.

## Claim ceiling

Allowed:

> In the only prospectively qualified independent pollination-mode comparison, the Palearctic accessibility response did not differ detectably between biotically and wind-pollinated plants.

Blocked:

- pollinator loss causes H3;
- wind and biotic plants are universally equivalent;
- pollinators are irrelevant to island floral assembly;
- reopening the failed N1/N2 chain;
- using the null H5c result to negate local mechanistic tests in `izu-core`.
