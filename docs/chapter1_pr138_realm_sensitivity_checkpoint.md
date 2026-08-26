# PR138 formal biogeographic-realm sensitivity checkpoint

Status: **provisional / not frozen**.

This checkpoint tests whether the broad northern/tropical/southern PR138 result is an artefact of the latitude-based context labels.

## External realm layer

- RESOLVE Ecoregions 2017 / terrestrial biogeographic realms
- source: `https://storage.googleapis.com/teow2016/Ecoregions2017.zip`
- license: CC-BY 4.0
- island assignment: point-in-polygon only
- no nearest-realm imputation

Run: `32927060210`
Artifact: `pr138-realm-sensitivity-32927060210`
Digest: `sha256:6556a7938a52547ddf7b630a822d743b2ba5ca479589716919441b96dcfa0521`

## Realm assignment

Of 8,265 islands, 5,583 receive a formal terrestrial-realm assignment without imputation.

- Nearctic: 1,456
- Neotropical: 1,045
- Palearctic: 947
- Indomalayan: 896
- Australasian: 843
- Afrotropical: 205
- Oceanian: 158
- Antarctic: 30

Realm assignment is much broader than current syndrome evidence support. Confirmatory syndrome support is presently concentrated in Palearctic and Neotropical islands; Oceanian is pilot-only. The strict source-backed Nearctic rerun (`32938637246`; artifact `pr138-nearctic-status-pilot-32938637246`) resolves 277/355 San Nicolas status rows but admits no binary Cedros statuses from its qualified `presumably introduced` marks. Nearctic all-native support is consequently 29 islands for `butterfly_like`, 28 for `generalized_accessible` and `large_bee_like`, and 27 for the other syndromes. It remains below the 30-island pilot gate despite many realm-assigned islands.

## Palearctic

All-native confirmatory slopes:

- `large_bee_like`: -0.1166, p=0.00145
- `generalized_accessible`: +0.1086, p=1.46e-7
- `selfing_core`: +0.0784, p=1.47e-6
- `selfing_syndrome`: +0.0828, p=8.79e-8
- `bird_like`: -0.1376, p=0.00181
- `butterfly_like`: -0.0825, p=0.00544

Predeclared northern classic projection:

- all native: +0.1027, FDR q=3.46e-6
- native non-endemic: +0.1017, FDR q=9.35e-7

Thus the northern classic syndrome is strongly reproduced inside the Palearctic and persists after endemic species are removed.

## Neotropical

All-native confirmatory data show a weaker classic-direction projection (+0.0325, q=0.0194), driven mainly by positive generalized accessibility and reproductive-assurance axes rather than loss of large-bee-like or gain of bird/butterfly-like scores.

Crucially, this does **not** persist in native non-endemics:

- native-nonendemic classic projection: +0.0175, q=0.0724, not supported
- native-nonendemic syndrome-vector omnibus: q=0.0573, not supported

The Palearctic-vs-Neotropical syndrome vector differs directly in native non-endemics (q=0.00401). In all natives the direct vector difference is not supported (q=0.192), indicating that endemic/floristic composition contributes materially to the Neotropical all-native pattern.

## Oceanian

Oceanian support is pilot-only (~34-40 islands by syndrome). Its northern-classic projection is negative and unsupported. Several individual slopes are heterogeneous; these are retained as pilot signals only.

## What this does and does not establish

Supported:

1. The northern classic syndrome is not only a broad-latitude artefact: it is strongly reproduced within the Palearctic.
2. The signal persists in native non-endemics in the Palearctic.
3. A weak Neotropical all-native classic-direction signal disappears after endemic species are removed, exposing a direct Palearctic-vs-Neotropical difference.
4. Realm-level coverage is currently insufficient to claim replication across both major northern realms: Nearctic is not yet pilot-testable.

Not supported / not identified:

- universal northern-realm replication;
- causal Bombus loss;
- causal butterfly or bird replacement;
- a universal tropical counter-syndrome across all tropical realms.

This sensitivity therefore strengthens the **biogeographic contingency** interpretation while also showing that floristic status matters and that formal-realm evidence coverage remains uneven.
