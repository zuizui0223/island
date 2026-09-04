# Chapter 1 area/capacity moderation checkpoint

Updated: 2026-08-27

## Result first

The frozen PR138-139 data support one stable area-conditioned source-pool pattern:

> In the broad Palearctic flora, isolation-associated enrichment in which
> source-available genera are represented is stronger on smaller islands.

This result holds in **8/8** prespecified source-mode × floristic-stratum
scenarios. It is not a genetic founder-effect estimate and it does not identify
a pollinator mechanism.

The plant-response evidence is narrower. Reproductive assurance is more strongly
associated with distance on smaller Neotropical islands in both primary floristic
strata. The direct Palearctic-versus-Neotropical reproductive-assurance moderation
axis is not FDR-supported, so this cannot be labelled a uniquely Neotropical Baker
effect. The broad northern-midlatitude and tropical plant-response vectors do not
show supported within-context area amplification.

## Frozen analysis contract

The model is predeclared as:

```text
response ~ distance + area + distance:area + climate_pc1..4
```

Inference uses spatial-block cluster-robust covariance. Area and distance are
continuous; no threshold was searched. A joint two-response vector must pass before
an individual response is classified. Regional moderation is tested directly rather
than inferred by comparing significance labels.

The two plant responses are:

- `accessibility_generalization` for the attraction/access branch;
- `reproductive_assurance` for the Baker-compatible reproductive branch.

The two source-lineage responses are:

- `entry_enrichment`: which source-available genera are represented, with one vote
  per genus;
- `loading_increment`: extra species weighting after a genus is represented.

Source matching fixes source-region prevalence and source-species-richness classes.
Thus `entry_enrichment` is a functional-composition response among available genera,
not raw genus richness or historical source ancestry.

## 1. Plant source-pool filter × island area

### Broad Palearctic genus entry: supported

Across all four frozen source definitions and both primary floristic strata:

- the distance coefficient at mean area is positive and FDR-supported;
- the distance × area coefficient is negative and FDR-supported;
- the joint entry/loading vector is supported;
- the classified state is `distance_effect_stronger_on_smaller_islands` in **8/8**
  scenarios.

The mean-area distance coefficient ranges from `+0.0672` to `+0.0965`; the
distance × area coefficient ranges from `-0.0444` to `-0.0333`. Therefore the
positive source-matched genus-entry enrichment associated with isolation is
consistently steeper on smaller Palearctic islands.

This advances the source-pool branch from “distance is associated with lineage
representation” to:

```text
distance/source accessibility × island area
  -> non-random representation of source-available genus functional positions
```

It does not show that fewer genera arrived, that a particular founding population
was small, or that drift generated the trait distribution.

### Regional contingency: supported, with different strength by contrast

For broad genus entry, the north-versus-tropical distance × area difference is
axis-supported in **8/8** source-mode × stratum scenarios. Northern interactions are
negative while tropical interactions are positive: the two regional entry branches
converge toward weaker slopes on larger islands.

The Palearctic-versus-Neotropical entry-moderation difference is supported in **5/8**
scenarios. Palearctic within-context small-island amplification is 8/8; Neotropical
within-context vectors are unsupported in 8/8. The formal realm difference is thus
substantial but not fully source-definition invariant.

### Species loading: not stable

Broad Palearctic `loading_increment` is classified as stronger on larger islands in
only **3/8** scenarios and lacks that axiswise classification in 5/8. The exact-SI
sensitivity recovers larger-island loading amplification in only 2/8 northern and
2/8 Palearctic scenarios, with the remaining source definitions unsupported or
vector-only.

Therefore Chapter 1 should not claim a general area law for species loading within
represented genera.

## 2. Baker-compatible founder filter

The reproductive-assurance branch is the only present assemblage-level proxy for a
Baker-type establishment filter. In the Neotropical realm it shows smaller-island
amplification in both primary strata:

| Stratum | Distance at mean area | Distance × area | Slope at area z = -1 | Slope at area z = +1 |
|---|---:|---:|---:|---:|
| all native | +0.1264 | -0.0672 | +0.1936 | +0.0592 |
| native non-endemic | +0.1088 | -0.0584 | +0.1673 | +0.0504 |

Both interaction axes pass FDR after the joint-vector gate. However:

- the Neotropical estimate is based on 59 response-supported islands in 14 spatial
  blocks, so small-cluster uncertainty remains a limitation;
- the broad northern-midlatitude and tropical within-context area vectors are not
  supported;
- the direct Palearctic-versus-Neotropical reproductive-assurance interaction axis
  is not supported;
- founding events, founder number, autonomous seed set at colonisation, and genetic
  drift are not observed.

The correct statement is therefore **Baker-compatible assemblage evidence in the
Neotropical sample**, not proof of Baker's law or a genetic founder effect.

## 3. Pollination-channel filter

The attraction/accessibility branch does not yield a simple “small islands lose
pollination service” result.

- neither northern-midlatitude nor tropical within-context plant vector supports
  axiswise area amplification;
- nevertheless, the direct tropical-minus-northern accessibility interaction is
  supported in both strata (`+0.1012`, `q = 0.00780`; `+0.0879`, `q = 0.0134`);
- Palearctic and Neotropical plant vectors are jointly area-sensitive, but no
  accessibility interaction axis or direct realm accessibility difference survives
  the hierarchical FDR gates.

Thus island area conditions the **regional contrast** in attraction/access geography,
but present data do not identify pollinator arrival, establishment, population
persistence, visitation, per-visit effectiveness, or local extinction. Pollinator
mobility and guild-labelled syndromes remain Discussion-level compatibility checks.

## 4. What island area means here

Area is a shared environmental-capacity proxy that can combine target size, habitat
and resource diversity, population persistence, and local-extinction risk. It is not
a direct measure of any one of those mechanisms. Connectivity is represented only by
the frozen mainland-distance/source-assignment variables; within-archipelago stepping
stone connectivity is not separately measured.

## Revised Chapter 1 evidence graph

```text
mainland distance/source accessibility × island area
|
|-- plant source-pool filter
|     |-- source-available genus functional-position entry: supported
|     |     `-- strongest stable result: small-island Palearctic amplification
|     |-- additional within-genus species loading: area law not stable
|     |-- Baker-compatible reproductive assurance: bounded Neotropical support
|     `-- genetic founder effect / drift: not measured
|
`-- pollination-channel filter
      |-- area-dependent regional accessibility contrast: supported
      `-- pollinator arrival, persistence and effective service: not measured
```

## Claim ceiling

Supported:

> Distance/source accessibility and island area jointly structure source-matched
> lineage representation, most consistently through genus entry in the broad
> Palearctic flora. Area also conditions some regional plant-response contrasts.

Not established:

- genetic founder effects or drift;
- a historical colonisation probability;
- causal Baker filtering at founding;
- pollinator mobility, loss, replacement or effective service;
- a single universal area-dependent island floral syndrome.

## Validation state

The dedicated frozen-input workflow is green:

- workflow: `Run Chapter 1 area capacity moderation`;
- submission-candidate run: `33067239884`;
- validated head: `6c2d7e267f872dc24983969cb6475ca2a8051975`;
- artifact: `chapter1-area-capacity-moderation-33067239884`;
- artifact id: `9644240178`;
- digest: `sha256:8cbc18dbe8369e1bf04d7e7b8d1cbe88d83ce478790e1d15dedcc2fb7a42d486`.

This later run supersedes the validation-only `33066985270` reference without changing the
frozen numerical estimands. It is the artifact bound by
`config/chapter1_submission_freeze_lock.json`.

The frozen-input local run produced 248 coefficient rows, 144 within-context vector
rows, 104 between-context coefficient rows and 72 between-context vector rows. A
clean second run reproduced all five outputs byte-for-byte. Dedicated tests pass
(`5 passed`), Ruff passes, key uniqueness and hierarchical classification gates pass,
and all q-values lie in `[0, 1]`.

Independent comparison of the downloaded Actions artifact with the local run matches
all non-numeric fields and row keys. Maximum cross-platform numerical differences are
`1.87e-12` or smaller, and the JSON manifest is semantically identical. Tracked config
hashes are newline-normalized so Windows and Linux verify the same frozen contract.

The dedicated Actions workflow reuses immutable PR138-139 artifact
`pr138-lineage-representation-bridge-33064250362` and verifies the frozen input hashes,
row counts, headline coefficients, state classifications and regional contrast counts.
