# Chapter 1 H5 — GloPL floral-architecture moderation result

## Question

Does generally accessible floral architecture buffer the increase in experimental pollen limitation with distance from major continental landmasses?

This is a species-level functional association test. It does not treat floral phenotype as observed pollinator identity and it does not test causal mediation.

## Frozen design

The design was frozen at commit `4f478c76b4ac0ed963a0b7d0ae92a8b143e72283`, before the GloPL effect column was joined to the floral-architecture overlap table.

Three predeclared atomic architecture contrasts were tested separately:

1. `generalized_form`: open-radial / brush-puff / composite-head versus specialized forms;
2. `actinomorphic_symmetry`: actinomorphic versus zygomorphic;
3. `shallow_open_tube`: absent-or-open / shallow versus deep.

The prediction was a negative distance × accessible-architecture interaction: the increase in pollen limitation with distance should be weaker in the generally accessible state.

The analysis retained the full-GloPL parent conventions: publication total weight = 1, publication-cluster robust uncertainty, measurement fixed effects, and the frozen parent distance standardization. Supplemental-only and no-zero-constant fits were frozen sensitivities. Route B required at least two of the three atomic contrasts to satisfy the buffering rule.

## TDD and provenance

- RED run: `35094224454` — expected module-absent collection failure.
- Canonical GREEN run: `35094588521`, job `104788690975`.
- Validation: 6 tests passed; Ruff passed.
- Artifact: `10445189257`.
- Artifact digest: `sha256:ef81551d41cc747b07e927b9a5c4f13eaf0fc013a5cf8ac3fb1a5491aadec8c4`.

## Outcome-blind overlap

The preflight used all 2,969 GloPL rows but did not materialize effect columns.

- total exact-matched architecture overlap: 624 species;
- `generalized_form`: 246 matched species, global gate passed;
- `actinomorphic_symmetry`: 582 matched species, global gate passed;
- `shallow_open_tube`: 52 matched species, global gate failed.

### Generalized form

Restricted architecture:
- 169 species;
- 245 sites;
- 174 publications;
- 56 offshore sites.

Generally accessible architecture:
- 77 species;
- 114 sites;
- 79 publications;
- 24 offshore sites.

The North–Tropical secondary diagnostic was not admitted because tropical accessible support was only 5 species, 8 sites, 5 publications and 1 offshore site.

### Actinomorphic symmetry

Zygomorphic:
- 213 species;
- 263 sites;
- 180 publications;
- 89 offshore sites.

Actinomorphic:
- 369 species;
- 434 sites;
- 315 publications;
- 117 offshore sites.

Both northern-midlatitude and tropical context gates passed.

### Shallow/open tube

Deep:
- 34 species;
- 41 sites;
- 33 publications;
- 4 offshore sites.

Open/shallow:
- 18 species;
- 22 sites;
- 15 publications;
- 7 offshore sites.

This contrast failed the frozen global support gate and is not a biological negative result.

## Results

### Generalized form

Primary:
- restricted distance slope = `+0.07409`;
- accessible distance slope = `+0.00735`;
- distance × accessible interaction = `-0.06674`;
- one-sided negative p = `0.2441`.

The point estimate is in the predicted buffering direction, but it is imprecise. It also reverses direction in the supplemental-only sensitivity (`+0.06455`), while the no-zero-constant sensitivity remains negative (`-0.07565`).

**Decision: unsupported; sensitivity direction not retained.**

### Actinomorphic symmetry

Primary:
- zygomorphic distance slope = `+0.08518`;
- actinomorphic distance slope = `+0.12013`;
- interaction = `+0.03496`;
- one-sided negative p = `0.6690`.

This is opposite the buffering prediction in the primary fit. Supplemental-only is slightly negative (`-0.01374`), but no-zero-constant is again positive (`+0.03989`).

North–Tropical diagnostic:
- northern interaction = `+0.04181`;
- tropical interaction = `+0.12071`;
- tropical-minus-northern = `+0.07891`, p = `0.7130`.

**Decision: unsupported; sensitivity direction not retained.**

### Shallow/open tube

Not evaluated because the frozen outcome-blind support gate failed.

## Route-B decision

Supported atomic contrasts: **0/3**.

Classification: **`floral_architecture_buffering_not_supported`**.

The result does not say that floral architecture is irrelevant. It says that, in the current exact species-matched GloPL sample, the broad distance-associated increase in experimental pollen limitation is not robustly weaker in species carrying the three predeclared generally accessible architecture states.

## H5 implication

Together with the reproductive-assurance moderation test, the result narrows the H5 story:

```text
isolation -> experimental pollen limitation
    supported globally
        |
        +-- reproductive-assurance buffering: not supported
        |
        +-- atomic accessible-architecture buffering: not supported
        |
        +-- North-Tropical service branching: not supported
```

Thus pollen limitation is supported as a broad correlate/stressor of geographic isolation, but current species-level trait moderation does not provide the missing functional bridge from that stressor to the context-dependent six-atomic Chapter 1 H2 plant response.

## Claim ceiling

This result does not establish absence of local pollinator-mediated selection, absence of floral adaptation, trait evolution caused by isolation, causal selection by pollen limitation, H2 mediation, or named pollinator identity.
