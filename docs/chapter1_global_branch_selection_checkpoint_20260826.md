# Chapter 1 global branch observation-selection checkpoint

Date: 2026-08-26

## Result first

The broad conclusion that northern-midlatitude and tropical island floras follow different
plant-response vectors survives outcome-blind observation-selection weighting. The narrower
claim that the northern-midlatitude band independently shows a supported accessibility-
generalization vector does not survive.

The sensitivity uses the full 8,265-island geography universe to model whether each primary
axis score is observed. Each axis and floristic stratum gets a separate selection model using
only island area, climate, mainland distance, broad context and distance-by-context terms.
Branch-score values, ecological effect estimates and p-values are prohibited inputs.

## Headline results

### Analysis-regime layer

For all native species, the direct northern-midlatitude versus tropical two-axis vector
difference is retained under every mode:

| Selection mode | Direct vector-difference q | North state | Tropical state |
|---|---:|---|---|
| unweighted | `3.20e-6` | accessibility generalization | reproductive assurance |
| stabilized IPW, propensity clipped at 0.05 | `4.23e-6` | no supported context vector | reproductive assurance |
| stabilized IPW, propensity clipped at 0.025 | `1.77e-6` | no supported context vector | reproductive assurance |

Native non-endemics retain the direct regional vector difference as well:

- unweighted: `q = 3.63e-7`;
- IPW 0.05: `q = 1.60e-5`;
- IPW 0.025: `q = 1.64e-5`.

The northern all-native accessibility slope remains positive but loses axis and vector
support after weighting:

| Selection mode | Accessibility slope | Axis q | Context-vector q |
|---|---:|---:|---:|
| unweighted | `+0.04755` | `0.0282` | `0.0380` |
| IPW 0.05 | `+0.04043` | `0.174` | `0.291` |
| IPW 0.025 | `+0.03399` | `0.404` | `0.523` |

Therefore the broad-regime result should be stated as a **supported regional vector
difference**, not as an observation-robust independent northern accessibility branch.

### Biogeographic-realm layer

Palearctic remains strongly supported as accessibility generalization plus reproductive
assurance in both floristic strata and both IPW modes. For all native species:

- IPW 0.05: accessibility `+0.12365`, `q = 2.76e-9`; reproductive assurance
  `+0.05354`, `q = 2.52e-4`;
- IPW 0.025: accessibility `+0.11284`, `q = 1.76e-6`; reproductive assurance
  `+0.04881`, `q = 1.24e-4`.

Neotropical resolves to reproductive assurance only after weighting. The direct
Palearctic-versus-Neotropical two-axis comparison, marginal in the unweighted all-native
analysis (`q = 0.0574`), becomes supported under IPW (`q = 0.00109` and `0.00105`). Native
non-endemic direct differences remain supported in all modes.

This strengthens the realm-specific WHERE result, but IPW is a measured-selection
sensitivity and must not be described as source-pool standardization.

## Positivity and information loss

The analysis predeclares two propensity lower bounds and caps weights at 20. The observed
maximum weight is `11.95`; no row reaches the cap. The lowest effective-sample fraction is
`0.388`, in the realm/native-nonendemic 0.025-clipping analysis. Headline examples:

- analysis-regime all-native axes retain approximately 43-56% of the nominal observed
  information under IPW;
- realm all-native axes retain approximately 39-49%;
- 24-67 headline observations have propensities clipped, depending on axis, layer and mode.

This is adequate for a declared sensitivity but not strong overlap. The loss of the broad
northern context-vector support is therefore material and must be retained as a claim
boundary rather than repaired by post-hoc weighting choices.

The analysis-regime endemic selection models did not converge. Their IPW results fail
closed to `not_testable`; unweighted endemic results remain separate. Endemic IPW is not
used in any headline statement.

## Updated Chapter 1 claim

Supported after selection sensitivity:

> The distance-associated accessibility/reproductive-assurance response vector differs
> between northern-midlatitude and tropical island floras. At the realm level, Palearctic
> floras robustly show increasing accessibility generalization plus reproductive assurance,
> whereas Neotropical floras resolve mainly to reproductive assurance after measured
> observation-selection adjustment.

Not supported as an observation-robust general statement:

> Northern-midlatitude island floras as a whole independently exhibit a supported
> accessibility-generalization branch.

Still unidentified:

- source-region pollination-channel availability;
- source-pool-standardized assembly;
- pollinator mobility, establishment or loss;
- effective pollination service;
- attraction relaxation, reinforcement or functional replacement.

## Reproducibility and claim boundary

The locked sensitivity contract is `config/chapter1_global_branch_selection.yml`. The
executable is `src/island_v2/chapter1_global_branch_selection.py`. It outputs selection
model fits, stabilized weights, effective-sample diagnostics, weighted slopes, context
vectors, direct context contrasts and fail-closed branch states.

Persistence under IPW means the main regional vector difference is not readily explained
by the measured predictors of score availability. It does not establish missing at random
for unmeasured observation processes and does not raise the result to pollination causation.
