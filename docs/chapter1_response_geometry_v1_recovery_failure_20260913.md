# Chapter 1 response-geometry V1 recovery audit — frozen failure

## Decision

The first prospective response-geometry design audit ran successfully on the frozen final PR142 Chapter 1 artifact and **did not qualify any cell for observed nonlinear fitting**.

This is a method-identifiability failure, not evidence that the biological response is linear and not evidence that a breakpoint or reversal is absent.

The observed nonlinear Chapter 1 response remains unopened.

## Frozen run

- contract: `chapter1_response_geometry_identifiability_v1`
- branch: `ch1/el-hierarchical-depth-audit`
- workflow run: `34764480158`
- workflow commit: `52498fb51182803b4d53536d58d4d150e1368aa7`
- input PR142 run: `34232450884`
- input PR142 artifact: `chapter1-progressive-analysis-34232450884`
- input PR142 artifact ID: `10058653212`
- input PR142 digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`
- recovery artifact: `chapter1-response-geometry-recovery-34764480158`
- recovery artifact ID: `10319724218`
- recovery artifact digest: `sha256:d0772216dc47495b923617cb699ed1441d419906abe9d13b0a7d1ee99f884928`

The workflow completed successfully. All tests, lint checks, input artifact digest checks, recovery simulation, fail-closed output checks and artifact upload passed.

## Design support

All 24 requested evidence-scope × response × context × stratum design cells passed the predeclared >=50-island support floor.

The 12 cross-evidence cells were:

- 3 measurement-domain representatives (`plain_colour`, `generalized_form`, `self_compatibility`);
- northern-midlatitude and tropical contexts;
- all-native and native-nonendemic strata.

Each was evaluated separately in all-analysis-eligible and direct-only evidence.

## Why V1 failed

At the frozen target signal (`0.50` logit SD) under the declared spatial-block random-intercept stress (`SD = 0.25`), naive AICc shape selection too often converted a true smooth cline into a nonlinear response.

Across the 24 scope-specific cells:

- false selection of a nonlinear shape under a true cline ranged from **0.2083 to 0.5729**;
- false step-or-hinge selection under a true cline ranged from **0.2083 to 0.5729**;
- the predeclared ceiling was `0.10`;
- therefore **0 / 24 scope-specific cells** qualified;
- consequently **0 / 12 cross-evidence cells** qualified for a headline geometry analysis.

The failure is not driven by inability to recover a true step. At the same target signal, step recovery was approximately **0.958 to 1.000** across scope-specific cells. The principal V1 problem is false nonlinear selection when the truth is a smooth cline.

This is exactly the failure mode that a pre-analysis recovery audit was intended to expose.

## Representative diagnostics

| evidence | response | context | stratum | cline recovery | step recovery | hinge recovery | reversal recovery | false nonlinear under cline |
|---|---|---|---|---:|---:|---:|---:|---:|
| all-analysis | plain colour | tropical | all native | 0.427 | 0.972 | 0.920 | 0.896 | 0.573 |
| all-analysis | generalized form | tropical | native non-endemic | 0.781 | 0.979 | 0.781 | 0.802 | 0.219 |
| all-analysis | self compatibility | tropical | native non-endemic | 0.792 | 0.979 | 0.819 | 0.792 | 0.208 |
| direct-only | plain colour | northern midlatitude | all native | 0.500 | 1.000 | 0.962 | 0.844 | 0.500 |
| direct-only | generalized form | tropical | native non-endemic | 0.719 | 0.958 | 0.618 | 0.688 | 0.281 |

Breakpoint-location recovery for true steps was generally very accurate, whereas hinge breakpoint error was larger, especially in tropical cells. Thus breakpoint localization is not the main V1 problem for a true step; discrimination of cline versus nonlinear alternatives is.

## Frozen consequence

V1 must not be used to inspect or report an observed breakpoint.

In particular, the following remain prohibited:

- reporting the observed best step / hinge / reversal;
- selecting a breakpoint range from the observed response;
- lowering the false-selection ceiling after seeing an observed result;
- treating a selected step as evidence of pollinator loss;
- choosing direct-only evidence when all-analysis evidence disagrees.

## V2 repair principle

V2 may change the **calibration method**, but not because an observed biological breakpoint was seen: no observed nonlinear fit has been opened.

The repair must address the diagnosed failure directly:

1. retain the same five candidate geometries and frozen design cells;
2. use simulations under the monotonic null family with the actual spatial-block structure to calibrate a per-cell nonlinear-evidence threshold;
3. estimate that threshold on calibration simulations and assess error/recovery on independent validation simulations;
4. control false nonlinear selection at <=10% before any observed geometry is opened;
5. qualify step, hinge and reversal separately, because V1 shows that their recoverability differs materially;
6. require all-analysis and direct-only evidence to qualify before a nonlinear shape can be used as a headline result.

The V1 failure remains part of the audit trail and is not overwritten by V2.
