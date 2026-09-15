# P0 immutable claim-to-result ledger — 2026-09-15

Status: executable P0 inventory for the island-first NEE development plan.

This ledger does not replace any frozen result. It records what is already supported, the exact immutable source used for the current manuscript, and which defenses remain open for P1.

## Immutable source backbone

Primary biological analysis:

- workflow run: `34232450884`
- artifact ID: `10058653212`
- artifact name: `chapter1-progressive-analysis-34232450884`
- digest: `sha256:b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999`

Canonical submission figure locks:

- Figure 2: `config/chapter1_v8_figure2_result_lock.json`
- Figure 3: `config/chapter1_v8_figure3_submission_result_lock.json`
- Figure 4: `config/chapter1_v8_figure4_result_lock.json`

Canonical manuscript surface:

- `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`
- `docs/chapter1_submission_freeze_20260915.md`

## Claim ledger

| ID | Island-biological statement | Frozen estimand / result | Scope / denominator | Immutable source | Current status | Claim ceiling / P1 need |
|---|---|---|---|---|---|---|
| C1 | There is no single global floral island syndrome. | Primary H2 response is the joint accessibility/generalization + reproductive-assurance response; context branching is the manuscript-level interpretation of the frozen primary analysis. | Fixed 8,265-island universe; 106,295 accepted angiosperm denominator; all-analysis and direct-only evidence kept separate. | Primary run `34232450884`, artifact `10058653212`; Figure 2 lock run `34810404509`, artifact `10334412842`. | **SUPPORTED, subject to P2 reconciliation.** | Figure 2 is descriptive visualization of frozen H2 estimates; separate x/y intervals are not a new between-context multivariate test. P2 must keep direct between-context vector inference as the formal basis. |
| C2 | Palearctic accessibility/generalization and reproductive assurance both increase with source separation, while tropical accessibility can move in the opposite direction. | Palearctic all-analysis all-native: accessibility `0.0795`, q=`0.00463`; reproductive assurance `0.0534`, q=`0.00463`. Palearctic direct-only all-native: accessibility `0.0616`, q=`0.0389`; reproductive assurance `0.0966`, q=`8.48e-14`. Tropical direct-only all-native: accessibility `-0.1005`, q=`0.00275`; reproductive assurance `0.1360`, q=`0.00428`. | Same primary response definitions; evidence scopes reported separately. | Primary run `34232450884`; Figure 2 lock `chapter1_v8_figure2_result_lock_v1`. | **SUPPORTED as frozen component estimates.** | Do not infer context difference merely from separate significance. P2 must audit common-island support, covariance and direct vector contrast. Tropical single-axis accessibility is less observation-robust than the Palearctic branch. |
| C3 | The strongest Palearctic island syndrome is retained after family adjustment but not after source-matched genus adjustment. | Support ladder `4/4 -> 4/4 -> 0/4`. Family attenuation `19.6–33.4%`; genus attenuation `78.8–85.9%`; conditional family→genus attenuation `70.6–79.1%`. | Four primary evidence-scope × floristic-stratum cells; 16 frozen attenuation profiles across source modes. | Figure 3 run `34816668470`, artifact `10336875587`, digest `sha256:853fd1c590d86341d6302c1da8e998af93b0b5fa3b35c464dd3762378385bd4d`; source effect-fingerprint run `34796763611`, artifact `10330230959`; source primary run `34232450884`. | **SUPPORTED as taxonomic attenuation; CENTRAL P1 CLAIM remains OPEN to stronger defense.** | Existing result localizes hierarchical expression but does not yet prove that genus attenuation exceeds matched-complexity arbitrary grouping, is independent of all sample-loss paths, or excludes model-flexibility artifacts. P1 must defend those points before upgrading to a strong `measurable assembly depth` claim. |
| C4 | The specified species-list observation bias is a poor explanation for the positive Palearctic accessibility branch, while the tropical branch is more fragile. | V6 survival: Palearctic `99/100`; North–Tropical vector `70/75`; tropical `35/75`. In the biologically concerning remote-under-survey direction, Palearctic survives `80/80`. | Baseline-supported frozen V6 surfaces; original observed information weights preserved. | V6 run `34800498716`, artifact `10331282464`, digest `sha256:095506a214143a2312bcb1d115644ab21d2d130b7813261019bc386206fd3849`; Figure 4 lock. | **SUPPORTED within the frozen sensitivity grid.** | The grid is an assumption set, not an estimate/posterior of true completeness. Arbitrary taxon-dependent omission is not ruled out. P3 may add only prospectively labelled joint sensitivity. |
| C5 | No common nonlinear global assemblage threshold is promoted. | `0/12` observed broad atomic cells pass the frozen all-analysis + direct-only nonlinear promotion rule. | Calibrated design cells with held-out false-nonlinear control; common composite distance exposure. | Geometry run `34765070227`, artifact `10320630085`, digest `sha256:5ac34231b0cfa4acc1d536ab798a574a3cbcdd9423b6f0f2995195aab78a2f23`; Figure 4 lock. | **SUPPORTED negative promotion result.** | Does not imply local/lineage thresholds are absent. Cannot be converted into support for smooth mechanism. P4 cannot reopen this gate. |
| C6 | The global data do not identify a pollinator-specific upstream mechanism. | H5c distance×biotic `+0.0649525`, 95% CI `[-0.0902968, 0.2202019]`, p=`0.4122068` in the sole prospectively qualified cell. N1 and sampled source breadth also do not promote a global channel mechanism. | H5c: direct-only × native-nonendemic × Palearctic qualified cell only. | H5c run `34803837463`, artifact `10332871066`, digest `sha256:fb8927bcf2b6e1c762b3ecbe59c42360faac679cf37b1ad0ed31386631a0f405`; Figure 4 lock. | **SUPPORTED as non-identification, not evidence of no pollinator effect.** | N1 remains unsupported and N2 remains unopened. H5c p=.412 must never be described as evidence that pollinators are irrelevant. |
| C7 | Distributed lineage thresholds are not identifiable against heterogeneous smooth clines in the present global design. | H5d qualified `0/8`; frozen simulation showed insufficient classification accuracy and excessive false distributed-threshold selection under smooth alternatives. | Outcome-closed source-mode × stratum feasibility cells. | H5d run `34803545574`, artifact `10332671123`, digest `sha256:13ea1c3cd0066950d37fa3399f771d41ea2755c9746d02ac18594e224d895539`; Figure 4 lock. | **SUPPORTED identifiability boundary.** | Observed genus threshold distributions remain closed. This motivates local `izu-core` resolution but cannot be inverted into support for either threshold or cline generators. |

## Reconciliation checks completed in P0

1. The canonical Figure 2, 3 and 4 locks all point back to immutable artifacts with explicit digests.
2. The current submission freeze and manuscript use the same `4/4 -> 4/4 -> 0/4`, `78.8–85.9%`, `99/100`, `35/75`, `70/75`, `0/12`, H5c `p=0.41221` and H5d `0/8` anchors.
3. The manuscript value for the source-trained shared floral-architecture factor was already synchronized to `86.85%` all-analysis and `86.44%` direct-only before this P0 ledger.
4. The design-provenance identifier `e4a73796a` supplied outside the repository is not resolvable as a current repository commit and is therefore not used as immutable evidence.

## P0 findings that constrain P1

The main unresolved issue is not whether the frozen genus attenuation exists. It does.

The unresolved issue is whether the magnitude and localization of that attenuation can withstand the strongest artifact explanations:

- changing species/island support across observed → family → genus stages;
- genus-imputed/Validated-Low evidence;
- model flexibility from fine grouping;
- ordinary trait taxonomic autocorrelation that would be reproduced by arbitrary matched partitions;
- spatial clustering not already covered by frozen safeguards.

Therefore P1 must not search for a stronger endpoint. It must attack these five failure modes on the same primary island result.

## P0 decision

**P0 passes as an inventory/reconciliation gate.**

No headline manuscript discrepancy currently requires reopening a frozen empirical result. P1 may proceed, but the strong wording `the floral island syndrome has a measurable assembly depth` remains conditional on the P1 defenses above.
