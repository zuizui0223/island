from __future__ import annotations

import argparse
from pathlib import Path

BASE = Path("docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md")
OUTPUT = Path("docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md")


def _replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"expected exactly one replacement target, found {count}: {old[:80]!r}")
    return text.replace(old, new, 1)


def _replace_section(text: str, heading: str, next_heading: str, body: str) -> str:
    start = text.find(heading)
    if start < 0:
        raise ValueError(f"missing section heading: {heading}")
    end = text.find(next_heading, start + len(heading))
    if end < 0:
        raise ValueError(f"missing following heading: {next_heading}")
    return text[:start] + heading + "\n\n" + body.strip() + "\n\n" + text[end:]


def promote(text: str) -> str:
    text = _replace_once(
        text,
        "## Full working manuscript v10 — island-first, P1/P2-defended draft — 2026-09-15",
        "## Full working manuscript v11 — island-first, P1/P2/P3-defended draft — 2026-09-15",
    )

    old_abstract = (
        "The Palearctic result survived a finite MNAR trait-missingness stress test and a separate "
        "species-list detection tipping analysis: 99/100 baseline-supported Palearctic accessibility "
        "surfaces remained supported, including all 80/80 scenarios in which remote islands were "
        "less complete. By contrast, tropical accessibility was more sensitive. Neither continuous "
        "island area, a common nonlinear breakpoint, coarse pollinator-channel heterogeneity, sampled "
        "source interaction breadth, nor an independently defined biotic-versus-wind pollination "
        "contrast supplied a promoted global mechanism. A prospective distributed-threshold audit "
        "further showed that the present macroecological design cannot reliably distinguish "
        "lineage-specific thresholds from heterogeneous smooth clines."
    )
    new_abstract = (
        "After reproducing V5 trait-resolution MNAR and V6 species-list detection separately on the "
        "pinned input, a prospectively frozen joint observation-bias analysis crossed both processes. "
        "The direct-only native-nonendemic North–Tropical vector difference survived 1,541/1,575 "
        "finite-grid cells, while Palearctic accessibility survived 1,575/1,575 in both evidence "
        "scopes. Tropical accessibility was less stable (1,161/1,575 all-analysis; 1,269/1,575 "
        "direct-only). Deterministic partial-identification bounds preserved the positive Palearctic "
        "accessibility sign, but did not identify support for the formal context contrast across all "
        "corners; tropical accessibility crossed zero. Thus the Palearctic branch is the strongest "
        "observation-defended component, whereas the context contrast is finite-domain robust but only "
        "partially identified and the tropical single axis is observation-fragile. Neither continuous "
        "island area, a common nonlinear breakpoint, nor independent pollination tests supplied a "
        "promoted global mechanism."
    )
    text = _replace_once(text, old_abstract, new_abstract)

    methods_body = """
V5 and V6 attacked different missing-data processes and were first reproduced separately on the same pinned PR142 input before any joint extension was opened. **V5** varied trait-state-dependent resolution of `selfing_core` among recorded flora using the prospectively frozen MNAR odds-ratio grid. **V6** varied distance-dependent flora-list completeness together with state-dependent recording of `generalized_accessible`, while retaining the original observed information weights so hypothetical species never increased regression precision.

**P3** then combined those two frozen processes without re-estimating either from the outcomes. The primary joint surface crossed nine V5 trait-resolution odds ratios with five median-completeness levels, five distance-completeness odds ratios and seven state-recording odds ratios, yielding 1,575 assumption surfaces per evidence scope. V5 acted on reproductive assurance and V6 on accessibility/generalization; unresolved or hypothetical species were not treated as biological zeros. Grid-cell fractions summarize the geometry of this predeclared sensitivity domain and are not posterior probabilities over the true observation process.

A second P3 layer used deterministic partial-identification bounds rather than a fitted latent observation model. Six predeclared V5 extreme assignments were crossed with eight prespecified V6 corner settings, yielding 48 bound surfaces per evidence scope. If a sign crossed zero or formal support was lost anywhere in this envelope, the corresponding claim was labelled observation-fragile or only partially identified rather than rescued by a Bayesian or occupancy-style latent fit.

P3 cannot estimate true island-flora completeness, establish missing-at-random assumptions, or identify a latent biological truth under arbitrary missingness. Its role is narrower: to locate which P1/P2 conclusions survive a prospectively defined joint observation-bias domain and which conclusions depend on stronger assumptions.
"""
    text = _replace_section(
        text,
        "### Missingness and observation-process stress tests",
        "### H5: independent mechanistic tests and claim ceiling",
        methods_body,
    )

    results_body = """
V5 and V6 reproduced their frozen no-differential baselines before the joint analysis was opened; the maximum score discrepancy at the no-differential point was `3.33 × 10^-16`. The P3 finite surface then evaluated 1,575 joint observation assumptions per evidence scope without increasing information weights for hypothetical or imputed species (Fig. 4A,B).

For the formal P2 comparison, the primary direct-only native-nonendemic North–Tropical vector difference remained supported in `1,541/1,575` joint cells (`97.8%`). The 34 failures were concentrated under strongly adversarial combinations, especially the strongest distance-dependent incompleteness (`OR_C=0.25`) together with severe positive-state under-recording and low trait-resolution odds. The formal context contrast is therefore highly stable over the fixed finite domain, but not invariant to every allowed observation mechanism.

The Palearctic accessibility branch was the strongest observation-defended scalar result. In native non-endemics it survived `1,575/1,575` joint cells in both all-analysis and direct-only evidence. Across the deterministic 48-corner envelope, the all-analysis Palearctic slope remained positive and supported (`0.0402` to `0.1381`). The direct-only envelope also remained strictly positive (`0.0392` to `0.1428`), although FDR support was not retained at every extreme corner. Thus the direction is identified across the frozen envelope in both evidence scopes, with full support identification only in the broader ledger (Fig. 4B,C).

Tropical accessibility was markedly more fragile. In native non-endemics it survived `1,161/1,575` all-analysis cells (`73.7%`) and `1,269/1,575` direct-only cells (`80.6%`). Its partial-identification envelopes crossed zero in both scopes (all-analysis `-0.0856` to `0.0399`; direct-only `-0.0942` to `0.0391`) and formal support was not retained across all corners. The tropical single-axis direction is therefore not partially identified under the frozen deterministic envelope.

The strict envelope also narrows the formal context claim. Although the direct-only native-nonendemic North–Tropical vector difference is robust in 97.8% of the finite grid, support is not retained across every deterministic corner. We therefore describe it as **finite-domain robust but partially identified**, not as observation-proof. Grid fractions quantify coverage of the declared assumption set; they are not probabilities that the corresponding biological claim is true.
"""
    text = _replace_section(
        text,
        "### V5 and V6 constrained two different observation-bias explanations",
        "### H5: independent global pollination tests did not identify the upstream mechanism",
        results_body,
    )

    discussion_body = """
P3 makes the asymmetry among contexts explicit rather than smoothing it away. The Palearctic accessibility direction is the observation-robust core: the native-nonendemic result survives every finite joint cell in both evidence scopes, and its sign remains positive across every deterministic bound corner. The formal same-layer North–Tropical vector difference is also highly stable over the finite joint domain (`1,541/1,575` direct-only native-nonendemic cells), but deterministic extremes can remove formal support. It should therefore be described as finite-domain robust but only partially identified under unrestricted corner assignments.

The tropical accessibility component occupies a different category. A substantial fraction of the finite grid breaks, and the deterministic envelope crosses zero. This does not invalidate the broader H2 evidence, because the formal multivariate context comparison contains both primary axes and has separate P2 common-island and common-species defenses. It does mean that the negative tropical accessibility coefficient cannot carry the same rhetorical weight as the positive Palearctic branch.

This separation is useful. Sensitivity analysis should identify where inference depends on unmeasured observation processes rather than convert every result into an equally reassuring robustness statement. We therefore distinguish three labels in the manuscript and Fig. 4: **observation-robust core** for Palearctic accessibility, **finite-domain robust / partially identified** for the formal North–Tropical vector difference, and **observation-fragile** for tropical accessibility. None of those labels estimates the probability of a bias mechanism or a latent true flora.
"""
    text = _replace_section(
        text,
        "### Robustness is asymmetric across contexts",
        "### Why the global pollinator mechanism remains unidentified",
        discussion_body,
    )

    old_conclusion = (
        "The result also survives substantial trait-missingness and species-detection challenges, "
        "while several simpler explanations fail prospective tests."
    )
    new_conclusion = (
        "A prospectively frozen joint observation-bias analysis further localizes that robustness: "
        "Palearctic accessibility is the strongest observation-defended component, the formal "
        "North–Tropical vector difference is highly stable over the finite joint domain but not "
        "identified across all deterministic bounds, and tropical accessibility remains "
        "observation-fragile. Several simpler explanations also fail prospective tests."
    )
    text = _replace_once(text, old_conclusion, new_conclusion)

    old_fig4 = (
        "**Figure 4 | Cross-examination narrows the interpretation without erasing the plant-side "
        "result.** Species-list detection sensitivity is asymmetric across contexts (A): the "
        "Palearctic accessibility branch survives `99/100` baseline-supported bias surfaces, whereas "
        "the tropical accessibility branch survives `35/75`; the direct North–Tropical multivariate "
        "contrast is more robust (`70/75`). Calibrated response-geometry tests require nonlinear "
        "evidence to exceed a design-specific threshold in both all-analysis and direct-only scopes; "
        "none of the 12 observed cells meets that promotion rule (B). An independent biotic-versus-"
        "wind specificity test in the sole prospectively qualified cell yields a positive but "
        "imprecise interaction (`+0.06495`, 95% CI `-0.09030` to `0.22020`, `p=0.41221`) and does not "
        "identify a pollinator-specific global mechanism (C). Finally, distributed lineage thresholds "
        "cannot be distinguished reliably from heterogeneous smooth clines in the realized global "
        "design (`0/8` qualified; D). Negative or non-identified results are not inverted into evidence "
        "that pollinators, local thresholds or other mechanisms are irrelevant."
    )
    new_fig4 = (
        "**Figure 4 | Joint observation bias separates robust pattern from fragile identification.** "
        "The prospectively frozen P3 surface combines V5 trait-resolution MNAR with V6 distance-"
        "dependent list incompleteness and state-dependent species recording without increasing "
        "precision for hypothetical species. Robust/fragile regions are shown over the fixed "
        "assumption domain; grid fractions are not probabilities (A,B). In the primary direct-only "
        "native-nonendemic profile, the formal North–Tropical vector difference survives `1541/1575` "
        "joint cells, Palearctic accessibility survives `1575/1575`, and tropical accessibility "
        "survives `1269/1575`. Deterministic partial-identification bounds preserve the positive "
        "Palearctic accessibility sign but do not preserve formal support for every direct-only corner; "
        "tropical accessibility crosses zero (C). Existing response-geometry and H5 boundaries remain "
        "closed: `0/12` nonlinear shapes promoted, H5c `p=0.412`, and `0/8` distributed-threshold cells "
        "qualified (D). P3 localizes observation robustness; it does not estimate true completeness or "
        "promote a pollination mechanism."
    )
    text = _replace_once(text, old_fig4, new_fig4)

    text = _replace_once(text, "## Claim ceiling for v10", "## Claim ceiling for v11")
    old_ceiling = (
        "It must also state that the exact incremental family-to-genus attenuation is spatially "
        "imprecise, that Palearctic and tropical are not the two labels of one formal direct context "
        "contrast, and that the stronger P2 vector non-collinearity geometry is not established."
    )
    new_ceiling = (
        "It must also state that the exact incremental family-to-genus attenuation is spatially "
        "imprecise, that Palearctic and tropical are not the two labels of one formal direct context "
        "contrast, and that the stronger P2 vector non-collinearity geometry is not established. P3 "
        "further requires that the formal North–Tropical contrast be described as highly finite-domain "
        "robust but not identified across all deterministic joint bounds, and that tropical "
        "accessibility be labelled observation-fragile."
    )
    text = _replace_once(text, old_ceiling, new_ceiling)

    required = [
        "1,541/1,575",
        "1,575/1,575",
        "1,269/1,575",
        "finite-domain robust but partially identified",
        "observation-fragile",
        "## Claim ceiling for v11",
    ]
    for phrase in required:
        if phrase not in text:
            raise ValueError(f"missing required P3 manuscript phrase: {phrase}")
    return text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    promoted = promote(args.base.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(promoted, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
