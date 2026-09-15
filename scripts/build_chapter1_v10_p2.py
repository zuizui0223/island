from __future__ import annotations

import re
from pathlib import Path

SOURCE = Path("docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md")
TARGET = Path("docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md")


def sub_once(text: str, pattern: str, replacement: str, *, flags: int = 0) -> str:
    out, n = re.subn(pattern, replacement, text, count=1, flags=flags)
    if n != 1:
        raise RuntimeError(f"expected one replacement, got {n}: {pattern[:80]}")
    return out


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = text.replace(
        "## Full working manuscript v9 — island-first, P1-defended draft — 2026-09-15",
        "## Full working manuscript v10 — island-first, P1/P2-defended draft — 2026-09-15",
        1,
    )

    abstract_block = """A universal floral island syndrome was not recovered. The formal same-layer H2 comparison is the **northern-midlatitude versus tropical** contrast within the frozen analysis-regime layer. In the post-baseline P2 audit, the direct-only native-nonendemic joint vector difference survived forcing accessibility/generalization and reproductive assurance onto exactly the same islands (`p=0.000947`; 348 islands, 78 spatial blocks) and remained supported when both axes were rebuilt from the same 853 co-observed species (`p=0.00517`). The stronger geometric claim was not established: the paired-block 95% interval for the 2D determinant included zero. Palearctic is a separate biogeographic-realm context; it provides the clearest within-context branch and the focal H3/P1 genus-structure result, but is not relabelled as the northern side of the formal North–Tropical direct contrast.

The broad Palearctic response persisted in native non-endemics and through family adjustment, but the predeclared multivariate gate failed after source-matched genus adjustment in all four evidence-scope-by-floristic-stratum combinations (`4/4 -> 4/4 -> 0/4`). Genus adjustment attenuated 78.8–85.9% of the observed response-vector magnitude. A post-baseline, outcome-blind matched-complexity test showed that true genus boundaries produced stronger conditional attenuation than arbitrary within-family partitions preserving the exact genus-count and genus-size structure (observed median 0.721; null median 0.224; one-sided randomization `p=0.029`, 2,000 permutations). However, paired spatial-block bootstrap intervals for the *additional* family-to-genus attenuation included zero in all eight direct-only primary profiles, so the exact incremental depth is not precisely localized.
"""
    text = sub_once(
        text,
        r"A universal floral island syndrome was not recovered\..*?The Palearctic result survived",
        abstract_block + "\nThe Palearctic result survived",
        flags=re.S,
    )

    p2_methods = """### P2 post-baseline component non-concordance audit

P2 did not redefine H2. It tested whether the frozen primary two-axis response difference could be an artefact of mismatched observation support. The formal direct comparison was fixed as `northern_midlatitude` versus `tropical` within the `analysis_regime` layer. `Palearctic` belongs to the separate `biogeographic_realm` layer and was prohibited from being substituted as the northern label of this direct contrast.

**P2a** restricted both primary axes to the intersection of islands with finite scores for both axes within each evidence scope and floristic stratum, then refit the frozen distance, climate, area and spatial-block model. A 2,000-draw paired spatial-block bootstrap propagated uncertainty in component slopes, the between-context difference, vector angle and the 2D determinant. **P2b** imposed a stronger denominator restriction by rebuilding both axes from only species with finite species-level scores for both `generalized_accessible` and `selfing_core`; this changed the estimand and was treated as sensitivity analysis rather than replacement H2. **P2c** retained the frozen Palearctic–Neotropical direct tests only as a claim boundary. Strong geometric non-concordance required the determinant interval to exclude zero; a significant joint vector difference alone was not sufficient.

"""
    marker = "### H3: source-matched taxonomic structure and P1 defense"
    if marker not in text:
        raise RuntimeError("H3 marker missing")
    text = text.replace(marker, p2_methods + marker, 1)

    h2_results = """### H2/P2: the same-layer North–Tropical joint response difference survives common support

The formal H2 direct comparison is northern-midlatitude versus tropical within the frozen analysis-regime layer, not Palearctic versus tropical. P2a forced accessibility/generalization and reproductive assurance onto exactly the same islands before refitting. In the primary direct-only native-nonendemic profile, the northern-midlatitude vector was `(0.0254, 0.0374)` and the tropical vector was `(-0.1039, 0.1747)` for `(accessibility/generalization, reproductive assurance)`. The direct joint vector difference remained supported (`p=0.000947`; 348 islands, 78 spatial blocks). The direct-only all-native profile also retained support (`p=0.0172`), as did the all-analysis native-nonendemic profile (`p=0.000229`); the all-analysis all-native profile was weaker (`p=0.177`). Thus the primary direct comparison is not generated simply by using different islands for the two components (Fig. 2A).

The stronger geometric claim did not pass. In the direct-only native-nonendemic profile the 2D determinant was `0.00746`, but its paired spatial-block 95% interval was `[-0.00319, 0.01599]`. The estimated vector angle was `67.1°`, with a wide 95% interval of `12.3–169.2°`. All four determinant intervals included zero. We therefore retain a **joint North–Tropical vector difference**, but do not claim that the two context vectors are precisely demonstrated to be non-collinear rather than noisy scaled or rotated alternatives (Fig. 2B).

P2b imposed the stronger common-species denominator. Among direct-only evidence, 853 species had finite scores for both primary source axes. Rebuilding both island responses from those same species retained the North–Tropical joint vector difference in native non-endemics (`p=0.00517`; 280 islands) and all natives (`p=0.00395`; 284 islands). Under this changed estimand, the direct-only native-nonendemic northern vector was `(0.0367, -0.0183)` and the tropical vector was `(-0.0668, 0.2073)`. Differing species denominators therefore do not explain away the joint context difference, but the sensitivity analysis also shows why individual within-context component slopes should not be treated as invariant constants across denominator definitions (Fig. 2C).

The Palearctic remains biologically important but plays a different inferential role. It is the clearest within-context realm branch and the focal H3/P1 genus-structure system. Its accessibility/generalization and reproductive-assurance slopes are positive in the frozen primary analysis, and the Palearctic attraction/access component remains positive after conditioning on `selfing_core` across four fixed source definitions. However, the frozen **Palearctic–Neotropical direct vector tests are unsupported** in all four primary scope-by-stratum profiles (`p=0.078–0.396`). The manuscript therefore does not use “Palearctic versus tropical” as a formal direct H2 contrast (Fig. 2D).

Secondary large-bee-like, butterfly-like, and bird-like templates differed among regions, but a source-trained common factor explained 86.85% of their variance in all-analysis evidence and 86.44% in direct-only evidence. These named templates therefore mainly captured overlapping plant architecture. They were not treated as evidence for realized visitor identity.

"""
    text = sub_once(
        text,
        r"### H2: reproductive assurance and floral architecture branched by biogeographic context\n\n.*?### H3: the Palearctic floral-island response is genus-structured beyond matched grouping complexity",
        h2_results + "### H3: the Palearctic floral-island response is genus-structured beyond matched grouping complexity",
        flags=re.S,
    )

    text = text.replace(
        "First, the response is not globally uniform: Palearctic and tropical assemblages can combine reproductive and floral changes differently.",
        "First, the formal North–Tropical response vectors differ directly within the same analysis-regime layer, and that difference survives common-island and common-species restrictions; the stronger non-collinearity geometry remains imprecise.",
    )
    text = text.replace(
        "The decoupling of reproductive assurance and floral accessibility is equally important. In the Palearctic, both components increase with separation. In tropical assemblages, reproductive assurance can increase while floral accessibility decreases. This means the classical verbal package—less reliable pollination, more selfing, simpler flowers—cannot be assumed to operate as one serial pathway globally.",
        "The component comparison is equally important, but it must respect context layers. The formal same-layer North–Tropical joint vector difference survives common support, while its exact 2D non-collinearity remains uncertain. Separately, the frozen Palearctic within-context branch has positive accessibility/generalization and reproductive-assurance slopes, whereas tropical within-context estimates can combine increasing reproductive assurance with decreasing accessibility. These results reject a requirement that the classical verbal package—less reliable pollination, more selfing, simpler flowers—must recur as one fixed scalar response everywhere, without claiming a precisely estimated universal geometry.",
    )
    text = text.replace(
        "The tropical result remains useful because the multivariate North–Tropical contrast is more robust than the single axis, but the two branches should not be presented as equally defended.",
        "The tropical result remains useful because the formal multivariate North–Tropical contrast is more robust than the tropical single axis, but that direct contrast must not be conflated with the separate Palearctic realm result.",
    )

    fig2 = """**Figure 2 | The formal North–Tropical floral/reproductive response difference survives common observation support.** The P2 post-baseline audit uses the frozen `analysis_regime` contrast, northern-midlatitude versus tropical, and never substitutes the separate Palearctic realm result as one side of that test. Common-island refits show the two-axis response vectors across evidence scopes and floristic strata (A); in the primary direct-only native-nonendemic profile the joint vector difference remains supported (`p=0.000947`). Paired spatial-block uncertainty shows that the stronger 2D non-collinearity claim is not established: the determinant interval includes zero and the vector-angle interval is broad (B). Rebuilding both axes from the same 853 co-observed direct-only species retains the native-nonendemic North–Tropical joint difference (`p=0.00517`; C), showing that differing species denominators are not sufficient to explain the context contrast. Frozen Palearctic–Neotropical direct tests remain unsupported (`p=0.078–0.396`; D), so Palearctic is treated as a separate within-context branch and the focal H3/P1 genus-structure system rather than as the northern label of the formal H2 contrast.
"""
    text = sub_once(
        text,
        r"\*\*Figure 2 \|.*?(?=\n\n\*\*Figure 3 \|)",
        fig2.rstrip(),
        flags=re.S,
    )

    text = text.replace("## Claim ceiling for v9", "## Claim ceiling for v10", 1)
    old_ceiling = "The manuscript may state that the strongest Palearctic floral/reproductive syndrome is strongly structured by source-matched genus composition, that true genus boundaries attenuate the response more than matched arbitrary within-family partitions of identical grouping complexity, and that trait-syndrome direction and taxonomic representation depend on biogeographic context. It may state that several simple alternatives fail explicit falsification tests."
    new_ceiling = "The manuscript may state that the formal same-layer North–Tropical joint floral/reproductive response difference survives common-island and common-species denominator safeguards, while strong vector non-collinearity is not precisely established. It may separately state that the strongest Palearctic within-context floral/reproductive syndrome is strongly structured by source-matched genus composition and that true genus boundaries attenuate the response more than matched arbitrary within-family partitions of identical grouping complexity. It may state that several simple alternatives fail explicit falsification tests."
    if old_ceiling not in text:
        raise RuntimeError("claim ceiling paragraph missing")
    text = text.replace(old_ceiling, new_ceiling, 1)
    text = text.replace(
        "It must also state that the exact incremental family-to-genus attenuation is spatially imprecise.",
        "It must also state that the exact incremental family-to-genus attenuation is spatially imprecise, that Palearctic and tropical are not the two labels of one formal direct context contrast, and that the stronger P2 vector non-collinearity geometry is not established.",
        1,
    )

    TARGET.write_text(text, encoding="utf-8")
    print(TARGET)


if __name__ == "__main__":
    main()
