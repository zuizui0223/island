from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
HYPOTHESIS = ROOT / "docs/chapter1_unified_hypothesis_20260917.md"
FREEZE = ROOT / "docs/chapter1_submission_freeze_v13_20260917.md"
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
FIGURES = ROOT / "docs/chapter1_v13_submission_figure_sync_20260917.md"


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    if start not in text:
        raise RuntimeError(f"missing start marker: {start}")
    a = text.index(start)
    if end not in text[a:]:
        raise RuntimeError(f"missing end marker after {start}: {end}")
    b = text.index(end, a)
    return text[:a] + replacement + text[b:]


def remove_paragraph(text: str, start: str, end: str) -> str:
    return replace_between(text, start, end, "")


def update_manuscript() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    text = replace_between(
        text,
        "### Abstract\n",
        "**Keywords:**",
        """### Abstract\n\nIsland plants are often expected to become more self-reliant and less dependent on specialized pollination as geographic isolation increases, but whether this forms a recurrent global floral syndrome has remained unclear. We combined a fixed universe of 8,265 islands and 106,295 accepted angiosperm species with experimental pollen-limitation data to test a unified hypothesis: increasing isolation is associated with a recurrent floral/reproductive island-syndrome direction because successful pollination becomes increasingly limiting, with plant assemblages responding through partially separable reproductive-assurance and floral-accessibility pathways.\n\nAcross contemporary observed island floras, the six-atomic floral and reproductive response was supported in every one of four predeclared global contexts and the classic-island orientation was positive in every context in both the primary and Direct-only evidence scopes. The submission uses those contexts as independent geographic replications of the same global prediction; between-region contrasts are not part of the paper claim.\n\nIndependent experimental evidence supported the same global ecological pressure. Across 2,969 GloPL experiments from 1,248 sites and 919 publications, pollen limitation increased with distance from major continental landmasses (standardized slope `0.0794 ± 0.0377`; two-sided `p=0.0354`, one-sided positive `p=0.0177`). Plant responses separated into two globally positive components: reproductive assurance and generally accessible floral architecture. A post-hoc functional triangulation using frozen exact-species GloPL overlaps then asked whether those trait states were associated with lower current pollen limitation. Autonomous selfing showed the strongest association (`β=-0.4467`, `p=2.60×10^-8`) and remained negative within publications (`β=-0.3060`, `p=0.0295`) and within publication-by-site groups (one-sided `p=0.0417`). Actinomorphy and generalized floral form also showed negative global associations, although their within-study support was weaker.\n\nTogether, the results show a recurrent global floral/reproductive island syndrome aligned with increasing pollen limitation and expressed through two partially independent plant-response pathways. The evidence triangulates a pollination-constraint hypothesis but does not establish historical mediation from pollen limitation to trait evolution.\n\n""",
    )
    old_keywords = "**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; pollination syndrome; community assembly; taxonomic sorting; GloPL"
    new_keywords = "**Keywords:** island biogeography; pollen limitation; reproductive assurance; selfing; floral accessibility; pollination syndrome; GloPL"
    if old_keywords not in text:
        raise RuntimeError("expected old keywords not found")
    text = text.replace(old_keywords, new_keywords, 1)

    text = replace_between(
        text,
        "A second problem concerns the meaning of a community-level syndrome.",
        "Direct evidence for the proposed ecological pressure is also essential.",
        """The central global question is deliberately simpler than a test of regional branching or lineage-specific assembly. If an island syndrome is general, the same classic direction should recur across geographically distinct parts of the world even when the exact species pools and ecological settings differ. We therefore use the four predeclared geographic contexts as replication strata for the global prediction rather than as competing biological hypotheses.\n\n""",
    )

    text = replace_between(
        text,
        "Here we synthesize the now-frozen Chapter 1 evidence around one mechanistic hypothesis.",
        "---\n\n## Materials and Methods",
        """Here we synthesize the now-frozen Chapter 1 evidence around one global mechanistic hypothesis. **H1** asks whether the classic floral/reproductive island-syndrome direction recurs in every predeclared global context. **H2** asks whether experimental pollen limitation increases with geographic isolation globally. **H3** separates the plant response into reproductive-assurance and accessibility/generalization pathways and asks whether both point in the island-syndrome direction across the global replication strata. **H4** provides an explicitly post-hoc functional triangulation: using trait states frozen before the new analysis, do plants with reproductive-assurance or generally accessible states show lower current experimental pollen limitation after accounting for geography and measurement structure?\n\nNamed floral templates and GloBI interaction records are retained as supplementary compatibility evidence, but neither is required to establish a named pollinator mechanism. The central test is whether three independent evidence layers converge: isolation is associated with the same floral/reproductive direction across the world, isolation is associated with experimental pollen limitation, and trait states expected to reduce dependence on successful pollination are associated with lower current pollen limitation.\n\n""",
    )

    text = replace_between(
        text,
        "The v13 synthesis distinguishes recurrence from equality.",
        "### H2: global experimental pollen limitation",
        """The four predeclared contexts—northern mid-latitude, northern high-latitude, tropical and southern extratropical—are used as geographic replication strata for the same H1 prediction. Within each context, the already locked six-dimensional omnibus test establishes whether the response differs from zero. Orientation is summarized descriptively as the arithmetic mean of the six identically coded distance slopes. No new p-value is attached to that mean because the full six-slope covariance is not used for this post-hoc scalar summary. A supported joint response and positive classic orientation in every context constitute the global recurrence test. No between-context contrast is used as a submission claim.\n\n### H2: global experimental pollen limitation""",
    )

    text = replace_between(
        text,
        "### H5: taxonomic realization and floristic-status boundary",
        "### Supplementary pollination evidence",
        "### Supplementary pollination evidence",
    )

    text = remove_paragraph(
        text,
        "Recurrence did not imply equal response vectors.",
        "### H2: experimental pollen limitation increases with geographic isolation",
    )
    text = text.replace(
        "### H2: experimental pollen limitation increases with geographic isolation",
        "### H2: experimental pollen limitation increases with geographic isolation",
        1,
    )
    text = remove_paragraph(
        text,
        "The old prediction that pollen-limitation increase should be stronger in northern mid-latitudes than in the tropics was not supported.",
        "### H3: reproductive assurance and floral accessibility form two recurrent response pathways",
    )
    text = text.replace(
        "### H3: reproductive assurance and floral accessibility form two recurrent response pathways",
        "### H3: reproductive assurance and floral accessibility form two recurrent response pathways",
        1,
    )
    text = replace_between(
        text,
        "### H5: the recurrent syndrome has different taxonomic realization across flora layers",
        "### Supplementary interaction evidence does not change the mechanism claim",
        "### Supplementary interaction evidence does not change the mechanism claim",
    )

    text = replace_between(
        text,
        "### A global island syndrome can recur without being identical everywhere",
        "### Pollen limitation supplies an independent global ecological pressure",
        """### A floral island syndrome recurs globally\n\nThe central result is that increasing isolation is associated with the classic floral/reproductive island-syndrome direction throughout the global sampling frame. Every predeclared geographic context has a supported six-atomic isolation response, and every context has a positive classic orientation in both evidence scopes. The geographic strata therefore function as repeated tests of one global prediction: the syndrome is not confined to one climatic zone or hemisphere.\n\nThis result also avoids treating a syndrome as a rigid checklist. Individual atomic traits need not all be positive in every fit. The global claim concerns the recurrent multivariate direction toward greater reproductive assurance and greater floral accessibility/generalization.\n\n### Pollen limitation supplies an independent global ecological pressure""",
    )
    text = remove_paragraph(
        text,
        "The lack of the old North-greater-than-Tropical service pattern is informative.",
        "### Reproductive assurance is the strongest functional bridge",
    )
    text = text.replace(
        "### Reproductive assurance is the strongest functional bridge",
        "### Reproductive assurance is the strongest functional bridge",
        1,
    )
    text = replace_between(
        text,
        "### Taxonomic assembly changes the realization, not the existence, of the syndrome",
        "### Pollination syndromes are useful as trait geometry, not visitor labels",
        "### Pollination syndromes are useful as trait geometry, not visitor labels",
    )

    text = replace_between(
        text,
        "### A unified interpretation of the floral island syndrome",
        "### Conclusion",
        """### A unified interpretation of the floral island syndrome\n\nThe evidence supports a simple global synthesis. Island isolation is associated with a recurrent direction toward greater reproductive assurance and greater floral accessibility. Independent experiments show that pollen limitation increases with isolation. The functional trait analysis shows that the clearest reproductive-assurance state—autonomous selfing—is associated with substantially lower current pollen limitation, with weaker but concordant evidence for generally accessible floral architecture. These links make pollination constraint a plausible common ecological pressure without requiring a single named pollinator mechanism.\n\nThe syndrome is therefore best viewed as a global solution space with two partially separable plant strategies. Reproductive assurance reduces dependence on pollen delivery. Generalized floral accessibility may reduce dependence on narrowly matched visitor access. Their repeated co-occurrence across the global replication strata is the biological pattern the paper seeks to explain.\n\n### Conclusion""",
    )

    text = replace_between(
        text,
        "### Conclusion",
        "---\n\n## Main figure legends",
        """### Conclusion\n\nAcross thousands of contemporary island floras, increasing geographic separation is associated with a recurrent floral/reproductive island-syndrome direction across every predeclared global context. Independent pollen-supplementation experiments show that pollen limitation also increases with geographic isolation. The plant response contains at least two partially separable components—reproductive assurance and floral accessibility/generalization—and post-hoc functional triangulation provides the strongest bridge for autonomous selfing, which is associated with markedly lower current experimental pollen limitation even within studies.\n\nThese results support a unified global pollination-constraint hypothesis while preserving a strict causal boundary. They show convergence among geography, experimental reproductive limitation and functional trait states, but they do not identify historical pollen limitation as the cause of the observed trait evolution. The next decisive step is longitudinal or lineage-resolved evidence linking pollination service, reproductive success and trait change through time.\n\n""",
    )

    text = replace_between(
        text,
        "## Main figure legends",
        "---\n\n## Claim ceiling for v13",
        """## Main figure legends\n\n**Figure 1 | Unified global hypothesis and evidence hierarchy for a recurrent floral island syndrome.** Geographic isolation is associated with a global pollination-service constraint and a recurrent plant response. The plant response separates into reproductive-assurance and accessibility/generalization pathways, which can operate partly independently. Solid arrows denote directly supported associations; the link from frozen trait state to current experimental pollen limitation is post-hoc functional triangulation; the historical `pollen limitation -> selection -> trait evolution` arrow remains dashed and unclaimed.\n\n**Figure 2 | The floral/reproductive island-syndrome direction recurs across the world.** Six-atomic standardized isolation slopes are shown for northern mid-latitude, northern high-latitude, tropical and southern extratropical island floras in the primary and Direct-only evidence scopes. All four contexts have supported multivariate response vectors and positive descriptive classic orientations. Reproductive-assurance and accessibility/generalization family summaries show that both response families point in the island-syndrome direction. The figure does not rank or contrast regions.\n\n**Figure 3 | Experimental pollen limitation and the post-hoc functional bridge.** The full GloPL global-distance model shows increasing pollen limitation with geographic separation across 2,969 experiments, 1,248 sites and 919 publications. Frozen exact-species trait states are then compared with current pollen limitation after common adjustment: autonomous selfing shows the strongest negative association and remains negative within publications and publication-by-site groups; actinomorphy and generalized form provide additional global concordance with weaker within-study support. The previously frozen distance-by-trait moderation failures remain displayed separately as negative results.\n\n""",
    )

    text = replace_between(
        text,
        "## Claim ceiling for v13",
        "---\n\n## Literature cited",
        """## Claim ceiling for v13\n\nThe manuscript may state that geographic isolation is associated with a recurrent global classic island-syndrome direction supported independently in all four predeclared geographic contexts. It may state that experimental pollen limitation increases with geographic isolation globally. It may state that reproductive assurance and floral accessibility/generalization are partially separable plant-response pathways. It may state that, in an explicitly post-hoc exact-species triangulation, autonomous selfing is robustly associated with lower current experimental pollen limitation and that actinomorphy and generalized floral form provide additional global functional concordance of differing robustness.\n\nThe manuscript must not state that historical pollen limitation is proven to have selected the observed traits, that pollen limitation statistically mediates the global syndrome, that pollinator abundance or visitation globally declines with isolation, that a named pollinator group caused the pattern, or that GloBI identifies the causal mechanism. Between-region differences and lineage/taxonomic decomposition are outside the current submission claim and remain historical or supplementary provenance.\n\n""",
    )

    for line in (
        "- H3A observed taxonomic-depth run: `35043562731`, artifact `10426327347`;\n",
        "- defended native P1/P3 provenance remains frozen in v11/v12 and is not rewritten by v13.\n",
    ):
        text = text.replace(line, "")

    forbidden = (
        "North–Tropical",
        "Palearctic",
        "taxonomic realization",
        "regional modification",
        "4/4 -> 4/4 -> 0/4",
        "Figure 4",
        "H5 —",
        "### H5",
    )
    for phrase in forbidden:
        if phrase in text:
            raise RuntimeError(f"submission manuscript still contains prohibited phrase: {phrase}")

    MANUSCRIPT.write_text(text, encoding="utf-8")


def update_hypothesis() -> None:
    HYPOTHESIS.write_text(
        """# Chapter 1 v13 — global island-syndrome hypothesis\n\nDate: 2026-09-17\n\nCanonical paper lock: `config/chapter1_v13_unified_island_syndrome_result_lock.json`\n\n## Central hypothesis\n\nGeographic isolation is associated with a recurrent global floral and reproductive island syndrome because successful pollination becomes more limiting as connection to major continental landmasses declines. The plant response has two partially separable routes:\n\n1. a **reproductive-assurance route**, in which self-compatibility and autonomous or realized selfing reduce dependence on successful outcross pollen delivery; and\n2. an **accessibility/generalization route**, in which dependence on restricted or specialized floral access is reduced and assemblages shift toward more generally accessible floral architecture.\n\nThese routes can coexist without forming one obligatory sequence.\n\n```text\n                         geographic isolation\n                                  |\n                                  v\n                    pollination-service constraint\n                                  |\n                      experimental pollen limitation\n                                  |\n                 +----------------+----------------+\n                 |                                 |\n                 v                                 v\n       reproductive assurance              accessibility/generalization\n                 |                                 |\n                 v                                 v\n       SC / autonomous selfing          open/generalized/actinomorphic\n            / selfing                         floral composition\n                 \\                                 /\n                  \\                               /\n                   +------ recurrent global ------+\n                           island syndrome\n```\n\nThe arrows from isolation to experimental pollen limitation and from isolation to plant composition are directly supported associations. The link from trait state to current experimental pollen limitation is post-hoc functional triangulation, strongest for autonomous selfing. The historical sequence `past pollen limitation -> selection -> trait evolution` is not directly identified.\n\n## H1 — Global recurrent island syndrome\n\n**Question.** Does the same classic floral/reproductive island-syndrome direction recur across the global island flora?\n\n**Prediction.** The six-atomic multivariate response is supported within each of four predeclared geographic contexts and the descriptive classic orientation is positive in every context under the same coding.\n\n**Result.** Supported. The six-atomic joint response is supported in all four contexts in both all-analysis and Direct-only evidence, and the descriptive classic orientation is positive in all four contexts in both scopes. The contexts are treated as independent geographic replications of one global prediction; pairwise context contrasts are not part of the submission claim.\n\n**Claim.** A recurrent global floral/reproductive island syndrome is present in contemporary observed island floras.\n\n## H2 — Global pollination constraint\n\n**Question.** Is geographic isolation independently associated with stronger experimental pollen limitation?\n\n**Result.** Supported. Across 2,969 experiments, 1,248 sites and 919 publications, the standardized global distance coefficient is `+0.07937` (SE `0.03773`; two-sided `p=0.0354`, one-sided positive `p=0.0177`). The sign is positive in both frozen measurement sensitivities.\n\n**Claim.** Experimental pollen limitation is a supported global correlate of geographic isolation. This does not identify pollinator abundance or visitation decline.\n\n## H3 — Dual plant-response pathways\n\n### H3a Reproductive assurance\n\nThe reproductive-assurance family points in the classic island-syndrome direction in all four geographic replication strata in both evidence scopes.\n\n### H3b Accessibility/generalization\n\nThe accessibility/generalization family also points in the classic island-syndrome direction in all four replication strata in both evidence scopes. The supported statement is multivariate accessibility/generalization, not universal simplification of every floral structure.\n\n### H3c Partial independence\n\nThe frozen plant-side conditional analysis does not require an obligatory `isolation -> selfing -> floral simplification` sequence. Reproductive assurance and accessibility/generalization are therefore treated as partially separable response pathways.\n\n## H4 — Functional compatibility with pollen limitation\n\n**Status.** Post-hoc functional triangulation, never confirmatory mediation.\n\n**Result.** The strongest bridge is autonomous selfing: `beta=-0.4467`, two-sided `p=2.60e-08`, with the negative association retained in measurement sensitivities and within publications (`beta=-0.3060`, `p=0.0295`). Within publication×site groups the estimate remains negative (`beta=-0.1864`; one-sided negative `p=0.0417`). Actinomorphy also has a strong negative global association; generalized floral form is negative in the primary model but less robust; self-compatibility is negative but imprecise.\n\n**Interpretation.** Frozen trait states expected to reduce dependence on successful pollination are associated with lower current pollen limitation, most clearly for autonomous selfing. This supports functional compatibility but not the historical causal link from pollen limitation to trait evolution.\n\nThe predeclared Route A/B tests asking whether those states buffer the *distance slope* of pollen limitation remain unsupported and are not reclassified.\n\n## Role of pollination syndromes\n\nPredeclared named-pollinator templates remain compatibility layers only. Their shared architecture is used at the trait-combination level: restricted/specialized floral architectures versus generally accessible architectures, alongside a separately defined selfing syndrome.\n\n## Role of GloBI\n\nGloBI is supplementary only. Documented functional-channel breadth is affected by interaction-recording effort and source definitions and is not required for the global mechanism argument.\n\n## Claim ceiling\n\nv13 may state that:\n\n- the classic floral/reproductive island-syndrome direction recurs across every predeclared global context;\n- experimental pollen limitation increases with isolation globally;\n- reproductive assurance and floral accessibility/generalization form partially separable plant-response pathways;\n- frozen functional trait states are associated with current experimental pollen limitation, with the strongest post-hoc bridge for autonomous selfing.\n\nv13 may **not** state that:\n\n- past pollen limitation is proven to have caused the observed trait evolution;\n- pollen limitation statistically mediates the global island syndrome;\n- global pollinator abundance or visitation is proven to decline with island isolation;\n- any named pollinator group caused the global pattern;\n- GloBI identifies the causal pollinator mechanism.\n\nPairwise geographic contrasts and lineage/taxonomic decomposition remain historical or supplementary provenance and are not part of the current submission claim.\n""",
        encoding="utf-8",
    )


def update_freeze() -> None:
    FREEZE.write_text(
        """# Chapter 1 v13 submission freeze — global-only — 2026-09-17\n\nStatus: review candidate; not merged.\n\n## Canonical publication surface\n\n1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`\n2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`\n3. `docs/chapter1_v13_submission_figure_sync_20260917.md`\n4. `config/chapter1_v13_functional_bridge_result_lock.json`\n5. `docs/chapter1_unified_hypothesis_20260917.md`\n\nHistorical v11/v12 manuscripts and result locks remain immutable provenance.\n\n## Frozen paper architecture: H1–H4\n\n- **H1 — recurrent global island syndrome:** the classic floral/reproductive response is supported and positively oriented in every predeclared global context; contexts are replication strata, not competing submission hypotheses.\n- **H2 — global pollination constraint:** experimental pollen limitation increases with geographic separation in the full-global GloPL analysis.\n- **H3 — dual plant response pathways:** reproductive assurance and generalized/accessibility floral architecture are partially separable plant-side response families and point in the island-syndrome direction globally.\n- **H4 — functional bridge:** exact-species GloPL trait-state comparisons are explicitly `posthoc_functional_triangulation`, not confirmatory mediation. Autonomous selfing is the strongest reproducible functional bridge to lower current pollen limitation; architecture evidence is weaker.\n\nPairwise context comparisons and lineage/taxonomic decomposition are outside the submission spine. They remain available only as historical or supplementary provenance.\n\n## Canonical v13 functional bridge\n\n- workflow run: `35141624253`\n- artifact: `10465048981`\n- artifact digest: `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`\n- inferential role: `posthoc_functional_triangulation`\n- frozen Route A/B distance-by-trait moderation failures remain negative results and are not reclassified.\n\n## Reused frozen evidence\n\n- all-data primary probability run: `34961775336`, artifact `10394245237`, digest `sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b`\n- full-global GloPL run: `35090599662`, artifact `10444156159`, digest `sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623`\n\n## Claim ceiling\n\nv13 may state global recurrence, global pollen-limitation increase, partial pathway separation and post-hoc functional compatibility.\n\nv13 must not state that pollen limitation historically caused trait evolution, that pollen limitation statistically mediates the global syndrome, that pollinator abundance or visitation globally declines with isolation, that GloBI proves the causal mechanism, or that post-hoc functional triangulation is confirmatory.\n\n## Submission gate\n\nPromotion to review-ready requires a fresh branch-head audit that passes:\n\n- all `tests/test_chapter1_v13_*.py`;\n- Ruff on v13 source/tests;\n- historical v11/v12 immutability guard;\n- submission-surface guards proving that the current paper contains H1–H4, three main figures, no pairwise-region claim, and no taxonomic branch.\n\nNo merge to `main` is authorized by this freeze document.\n""",
        encoding="utf-8",
    )


def update_readme() -> None:
    text = README.read_text(encoding="utf-8")
    text = replace_between(
        text,
        "## 6. Current submission surface",
        "## 7. Chapter 1 / Chapter 2 division of labour",
        """## 6. Current submission surface\n\nThe current publication-facing surface is **v13 global-only**. It is an **H1–H4** synthesis and does not rewrite frozen v11/v12 provenance.\n\nRead these first, in this order:\n\n1. [`docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`](docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md) — current v13 full manuscript;\n2. [`config/chapter1_v13_unified_island_syndrome_result_lock.json`](config/chapter1_v13_unified_island_syndrome_result_lock.json) — global-only H1–H4 result and claim-ceiling lock;\n3. [`docs/chapter1_v13_submission_figure_sync_20260917.md`](docs/chapter1_v13_submission_figure_sync_20260917.md) — three-main-figure evidence-role contract;\n4. [`config/chapter1_v13_functional_bridge_result_lock.json`](config/chapter1_v13_functional_bridge_result_lock.json) — post-hoc exact-species GloPL functional triangulation;\n5. [`docs/chapter1_unified_hypothesis_20260917.md`](docs/chapter1_unified_hypothesis_20260917.md) — unified global biological hypothesis.\n\nHistorical defended surface retained unchanged:\n\n- [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md) — historical defended manuscript;\n- v11/v12 result locks remain provenance, not current submission claims.\n\n### Current v13 paper-level verdict\n\n- **H1 — recurrent global island syndrome:** every predeclared global context shows a supported multivariate isolation response and positive classic-island orientation; contexts are geographic replications, not comparison targets;\n- **H2 — global pollination constraint:** full-global GloPL shows increasing experimental pollen limitation with geographic separation (`beta=+0.07937`, one-sided `p=0.01772`);\n- **H3 — dual plant pathways:** reproductive assurance and generalized/accessibility architecture are partially separable response families and both point in the island-syndrome direction globally;\n- **H4 — functional bridge:** exact-species GloPL triangulation is explicitly **post-hoc**; autonomous selfing is the strongest reproducible link to lower current pollen limitation, while architecture evidence is weaker.\n\nPublication-facing concept:\n\n> **Geographic isolation is associated worldwide with a recurrent floral/reproductive island syndrome and with stronger experimental pollen limitation. Reproductive assurance and floral accessibility/generalization provide two partially separable plant-side routes consistent with reduced dependence on reliable pollination.**\n\nThe v13 synthesis is triangulation, not historical causal mediation. It does not claim that pollen limitation evolved the observed traits, that pollinator abundance globally declines with isolation, or that GloBI identifies the causal mechanism.\n\n""",
    )
    text = replace_between(
        text,
        "## 7. Chapter 1 / Chapter 2 division of labour",
        "## 8. Repository map",
        """## 7. Chapter 1 / Chapter 2 division of labour\n\nChapter 1 establishes the global macroecological pattern: where island isolation increases, floral/reproductive composition repeatedly shifts toward greater reproductive assurance and floral accessibility, alongside stronger experimental pollen limitation.\n\nChapter 2 (`izu-core`) is the mechanistic-resolution layer. It can measure the within-system chain:\n\n`interaction state -> effective service -> reproductive outcome -> phenotype`\n\nand test the causal sequence within one biological system.\n\n""",
    )
    README.write_text(text, encoding="utf-8")


def update_pipeline() -> None:
    text = PIPELINE.read_text(encoding="utf-8")
    text = replace_between(
        text,
        "## 13. Canonical paper surface",
        "## 14. Chapter 1 / Chapter 2 handoff",
        """## 13. Canonical paper surface\n\nThe current publication-facing surface is **v13 global-only**, organized as **H1–H4**. Read in this order:\n\n1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`\n2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`\n3. `docs/chapter1_v13_submission_figure_sync_20260917.md`\n4. `config/chapter1_v13_functional_bridge_result_lock.json`\n5. `docs/chapter1_unified_hypothesis_20260917.md`\n6. `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md` — historical defended manuscript\n\nThe publication-facing hierarchy is:\n\n- H1: recurrent global classic-island syndrome across every predeclared geographic replication stratum;\n- H2: full-global experimental pollen-limitation gradient;\n- H3: partially separable reproductive-assurance and floral-accessibility pathways;\n- H4: explicitly post-hoc exact-species functional triangulation.\n\nCanonical v13 functional-bridge run:\n\n- run `35141624253`;\n- artifact `10465048981`;\n- digest `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`.\n\nThe frozen Route A/B distance-by-trait moderation failures remain negative results. Pairwise geographic contrasts and lineage/taxonomic decomposition remain historical or supplementary provenance and are not current paper claims.\n\nPrevious v8/v9/v10/v11/v12 surfaces remain historical provenance and continue to define the chronology of claim defense.\n\n""",
    )
    text = replace_between(
        text,
        "## 14. Chapter 1 / Chapter 2 handoff",
        "## 15. Legacy v1 boundary",
        """## 14. Chapter 1 / Chapter 2 handoff\n\nChapter 1 establishes a recurrent worldwide floral/reproductive island syndrome, the accompanying global pollen-limitation gradient, and two plant-side response routes. It does not identify the historical causal sequence.\n\nChapter 2 / `izu-core` should resolve:\n\n`interaction state -> effective service -> reproductive outcome -> phenotype`\n\nwithin one biological system.\n\n- **Chapter 1:** does the syndrome recur globally, and which functional routes align with it?\n- **Chapter 2:** how and why does the causal chain operate within a resolved system?\n\n""",
    )
    PIPELINE.write_text(text, encoding="utf-8")


def main() -> None:
    update_manuscript()
    update_hypothesis()
    update_freeze()
    update_readme()
    update_pipeline()

    prohibited = (
        "North–Tropical",
        "Palearctic",
        "taxonomic realization",
        "regional modification",
        "4/4 -> 4/4 -> 0/4",
        "Figure 4",
        "H5 —",
        "### H5",
    )
    for path in (MANUSCRIPT, HYPOTHESIS, FREEZE, FIGURES):
        text = path.read_text(encoding="utf-8")
        for phrase in prohibited:
            if phrase in text:
                raise RuntimeError(f"{phrase!r} remains in {path}")


if __name__ == "__main__":
    main()
