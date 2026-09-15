from __future__ import annotations

from pathlib import Path


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"expected exactly one occurrence: {old[:80]!r}, found {text.count(old)}")
    return text.replace(old, new, 1)


def update_readme() -> None:
    path = Path("README.md")
    text = path.read_text(encoding="utf-8")
    text = replace_once(text, "  -> P0 provenance + P1 assembly-depth defense\n", "  -> P0 provenance + P1 assembly defense + P2 component-support defense\n")
    text = replace_once(text, "  -> canonical island-first v9 manuscript", "  -> canonical island-first v10 manuscript")
    text = replace_once(text, "## 5. P0/P1 defense of the assembly inference", "## 5. P0/P1/P2 defense of the island-syndrome inference")
    anchor = "- [`config/chapter1_p1_final_decision_result_lock.json`](config/chapter1_p1_final_decision_result_lock.json)\n\n"
    p2 = """- [`config/chapter1_p1_final_decision_result_lock.json`](config/chapter1_p1_final_decision_result_lock.json)\n\n### P2 — common-support defense of the H2 component contrast\n\nThe formal direct H2 comparison is `northern_midlatitude` versus `tropical` within the `analysis_regime` layer. `Palearctic` remains a separate within-context realm result and the focal H3/P1 genus-structure system.\n\n- run **34945548775**;\n- artifact **10387261614**;\n- direct-only native-nonendemic common-island North–Tropical vector difference: **p=0.000947**;\n- 2,000 paired spatial-block draws: determinant CI includes zero, so strong vector non-collinearity is **not established**;\n- common-species sensitivity: **853** direct-only co-observed species, native-nonendemic joint vector difference **p=0.00517**;\n- frozen Palearctic–Neotropical direct tests remain unsupported (`p=0.078–0.396`).\n\nCanonical lock:\n\n- [`config/chapter1_p2_component_nonconcordance_result_lock.json`](config/chapter1_p2_component_nonconcordance_result_lock.json)\n\n"""
    text = replace_once(text, anchor, p2)
    text = replace_once(text, "1. [`docs/chapter1_submission_freeze_20260915_p1_defended.md`](docs/chapter1_submission_freeze_20260915_p1_defended.md) — current v9 submission state and claim ceiling;", "1. [`docs/chapter1_submission_freeze_20260915_p1_p2_defended.md`](docs/chapter1_submission_freeze_20260915_p1_p2_defended.md) — current v10 submission state and claim ceiling;")
    text = replace_once(text, "2. [`docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md`](docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md) — canonical island-first v9 manuscript;", "2. [`docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`](docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md) — canonical island-first v10 manuscript;")
    text = replace_once(text, "3. [`docs/chapter1_v9_submission_figure_sync_20260915.md`](docs/chapter1_v9_submission_figure_sync_20260915.md) — final panel mapping and figure-reference contract;", "3. [`docs/chapter1_v10_submission_figure_sync_20260915.md`](docs/chapter1_v10_submission_figure_sync_20260915.md) — final panel mapping and figure-reference contract;")
    text = replace_once(text, "4. [`docs/chapter1_p1_final_decision_20260915.md`](docs/chapter1_p1_final_decision_20260915.md) — integrated P1a/P1c/P1d decision;", "4. [`docs/chapter1_p2_component_nonconcordance_result_20260915.md`](docs/chapter1_p2_component_nonconcordance_result_20260915.md) — P2 common-support and claim-boundary result;\n5. [`docs/chapter1_p1_final_decision_20260915.md`](docs/chapter1_p1_final_decision_20260915.md) — integrated P1a/P1c/P1d decision;")
    text = text.replace("5. [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`]", "6. [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`]", 1)
    text = text.replace("6. [`docs/chapter1_literature_positioning_20260909.md`]", "7. [`docs/chapter1_literature_positioning_20260909.md`]", 1)
    text = replace_once(text, "- [`docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md`](docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md)", "- [`docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`](docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md)")
    text = replace_once(text, "v8 remains available as the previous submission surface and should not be used to quote the family→genus increment without the P1d caveat.", "v9 remains available as the previous P1-defended surface. v10 additionally separates the formal North–Tropical H2 contrast from the Palearctic H3/P1 branch and must be used for current quoting.")
    text = replace_once(text, "- **Figure 2:** `config/chapter1_v8_figure2_result_lock.json` — primary biogeographic branching;", "- **Figure 2:** `config/chapter1_v10_figure2_p2_result_lock.json` — formal same-layer North–Tropical contrast plus P2 common-support defense;")
    text = replace_once(text, "- **H2:** source-separation responses branch among biogeographic contexts and trait components can decouple;", "- **H2/P2:** the formal same-layer North–Tropical joint response difference survives common-island and common-species restrictions; strong geometric non-collinearity is not precisely established;")
    text = replace_once(text, "`assembly depth` remains useful as a localization concept, but v9 does not present the family→genus increment as a perfectly sharp taxonomic breakpoint.", "`assembly depth` remains useful as a localization concept, but v10 neither presents the family→genus increment as a perfectly sharp taxonomic breakpoint nor treats Palearctic and tropical as labels of one formal direct contrast.")
    text = replace_once(text, "docs/PAPER_PIPELINE.md\n  <- database -> H1-H5 -> robustness/falsification -> P0/P1 -> v9", "docs/PAPER_PIPELINE.md\n  <- database -> H1-H5 -> robustness/falsification -> P0/P1/P2 -> v10")
    text = replace_once(text, "config/chapter1_p1_final_decision_result_lock.json\n", "config/chapter1_p1_final_decision_result_lock.json\nconfig/chapter1_p2_component_nonconcordance_result_lock.json\n")
    text = replace_once(text, "config/chapter1_v8_figure2_result_lock.json\n", "config/chapter1_v10_figure2_p2_result_lock.json\n")
    text = replace_once(text, "docs/chapter1_submission_freeze_20260915_p1_defended.md\ndocs/chapter1_v9_submission_figure_sync_20260915.md\ndocs/chapter1_p1_final_decision_20260915.md\ndocs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md", "docs/chapter1_submission_freeze_20260915_p1_p2_defended.md\ndocs/chapter1_v10_submission_figure_sync_20260915.md\ndocs/chapter1_p2_component_nonconcordance_result_20260915.md\ndocs/chapter1_p1_final_decision_20260915.md\ndocs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md")
    path.write_text(text, encoding="utf-8")


def update_pipeline() -> None:
    path = Path("docs/PAPER_PIPELINE.md")
    text = path.read_text(encoding="utf-8")
    text = text.replace("current v9 paper", "current v10 paper", 1)
    text = text.replace("court-style evidence ledger\n                    |\n                    v\n[10] canonical v9 manuscript", "P0/P1/P2 defense\n                    |\n                    v\n[10] canonical v10 manuscript", 1)
    old = "The manuscript-facing evidence table is:\n\n- `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`"
    new = """The submission-facing defense now includes:\n\n- `docs/chapter1_p0_claim_ledger_20260915.md`\n- `config/chapter1_p1_final_decision_result_lock.json`\n- `config/chapter1_p2_component_nonconcordance_result_lock.json`\n- `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`"""
    text = replace_once(text, old, new)
    text = text.replace("The strongest claim surviving cross-examination is:\n\n> geographic isolation is associated with reproducible but biogeographically contingent floral/reproductive assemblage change, and the strongest Palearctic syndrome is principally expressed at the family-to-genus assembly transition.", "The strongest claim surviving cross-examination is:\n\n> the formal North–Tropical joint plant response differs under common observation support, while the strongest Palearctic within-context syndrome is strongly structured by real genus composition beyond matched arbitrary grouping complexity. The exact vector geometry and the exact family→genus increment remain imprecise.", 1)
    text = replace_once(text, "## 9. Canonical paper surface", "## 9. Canonical paper surface")
    text = text.replace("1. `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`\n2. `docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md`\n3. `docs/chapter1_literature_positioning_20260909.md`\n4. `docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md`", "1. `docs/chapter1_submission_freeze_20260915_p1_p2_defended.md`\n2. `docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`\n3. `docs/chapter1_v10_submission_figure_sync_20260915.md`\n4. `docs/chapter1_p2_component_nonconcordance_result_20260915.md`\n5. `docs/chapter1_p1_final_decision_20260915.md`", 1)
    text = text.replace("Previous v7 and 2026-09-09 framing documents remain historical development surfaces and git provenance; they no longer define the canonical narrative.", "Previous v8/v9 and 2026-09-09 framing documents remain historical development surfaces and git provenance; they no longer define the canonical narrative.", 1)
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    update_readme()
    update_pipeline()
