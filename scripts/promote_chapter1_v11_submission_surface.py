from __future__ import annotations

import json
from pathlib import Path

README = Path("README.md")
PIPELINE = Path("docs/PAPER_PIPELINE.md")
MANUSCRIPT = Path(
    "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md"
)
P3_LOCK = Path("config/chapter1_p3_joint_observation_bias_result_lock.json")
FIG4_LOCK = Path("config/chapter1_v11_figure4_result_lock.json")
FIG_SYNC = Path("docs/chapter1_v11_submission_figure_sync_20260915.md")
FREEZE = Path("docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md")


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"expected one replacement, found {count}: {old[:90]!r}")
    return text.replace(old, new, 1)


def require_inputs() -> tuple[dict, dict]:
    for path in (MANUSCRIPT, P3_LOCK, FIG4_LOCK):
        if not path.is_file():
            raise FileNotFoundError(path)
    p3 = json.loads(P3_LOCK.read_text(encoding="utf-8"))
    fig4 = json.loads(FIG4_LOCK.read_text(encoding="utf-8"))
    if p3.get("contract") != "chapter1_p3_joint_observation_bias_result_lock_v1":
        raise ValueError("unexpected P3 lock")
    if fig4.get("contract") != "chapter1_v11_figure4_result_lock_v1":
        raise ValueError("unexpected Figure 4 lock")
    if not fig4.get("visual_review", {}).get("performed"):
        raise ValueError("Figure 4 must be visually reviewed before surface promotion")
    return p3, fig4


def update_readme(text: str) -> str:
    text = replace_once(
        text,
        "  -> P0 provenance + P1 assembly defense + P2 component-support defense\n"
        "  -> court-style evidence ledger\n"
        "  -> locked Figures 1-4\n"
        "  -> canonical island-first v10 manuscript",
        "  -> P0 provenance + P1 assembly defense + P2 component-support defense\n"
        "  -> P3 joint observation-bias / partial-identification defense\n"
        "  -> court-style evidence ledger\n"
        "  -> locked Figures 1-4\n"
        "  -> canonical island-first v11 manuscript",
    )
    insertion = """
### P3 — joint observation-bias and partial-identification defense

V5 trait-resolution MNAR and V6 species-list detection were first reproduced separately on the pinned PR142 input. Only then were the two processes crossed under a prospectively frozen joint contract.

- canonical P3 run: **34949880409**;
- artifact: **10389197309**;
- finite joint surfaces: **1,575 per evidence scope**;
- primary direct-only native-nonendemic North–Tropical vector: **1,541/1,575** robust cells;
- Palearctic accessibility, native non-endemics: **1,575/1,575** robust in both evidence scopes;
- tropical accessibility, native non-endemics: **1,161/1,575** all-analysis and **1,269/1,575** direct-only;
- deterministic partial-identification envelope preserves the positive Palearctic accessibility sign but does not identify formal North–Tropical support across every corner; tropical accessibility crosses zero.

Grid-cell fractions describe the declared sensitivity domain and are **not probabilities**. P3 estimates neither true flora completeness nor an arbitrary-MNAR latent truth.

Canonical lock:

- [`config/chapter1_p3_joint_observation_bias_result_lock.json`](config/chapter1_p3_joint_observation_bias_result_lock.json)

"""
    marker = "## 6. Current submission surface\n"
    if marker not in text:
        raise ValueError("README current-submission marker missing")
    text = text.replace(marker, insertion + marker, 1)

    old_surface = """1. [`docs/chapter1_submission_freeze_20260915_p1_p2_defended.md`](docs/chapter1_submission_freeze_20260915_p1_p2_defended.md) — current v10 submission state and claim ceiling;
2. [`docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`](docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md) — canonical island-first v10 manuscript;
3. [`docs/chapter1_v10_submission_figure_sync_20260915.md`](docs/chapter1_v10_submission_figure_sync_20260915.md) — final panel mapping and figure-reference contract;
4. [`docs/chapter1_p2_component_nonconcordance_result_20260915.md`](docs/chapter1_p2_component_nonconcordance_result_20260915.md) — P2 common-support and claim-boundary result;
5. [`docs/chapter1_p1_final_decision_20260915.md`](docs/chapter1_p1_final_decision_20260915.md) — integrated P1a/P1c/P1d decision;
6. [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`](docs/chapter1_court_evidence_and_theory_synthesis_20260914.md) — historical court-style evidence ledger;
7. [`docs/chapter1_literature_positioning_20260909.md`](docs/chapter1_literature_positioning_20260909.md) — frozen literature-positioning note."""
    new_surface = """1. [`docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md`](docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md) — current v11 submission state and claim ceiling;
2. [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md) — canonical island-first v11 manuscript;
3. [`docs/chapter1_v11_submission_figure_sync_20260915.md`](docs/chapter1_v11_submission_figure_sync_20260915.md) — final v11 panel mapping and figure-reference contract;
4. [`docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md`](docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md) — P3 finite-grid and partial-identification result;
5. [`docs/chapter1_p2_component_nonconcordance_result_20260915.md`](docs/chapter1_p2_component_nonconcordance_result_20260915.md) — P2 common-support and claim-boundary result;
6. [`docs/chapter1_p1_final_decision_20260915.md`](docs/chapter1_p1_final_decision_20260915.md) — integrated P1a/P1c/P1d decision;
7. [`docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`](docs/chapter1_court_evidence_and_theory_synthesis_20260914.md) — historical court-style evidence ledger;
8. [`docs/chapter1_literature_positioning_20260909.md`](docs/chapter1_literature_positioning_20260909.md) — frozen literature-positioning note."""
    text = replace_once(text, old_surface, new_surface)
    text = replace_once(
        text,
        "- [`docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`](docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md)",
        "- [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md)",
    )
    text = replace_once(
        text,
        "v9 remains available as the previous P1-defended surface. v10 additionally separates the formal North–Tropical H2 contrast from the Palearctic H3/P1 branch and must be used for current quoting.",
        "v10 remains available as the previous P1/P2-defended surface. v11 additionally localizes joint observation robustness through P3 and must be used for current quoting.",
    )
    text = replace_once(
        text,
        "- **Figure 4:** `config/chapter1_v8_figure4_result_lock.json` — falsification and claim boundaries.",
        "- **Figure 4:** `config/chapter1_v11_figure4_result_lock.json` — joint observation-bias robustness, partial identification, and retained mechanism boundaries.",
    )
    text = replace_once(
        text,
        "- **V5/V6:** the Palearctic core survives strong trait-missingness and specified species-detection challenges, while tropical accessibility is less robust;",
        "- **P3:** the Palearctic accessibility branch is the observation-robust core; the formal North–Tropical vector is highly finite-domain robust but only partially identified, while tropical accessibility is observation-fragile;",
    )
    text = replace_once(text, "`assembly depth` remains useful as a localization concept, but v10", "`assembly depth` remains useful as a localization concept, but v11")
    text = replace_once(
        text,
        "  <- database -> H1-H5 -> robustness/falsification -> P0/P1/P2 -> v10",
        "  <- database -> H1-H5 -> robustness/falsification -> P0/P1/P2/P3 -> v11",
    )
    text = replace_once(
        text,
        "config/chapter1_p2_component_nonconcordance_result_lock.json\nconfig/chapter1_v8_figure1_result_lock.json",
        "config/chapter1_p2_component_nonconcordance_result_lock.json\nconfig/chapter1_p3_joint_observation_bias_result_lock.json\nconfig/chapter1_v8_figure1_result_lock.json",
    )
    text = replace_once(text, "config/chapter1_v8_figure4_result_lock.json", "config/chapter1_v11_figure4_result_lock.json")
    text = replace_once(
        text,
        "docs/chapter1_submission_freeze_20260915_p1_p2_defended.md\ndocs/chapter1_v10_submission_figure_sync_20260915.md",
        "docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md\ndocs/chapter1_v11_submission_figure_sync_20260915.md",
    )
    text = replace_once(
        text,
        "docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md",
        "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md",
    )
    return text


def update_pipeline(text: str) -> str:
    text = replace_once(
        text,
        "current **island-first, P1/P2-defended v10 paper**",
        "current **island-first, P1/P2/P3-defended v11 paper**",
    )
    text = replace_once(
        text,
        "[11] P2 component-support defense\n     same islands -> paired blocks -> same species denominator\n                    |\n                    v\n[12] canonical island-first v10 manuscript",
        "[11] P2 component-support defense\n     same islands -> paired blocks -> same species denominator\n                    |\n                    v\n[12] P3 joint observation-bias defense\n     V5 reproduction -> V6 reproduction -> joint surface -> partial identification\n                    |\n                    v\n[13] canonical island-first v11 manuscript",
    )
    p3 = """
## 11. P3 — jointly bound observation bias

Canonical result:

- `config/chapter1_p3_joint_observation_bias_result_lock.json`
- `docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md`
- Run `34949880409` / artifact `10389197309`.

Execution order was fail-closed: V5 reproduced first, V6 reproduced second, both were reconciled to the frozen PR142 baseline, and only then was the joint surface opened. The finite domain contains 1,575 V5×V6 assumption surfaces per evidence scope and never increases regression precision for hypothetical species.

Primary direct-only native-nonendemic results:

- North–Tropical vector difference: `1541/1575` robust cells;
- Palearctic accessibility: `1575/1575` robust cells;
- tropical accessibility: `1269/1575` robust cells.

The deterministic 48-corner partial-identification envelope preserves the positive Palearctic accessibility sign, but support for the formal North–Tropical vector is not identified across every corner and tropical accessibility crosses zero. The resulting labels are **observation-robust core**, **finite-domain robust / partially identified**, and **observation-fragile**, respectively.

Grid fractions are assumption-domain coverage, not probabilities. P3 does not estimate true species-list completeness or arbitrary-MNAR latent truth.

"""
    marker = "## 11. Locked figures\n"
    if marker not in text:
        raise ValueError("pipeline locked-figures marker missing")
    text = text.replace(marker, p3 + "## 12. Locked figures\n", 1)
    text = replace_once(text, "## 12. Canonical paper surface", "## 13. Canonical paper surface")
    text = replace_once(text, "## 13. Chapter 1 / Chapter 2 handoff", "## 14. Chapter 1 / Chapter 2 handoff")
    text = replace_once(text, "## 14. Legacy v1 boundary", "## 15. Legacy v1 boundary")
    text = replace_once(
        text,
        "- Figure 4: `config/chapter1_v8_figure4_result_lock.json`",
        "- **Figure 4: `config/chapter1_v11_figure4_result_lock.json`**",
    )
    old_list = """1. `docs/chapter1_submission_freeze_20260915_p1_p2_defended.md`
2. `docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md`
3. `docs/chapter1_v10_submission_figure_sync_20260915.md`
4. `docs/chapter1_p2_component_nonconcordance_result_20260915.md`
5. `docs/chapter1_p1_final_decision_20260915.md`
6. `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`
7. `docs/chapter1_literature_positioning_20260909.md`"""
    new_list = """1. `docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md`
2. `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`
3. `docs/chapter1_v11_submission_figure_sync_20260915.md`
4. `docs/chapter1_p3_joint_observation_bias_checkpoint_20260915.md`
5. `docs/chapter1_p2_component_nonconcordance_result_20260915.md`
6. `docs/chapter1_p1_final_decision_20260915.md`
7. `docs/chapter1_court_evidence_and_theory_synthesis_20260914.md`
8. `docs/chapter1_literature_positioning_20260909.md`"""
    text = replace_once(text, old_list, new_list)
    text = replace_once(
        text,
        "Previous v8/v9 surfaces remain historical provenance; they do not define the current P1/P2-defended claim ceiling.",
        "Previous v8/v9/v10 surfaces remain historical provenance; they do not define the current P1/P2/P3-defended claim ceiling.",
    )
    return text


def write_sync_and_freeze(p3: dict, fig4: dict) -> None:
    figure_sync = f"""# Chapter 1 v11 submission figure sync — 2026-09-15

Canonical manuscript: `{MANUSCRIPT}`

Canonical P3 result lock: `{P3_LOCK}`

## Figure 4 — P3 joint observation-bias boundaries

Canonical lock: `{FIG4_LOCK}`

- run `{fig4['workflow_run_id']}`;
- artifact `{fig4['artifact_id']}`;
- digest `{fig4['artifact_digest']}`;
- visual review: passed.

Final Figure 4 role:

- **A:** formal direct-only native-nonendemic North–Tropical vector joint-bias surface (`1541/1575` robust cells);
- **B1:** Palearctic accessibility (`1575/1575` robust cells);
- **B2:** tropical accessibility (`1269/1575` robust cells);
- **C:** deterministic partial-identification envelopes separating sign identification from support identification;
- **D:** retained geometry/H5 claim boundaries (`0/12` nonlinear promotions; H5c `p=0.412`; H5d `0/8` qualified).

Grid fractions are predeclared sensitivity-domain coverage, not probabilities.

Figures 1–3 remain the P1/P2-defended canonical figures. Figure 4 is replaced by the P3-defended rendering; no P1/P2 estimate is refit.

## Final Figure 4 legend

**Figure 4 | Joint observation bias separates robust pattern from fragile identification.** The prospectively frozen P3 surface combines V5 trait-resolution MNAR with V6 distance-dependent list incompleteness and state-dependent species recording without increasing precision for hypothetical species. Robust/fragile regions are shown over the fixed assumption domain; grid fractions are not probabilities (A,B). In the primary direct-only native-nonendemic profile, the formal North–Tropical vector difference survives `1541/1575` joint cells, Palearctic accessibility survives `1575/1575`, and tropical accessibility survives `1269/1575`. Deterministic partial-identification bounds preserve the positive Palearctic accessibility sign but do not preserve formal support for every direct-only corner; tropical accessibility crosses zero (C). Existing response-geometry and H5 boundaries remain closed: `0/12` nonlinear shapes promoted, H5c `p=0.412`, and `0/8` distributed-threshold cells qualified (D). P3 localizes observation robustness; it does not estimate true completeness or promote a pollination mechanism.
"""
    FIG_SYNC.write_text(figure_sync, encoding="utf-8")

    freeze = f"""# Chapter 1 submission freeze — P1/P2/P3 defended — 2026-09-15

Canonical manuscript:

- `{MANUSCRIPT}`

Canonical P3 result:

- `{P3_LOCK}`
- run `{p3['workflow_run_id']}` / artifact `{p3['artifact_id']}`;
- digest `{p3['artifact_digest']}`.

Canonical Figure 4:

- `{FIG4_LOCK}`
- run `{fig4['workflow_run_id']}` / artifact `{fig4['artifact_id']}`;
- digest `{fig4['artifact_digest']}`;
- visual review passed.

## Frozen manuscript-level inference

1. **H1:** no single universal floral/reproductive island syndrome is recovered.
2. **H2/P2:** the formal same-layer North–Tropical joint response difference survives common-island and common-species safeguards; strong vector non-collinearity is not precisely established.
3. **H3/P1:** the strongest Palearctic response is strongly structured by true genus membership beyond matched arbitrary grouping complexity, but the exact additional family-to-genus attenuation remains spatially imprecise.
4. **P3 observation boundary:** Palearctic accessibility is the observation-robust core; the formal North–Tropical vector is highly finite-domain robust but only partially identified under deterministic joint bounds; tropical accessibility is observation-fragile.
5. **H5/geometry:** no global pollinator-specific mechanism, common nonlinear assemblage threshold, or distributed-threshold generator is promoted.

## Claim ceiling

The manuscript may say that the primary direct-only native-nonendemic North–Tropical vector difference survives `1541/1575` predeclared finite joint observation-bias cells, and that Palearctic accessibility survives `1575/1575` in both evidence scopes. It must simultaneously state that deterministic partial-identification bounds do not retain formal North–Tropical support across every corner, that direct-only Palearctic support is not identified across every extreme corner despite the sign remaining positive, and that tropical accessibility crosses zero and is observation-fragile.

P3 does not estimate the probability of any missingness mechanism, true flora completeness, arbitrary-MNAR latent truth, pollinator loss, or a causal assembly process. Grid fractions are not probabilities.

No P1/P2 estimate, H1–H5 threshold, primary data, failed mechanism gate, or historical result lock is reopened by this freeze.
"""
    FREEZE.write_text(freeze, encoding="utf-8")


def main() -> None:
    p3, fig4 = require_inputs()
    README.write_text(update_readme(README.read_text(encoding="utf-8")), encoding="utf-8")
    PIPELINE.write_text(
        update_pipeline(PIPELINE.read_text(encoding="utf-8")), encoding="utf-8"
    )
    write_sync_and_freeze(p3, fig4)

    for path in (README, PIPELINE, FIG_SYNC, FREEZE):
        text = path.read_text(encoding="utf-8")
        for required in (
            "chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md",
            "chapter1_p3_joint_observation_bias_result_lock.json",
            "chapter1_v11_figure4_result_lock.json",
        ):
            if required not in text:
                raise ValueError(f"{path} missing canonical v11 pointer: {required}")
        if "1541/1575" not in text.replace(",", ""):
            raise ValueError(f"{path} missing primary P3 finite-grid result")


if __name__ == "__main__":
    main()
