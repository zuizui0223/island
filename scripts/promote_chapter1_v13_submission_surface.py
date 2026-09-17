from __future__ import annotations

import json
from pathlib import Path

README = Path("README.md")
PIPELINE = Path("docs/PAPER_PIPELINE.md")
MANUSCRIPT = Path("docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md")
FIGURE_SYNC = Path("docs/chapter1_v13_submission_figure_sync_20260917.md")
HYPOTHESIS = Path("docs/chapter1_unified_hypothesis_20260917.md")
FREEZE = Path("docs/chapter1_submission_freeze_v13_20260917.md")
LOCK = Path("config/chapter1_v13_unified_island_syndrome_result_lock.json")
V11_MANUSCRIPT = Path("docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md")


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError(f"section markers not unique: {start!r}, {end!r}")
    prefix, tail = text.split(start, 1)
    _, suffix = tail.split(end, 1)
    return prefix + start + "\n" + replacement.strip() + "\n\n" + end + suffix


def require_inputs() -> dict:
    for path in (
        README,
        PIPELINE,
        MANUSCRIPT,
        FIGURE_SYNC,
        HYPOTHESIS,
        FREEZE,
        LOCK,
        V11_MANUSCRIPT,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if lock.get("contract") != "chapter1_v13_unified_island_syndrome_result_lock_v2":
        raise ValueError("unexpected v13 paper lock")
    if lock.get("status") != "submission_candidate_global_only_frozen":
        raise ValueError("v13 global-only paper lock is not frozen")
    if set(lock.get("architecture", {})) != {"H1", "H2", "H3", "H4"}:
        raise ValueError("v13 publication surface must remain H1-H4")
    if lock.get("H4_functional_bridge", {}).get("inferential_role") != "posthoc_functional_triangulation":
        raise ValueError("v13 functional bridge must remain post-hoc")
    return lock


def update_readme(text: str) -> str:
    section = """
The current publication-facing surface is **v13 global-only**, organized as **H1–H4**. This revision keeps the merged v13 evidence but removes assembly-depth and between-stratum branches from the submission spine. Frozen v11/v12 files remain historical provenance.

Read these first, in this order:

1. [`docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`](docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md) — current global-only v13 manuscript;
2. [`config/chapter1_v13_unified_island_syndrome_result_lock.json`](config/chapter1_v13_unified_island_syndrome_result_lock.json) — H1–H4 result and claim-ceiling lock;
3. [`docs/chapter1_v13_submission_figure_sync_20260917.md`](docs/chapter1_v13_submission_figure_sync_20260917.md) — three-main-figure evidence-role contract;
4. [`config/chapter1_v13_functional_bridge_result_lock.json`](config/chapter1_v13_functional_bridge_result_lock.json) — post-hoc exact-species GloPL functional triangulation;
5. [`docs/chapter1_unified_hypothesis_20260917.md`](docs/chapter1_unified_hypothesis_20260917.md) — global-only biological hypothesis and evidence map.

Historical defended surface retained unchanged:

- [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md) — historical defended manuscript;
- `config/chapter1_v12_two_panel_result_lock.json` and `config/chapter1_v12_h5_glopl_extension_result_lock.json` — historical integration locks.

### Current v13 paper-level verdict

- **H1 — recurrent global island syndrome:** all four geographic replication strata show supported multivariate isolation responses and a positive descriptive classic-island direction in both evidence scopes; no between-stratum comparison is promoted in the current submission;
- **H2 — global pollination constraint:** full-global GloPL shows increasing experimental pollen limitation with geographic separation (`beta=+0.07937`, one-sided `p=0.01772`);
- **H3 — dual plant pathways:** reproductive assurance and generalized/accessibility architecture are partially separable response families rather than one obligatory selfing-to-simplification sequence;
- **H4 — functional bridge:** exact-species GloPL triangulation is explicitly **post-hoc**; autonomous selfing is the strongest reproducible link to lower current pollen limitation, while architecture evidence is weaker.

Publication-facing concept:

> **Geographic isolation is associated globally with stronger experimental pollen limitation and with a recurrent floral/reproductive island-syndrome direction expressed through partially separable reproductive-assurance and floral-accessibility pathways.**

The v13 synthesis is triangulation, not historical causal mediation. It does not claim that pollen limitation caused the observed trait evolution, that pollinator abundance globally declines with isolation, or that GloBI identifies the causal mechanism.
"""
    return replace_section(
        text,
        "## 6. Current submission surface\n",
        "## 7. Chapter 1 / Chapter 2 division of labour\n",
        section,
    )


def update_pipeline(text: str) -> str:
    section = """
The current publication-facing surface is **v13 global-only**, with an **H1–H4** scientific spine. Historical v11/v12 surfaces remain immutable provenance.

Read in this order:

1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
3. `docs/chapter1_v13_submission_figure_sync_20260917.md`
4. `config/chapter1_v13_functional_bridge_result_lock.json`
5. `docs/chapter1_unified_hypothesis_20260917.md`
6. `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md` — historical defended manuscript

The current publication-facing hierarchy is:

- H1: recurrent global classic-island direction across four geographic replication strata;
- H2: full-global experimental pollen-limitation gradient;
- H3: partially separable reproductive-assurance and floral-accessibility pathways;
- H4: explicitly post-hoc exact-species functional triangulation.

The four geographic strata are used to demonstrate recurrence of the global direction. Between-stratum differences are not part of the current submission spine.

Canonical v13 functional-bridge run:

- run `35141624253`;
- artifact `10465048981`;
- digest `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`.

The frozen Route A/B distance-by-trait moderation failures remain negative results. v13 does not relabel them as support. GloBI remains supplementary sampling-sensitive evidence.

Previous v8/v9/v10/v11/v12 surfaces remain historical provenance and continue to document the chronology of claim defense.
"""
    return replace_section(
        text,
        "## 13. Canonical paper surface\n",
        "## 14. Chapter 1 / Chapter 2 handoff\n",
        section,
    )


def main() -> None:
    require_inputs()
    README.write_text(update_readme(README.read_text(encoding="utf-8")), encoding="utf-8")
    PIPELINE.write_text(update_pipeline(PIPELINE.read_text(encoding="utf-8")), encoding="utf-8")


if __name__ == "__main__":
    main()
