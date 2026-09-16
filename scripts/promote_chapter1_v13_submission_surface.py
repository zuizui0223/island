from __future__ import annotations

import json
from pathlib import Path

README = Path("README.md")
PIPELINE = Path("docs/PAPER_PIPELINE.md")
MANUSCRIPT = Path("docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md")
FIGURE_SYNC = Path("docs/chapter1_v13_submission_figure_sync_20260917.md")
LOCK = Path("config/chapter1_v13_unified_island_syndrome_result_lock.json")
V11_MANUSCRIPT = Path("docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md")


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"expected one replacement, found {count}: {old[:100]!r}")
    return text.replace(old, new, 1)


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError(f"section markers not unique: {start!r}, {end!r}")
    prefix, tail = text.split(start, 1)
    _, suffix = tail.split(end, 1)
    return prefix + start + replacement.rstrip() + "\n\n" + end + suffix


def require_inputs() -> dict:
    for path in (README, PIPELINE, MANUSCRIPT, FIGURE_SYNC, LOCK, V11_MANUSCRIPT):
        if not path.is_file():
            raise FileNotFoundError(path)
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if lock.get("contract") != "chapter1_v13_unified_island_syndrome_result_lock_v1":
        raise ValueError("unexpected v13 paper lock")
    if lock.get("paper_status") != "submission_surface_candidate":
        raise ValueError("v13 paper lock is not a submission-surface candidate")
    if lock.get("H4_functional_bridge", {}).get("inferential_role") != "posthoc_functional_triangulation":
        raise ValueError("v13 functional bridge must remain post-hoc")
    return lock


def update_readme(text: str) -> str:
    text = replace_once(
        text,
        "  -> canonical island-first v11 manuscript",
        "  -> historical v11 defended manuscript\n  -> v13 unified global island-syndrome manuscript",
    )
    section = """

The current publication-facing surface is **v13**. It is a new synthesis layer and does not rewrite the frozen v11/v12 provenance.

Read these first, in this order:

1. [`docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`](docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md) — current v13 full manuscript;
2. [`config/chapter1_v13_unified_island_syndrome_result_lock.json`](config/chapter1_v13_unified_island_syndrome_result_lock.json) — paper-level H1–H5 result and claim-ceiling lock;
3. [`docs/chapter1_v13_submission_figure_sync_20260917.md`](docs/chapter1_v13_submission_figure_sync_20260917.md) — v13 Figure 1–4 evidence-role contract;
4. [`config/chapter1_v13_functional_bridge_result_lock.json`](config/chapter1_v13_functional_bridge_result_lock.json) — post-hoc exact-species GloPL functional triangulation;
5. [`docs/chapter1_unified_hypothesis_20260917.md`](docs/chapter1_unified_hypothesis_20260917.md) — unified biological hypothesis and evidence map.

Historical defended surface retained unchanged:

- [`docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md`](docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md) — historical v11 island-first P1/P2/P3-defended manuscript;
- `config/chapter1_v12_two_panel_result_lock.json` and `config/chapter1_v12_h5_glopl_extension_result_lock.json` — historical v12 integration locks.

### Current v13 paper-level verdict

- **H1 — recurrent global island syndrome:** all four broad contemporary-flora contexts show supported multivariate isolation responses and a positive descriptive classic-island direction; exact regional response vectors are not assumed equal;
- **H2 — global pollination constraint:** full-global GloPL shows increasing experimental pollen limitation with geographic separation (`beta=+0.07937`, one-sided `p=0.01772`);
- **H3 — dual plant pathways:** reproductive assurance and generalized/accessibility architecture are treated as partially separable response families rather than one obligatory selfing-to-simplification sequence;
- **H4 — functional bridge:** exact-species GloPL triangulation is explicitly **post-hoc**; autonomous selfing is the strongest reproducible link to lower current pollen limitation, while architecture evidence is weaker and heterogeneous;
- **H5 — taxonomic realization:** the broad all-observed response survives source-free genus residualization, whereas the defended native Palearctic response is strongly genus-structured (`4/4 -> 4/4 -> 0/4`; matched pseudo-genus `p=0.02899`).

Publication-facing concept:

> **Geographic isolation is associated with a recurrent floral/reproductive island-syndrome direction and with stronger experimental pollen limitation. The syndrome is expressed through partially separable reproductive-assurance and floral-accessibility pathways, while taxonomic assembly determines how that shared pattern is realized across flora layers.**

The v13 synthesis is triangulation, not historical causal mediation. It does not claim that pollen limitation evolved the observed traits, that pollinator abundance globally declines with isolation, or that GloBI identifies the causal mechanism.
"""
    return replace_section(
        text,
        "## 6. Current submission surface\n",
        "## 7. Chapter 1 / Chapter 2 division of labour\n",
        section,
    )


def update_pipeline(text: str) -> str:
    text = replace_once(
        text,
        "current **island-first, P1/P2/P3-defended v11 paper**",
        "current **v13 unified global island-syndrome paper**, with v11/v12 retained as historical defended provenance",
    )
    text = replace_once(
        text,
        "[13] canonical island-first v11 manuscript",
        "[13] historical island-first v11 defended manuscript\n                    |\n                    v\n[14] v13 global syndrome + GloPL + functional bridge synthesis",
    )
    section = """

The current publication-facing surface is **v13**. Read in this order:

1. `docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md`
2. `config/chapter1_v13_unified_island_syndrome_result_lock.json`
3. `docs/chapter1_v13_submission_figure_sync_20260917.md`
4. `config/chapter1_v13_functional_bridge_result_lock.json`
5. `docs/chapter1_unified_hypothesis_20260917.md`
6. `docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md` — historical defended manuscript

The v13 paper does **not** retroactively redefine the original progressive H1–H5 contract. Instead it synthesizes already-frozen and newly versioned evidence into a publication-facing hierarchy:

- H1: recurrent global classic-island direction with regional modification;
- H2: full-global experimental pollen-limitation gradient;
- H3: partially separable reproductive-assurance and floral-accessibility pathways;
- H4: explicitly post-hoc exact-species functional triangulation;
- H5: flora-layer-specific taxonomic realization.

Canonical new v13 functional-bridge run:

- run `35141624253`;
- artifact `10465048981`;
- digest `sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e`.

The frozen Route A/B distance-by-trait moderation failures remain negative results. v13 does not relabel them as support; it asks a distinct post-hoc functional question about current pollen-limitation levels under the already-frozen exact-species trait states.

Previous v8/v9/v10/v11/v12 surfaces remain historical provenance and continue to define the chronology of claim defense.
"""
    return replace_section(
        text,
        "## 13. Canonical paper surface\n",
        "## 14. Chapter 1 / Chapter 2 handoff\n",
        section,
    )


def update_manuscript(text: str) -> str:
    return replace_once(
        text,
        "This does not establish that historical pollen limitation caused the observed trait evolution.",
        "This does not identify historical pollen limitation as the cause of the observed trait evolution.",
    )


def main() -> None:
    require_inputs()
    README.write_text(update_readme(README.read_text(encoding="utf-8")), encoding="utf-8")
    PIPELINE.write_text(update_pipeline(PIPELINE.read_text(encoding="utf-8")), encoding="utf-8")
    MANUSCRIPT.write_text(update_manuscript(MANUSCRIPT.read_text(encoding="utf-8")), encoding="utf-8")


if __name__ == "__main__":
    main()
