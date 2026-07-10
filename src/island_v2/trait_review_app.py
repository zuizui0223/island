"""Streamlit helper for one-row-at-a-time trait candidate review.

The app is only an adjudication interface. It writes review decisions back to
the staging queue and reuses the queue validator before saving, so it cannot
promote trait values by itself.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from island_v2 import trait_candidate_review_queue as queue_tools

DEFAULT_QUEUE_CSV = Path(
    "data/v2/staging/traits/validation_pilot/candidate_review_queue/"
    "trait_candidate_review_queue.csv"
)
DEFAULT_ONTOLOGY = Path("config/trait_ontology.yml")
DEFAULT_GROUPS = [
    "P1_species_direct_reported",
    "P1_reported_ecology_llm_source_check",
    "P1_visual_trait_candidate",
    "P2_species_indirect_source_check",
    "P3_proxy_candidate",
    "P4_descriptive_scout_source_check",
]

DECISIONS = ["accepted", "rejected", "needs_source_check"]

COLOUR_VALUE_HINTS = {
    "white": ["white"],
    "cream": ["white", "other_described"],
    "yellow": ["yellow_orange"],
    "orange": ["yellow_orange"],
    "red": ["red_pink"],
    "pink": ["red_pink"],
    "purple": ["blue_purple"],
    "blue": ["blue_purple"],
    "green": ["green_brown_inconspicuous"],
    "brown": ["green_brown_inconspicuous"],
    "black": ["green_brown_inconspicuous"],
}


def load_ontology(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        return {}
    return data


def read_queue(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("")


def pending_indices(table: pd.DataFrame, groups: list[str]) -> list[int]:
    if table.empty:
        return []
    subset = table.loc[
        table["adjudication_decision"].astype(str).str.strip().eq("")
        & table["review_group"].astype(str).isin(groups)
    ].copy()
    if subset.empty:
        return []
    subset["_priority_sort"] = subset["review_priority"].astype(int)
    subset = subset.sort_values(
        ["_priority_sort", "accepted_species", "trait_name", "candidate_value"],
        kind="stable",
    )
    return [int(index) for index in subset.index]


def next_pending_index(table: pd.DataFrame, groups: list[str]) -> int | None:
    indices = pending_indices(table, groups)
    return indices[0] if indices else None


def final_value_options(row: pd.Series, ontology: dict[str, Any]) -> list[str]:
    trait_name = str(row.get("trait_name", "")).strip()
    allowed = queue_tools._allowed_final_values(trait_name, ontology)  # noqa: SLF001
    if allowed is None:
        value = str(row.get("candidate_value", "")).strip()
        return [value] if value else [""]

    allowed_list = sorted(allowed)
    candidate = str(row.get("candidate_value", "")).strip()
    preferred: list[str] = []
    if candidate in allowed:
        preferred.append(candidate)
    if trait_name == "flower_primary_color":
        preferred.extend(value for value in COLOUR_VALUE_HINTS.get(candidate, []) if value in allowed)

    seen: set[str] = set()
    ordered = [value for value in preferred + allowed_list if not (value in seen or seen.add(value))]
    return ordered


def default_reason(row: pd.Series, decision: str) -> str:
    if decision == "accepted":
        if str(row.get("source_lane", "")) == "visual_trait_candidate":
            return "Image visibly supports the floral trait candidate; no reproductive or pollinator inference made."
        if str(row.get("source_lane", "")) == "llm_reported_ecology":
            return "Source excerpt explicitly supports the reported reproductive, selfing, pollen-vector, or visitor evidence."
        return "Source explicitly describes the floral trait for the accepted species."
    if decision == "rejected":
        return "Source excerpt does not explicitly support this trait value for the accepted species."
    return "Needs source-page inspection before accepting or rejecting."


def evidence_image_url(row: pd.Series) -> str:
    source_lane = str(row.get("source_lane", "")).strip()
    source_type = str(row.get("source_type", "")).strip().lower()
    source_url = str(row.get("source_url", "")).strip()
    if source_lane == "visual_trait_candidate" and source_url:
        return source_url
    if source_type in {"image", "specimen_image", "gbif_media"} and source_url:
        return source_url
    return ""


def apply_decision(
    table: pd.DataFrame,
    row_index: int,
    *,
    decision: str,
    final_value: str,
    reason: str,
    reviewer: str,
    review_date: str,
) -> pd.DataFrame:
    updated = table.copy()
    updated.loc[row_index, "adjudication_decision"] = decision
    updated.loc[row_index, "final_value"] = final_value if decision == "accepted" else ""
    updated.loc[row_index, "decision_reason"] = reason
    updated.loc[row_index, "reviewer"] = reviewer
    updated.loc[row_index, "review_date"] = review_date
    return updated


def _render_app() -> None:
    try:
        import streamlit as st
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Streamlit is not installed. Install it with: pip install -e ."
        ) from exc

    st.set_page_config(page_title="Trait Candidate Review", layout="wide")
    st.title("Trait Candidate Review")

    queue_path = Path(
        st.sidebar.text_input("Queue CSV", value=str(DEFAULT_QUEUE_CSV).replace("\\", "/"))
    )
    ontology_path = Path(
        st.sidebar.text_input("Trait ontology", value=str(DEFAULT_ONTOLOGY).replace("\\", "/"))
    )
    reviewer = st.sidebar.text_input("Reviewer", value="")
    selected_groups = st.sidebar.multiselect(
        "Review groups",
        DEFAULT_GROUPS,
        default=[group for group in DEFAULT_GROUPS if group.startswith("P1")],
    )

    if not queue_path.exists():
        st.error(f"Queue CSV not found: {queue_path}")
        return
    if not ontology_path.exists():
        st.error(f"Ontology not found: {ontology_path}")
        return

    table = read_queue(queue_path)
    ontology = load_ontology(ontology_path)
    summary = queue_tools.review_queue_summary(table)
    c1, c2, c3 = st.columns(3)
    c1.metric("Rows", summary["n_review_candidates"])
    c2.metric("Decisions", summary["n_review_decisions_recorded"])
    c3.metric("Accepted ready", summary["n_accepted_ready_for_curation"])

    pending = pending_indices(table, selected_groups)
    if not pending:
        st.success("No pending rows in the selected group(s).")
        return

    row_index = st.selectbox(
        "Pending row",
        pending,
        format_func=lambda idx: (
            f"{table.loc[idx, 'review_priority']} | {table.loc[idx, 'accepted_species']} | "
            f"{table.loc[idx, 'trait_name']}={table.loc[idx, 'candidate_value']}"
        ),
    )
    row = table.loc[row_index]

    left, right = st.columns([2, 1])
    with left:
        st.subheader(str(row["accepted_species"]))
        st.write(f"Trait: `{row['trait_name']}`")
        st.write(f"Candidate: `{row['candidate_value']}`")
        st.write(f"Evidence scope: `{row['evidence_scope']}`")
        if str(row.get("source_url", "")).strip():
            st.markdown(f"[Open source]({row['source_url']})")
        image_url = evidence_image_url(row)
        if image_url:
            st.image(image_url, caption="Image evidence", use_container_width=True)
        st.text_area("Source excerpt", value=str(row.get("raw_description", "")), height=180)
        st.text_area("Review boundary", value=str(row.get("review_action", "")), height=90)

    with right:
        decision = st.radio("Decision", DECISIONS, horizontal=False)
        options = final_value_options(row, ontology)
        final_value = ""
        if decision == "accepted":
            final_value = st.selectbox("Final value", options)
        reason = st.text_area("Decision reason", value=default_reason(row, decision), height=110)
        review_date = st.date_input("Review date", value=date.today()).isoformat()

        if st.button("Save decision and move on", type="primary"):
            updated = apply_decision(
                table,
                int(row_index),
                decision=decision,
                final_value=final_value,
                reason=reason,
                reviewer=reviewer,
                review_date=review_date,
            )
            errors = queue_tools.validate_review_queue(updated, ontology=ontology)
            if errors:
                st.error("Not saved. Fix validation errors first.")
                for error in errors[:20]:
                    st.write(f"- {error}")
                return
            updated.to_csv(queue_path, index=False)
            queue_tools.write_review_outputs(updated, queue_path.parent, ontology=ontology)
            st.success("Saved and validated.")
            st.rerun()

    st.divider()
    st.dataframe(
        table.loc[pending, ["review_priority", "review_group", "accepted_species", "trait_name", "candidate_value"]],
        use_container_width=True,
    )


def main() -> None:
    _render_app()


if __name__ == "__main__":
    main()
