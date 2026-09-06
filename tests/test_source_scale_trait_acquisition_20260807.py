from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FETCH = load_script(
    "fetch_europe_pmc_reproduction_scale_20260807",
    ROOT / "analysis/fetch_europe_pmc_reproduction_scale_20260807.py",
)
INTEGRATE = load_script(
    "integrate_europe_pmc_reproduction_scale_20260807",
    ROOT / "analysis/integrate_europe_pmc_reproduction_scale_20260807.py",
)
VALUE_PATTERNS = FETCH.VALUE_PATTERNS
BATCH2_CORRECTED_PROMOTIONS = INTEGRATE.BATCH2_CORRECTED_PROMOTIONS
BATCH2_CORROBORATING_NOT_COUNTED = INTEGRATE.BATCH2_CORROBORATING_NOT_COUNTED
BATCH2_INVALIDATED_PRIOR = INTEGRATE.BATCH2_INVALIDATED_PRIOR
BATCH2_PROMOTIONS = INTEGRATE.BATCH2_PROMOTIONS
BATCH2_REJECTIONS = INTEGRATE.BATCH2_REJECTIONS
INVALIDATED_PRIOR = INTEGRATE.INVALIDATED_PRIOR
PROMOTIONS = INTEGRATE.PROMOTIONS
common_evidence = INTEGRATE.common_evidence
review_candidates = INTEGRATE.review_candidates

CANDIDATES = (
    ROOT
    / "data/v2/staging/traits/direct_llm_pilot"
    / "20260807_europe_pmc_reproduction_scale"
    / "fulltext_reproductive_candidates.csv.gz"
)


def test_every_fulltext_candidate_has_an_explicit_review_decision() -> None:
    candidates = pd.read_csv(CANDIDATES, dtype=str).fillna("")
    review = review_candidates(candidates)

    assert len(review) == 43
    assert review["candidate_id"].nunique() == 43
    assert int(review["accepted_correct"].sum()) == 17
    assert int(review["promotion_accepted"].sum()) == 8
    assert review["review_decision"].value_counts().to_dict() == {
        "rejected": 26,
        "corroborating_not_counted": 9,
        "accepted": 8,
    }


def test_promoted_rows_are_unique_direct_species_traits() -> None:
    candidates = pd.read_csv(CANDIDATES, dtype=str).fillna("")
    common = common_evidence(review_candidates(candidates))

    assert len(common) == len(PROMOTIONS) == 8
    assert not common.duplicated(["accepted_species", "trait_name"]).any()
    assert set(common["quality"]) == {"high", "medium"}
    assert set(common["evidence_scope"]) == {"species_direct"}
    assert common["source_excerpt"].str.len().gt(0).all()
    assert common["source_lineage"].str.startswith(("doi:", "pmcid:")).all()


def test_known_panicum_comparison_misattribution_is_invalidated() -> None:
    assert (
        "Panicum hallii",
        "self_incompatibility",
        "SI",
        "url:https://europepmc.org/article/MED/28449656",
    ) in INVALIDATED_PRIOR


def test_self_fertile_is_not_an_autonomous_selfing_pattern() -> None:
    autonomous_patterns = VALUE_PATTERNS["autonomous_selfing_capacity"]
    assert not any(
        pattern.search("The plant is self-fertile.") for _, pattern in autonomous_patterns
    )


def test_batch2_review_is_complete_and_preserves_correction_audit() -> None:
    candidate_ids = sorted(
        set(BATCH2_PROMOTIONS)
        | set(BATCH2_CORRECTED_PROMOTIONS)
        | set(BATCH2_CORROBORATING_NOT_COUNTED)
        | set(BATCH2_REJECTIONS)
    )
    candidates = pd.DataFrame(
        {
            "candidate_id": candidate_ids,
            "accepted_species": [f"Species {index}" for index in range(len(candidate_ids))],
            "trait_name": "self_incompatibility",
            "axis": "reproductive_assurance",
            "raw_value": "SI",
            "normalized_value": "SI",
            "pmcid": [f"PMC{index}" for index in range(len(candidate_ids))],
        }
    )

    review = review_candidates(candidates, "20260808_batch2")

    assert len(review) == 53
    assert int(review["accepted_correct"].sum()) == 37
    assert int(review["promotion_accepted"].sum()) == 12
    assert review["review_decision"].value_counts().to_dict() == {
        "corroborating_not_counted": 27,
        "rejected": 14,
        "accepted": 10,
        "corrected_and_accepted": 2,
    }
    corrected = review.loc[review["review_decision"].eq("corrected_and_accepted")]
    assert set(corrected["extracted_trait_name"]) == {"self_incompatibility"}
    assert set(corrected["trait_name"]) == {"mating_system"}
    assert not corrected["accepted_correct"].any()
    assert len(BATCH2_INVALIDATED_PRIOR) == 5


def test_completed_query_logs_are_excluded_before_next_batch(tmp_path: Path) -> None:
    queue = pd.DataFrame(
        [
            {
                "genus": "Done",
                "trait_name": "self_incompatibility",
                "current_support": "2",
            },
            {
                "genus": "Next",
                "trait_name": "self_incompatibility",
                "current_support": "2",
            },
        ]
    )
    queue_path = tmp_path / "queue.csv"
    queue.to_csv(queue_path, index=False)
    completed_path = tmp_path / "completed.csv"
    pd.DataFrame([{"genus": "Done", "trait_name": "self_incompatibility"}]).to_csv(
        completed_path, index=False
    )

    tasks = FETCH.build_tasks(
        queue_path,
        top_n=10,
        completed_query_logs=[completed_path],
    )

    assert tasks[["genus", "trait_name"]].to_dict("records") == [
        {"genus": "Next", "trait_name": "self_incompatibility"}
    ]


def test_latest_direct_support_and_agreement_control_the_next_queue(
    tmp_path: Path,
) -> None:
    queue_path = tmp_path / "queue.csv"
    pd.DataFrame(
        [
            {
                "genus": "Agree",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "current_support": "1",
            },
            {
                "genus": "Conflict",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "current_support": "2",
            },
        ]
    ).to_csv(queue_path, index=False)
    direct_path = tmp_path / "direct.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Agree one",
                "genus": "Agree",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "state_set": '["SI"]',
            },
            {
                "accepted_species": "Agree two",
                "genus": "Agree",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "state_set": '["SI"]',
            },
            {
                "accepted_species": "Conflict one",
                "genus": "Conflict",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "state_set": '["SI"]',
            },
            {
                "accepted_species": "Conflict two",
                "genus": "Conflict",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "state_set": '["SC"]',
            },
        ]
    ).to_csv(direct_path, index=False)

    tasks = FETCH.build_tasks(
        queue_path,
        top_n=10,
        current_direct_path=direct_path,
        require_agreement=True,
    )

    assert tasks[["genus", "trait_name"]].to_dict("records") == [
        {"genus": "Agree", "trait_name": "self_incompatibility"}
    ]


def test_support_one_is_only_included_when_explicitly_requested(tmp_path: Path) -> None:
    queue_path = tmp_path / "queue.csv"
    pd.DataFrame(
        [
            {
                "genus": "One",
                "trait_name": "self_incompatibility",
                "current_support": "1",
            },
            {
                "genus": "Two",
                "trait_name": "self_incompatibility",
                "current_support": "2",
            },
        ]
    ).to_csv(queue_path, index=False)

    default = FETCH.build_tasks(queue_path, top_n=10)
    expanded = FETCH.build_tasks(queue_path, top_n=10, support_levels=(1, 2))

    assert default["genus"].tolist() == ["Two"]
    assert expanded["genus"].tolist() == ["One", "Two"]


def test_exact_species_tasks_are_distinct_from_completed_genus_searches(
    tmp_path: Path,
) -> None:
    queue_path = tmp_path / "queue.csv"
    pd.DataFrame(
        [
            {
                "genus": "Examplea",
                "trait_name": "self_incompatibility",
                "current_support": "2",
                "query_scope": "exact_species",
                "target_species": "Examplea gamma",
            }
        ]
    ).to_csv(queue_path, index=False)
    completed_path = tmp_path / "completed.csv"
    pd.DataFrame(
        [{"genus": "Examplea", "trait_name": "self_incompatibility"}]
    ).to_csv(completed_path, index=False)

    tasks = FETCH.build_tasks(
        queue_path,
        top_n=10,
        completed_query_logs=[completed_path],
    )

    assert len(tasks) == 1
    assert tasks.iloc[0]["query_scope"] == "exact_species"
    assert tasks.iloc[0]["target_species"] == "Examplea gamma"
    assert '("Examplea gamma")' in tasks.iloc[0]["query"]


def test_article_resolves_established_binomial_abbreviation() -> None:
    raw = b"""<article><body>
      <p>Examplea alpha is an island herb.</p>
      <p>E. alpha is self-compatible, while E. beta is self-incompatible.</p>
    </body></article>"""
    tasks = pd.DataFrame([{"genus": "Examplea", "trait_name": "self_incompatibility"}])

    candidates = FETCH.article_candidates(
        raw,
        pmcid="PMC1",
        title="Example",
        doi="10.0000/example",
        task_rows=tasks,
        canonical={
            "examplea alpha": "Examplea alpha",
            "examplea beta": "Examplea beta",
        },
        by_genus={"Examplea": ["Examplea alpha", "Examplea beta"]},
        article_license="CC BY",
        retrieved_at="2026-08-28T00:00:00Z",
    )

    assert [(row["accepted_species"], row["normalized_value"]) for row in candidates] == [
        ("Examplea alpha", "SC")
    ]
    assert candidates[0]["matched_name"] == "E. alpha"
    assert (
        candidates[0]["name_match_method"]
        == "accepted_binomial_established_then_unambiguous_abbreviation"
    )


def test_article_rejects_cross_sentence_species_trait_bleed() -> None:
    raw = b"""<article><body><p>
      Examplea alpha is an island herb. Another plant is self-compatible.
    </p></body></article>"""
    tasks = pd.DataFrame([{"genus": "Examplea", "trait_name": "self_incompatibility"}])

    candidates = FETCH.article_candidates(
        raw,
        pmcid="PMC2",
        title="Example",
        doi="10.0000/example2",
        task_rows=tasks,
        canonical={"examplea alpha": "Examplea alpha"},
        by_genus={"Examplea": ["Examplea alpha"]},
        article_license="CC BY",
        retrieved_at="2026-08-28T00:00:00Z",
    )

    assert candidates == []


def test_article_retains_nonmaster_congener_for_backbone_review() -> None:
    raw = b"""<article><body><p>
      Examplea gamma is self-incompatible in hand-pollination experiments.
    </p></body></article>"""
    tasks = pd.DataFrame([{"genus": "Examplea", "trait_name": "self_incompatibility"}])

    candidates = FETCH.article_candidates(
        raw,
        pmcid="PMC3",
        title="Example",
        doi="10.0000/example3",
        task_rows=tasks,
        canonical={"examplea alpha": "Examplea alpha"},
        by_genus={"Examplea": ["Examplea alpha"]},
        article_license="CC BY",
        retrieved_at="2026-08-28T00:00:00Z",
    )

    assert len(candidates) == 1
    assert candidates[0]["accepted_species"] == "Examplea gamma"
    assert candidates[0]["source_species_name"] == "Examplea gamma"
    assert candidates[0]["in_fixed_master"] == "false"
    assert candidates[0]["rule_support_only"] == "true"
    assert candidates[0]["name_match_method"] == ("source_binomial_exact_pending_two_backbone")


def test_exact_species_task_does_not_rediscover_an_existing_congener() -> None:
    raw = b"""<article><body><p>
      Examplea alpha is self-compatible, while Examplea gamma is self-incompatible.
    </p></body></article>"""
    tasks = pd.DataFrame(
        [
            {
                "genus": "Examplea",
                "trait_name": "self_incompatibility",
                "query_scope": "exact_species",
                "target_species": "Examplea gamma",
            }
        ]
    )

    candidates = FETCH.article_candidates(
        raw,
        pmcid="PMC4",
        title="Example",
        doi="10.0000/example4",
        task_rows=tasks,
        canonical={
            "examplea alpha": "Examplea alpha",
            "examplea gamma": "Examplea gamma",
        },
        by_genus={"Examplea": ["Examplea alpha", "Examplea gamma"]},
        article_license="CC BY",
        retrieved_at="2026-08-28T00:00:00Z",
    )

    assert [(row["accepted_species"], row["normalized_value"]) for row in candidates] == [
        ("Examplea gamma", "SI")
    ]
