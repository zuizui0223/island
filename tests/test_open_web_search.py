from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import yaml

from island_v2.open_web_search import (
    FixtureSearchBackend,
    SearchHit,
    SearchTask,
    build_query_tasks,
    discover,
    load_strict_name_variants,
    load_strict_synonyms,
    write_csv,
)

ROOT = Path(__file__).resolve().parents[1]


def _config() -> dict:
    return yaml.safe_load(
        (ROOT / "config" / "open_web_trait_acquisition.yml").read_text(encoding="utf-8")
    )


def test_multilingual_query_generation_keeps_trait_and_name_kind() -> None:
    queue = [
        {
            "accepted_species": "Plantus alba",
            "trait_name": "self_incompatibility",
            "priority_score": 12,
            "potential_species_trait_unlocked": 8,
            "expected_information_gain": 3.5,
        }
    ]
    tasks = build_query_tasks(
        queue,
        config=_config(),
        synonyms={"Plantus alba": ["Plantus candidus"]},
        vernacular_names={"Plantus alba": [("白花草", "ja")]},
        include_pdf_query=False,
    )
    assert {task.name_kind for task in tasks} == {"accepted", "synonym", "vernacular"}
    assert {"en", "la", "ja"}.issubset({task.language for task in tasks})
    assert all(task.trait_name == "self_incompatibility" for task in tasks)
    assert any("自家不和合性" in task.query for task in tasks)
    assert not any("自動自殖" in task.query for task in tasks)
    vernacular = [task for task in tasks if task.name_kind == "vernacular"]
    assert {task.language for task in vernacular} == {"ja"}
    assert all("self-compatible" not in task.query for task in vernacular)


def test_self_fertile_query_is_never_relabelled_as_autonomous_selfing() -> None:
    config = _config()
    si = config["languages"]["en"]["query_labels"]["self_incompatibility"]
    autonomous = config["languages"]["en"]["query_labels"]["autonomous_selfing_capacity"]
    assert "self-fertile" in si
    assert "self-fertile" not in autonomous
    assert config["fail_closed"]["self_fertile_to_autonomous_selfing"] is False


def test_fixture_backend_is_discovery_only_and_reproducible() -> None:
    fixture = json.loads(
        (ROOT / "tests" / "fixtures" / "open_web" / "search_fixture.json").read_text(
            encoding="utf-8"
        )
    )
    task = SearchTask(
        task_id="task-1",
        accepted_species="Plantus alba",
        trait_name="flower_primary_color",
        searched_name="Plantus alba",
        name_kind="accepted",
        language="en",
        query='"Plantus alba" ("flower colour" OR "flower color" OR "corolla colour")',
        priority_score=1,
        potential_species_trait_unlocked=1,
        expected_information_gain=1,
    )
    hits, errors, report = discover(
        [task],
        backends=[FixtureSearchBackend(fixture)],
        max_queries=1,
        results_per_query=10,
    )
    assert not errors
    assert len(hits) == 1
    assert hits[0].result_url == "https://example.edu/flora/plantus-alba"
    assert hits[0].result_domain == "example.edu"
    assert hits[0].snippet_discovery_only == "The flowers are white."
    assert report["snippets_are_evidence"] is False
    stable = asdict(hits[0])
    assert stable["discovery_sha256"]


def test_query_budget_stops_before_second_task() -> None:
    task = SearchTask(
        task_id="one",
        accepted_species="Plantus alba",
        trait_name="flower_primary_color",
        searched_name="Plantus alba",
        name_kind="accepted",
        language="en",
        query="one",
        priority_score=1,
        potential_species_trait_unlocked=1,
        expected_information_gain=1,
    )
    backend = FixtureSearchBackend({"one": [], "two": []})
    second = SearchTask(**{**asdict(task), "task_id": "two", "query": "two"})
    _, _, report = discover([task, second], backends=[backend], max_queries=1)
    assert report["queries_used"] == 1
    assert report["attempted_task_ids"] == ["one"]


def test_empty_discovery_csv_retains_schema_for_resumable_workflow(tmp_path: Path) -> None:
    target = tmp_path / "search_hits.csv"
    write_csv(
        target,
        [],
        fieldnames=list(SearchHit.__dataclass_fields__),
    )
    assert target.read_text(encoding="utf-8").startswith(
        "task_id,accepted_species,trait_name,searched_name"
    )


def test_query_synonyms_require_two_backbone_snapshots(tmp_path: Path) -> None:
    target = tmp_path / "synonyms.csv"
    target.write_text(
        "accepted_species,synonym,backbone\n"
        "Plantus alba,Planta candida,GBIF\n"
        "Plantus alba,Planta candida,WFO\n"
        "Plantus alba,Planta dubia,GBIF\n",
        encoding="utf-8",
    )
    assert load_strict_synonyms(target) == {"Plantus alba": ["Planta candida"]}


def test_basionym_requires_two_backbones_and_keeps_name_kind(
    tmp_path: Path,
) -> None:
    target = tmp_path / "names.csv"
    target.write_text(
        "accepted_species,basionym,backbone\n"
        "Plantus alba,Basionyma alba,GBIF\n"
        "Plantus alba,Basionyma alba,WFO\n",
        encoding="utf-8",
    )
    variants = load_strict_name_variants(target)
    assert variants == {"Plantus alba": [("Basionyma alba", "basionym")]}
    tasks = build_query_tasks(
        [
            {
                "accepted_species": "Plantus alba",
                "trait_name": "floral_form",
            }
        ],
        config=_config(),
        synonyms=variants,
        include_pdf_query=False,
    )
    assert any(task.name_kind == "basionym" for task in tasks)
