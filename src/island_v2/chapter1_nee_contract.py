from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import typer
import yaml

app = typer.Typer(help="Validate the prospective Chapter 1 NEE double-filter contract.")

EXPECTED_CHANNELS = (
    "bombus",
    "non_bombus_bees",
    "lepidoptera",
    "flower_visiting_birds",
    "diptera",
)
EXPECTED_SOURCE_STATES = {"available", "structurally_absent", "unresolved"}
EXPECTED_ISLAND_STATES = {"retained", "disrupted", "unresolved", "structurally_absent"}
PROHIBITED_DEPENDENCY_FIELDS = {
    "flower_colour",
    "floral_form",
    "tube_depth",
    "floral_symmetry",
    "large_bee_like_score",
    "butterfly_like_score",
    "bird_like_score",
    "chapter1_focal_trait_response",
}
PROHIBITED_CHANNEL_INFERENCE = {
    "floral_colour_inference",
    "floral_form_inference",
    "pollination_syndrome_reverse_inference",
    "focal_trait_response_inference",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_true(mapping: dict[str, Any], key: str, section: str) -> None:
    _require(mapping.get(key) is True, f"{section}.{key} must remain true")


def validate_contract_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate the pre-outcome N0 lock and return a compact qualification receipt."""

    _require(
        payload.get("contract") == "chapter1_nee_double_filter_v1",
        "unexpected NEE contract id",
    )
    _require(
        payload.get("status") == "prospective_preoutcome_qualification",
        "contract must remain prospective_preoutcome_qualification",
    )

    parent = payload.get("parent_freeze", {})
    _require(
        parent.get("analysis_contract") == "chapter1_progressive_analysis_v1",
        "NEE challenge must remain attached to the frozen Chapter 1 parent contract",
    )
    _require(
        parent.get("frozen_results") == ["H1", "H2", "H3", "H4"],
        "H1-H4 frozen result set changed",
    )
    _require(
        parent.get("reopen_or_rescue_parent_results") is False,
        "NEE challenge may not reopen or rescue frozen H1-H4 results",
    )

    n0 = payload.get("N0_qualification", {})
    _require_true(n0, "locked_before_N1_outcome_inspection", "N0_qualification")
    _require_true(n0, "outcome_blind_channel_selection", "N0_qualification")

    ontology = n0.get("channel_ontology", {})
    channels = tuple(ontology.get("channels", []))
    _require(channels == EXPECTED_CHANNELS, "N0 channel ontology changed after freeze")
    _require_true(ontology, "no_posthoc_channel_addition_after_N1", "channel_ontology")
    _require_true(
        ontology,
        "exclusion_allowed_only_for_predeclared_data_quality_failure",
        "channel_ontology",
    )
    _require_true(
        ontology,
        "excluded_channels_remain_in_qualification_receipt",
        "channel_ontology",
    )
    _require_true(ontology, "floral_syndrome_scores_are_not_channels", "channel_ontology")

    source = n0.get("source_availability", {})
    _require(
        set(source.get("states", [])) == EXPECTED_SOURCE_STATES,
        "source availability states must remain available/structurally_absent/unresolved",
    )
    _require_true(
        source,
        "source_region_must_be_defined_without_focal_plant_outcomes",
        "source_availability",
    )
    _require_true(
        source,
        "structural_absence_requires_source_region_evidence",
        "source_availability",
    )
    _require_true(
        source,
        "absence_of_island_records_cannot_define_structural_absence",
        "source_availability",
    )

    island_state = n0.get("island_channel_state", {})
    _require(
        set(island_state.get("states", [])) == EXPECTED_ISLAND_STATES,
        "island channel states changed",
    )
    for key in (
        "retained_requires_source_available",
        "disrupted_requires_source_available",
        "disrupted_requires_adequate_observation_effort",
        "raw_presence_absence_is_insufficient",
        "climate_compatibility_alone_is_insufficient",
        "observation_model_must_separate_detection_from_non_detection",
    ):
        _require_true(island_state, key, "island_channel_state")

    channel_evidence = n0.get("channel_evidence", {})
    prohibited_channel = set(channel_evidence.get("prohibited", []))
    _require(
        PROHIBITED_CHANNEL_INFERENCE.issubset(prohibited_channel),
        "channel evidence no longer blocks focal floral/syndrome reverse inference",
    )

    dependency = n0.get("lineage_dependency", {})
    prohibited_dependency = set(dependency.get("prohibited_evidence", []))
    _require(
        PROHIBITED_DEPENDENCY_FIELDS.issubset(prohibited_dependency),
        "lineage dependency no longer blocks focal floral architecture leakage",
    )
    _require(
        dependency.get("reproductive_assurance_role")
        == "rival_Baker_filter_not_channel_dependency",
        "reproductive assurance must remain a Baker rival rather than channel dependency",
    )
    _require_true(
        dependency,
        "unresolved_dependency_is_not_zero",
        "lineage_dependency",
    )
    _require_true(
        dependency,
        "evidence_uncertainty_must_be_propagated",
        "lineage_dependency",
    )

    independence = n0.get("independence_checks", {})
    for key in (
        "channel_state_must_not_use_focal_plant_traits",
        "lineage_dependency_must_not_use_focal_plant_traits",
        "lineage_dependency_must_not_be_imputed_from_pollination_syndrome",
        "chapter1_pollination_architecture_concordance_is_interpretation_only",
    ):
        _require_true(independence, key, "independence_checks")

    model_families = n0.get("model_families", {})
    _require(model_families.get("plant_only", {}).get("id") == "P", "plant-only model id changed")
    _require(
        model_families.get("double_filter", {}).get("id") == "DP",
        "double-filter model id changed",
    )
    _require(
        model_families.get("double_filter", {}).get("extends") == "P",
        "DP must remain a strict extension of plant-only model P",
    )
    _require(
        model_families.get("double_filter", {}).get("additional_term")
        == "channel_retention_x_lineage_dependency",
        "DP additional term changed",
    )
    _require_true(model_families.get("baker_rival", {}), "required", "baker_rival")

    heldout = n0.get("heldout_prediction", {})
    _require(heldout.get("primary_unit") == "archipelago_id", "primary holdout must be archipelago")
    _require(
        heldout.get("split") == "leave_one_archipelago_out",
        "primary split must remain leave-one-archipelago-out",
    )
    _require(
        heldout.get("missing_primary_unit_policy") == "hard_stop",
        "missing archipelago id must remain a hard stop",
    )
    _require_true(heldout, "no_random_island_split_as_primary", "heldout_prediction")
    _require(
        heldout.get("primary_target") == "source_available_genus_entry",
        "N3 primary target must remain source-available genus entry",
    )
    _require(heldout.get("primary_metric") == "mean_log_loss", "N3 primary metric changed")
    _require_true(
        heldout,
        "trait_vector_is_secondary_downstream_endpoint",
        "heldout_prediction",
    )

    n1 = payload.get("N1_pollinator_filter", {})
    _require(
        n1.get("failure_action") == "stop_before_N2_and_keep_frozen_chapter1",
        "N1 failure must stop before N2",
    )
    n2 = payload.get("N2_channel_dependent_lineage_filter", {})
    _require(
        n2.get("failure_action") == "terminate_H5a_no_rescue",
        "N2 failure must terminate H5a without rescue",
    )
    _require(
        "Baker_colonization_assurance" in n2.get("required_rivals", []),
        "Baker/colonization assurance rival is mandatory",
    )

    n3 = payload.get("N3_heldout_double_filter_prediction", {})
    _require(n3.get("compare", {}).get("baseline") == "P", "N3 baseline must be P")
    _require(n3.get("compare", {}).get("candidate") == "DP", "N3 candidate must be DP")
    _require(
        n3.get("primary_endpoint") == "heldout_source_available_genus_entry",
        "N3 primary endpoint changed",
    )

    routing = payload.get("routing", {})
    _require(
        routing.get("N1_N2_N3_robust") == "NEE_first_submission",
        "success routing changed",
    )
    _require(
        routing.get("N2_fail") == "frozen_H3_genus_assembly_cause_unidentified",
        "N2 failure routing changed",
    )

    return {
        "contract": payload["contract"],
        "status": payload["status"],
        "parent_contract": parent["analysis_contract"],
        "frozen_results": parent["frozen_results"],
        "channels": list(channels),
        "primary_holdout": heldout["split"],
        "primary_target": heldout["primary_target"],
        "primary_metric": heldout["primary_metric"],
        "N1_failure_action": n1["failure_action"],
        "N2_failure_action": n2["failure_action"],
    }


def load_contract(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("contract YAML must contain a mapping")
    return payload


@app.command("validate")
def validate_command(
    config: Path = typer.Option(
        Path("config/chapter1_nee_double_filter.yml"),
        exists=True,
        dir_okay=False,
        help="Prospective NEE challenge YAML contract.",
    ),
    output_json: Path | None = typer.Option(None, help="Optional validation receipt path."),
) -> None:
    receipt = validate_contract_payload(load_contract(config))
    rendered = json.dumps(receipt, indent=2, sort_keys=True)
    if output_json is not None:
        output_json.parent.mkdir(parents=True, exist_ok=True)
        output_json.write_text(rendered + "\n", encoding="utf-8")
    typer.echo(rendered)


if __name__ == "__main__":
    app()
