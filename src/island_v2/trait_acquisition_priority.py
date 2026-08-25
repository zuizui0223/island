"""Evidence-audited extraction priority for the all-species trait campaign.

The priority is deliberately operational rather than biological.  It ranks
species for *direct evidence acquisition* using only reviewed geography,
current source-backed cell coverage, already-qualified genus inference
candidates, and frozen species-direct source text.  Missing inputs remain
explicitly unknown and are never converted to zero.
"""

from __future__ import annotations

import glob as globlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Iterable
from typing import Any

import pandas as pd
import typer

from island_v2 import acquisition_rate_report


TARGET_TRAITS: tuple[str, ...] = (
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
    "sex_system",
    "pollen_vector_mode",
    "flower_size_class",
    "tube_depth_class",
    "inflorescence_display",
    "reward_type",
    "dichogamy",
    "herkogamy",
)
TRAIT_PRIORITY_RANK = {trait: rank for rank, trait in enumerate(TARGET_TRAITS, start=1)}
TRAIT_PRIORITY_WEIGHT = {
    trait: len(TARGET_TRAITS) - rank + 1
    for trait, rank in TRAIT_PRIORITY_RANK.items()
}

DEFAULT_OCCURRENCE_CSV = Path(
    "data/v2/staging/gbif/collected/island_species_occurrences.csv.gz"
)
DEFAULT_ISLAND_REGISTRY_CSV = Path("data/v2/curation/pilot_island_registry.csv")
DEFAULT_CANONICAL_ISLAND_LATITUDES_CSV = Path(
    "data/v2/curation/island_canonical_latitudes.csv"
)
DEFAULT_ISLAND_ASSIGNMENTS_CSV = Path(
    "data/v2/curation/island_source_region_assignments.csv"
)
DEFAULT_ISLAND_ASSIGNMENT_EVIDENCE_CSV = Path(
    "data/v2/curation/island_assignment_evidence.csv"
)
DEFAULT_SOURCE_REGION_REGISTRY_CSV = Path(
    "data/v2/curation/source_region_bombus_registry.csv"
)
DEFAULT_SOURCE_REGION_EVIDENCE_CSV = Path(
    "data/v2/curation/source_region_evidence.csv"
)
DEFAULT_ACQUISITION_CONFIG = Path("config/acquisition_rate.yml")
DEFAULT_QUALIFIED_GENUS_GLOB = (
    "data/v2/staging/traits/**/genus_inference_candidates.csv*"
)
DEFAULT_FROZEN_PACKET_GLOB = (
    "data/v2/staging/traits/llm_evidence_extracted/*/packet_input.json"
)
TROPIC_OF_CANCER_LATITUDE = 23.4366
ARCTIC_CIRCLE_LATITUDE = 66.5634

GEOGRAPHY_AUDIT_COLUMNS = [
    "any_reviewed_NH_temperate",
    "any_NH_temperate_unreviewed",
    "any_NH_tropical",
    "any_SH",
    "crosses_regions",
    "scheduling_region_category",
    "priority_reviewed_NH_temperate_island_count",
    "priority_reviewed_NH_temperate_island_ids",
    "priority_NH_temperate_unreviewed_island_ids",
    "priority_NH_tropical_island_ids",
    "priority_SH_island_ids",
    "priority_unresolved_latitude_island_ids",
    "priority_observed_island_count",
    "priority_canonical_latitude_island_count",
    "priority_canonical_latitude_min",
    "priority_canonical_latitude_max",
    "priority_assignment_evidence_ids",
    "priority_source_urls",
    "priority_geography_evidence_basis",
    "priority_canonical_latitude_evidence_basis",
]
LEGACY_GEOGRAPHY_AUDIT_COLUMNS = [
    "any_reviewed_NH_temperate_bombus_applicable",
    "priority_nh_temperate_native_bombus",
    "priority_nh_temperate_native_bombus_island_count",
    "priority_island_ids",
    "priority_source_region_evidence_ids",
]
TRAIT_GAP_AUDIT_COLUMNS = [
    "priority_target_trait_cells",
    "source_backed_trait_cells",
    "unresolved_trait_cells",
    "highest_priority_unresolved_trait",
    "highest_priority_unresolved_rank",
    "weighted_unresolved_trait_score",
    "unresolved_priority_traits",
    "unresolved_trait_cells_state",
    "unresolved_trait_cells_basis",
]
GENUS_RESCUE_AUDIT_COLUMNS = [
    "potential_genus_rescue_cells",
    "potential_genus_rescue_state",
    "potential_genus_rescue_traits",
    "potential_genus_rescue_basis",
]
FROZEN_SOURCE_AUDIT_COLUMNS = [
    "frozen_species_direct_source",
    "frozen_species_direct_source_count",
    "frozen_species_direct_packet_ids",
    "frozen_species_direct_source_urls",
    "frozen_species_direct_evidence_scope",
    "frozen_species_direct_source_basis",
]
PRIORITY_AUDIT_COLUMNS = (
    GEOGRAPHY_AUDIT_COLUMNS
    + TRAIT_GAP_AUDIT_COLUMNS
    + GENUS_RESCUE_AUDIT_COLUMNS
    + FROZEN_SOURCE_AUDIT_COLUMNS
    + ["priority_reason"]
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _true(value: object) -> bool:
    return _text(value).casefold() in {"true", "1", "yes", "y"}


def _joined(values: Any) -> str:
    return "|".join(sorted({_text(value) for value in values if _text(value)}))


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {sorted(missing)}")


def _read_csv(path: Path, required: set[str], label: str) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    _require_columns(frame, required, label)
    return frame


def _unique_key(frame: pd.DataFrame, column: str, label: str) -> None:
    values = frame[column].map(_text)
    duplicated = values.ne("") & values.duplicated(keep=False)
    if duplicated.any():
        examples = sorted(values.loc[duplicated].unique())[:5]
        raise typer.BadParameter(f"{label} has duplicate {column}: {examples}")


def _contains_token(values: pd.Series, token: str) -> pd.Series:
    normalized = values.map(_text).str.casefold().str.replace(r"[^a-z]+", " ", regex=True)
    return normalized.str.contains(rf"(?:^| ){re.escape(token.casefold())}(?: |$)", regex=True)


def _empty_geography(master: pd.DataFrame, basis: str) -> pd.DataFrame:
    output = master[["accepted_species"]].copy()
    for column in (
        "any_reviewed_NH_temperate",
        "any_NH_temperate_unreviewed",
        "any_NH_tropical",
        "any_SH",
        "crosses_regions",
    ):
        output[column] = "unknown"
    output["scheduling_region_category"] = "unclassified"
    for column in (
        "priority_reviewed_NH_temperate_island_count",
        "priority_observed_island_count",
        "priority_canonical_latitude_island_count",
    ):
        output[column] = pd.Series(pd.NA, index=output.index, dtype="Int64")
    for column in ("priority_canonical_latitude_min", "priority_canonical_latitude_max"):
        output[column] = pd.Series(pd.NA, index=output.index, dtype="Float64")
    for column in (
        "priority_reviewed_NH_temperate_island_ids",
        "priority_NH_temperate_unreviewed_island_ids",
        "priority_NH_tropical_island_ids",
        "priority_SH_island_ids",
        "priority_unresolved_latitude_island_ids",
        "priority_assignment_evidence_ids",
        "priority_source_urls",
    ):
        output[column] = ""
    output["priority_geography_evidence_basis"] = basis
    output["priority_canonical_latitude_evidence_basis"] = basis
    return output[["accepted_species", *GEOGRAPHY_AUDIT_COLUMNS]]


def build_geography_audit(
    master: pd.DataFrame,
    *,
    occurrences: pd.DataFrame | None,
    canonical_island_latitudes: pd.DataFrame | None,
    island_registry: pd.DataFrame | None,
    island_assignments: pd.DataFrame | None,
    island_assignment_evidence: pd.DataFrame | None,
    source_region_registry: pd.DataFrame | None,
    source_region_evidence: pd.DataFrame | None,
    missing_basis: str = "unknown:reviewed geography inputs not provided",
) -> pd.DataFrame:
    """Classify geographic acquisition-priority strata.

    Hemisphere and latitude band are read only from the frozen polygon's
    canonical centroid.  The reviewed core additionally requires an explicit
    temperate biogeographic field and linked reviewed evidence.  Taxon, island,
    and region names are never parsed, and Bombus status is not a criterion.
    """

    # ``source_region_evidence`` is retained in the public signature for callers
    # that still pass the legacy Bombus-checklist table.  It has no role in this
    # geographic classification.
    del source_region_evidence
    if occurrences is None or canonical_island_latitudes is None:
        return _empty_geography(master, missing_basis)

    assert occurrences is not None
    assert canonical_island_latitudes is not None
    optional_review_inputs = {
        "island_registry": island_registry,
        "island_assignments": island_assignments,
        "island_assignment_evidence": island_assignment_evidence,
        "source_region_registry": source_region_registry,
    }
    missing_review_inputs = [
        label for label, frame in optional_review_inputs.items() if frame is None
    ]
    review_inputs_available = not missing_review_inputs

    island_registry = island_registry if island_registry is not None else pd.DataFrame(
        columns=[
            "island_id",
            "source_region_id",
            "centroid_lat",
            "identity_review_status",
            "identity_source_citation",
            "identity_source_url",
        ]
    )
    island_assignments = (
        island_assignments
        if island_assignments is not None
        else pd.DataFrame(
            columns=[
                "island_id",
                "source_region_id",
                "assignment_evidence_id",
                "review_status",
            ]
        )
    )
    island_assignment_evidence = (
        island_assignment_evidence
        if island_assignment_evidence is not None
        else pd.DataFrame(
            columns=[
                "assignment_evidence_id",
                "island_id",
                "source_citation",
                "source_url",
                "review_status",
            ]
        )
    )
    source_region_registry = (
        source_region_registry
        if source_region_registry is not None
        else pd.DataFrame(
            columns=["source_region_id", "biogeographic_region", "review_status"]
        )
    )

    occurrence_species = "species" if "species" in occurrences.columns else "accepted_species"
    _require_columns(occurrences, {occurrence_species, "island_id"}, "occurrence CSV")
    _require_columns(
        canonical_island_latitudes,
        {
            "island_id",
            "canonical_centroid_lat",
            "geometry_sha256",
            "latitude_method",
            "source_url",
            "source_file_sha256",
        },
        "canonical island latitude ledger",
    )
    _require_columns(
        island_registry,
        {
            "island_id",
            "source_region_id",
            "centroid_lat",
            "identity_review_status",
            "identity_source_citation",
            "identity_source_url",
        },
        "island registry",
    )
    _require_columns(
        island_assignments,
        {
            "island_id",
            "source_region_id",
            "assignment_evidence_id",
            "review_status",
        },
        "island assignments",
    )
    _require_columns(
        island_assignment_evidence,
        {
            "assignment_evidence_id",
            "island_id",
            "source_citation",
            "source_url",
            "review_status",
        },
        "island assignment evidence",
    )
    _require_columns(
        source_region_registry,
        {
            "source_region_id",
            "biogeographic_region",
            "review_status",
        },
        "source-region registry",
    )
    _unique_key(canonical_island_latitudes, "island_id", "canonical island latitude ledger")
    _unique_key(island_registry, "island_id", "island registry")
    _unique_key(island_assignments, "island_id", "island assignments")
    _unique_key(
        island_assignment_evidence,
        "assignment_evidence_id",
        "island assignment evidence",
    )
    _unique_key(source_region_registry, "source_region_id", "source-region registry")

    observed = occurrences[[occurrence_species, "island_id"]].rename(
        columns={occurrence_species: "accepted_species"}
    )
    observed["accepted_species"] = observed["accepted_species"].map(_text)
    observed["island_id"] = observed["island_id"].map(_text)
    observed = observed.loc[
        observed["accepted_species"].isin(set(master["accepted_species"].map(_text)))
        & observed["island_id"].ne("")
    ].drop_duplicates()

    latitudes = canonical_island_latitudes[
        [
            "island_id",
            "canonical_centroid_lat",
            "geometry_sha256",
            "latitude_method",
            "source_url",
            "source_file_sha256",
        ]
    ].rename(
        columns={
            "geometry_sha256": "latitude_geometry_sha256",
            "source_url": "latitude_source_url",
            "source_file_sha256": "latitude_source_file_sha256",
        }
    )

    islands = island_registry[
        [
            "island_id",
            "source_region_id",
            "centroid_lat",
            "identity_review_status",
            "identity_source_citation",
            "identity_source_url",
        ]
    ].rename(
        columns={
            "source_region_id": "island_source_region_id",
            "identity_review_status": "island_review_status",
            "identity_source_citation": "island_source_citation",
            "identity_source_url": "island_source_url",
        }
    )
    assignments = island_assignments[
        ["island_id", "source_region_id", "assignment_evidence_id", "review_status"]
    ].rename(columns={"review_status": "assignment_review_status"})
    assignment_evidence = island_assignment_evidence[
        [
            "assignment_evidence_id",
            "island_id",
            "source_citation",
            "source_url",
            "review_status",
        ]
    ].rename(
        columns={
            "island_id": "assignment_evidence_island_id",
            "source_citation": "assignment_source_citation",
            "source_url": "assignment_source_url",
            "review_status": "assignment_evidence_review_status",
        }
    )
    regions = source_region_registry[
        [
            "source_region_id",
            "biogeographic_region",
            "review_status",
        ]
    ].rename(columns={"review_status": "source_region_review_status"})

    joined = observed.merge(latitudes, on="island_id", how="left", validate="many_to_one")
    joined = joined.merge(islands, on="island_id", how="left", validate="many_to_one")
    joined = joined.merge(assignments, on="island_id", how="left", validate="many_to_one")
    joined = joined.merge(
        assignment_evidence,
        on="assignment_evidence_id",
        how="left",
        validate="many_to_one",
    )
    joined = joined.merge(regions, on="source_region_id", how="left", validate="many_to_one")

    latitude = pd.to_numeric(joined["canonical_centroid_lat"], errors="coerce")
    reviewed_latitude = pd.to_numeric(joined["centroid_lat"], errors="coerce")
    latitude_matches_review = latitude.sub(reviewed_latitude).abs().le(0.01)
    reviewed = (
        joined["island_review_status"].map(_text).str.casefold().eq("accepted")
        & joined["assignment_review_status"].map(_text).str.casefold().eq("accepted")
        & joined["assignment_evidence_review_status"]
        .map(_text)
        .str.casefold()
        .eq("accepted")
        & joined["source_region_review_status"].map(_text).str.casefold().eq("accepted")
        & joined["source_region_id"].map(_text).ne("")
        & joined["source_region_id"].map(_text).eq(joined["island_source_region_id"].map(_text))
        & joined["assignment_evidence_island_id"].map(_text).eq(joined["island_id"].map(_text))
        & joined["island_source_citation"].map(_text).ne("")
        & joined["island_source_url"].map(_text).ne("")
        & joined["assignment_source_citation"].map(_text).ne("")
        & joined["assignment_source_url"].map(_text).ne("")
        & latitude.notna()
        & latitude_matches_review
    )
    temperate = _contains_token(joined["biogeographic_region"], "temperate")
    nh_temperate = latitude.ge(TROPIC_OF_CANCER_LATITUDE) & latitude.lt(
        ARCTIC_CIRCLE_LATITUDE
    )
    reviewed_nh_temperate = (
        reviewed
        & nh_temperate
        & temperate
    )
    joined["_reviewed_nh_temperate"] = reviewed_nh_temperate
    joined["_latitude_known"] = latitude.notna()
    joined["_nh_temperate_unreviewed"] = nh_temperate & ~reviewed_nh_temperate
    joined["_nh_tropical"] = latitude.ge(0) & latitude.lt(TROPIC_OF_CANCER_LATITUDE)
    joined["_sh"] = latitude.lt(0)
    joined["_geographic_band_known"] = nh_temperate | joined["_nh_tropical"] | joined["_sh"]
    joined["_canonical_latitude"] = latitude

    # Keep the million-row occurrence join vectorized.  Building one pandas
    # sub-frame per species makes a wave refresh take minutes.
    counts = (
        joined.groupby("accepted_species", sort=False)
        .agg(
            _n_observed_islands=("island_id", "nunique"),
            _n_reviewed_nh_temperate=("_reviewed_nh_temperate", "sum"),
            _n_latitude_islands=("_latitude_known", "sum"),
            _n_classified_islands=("_geographic_band_known", "sum"),
            _n_nh_temperate_unreviewed=("_nh_temperate_unreviewed", "sum"),
            _n_nh_tropical=("_nh_tropical", "sum"),
            _n_sh=("_sh", "sum"),
            priority_canonical_latitude_min=("_canonical_latitude", "min"),
            priority_canonical_latitude_max=("_canonical_latitude", "max"),
        )
        .reset_index()
    )

    def any_state(count_column: str) -> pd.Series:
        state = pd.Series("unknown", index=counts.index, dtype="object")
        complete_classification = counts["_n_classified_islands"].eq(
            counts["_n_observed_islands"]
        )
        state.loc[complete_classification] = "false"
        state.loc[counts[count_column].gt(0)] = "true"
        return state

    counts["any_reviewed_NH_temperate"] = any_state("_n_reviewed_nh_temperate")
    if not review_inputs_available:
        counts["any_reviewed_NH_temperate"] = "unknown"
    counts["any_NH_temperate_unreviewed"] = any_state("_n_nh_temperate_unreviewed")
    counts["any_NH_tropical"] = any_state("_n_nh_tropical")
    counts["any_SH"] = any_state("_n_sh")
    any_nh_temperate = counts["_n_reviewed_nh_temperate"].add(
        counts["_n_nh_temperate_unreviewed"]
    ).gt(0)
    present_region_count = pd.concat(
        [
            any_nh_temperate,
            counts["_n_nh_tropical"].gt(0),
            counts["_n_sh"].gt(0),
        ],
        axis=1,
    ).sum(axis=1)
    counts["crosses_regions"] = "unknown"
    complete_classification = counts["_n_classified_islands"].eq(
        counts["_n_observed_islands"]
    )
    counts.loc[
        complete_classification & present_region_count.le(1), "crosses_regions"
    ] = "false"
    counts.loc[present_region_count.gt(1), "crosses_regions"] = "true"
    counts["scheduling_region_category"] = "unclassified"
    counts.loc[
        counts["any_NH_temperate_unreviewed"].eq("true"),
        "scheduling_region_category",
    ] = "nh_temperate_unreviewed"
    counts.loc[
        counts["any_reviewed_NH_temperate"].eq("true"),
        "scheduling_region_category",
    ] = "reviewed_nh_temperate"
    counts.loc[
        counts["any_NH_tropical"].eq("true"),
        "scheduling_region_category",
    ] = "nh_tropical"
    counts.loc[counts["any_SH"].eq("true"), "scheduling_region_category"] = (
        "southern_hemisphere"
    )
    counts.loc[~complete_classification, "scheduling_region_category"] = "unclassified"
    counts.loc[counts["crosses_regions"].eq("true"), "scheduling_region_category"] = (
        "cross_region"
    )

    def island_ids(mask: pd.Series, column: str) -> pd.DataFrame:
        selected = joined.loc[mask, ["accepted_species", "island_id"]]
        if selected.empty:
            return pd.DataFrame(columns=["accepted_species", column])
        return (
            selected.groupby("accepted_species", sort=False)["island_id"]
            .agg(_joined)
            .rename(column)
            .reset_index()
        )

    category_ids = island_ids(
        joined["_reviewed_nh_temperate"],
        "priority_reviewed_NH_temperate_island_ids",
    )
    for ids in (
        island_ids(
            joined["_nh_temperate_unreviewed"],
            "priority_NH_temperate_unreviewed_island_ids",
        ),
        island_ids(joined["_nh_tropical"], "priority_NH_tropical_island_ids"),
        island_ids(joined["_sh"], "priority_SH_island_ids"),
        island_ids(
            ~joined["_latitude_known"],
            "priority_unresolved_latitude_island_ids",
        ),
    ):
        category_ids = category_ids.merge(
            ids, on="accepted_species", how="outer", validate="one_to_one"
        )

    evidence_rows = joined.loc[joined["_reviewed_nh_temperate"]].copy()
    if evidence_rows.empty:
        evidence = pd.DataFrame(
            columns=[
                "accepted_species",
                "priority_assignment_evidence_ids",
                "priority_source_urls",
            ]
        )
    else:
        evidence = (
            evidence_rows.groupby("accepted_species", sort=False)
            .agg(
                priority_assignment_evidence_ids=("assignment_evidence_id", _joined),
            )
            .reset_index()
        )
        url_rows = evidence_rows.melt(
            id_vars=["accepted_species"],
            value_vars=[
                "latitude_source_url",
                "island_source_url",
                "assignment_source_url",
            ],
            value_name="source_url",
        )
        urls = (
            url_rows.groupby("accepted_species", sort=False)["source_url"]
            .agg(_joined)
            .rename("priority_source_urls")
            .reset_index()
        )
        evidence = evidence.merge(urls, on="accepted_species", how="left", validate="one_to_one")

    output = master[["accepted_species"]].copy()
    output["accepted_species"] = output["accepted_species"].map(_text)
    output = output.merge(counts, on="accepted_species", how="left", validate="one_to_one")
    for column in (
        "any_reviewed_NH_temperate",
        "any_NH_temperate_unreviewed",
        "any_NH_tropical",
        "any_SH",
        "crosses_regions",
    ):
        output[column] = output[column].fillna("unknown")
    output["scheduling_region_category"] = output["scheduling_region_category"].fillna(
        "unclassified"
    )
    if review_inputs_available:
        output["priority_reviewed_NH_temperate_island_count"] = output[
            "_n_reviewed_nh_temperate"
        ]
    else:
        output["priority_reviewed_NH_temperate_island_count"] = pd.NA
    output["priority_observed_island_count"] = output["_n_observed_islands"]
    output["priority_canonical_latitude_island_count"] = output["_n_latitude_islands"]
    output = output.merge(category_ids, on="accepted_species", how="left", validate="one_to_one")
    output = output.merge(evidence, on="accepted_species", how="left", validate="one_to_one")
    for column in (
        "priority_reviewed_NH_temperate_island_ids",
        "priority_NH_temperate_unreviewed_island_ids",
        "priority_NH_tropical_island_ids",
        "priority_SH_island_ids",
        "priority_unresolved_latitude_island_ids",
        "priority_assignment_evidence_ids",
        "priority_source_urls",
    ):
        output[column] = output[column].fillna("")
    output["priority_geography_evidence_basis"] = output["scheduling_region_category"].map(
        {
            "reviewed_nh_temperate": (
                "accepted linked island and source-region evidence explicitly identifies a "
                "temperate region between the Tropic of Cancer and Arctic Circle"
            ),
            "nh_temperate_unreviewed": (
                "canonical centroid latitude is north of the Tropic of Cancer but linked "
                "temperate geography evidence is not fully reviewed; Arctic latitudes are excluded"
            ),
            "nh_tropical": "canonical centroid latitude is in the Northern Hemisphere tropics",
            "southern_hemisphere": "canonical centroid latitude is in the Southern Hemisphere",
            "cross_region": "exact occurrence islands span more than one geographic latitude band",
            "unclassified": "one or more exact occurrence islands lack a complete geographic classification",
        }
    )
    no_occurrence = output["_n_observed_islands"].isna()
    output.loc[no_occurrence, "priority_geography_evidence_basis"] = (
        "unknown:no exact species-island occurrence"
    )
    if missing_review_inputs:
        review_missing_mask = output["priority_geography_evidence_basis"].ne(
            "unknown:no exact species-island occurrence"
        )
        output.loc[review_missing_mask, "priority_geography_evidence_basis"] = (
            output.loc[review_missing_mask, "priority_geography_evidence_basis"]
            + "; reviewed geography unavailable:missing="
            + ",".join(missing_review_inputs)
        )
    latitude_basis = (
        "exact occurrence island_id joined to frozen GSHHG 2.3.7 canonical centroid; "
        f"NH-temperate band={TROPIC_OF_CANCER_LATITUDE}<=latitude<"
        f"{ARCTIC_CIRCLE_LATITUDE}; methods="
        f"{_joined(canonical_island_latitudes['latitude_method'])}; source_urls="
        f"{_joined(canonical_island_latitudes['source_url'])}; source_file_sha256="
        f"{_joined(canonical_island_latitudes['source_file_sha256'])}"
    )
    output["priority_canonical_latitude_evidence_basis"] = latitude_basis
    output = output[["accepted_species", *GEOGRAPHY_AUDIT_COLUMNS]]
    for column in (
        "priority_reviewed_NH_temperate_island_count",
        "priority_observed_island_count",
        "priority_canonical_latitude_island_count",
    ):
        output[column] = pd.array(output[column], dtype="Int64")
    for column in ("priority_canonical_latitude_min", "priority_canonical_latitude_max"):
        output[column] = pd.array(output[column], dtype="Float64")
    return output


@dataclass(frozen=True)
class TraitCoverageSnapshot:
    """Compact current direct-coverage state for the priority target cells."""

    known_species: frozenset[str]
    applicable_species: frozenset[str]
    source_backed_pairs: frozenset[tuple[str, str]]
    target_traits: tuple[str, ...]
    basis: str

    def cell_is_unresolved(self, species: str, trait: str) -> bool:
        return (
            species in self.applicable_species
            and trait in self.target_traits
            and (species, trait) not in self.source_backed_pairs
        )


def trait_snapshot_from_acquisition_config(
    master: pd.DataFrame,
    config_path: Path,
) -> TraitCoverageSnapshot:
    """Build current direct-cell coverage from the repository's existing report contract."""

    config = acquisition_rate_report.load_config(config_path)
    master_path = Path(config["master_taxa_csv"])
    species_column = str(config["master_species_column"])
    configured_master = acquisition_rate_report.load_master_species(master_path, species_column)
    applicability = acquisition_rate_report.load_master_applicability(
        master_path,
        species_column,
        config,
    )
    if applicability is None:
        raise typer.BadParameter(
            "priority acquisition config must provide scope_config_path so non-angiosperms "
            "remain not_applicable"
        )

    broad_config = dict(config)
    broad_config["sources"] = {
        key: value
        for key, value in config["sources"].items()
        if bool(value.get("counts_as_broad", True))
    }
    contributions, _ = acquisition_rate_report.collect_contributions(
        broad_config, set(configured_master)
    )
    source_backed = contributions.loc[
        contributions["counts_as_broad"].astype(bool)
        & contributions["trait_name"].isin(TARGET_TRAITS),
        ["accepted_species", "trait_name"],
    ].drop_duplicates()

    master_species = set(master["accepted_species"].map(_text))
    applicability = applicability.copy()
    applicability["accepted_species"] = applicability["accepted_species"].map(_text)
    applicable_flag = applicability["angiosperm_analysis_eligible"].map(_true)
    known_species = frozenset(
        applicability.loc[applicability["accepted_species"].isin(master_species), "accepted_species"]
    )
    applicable_species = frozenset(
        applicability.loc[
            applicability["accepted_species"].isin(master_species) & applicable_flag,
            "accepted_species",
        ]
    )
    source_pairs = frozenset(
        (species, trait)
        for species, trait in source_backed.itertuples(index=False, name=None)
        if species in master_species
    )
    source_keys = ",".join(sorted(broad_config["sources"]))
    return TraitCoverageSnapshot(
        known_species=known_species,
        applicable_species=applicable_species,
        source_backed_pairs=source_pairs,
        target_traits=TARGET_TRAITS,
        basis=f"{config_path}:source-backed channels={source_keys}",
    )


def trait_snapshot_from_cells(
    master: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    basis: str,
) -> TraitCoverageSnapshot:
    """Normalize an explicit complete cell ledger for tests or external snapshots."""

    frame = cells.copy().fillna("")
    if {"accepted_species", "trait_name", "trait_applicable", "source_backed"}.issubset(
        frame.columns
    ):
        frame = frame[["accepted_species", "trait_name", "trait_applicable", "source_backed"]]
        frame["_applicable"] = frame["trait_applicable"].map(_true)
        frame["_source_backed"] = frame["source_backed"].map(_true)
    elif {"species", "field", "applicability", "selected_origin"}.issubset(frame.columns):
        frame = frame[["species", "field", "applicability", "selected_origin"]].rename(
            columns={"species": "accepted_species", "field": "trait_name"}
        )
        frame["_applicable"] = frame["applicability"].map(_text).str.casefold().eq("applicable")
        frame["_source_backed"] = frame["selected_origin"].map(_text).str.casefold().eq(
            "reported"
        )
    else:
        raise typer.BadParameter(
            "priority cell snapshot needs acquisition ledger columns or all-species value-channel "
            "columns"
        )

    frame["accepted_species"] = frame["accepted_species"].map(_text)
    frame["trait_name"] = frame["trait_name"].map(_text)
    frame = frame.loc[frame["trait_name"].isin(TARGET_TRAITS)].drop_duplicates(
        ["accepted_species", "trait_name"], keep="last"
    )
    master_species = set(master["accepted_species"].map(_text))
    expected = len(TARGET_TRAITS)
    counts = frame.loc[frame["accepted_species"].isin(master_species)].groupby(
        "accepted_species"
    )["trait_name"].nunique()
    known = set(counts.loc[counts.eq(expected)].index)
    applicable_by_species = frame.loc[frame["accepted_species"].isin(known)].groupby(
        "accepted_species"
    )["_applicable"].agg(["any", "all"])
    mixed = applicable_by_species[applicable_by_species["any"] != applicable_by_species["all"]]
    if not mixed.empty:
        raise typer.BadParameter(
            "priority target cells must share one angiosperm applicability state per species"
        )
    applicable = set(applicable_by_species.loc[applicable_by_species["all"]].index)
    direct = frame.loc[
        frame["accepted_species"].isin(applicable) & frame["_source_backed"]
    ]
    return TraitCoverageSnapshot(
        known_species=frozenset(known),
        applicable_species=frozenset(applicable),
        source_backed_pairs=frozenset(
            direct[["accepted_species", "trait_name"]].itertuples(index=False, name=None)
        ),
        target_traits=TARGET_TRAITS,
        basis=basis,
    )


def build_trait_gap_audit(
    master: pd.DataFrame,
    snapshot: TraitCoverageSnapshot | None,
    *,
    missing_basis: str = "unknown:current direct-cell input not provided",
) -> pd.DataFrame:
    output = master[["accepted_species"]].copy()
    rows: list[dict[str, Any]] = []
    for species in output["accepted_species"].map(_text):
        if snapshot is None or species not in snapshot.known_species:
            rows.append(
                {
                    "priority_target_trait_cells": pd.NA,
                    "source_backed_trait_cells": pd.NA,
                    "unresolved_trait_cells": pd.NA,
                    "highest_priority_unresolved_trait": "unknown",
                    "highest_priority_unresolved_rank": pd.NA,
                    "weighted_unresolved_trait_score": pd.NA,
                    "unresolved_priority_traits": "",
                    "unresolved_trait_cells_state": "unknown",
                    "unresolved_trait_cells_basis": missing_basis
                    if snapshot is None
                    else f"unknown:species absent from {snapshot.basis}",
                }
            )
        elif species not in snapshot.applicable_species:
            rows.append(
                {
                    "priority_target_trait_cells": pd.NA,
                    "source_backed_trait_cells": pd.NA,
                    "unresolved_trait_cells": pd.NA,
                    "highest_priority_unresolved_trait": "not_applicable",
                    "highest_priority_unresolved_rank": pd.NA,
                    "weighted_unresolved_trait_score": pd.NA,
                    "unresolved_priority_traits": "",
                    "unresolved_trait_cells_state": "not_applicable",
                    "unresolved_trait_cells_basis": f"non-angiosperm scope; {snapshot.basis}",
                }
            )
        else:
            unresolved = [
                trait
                for trait in snapshot.target_traits
                if (species, trait) not in snapshot.source_backed_pairs
            ]
            direct = len(snapshot.target_traits) - len(unresolved)
            highest = unresolved[0] if unresolved else "complete"
            highest_rank = TRAIT_PRIORITY_RANK.get(highest, len(TARGET_TRAITS) + 1)
            rows.append(
                {
                    "priority_target_trait_cells": len(snapshot.target_traits),
                    "source_backed_trait_cells": direct,
                    "unresolved_trait_cells": len(unresolved),
                    "highest_priority_unresolved_trait": highest,
                    "highest_priority_unresolved_rank": highest_rank,
                    "weighted_unresolved_trait_score": sum(
                        TRAIT_PRIORITY_WEIGHT.get(trait, 0) for trait in unresolved
                    ),
                    "unresolved_priority_traits": "|".join(unresolved),
                    "unresolved_trait_cells_state": "known",
                    "unresolved_trait_cells_basis": snapshot.basis,
                }
            )
    detail = pd.DataFrame(rows)
    for column in (
        "priority_target_trait_cells",
        "source_backed_trait_cells",
        "unresolved_trait_cells",
        "highest_priority_unresolved_rank",
        "weighted_unresolved_trait_score",
    ):
        detail[column] = pd.array(detail[column], dtype="Int64")
    return pd.concat([output.reset_index(drop=True), detail], axis=1)


def load_qualified_genus_candidates(pattern: str) -> tuple[pd.DataFrame | None, str]:
    paths = [Path(path) for path in sorted(globlib.glob(pattern, recursive=True))]
    if not paths:
        return None, f"unknown:no files matched {pattern}"
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        _require_columns(frame, {"genus", "trait_name"}, str(path))
        frame["_priority_input_file"] = str(path)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True), _joined(paths)


def build_genus_rescue_audit(
    master: pd.DataFrame,
    snapshot: TraitCoverageSnapshot | None,
    qualified_candidates: pd.DataFrame | None,
    *,
    genus_is_explicit: bool,
    input_basis: str,
) -> pd.DataFrame:
    output = master[["accepted_species"]].copy()
    output["potential_genus_rescue_cells"] = pd.Series(
        pd.NA, index=output.index, dtype="Int64"
    )
    output["potential_genus_rescue_state"] = "unknown"
    output["potential_genus_rescue_traits"] = ""
    if not genus_is_explicit:
        output["potential_genus_rescue_basis"] = (
            "unknown:master has no explicit genus column; species names were not parsed"
        )
        return output
    if snapshot is None:
        output["potential_genus_rescue_basis"] = "unknown:current direct-cell input missing"
        return output
    if qualified_candidates is None:
        output["potential_genus_rescue_basis"] = input_basis
        return output

    candidates = qualified_candidates.copy().fillna("")
    _require_columns(candidates, {"genus", "trait_name"}, "qualified genus candidates")
    if "evidence_scope" in candidates.columns:
        candidates = candidates.loc[
            candidates["evidence_scope"].map(_text).str.casefold().eq("genus_inference")
        ]
    if "candidate_kind" in candidates.columns:
        candidates = candidates.loc[
            candidates["candidate_kind"]
            .map(_text)
            .str.casefold()
            .eq("hierarchical_inference")
        ]
    if "review_status" in candidates.columns:
        candidates = candidates.loc[
            ~candidates["review_status"].map(_text).str.casefold().isin({"rejected", "invalid"})
        ]
    candidates["genus"] = candidates["genus"].map(_text)
    candidates["trait_name"] = candidates["trait_name"].map(_text)
    candidates = candidates.loc[candidates["trait_name"].isin(snapshot.target_traits)]

    master_lookup = master[["accepted_species", "genus"]].copy()
    master_lookup["accepted_species"] = master_lookup["accepted_species"].map(_text)
    master_lookup["genus"] = master_lookup["genus"].map(_text)
    _unique_key(master_lookup, "accepted_species", "priority master")
    genus_by_species = master_lookup.set_index("accepted_species")["genus"].to_dict()

    if "accepted_species" in candidates.columns:
        candidates["accepted_species"] = candidates["accepted_species"].map(_text)
        on_master = candidates["accepted_species"].isin(genus_by_species)
        mismatched = candidates.loc[
            on_master
            & candidates["genus"].ne(candidates["accepted_species"].map(genus_by_species))
        ]
        if not mismatched.empty:
            examples = mismatched["accepted_species"].drop_duplicates().head(5).tolist()
            raise typer.BadParameter(
                f"qualified genus candidates disagree with explicit master genus: {examples}"
            )
        candidate_cells = candidates.loc[
            on_master, ["accepted_species", "genus", "trait_name"]
        ].drop_duplicates()
    else:
        qualified_pairs = candidates[["genus", "trait_name"]].drop_duplicates()
        candidate_cells = master_lookup.merge(qualified_pairs, on="genus", how="inner")

    candidate_cells = candidate_cells.loc[
        candidate_cells.apply(
            lambda row: snapshot.cell_is_unresolved(
                _text(row["accepted_species"]), _text(row["trait_name"])
            ),
            axis=1,
        )
    ].drop_duplicates()
    grouped_cells: dict[tuple[str, str], set[str]] = {
        (genus, trait): set(group["accepted_species"].map(_text))
        for (genus, trait), group in candidate_cells.groupby(["genus", "trait_name"], sort=False)
    }

    counts: list[object] = []
    states: list[str] = []
    traits: list[str] = []
    bases: list[str] = []
    for record in master_lookup.to_dict("records"):
        species = _text(record["accepted_species"])
        genus = _text(record["genus"])
        if species not in snapshot.known_species:
            counts.append(pd.NA)
            states.append("unknown")
            traits.append("")
            bases.append(f"unknown:species absent from {snapshot.basis}")
            continue
        if species not in snapshot.applicable_species:
            counts.append(pd.NA)
            states.append("not_applicable")
            traits.append("")
            bases.append(f"non-angiosperm scope; {input_basis}")
            continue
        if not genus:
            counts.append(pd.NA)
            states.append("unknown")
            traits.append("")
            bases.append("unknown:blank explicit genus")
            continue

        by_trait: list[tuple[str, int]] = []
        for trait in snapshot.target_traits:
            if not snapshot.cell_is_unresolved(species, trait):
                continue
            rescued = grouped_cells.get((genus, trait), set()) - {species}
            if rescued:
                by_trait.append((trait, len(rescued)))
        counts.append(sum(value for _, value in by_trait))
        states.append("known")
        traits.append("|".join(f"{trait}:{value}" for trait, value in by_trait))
        bases.append(f"qualified genus candidates={input_basis}; exact explicit master genus")

    output["potential_genus_rescue_cells"] = pd.array(counts, dtype="Int64")
    output["potential_genus_rescue_state"] = states
    output["potential_genus_rescue_traits"] = traits
    output["potential_genus_rescue_basis"] = bases
    return output


def load_frozen_packet_sources(pattern: str) -> tuple[pd.DataFrame | None, str]:
    paths = [Path(path) for path in sorted(globlib.glob(pattern, recursive=True))]
    if not paths:
        return None, f"unknown:no files matched {pattern}"
    rows: list[dict[str, Any]] = []
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise typer.BadParameter(f"invalid frozen packet {path}: {exc}") from exc
        sources = payload.get("sources")
        if not isinstance(sources, list):
            raise typer.BadParameter(f"frozen packet {path} has no sources list")
        for source in sources:
            if not isinstance(source, dict):
                raise typer.BadParameter(f"frozen packet {path} contains a non-object source")
            if "accepted_species" not in source or "evidence_scope" not in source:
                raise typer.BadParameter(
                    f"frozen packet {path} source lacks accepted_species/evidence_scope"
                )
            rows.append(
                {
                    "accepted_species": _text(source.get("accepted_species")),
                    "evidence_scope": _text(source.get("evidence_scope")),
                    "source_chunk_id": _text(source.get("source_chunk_id")),
                    "source_url": _text(source.get("source_url")),
                    "source_citation": _text(source.get("source_citation")),
                    "packet_id": path.parent.name,
                }
            )
    columns = [
        "accepted_species",
        "evidence_scope",
        "source_chunk_id",
        "source_url",
        "source_citation",
        "packet_id",
    ]
    return pd.DataFrame(rows, columns=columns), _joined(paths)


def build_frozen_source_audit(
    master: pd.DataFrame,
    frozen_sources: pd.DataFrame | None,
    *,
    input_basis: str,
) -> pd.DataFrame:
    output = master[["accepted_species"]].copy()
    if frozen_sources is None:
        output["frozen_species_direct_source"] = "unknown"
        output["frozen_species_direct_source_count"] = pd.Series(
            pd.NA, index=output.index, dtype="Int64"
        )
        output["frozen_species_direct_packet_ids"] = ""
        output["frozen_species_direct_source_urls"] = ""
        output["frozen_species_direct_evidence_scope"] = ""
        output["frozen_species_direct_source_basis"] = input_basis
        return output

    sources = frozen_sources.copy().fillna("")
    _require_columns(sources, {"accepted_species", "evidence_scope"}, "frozen sources")
    sources["accepted_species"] = sources["accepted_species"].map(_text)
    direct = sources.loc[
        sources["evidence_scope"].map(_text).str.casefold().eq("species_direct")
    ].copy()
    if "source_chunk_id" not in direct.columns:
        direct["source_chunk_id"] = ""
    if "packet_id" not in direct.columns:
        direct["packet_id"] = ""
    if "source_url" not in direct.columns:
        direct["source_url"] = ""
    groups = {species: group for species, group in direct.groupby("accepted_species", sort=False)}
    states: list[str] = []
    counts: list[int] = []
    packets: list[str] = []
    urls: list[str] = []
    scopes: list[str] = []
    for species in output["accepted_species"].map(_text):
        group = groups.get(species)
        if group is None or group.empty:
            states.append("false")
            counts.append(0)
            packets.append("")
            urls.append("")
            scopes.append("")
        else:
            states.append("true")
            chunk_ids = {_text(value) for value in group["source_chunk_id"] if _text(value)}
            counts.append(len(chunk_ids) if chunk_ids else len(group))
            packets.append(_joined(group["packet_id"]))
            urls.append(_joined(group["source_url"]))
            scopes.append("species_direct")
    output["frozen_species_direct_source"] = states
    output["frozen_species_direct_source_count"] = pd.array(counts, dtype="Int64")
    output["frozen_species_direct_packet_ids"] = packets
    output["frozen_species_direct_source_urls"] = urls
    output["frozen_species_direct_evidence_scope"] = scopes
    output["frozen_species_direct_source_basis"] = input_basis
    return output


def _priority_reason(record: dict[str, Any]) -> str:
    def number(column: str, state_column: str) -> str:
        state = _text(record.get(state_column))
        value = record.get(column)
        if state == "known" and value is not None and not pd.isna(value):
            return str(int(value))
        return state or "unknown"

    return "; ".join(
        [
            "scheduling_region_category="
            + (_text(record.get("scheduling_region_category")) or "unclassified"),
            "reviewed_nh_temperate="
            + (_text(record.get("any_reviewed_NH_temperate")) or "unknown"),
            "nh_temperate_unreviewed="
            + (_text(record.get("any_NH_temperate_unreviewed")) or "unknown"),
            "nh_tropical=" + (_text(record.get("any_NH_tropical")) or "unknown"),
            "southern_hemisphere=" + (_text(record.get("any_SH")) or "unknown"),
            "unresolved_trait_cells="
            + number("unresolved_trait_cells", "unresolved_trait_cells_state"),
            "highest_priority_unresolved_trait="
            + (_text(record.get("highest_priority_unresolved_trait")) or "unknown"),
            "potential_genus_rescue_cells="
            + number("potential_genus_rescue_cells", "potential_genus_rescue_state"),
            "frozen_species_direct_source="
            + (_text(record.get("frozen_species_direct_source")) or "unknown"),
        ]
    )


def build_priority_audit(
    master: pd.DataFrame,
    *,
    geography: pd.DataFrame,
    trait_gaps: pd.DataFrame,
    genus_rescue: pd.DataFrame,
    frozen_sources: pd.DataFrame,
) -> pd.DataFrame:
    """Merge four one-row-per-species components and add a human-readable reason."""

    base = master[["accepted_species"]].copy()
    _unique_key(base, "accepted_species", "priority master")
    for label, component in (
        ("geography priority", geography),
        ("trait-gap priority", trait_gaps),
        ("genus-rescue priority", genus_rescue),
        ("frozen-source priority", frozen_sources),
    ):
        _require_columns(component, {"accepted_species"}, label)
        _unique_key(component, "accepted_species", label)
        base = base.merge(component, on="accepted_species", how="left", validate="one_to_one")
    base["priority_reason"] = [
        _priority_reason(record) for record in base.to_dict("records")
    ]
    missing = set(PRIORITY_AUDIT_COLUMNS).difference(base.columns)
    if missing:
        raise RuntimeError(f"priority audit construction omitted columns: {sorted(missing)}")
    return base[["accepted_species", *PRIORITY_AUDIT_COLUMNS]]


def build_repository_priority_audit(
    master: pd.DataFrame,
    *,
    master_genus_is_explicit: bool,
    occurrence_csv: Path = DEFAULT_OCCURRENCE_CSV,
    canonical_island_latitudes_csv: Path = DEFAULT_CANONICAL_ISLAND_LATITUDES_CSV,
    island_registry_csv: Path = DEFAULT_ISLAND_REGISTRY_CSV,
    island_assignments_csv: Path = DEFAULT_ISLAND_ASSIGNMENTS_CSV,
    island_assignment_evidence_csv: Path = DEFAULT_ISLAND_ASSIGNMENT_EVIDENCE_CSV,
    source_region_registry_csv: Path = DEFAULT_SOURCE_REGION_REGISTRY_CSV,
    source_region_evidence_csv: Path = DEFAULT_SOURCE_REGION_EVIDENCE_CSV,
    acquisition_config_path: Path = DEFAULT_ACQUISITION_CONFIG,
    current_cells_csv: Path | None = None,
    qualified_genus_glob: str = DEFAULT_QUALIFIED_GENUS_GLOB,
    frozen_packet_glob: str = DEFAULT_FROZEN_PACKET_GLOB,
) -> pd.DataFrame:
    """Load the repository's real inputs and produce a refreshable ledger audit."""

    required_geography_paths = {
        "occurrence": occurrence_csv,
        "canonical_island_latitudes": canonical_island_latitudes_csv,
    }
    missing_geography = [
        label for label, path in required_geography_paths.items() if not path.exists()
    ]
    if missing_geography:
        geography = build_geography_audit(
            master,
            occurrences=None,
            canonical_island_latitudes=None,
            island_registry=None,
            island_assignments=None,
            island_assignment_evidence=None,
            source_region_registry=None,
            source_region_evidence=None,
            missing_basis="unknown:missing geography inputs=" + ",".join(missing_geography),
        )
    else:
        geography = build_geography_audit(
            master,
            occurrences=_read_csv(
                occurrence_csv, {"species", "island_id"}, "occurrence CSV"
            ),
            canonical_island_latitudes=_read_csv(
                canonical_island_latitudes_csv,
                {"island_id", "canonical_centroid_lat"},
                "canonical island latitude ledger",
            ),
            island_registry=(
                _read_csv(island_registry_csv, {"island_id"}, "island registry")
                if island_registry_csv.exists()
                else None
            ),
            island_assignments=(
                _read_csv(island_assignments_csv, {"island_id"}, "island assignments")
                if island_assignments_csv.exists()
                else None
            ),
            island_assignment_evidence=(
                _read_csv(
                    island_assignment_evidence_csv,
                    {"assignment_evidence_id"},
                    "island assignment evidence",
                )
                if island_assignment_evidence_csv.exists()
                else None
            ),
            source_region_registry=(
                _read_csv(
                    source_region_registry_csv,
                    {"source_region_id"},
                    "source-region registry",
                )
                if source_region_registry_csv.exists()
                else None
            ),
            source_region_evidence=None,
        )

    snapshot: TraitCoverageSnapshot | None
    gap_missing_basis: str
    if current_cells_csv is not None:
        if current_cells_csv.exists():
            snapshot = trait_snapshot_from_cells(
                master,
                pd.read_csv(current_cells_csv, dtype=str).fillna(""),
                basis=str(current_cells_csv),
            )
            gap_missing_basis = ""
        else:
            snapshot = None
            gap_missing_basis = f"unknown:current cell file missing={current_cells_csv}"
    elif acquisition_config_path.exists():
        snapshot = trait_snapshot_from_acquisition_config(master, acquisition_config_path)
        gap_missing_basis = ""
    else:
        snapshot = None
        gap_missing_basis = f"unknown:acquisition config missing={acquisition_config_path}"
    trait_gaps = build_trait_gap_audit(
        master, snapshot, missing_basis=gap_missing_basis or "unknown:current cells unavailable"
    )

    qualified, qualified_basis = load_qualified_genus_candidates(qualified_genus_glob)
    genus_rescue = build_genus_rescue_audit(
        master,
        snapshot,
        qualified,
        genus_is_explicit=master_genus_is_explicit,
        input_basis=qualified_basis,
    )
    frozen, frozen_basis = load_frozen_packet_sources(frozen_packet_glob)
    frozen_audit = build_frozen_source_audit(master, frozen, input_basis=frozen_basis)
    return build_priority_audit(
        master,
        geography=geography,
        trait_gaps=trait_gaps,
        genus_rescue=genus_rescue,
        frozen_sources=frozen_audit,
    )


def merge_priority_audit(ledger: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    """Replace only priority audit columns while preserving campaign task state."""

    replaceable = [*PRIORITY_AUDIT_COLUMNS, *LEGACY_GEOGRAPHY_AUDIT_COLUMNS]
    old_audit = [column for column in replaceable if column in ledger.columns]
    result = ledger.drop(columns=old_audit).merge(
        audit,
        on="accepted_species",
        how="left",
        validate="one_to_one",
    )
    return result


def add_priority_sort_keys(
    frame: pd.DataFrame,
    *,
    target_traits: Iterable[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Add ascending sort keys, optionally restricted to one acquisition lane."""

    result = frame.copy()
    sort_columns: list[str] = []
    if "scheduling_region_category" in result.columns:
        result["_priority_geo"] = (
            result["scheduling_region_category"]
            .map(
                {
                    "reviewed_nh_temperate": 0,
                    "nh_temperate_unreviewed": 1,
                    "nh_tropical": 2,
                    "southern_hemisphere": 3,
                    "cross_region": 4,
                    "unclassified": 5,
                }
            )
            .fillna(5)
            .astype(int)
        )
        if "any_NH_temperate_unreviewed" in result.columns:
            result.loc[
                result["any_NH_temperate_unreviewed"].map(_text).eq("true"),
                "_priority_geo",
            ] = 1
        if "any_reviewed_NH_temperate" in result.columns:
            result.loc[
                result["any_reviewed_NH_temperate"].map(_text).eq("true"),
                "_priority_geo",
            ] = 0
        sort_columns.append("_priority_geo")
    trait_rank_column = "highest_priority_unresolved_rank"
    weighted_gap_column = "weighted_unresolved_trait_score"
    if target_traits is not None and "unresolved_priority_traits" in result.columns:
        target_set = {_text(value) for value in target_traits if _text(value)}
        ordered_targets = [trait for trait in TARGET_TRAITS if trait in target_set]
        if not ordered_targets:
            raise typer.BadParameter("target_traits contain no acquisition-priority traits")
        target_rank = {trait: TRAIT_PRIORITY_RANK[trait] for trait in ordered_targets}
        missing_by_row = result["unresolved_priority_traits"].map(
            lambda value: [
                trait
                for trait in _text(value).split("|")
                if trait in target_set
            ]
        )
        result["_target_highest_priority_unresolved_rank"] = missing_by_row.map(
            lambda traits: min((target_rank[trait] for trait in traits), default=len(TARGET_TRAITS) + 1)
        )
        result["_target_weighted_unresolved_trait_score"] = missing_by_row.map(
            lambda traits: sum(TRAIT_PRIORITY_WEIGHT[trait] for trait in traits)
        )
        trait_rank_column = "_target_highest_priority_unresolved_rank"
        weighted_gap_column = "_target_weighted_unresolved_trait_score"
    for column, key, descending in (
        (trait_rank_column, "_priority_trait_rank", False),
        (weighted_gap_column, "_priority_weighted_gap", True),
    ):
        if column not in result.columns:
            continue
        numeric = pd.to_numeric(result[column], errors="coerce")
        result[f"{key}_unknown"] = numeric.isna().astype(int)
        result[key] = -numeric.fillna(0) if descending else numeric.fillna(0)
        sort_columns.extend([f"{key}_unknown", key])
    for column, key in (
        ("unresolved_trait_cells", "_priority_unresolved"),
        ("potential_genus_rescue_cells", "_priority_genus_rescue"),
    ):
        if column not in result.columns:
            continue
        numeric = pd.to_numeric(result[column], errors="coerce")
        result[f"{key}_unknown"] = numeric.isna().astype(int)
        result[key] = -numeric.fillna(0)
        sort_columns.extend([f"{key}_unknown", key])
    if "frozen_species_direct_source" in result.columns:
        result["_priority_frozen"] = (
            result["frozen_species_direct_source"]
            .map({"true": 0, "unknown": 1, "false": 2})
            .fillna(1)
            .astype(int)
        )
        sort_columns.append("_priority_frozen")
    return result, sort_columns


def priority_summary(ledger: pd.DataFrame) -> dict[str, Any]:
    """Return JSON-safe operational counts for campaign and wave manifests."""

    summary: dict[str, Any] = {}
    categorical = (
        "any_reviewed_NH_temperate",
        "any_NH_temperate_unreviewed",
        "any_NH_tropical",
        "any_SH",
        "crosses_regions",
        "scheduling_region_category",
        "highest_priority_unresolved_trait",
        "unresolved_trait_cells_state",
        "potential_genus_rescue_state",
        "frozen_species_direct_source",
    )
    for column in categorical:
        if column in ledger.columns:
            summary[f"{column}_counts"] = {
                str(key): int(value)
                for key, value in ledger[column].fillna("unknown").value_counts().sort_index().items()
            }
    for column in (
        "unresolved_trait_cells",
        "weighted_unresolved_trait_score",
        "potential_genus_rescue_cells",
    ):
        if column in ledger.columns:
            numeric = pd.to_numeric(ledger[column], errors="coerce")
            summary[f"n_species_with_known_{column}"] = int(numeric.notna().sum())
            summary[f"sum_{column}"] = int(numeric.dropna().sum())
    return summary
