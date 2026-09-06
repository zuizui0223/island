"""Validate and materialize the frozen Wave48 reproductive evidence packet."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.wave37_europe_pmc_checkpoint import _sha256, _write_gzip_csv
from island_v2.wave45_promoted_reproductive_checkpoint import _validate_evidence

EXPECTED_SPECIES = 106_295
FORMAL_WAVE33_RUN_ID = 32_932_103_226
BASELINE_WAVE47_RUN_ID = 33_380_486_845
DIRECT_ROWS = 12
EXTERNAL_ROWS = 1
REVIEW_ROWS = 13
IDENTITY_ROWS = 10
REJECTED_ROWS = 9


def _normalise_text(value: str) -> str:
    text = " ".join(str(value).split())
    text = re.sub(r"\s+([,.;:)])", r"\1", text)
    text = re.sub(r"([(])\s+", r"\1", text)
    text = re.sub(r"\s*−\s*", "−", text)
    text = re.sub(r"μ\s+mol", "μmol", text)
    return text


def _fingerprint(value: str) -> str:
    return hashlib.sha256(_normalise_text(value).encode("utf-8")).hexdigest()


class _VisibleText(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


def _html_text(path: Path) -> str:
    parser = _VisibleText()
    parser.feed(path.read_text(encoding="utf-8"))
    return _normalise_text(" ".join(parser.parts))


def validate_packet(
    *,
    packet_dir: Path,
    target_coverage_csv: Path,
    retrieved_source_dir: Path,
    output_dir: Path,
    output_json: Path,
    expected_species: int = EXPECTED_SPECIES,
    packet_label: str = "wave48",
    baseline_formal_run_id: int = BASELINE_WAVE47_RUN_ID,
    expected_direct_rows: int | None = None,
    expected_external_rows: int | None = None,
    expected_review_rows: int | None = None,
    expected_identity_rows: int | None = None,
    expected_rejected_rows: int | None = None,
    contract: str = "wave48_multi_source_reproductive_checkpoint_v1",
    baseline_check_label: str = "immediate_wave47_baseline_pinned",
    allow_completed_direct_axis_enrichment: bool = False,
) -> dict[str, Any]:
    if not packet_label or any(token in packet_label for token in ("/", "\\")):
        raise ValueError("packet_label must be a plain filename label")
    expected_direct_rows = DIRECT_ROWS if expected_direct_rows is None else expected_direct_rows
    expected_external_rows = (
        EXTERNAL_ROWS if expected_external_rows is None else expected_external_rows
    )
    expected_review_rows = REVIEW_ROWS if expected_review_rows is None else expected_review_rows
    expected_identity_rows = (
        IDENTITY_ROWS if expected_identity_rows is None else expected_identity_rows
    )
    expected_rejected_rows = (
        REJECTED_ROWS if expected_rejected_rows is None else expected_rejected_rows
    )
    paths = {
        "manifest": packet_dir / "source_manifest.json",
        "direct": packet_dir / f"{packet_label}_reviewed_direct_evidence.csv",
        "external": packet_dir / f"{packet_label}_external_congener_evidence.csv",
        "review": packet_dir / f"{packet_label}_source_review_audit.csv",
        "identity": packet_dir / f"{packet_label}_identity_audit.csv",
        "rejected": packet_dir / f"{packet_label}_rejected_candidates.csv",
        "target_coverage": target_coverage_csv,
    }
    if missing := [str(path) for path in paths.values() if not path.is_file()]:
        raise ValueError(f"Wave48 packet is incomplete: {missing}")
    if not retrieved_source_dir.is_dir():
        raise ValueError("Wave48 retrieved-source directory is missing")

    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    target = pd.read_csv(target_coverage_csv, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "quality"}
    if missing := required_coverage.difference(target.columns):
        raise ValueError(f"Wave48 target coverage missing columns: {sorted(missing)}")
    if "trait_names" not in target.columns:
        if allow_completed_direct_axis_enrichment:
            raise ValueError(
                "Wave48 completed-axis enrichment requires baseline trait_names"
            )
        target["trait_names"] = ""
    target_species = set(target["accepted_species"])
    if len(target) != expected_species * 3 or len(target_species) != expected_species:
        raise ValueError("Wave48 target denominator mismatch")

    direct = pd.read_csv(paths["direct"], dtype=str).fillna("")
    external = pd.read_csv(paths["external"], dtype=str).fillna("")
    review = pd.read_csv(paths["review"], dtype=str).fillna("")
    identity = pd.read_csv(paths["identity"], dtype=str).fillna("")
    rejected = pd.read_csv(paths["rejected"], dtype=str).fillna("")
    _validate_evidence(direct, label="direct")
    _validate_evidence(external, label="external")
    if (
        len(direct) != expected_direct_rows
        or len(external) != expected_external_rows
        or len(review) != expected_review_rows
        or len(identity) != expected_identity_rows
        or len(rejected) != expected_rejected_rows
    ):
        raise ValueError("Wave48 frozen packet counts changed")

    if not direct["accepted_species"].isin(target_species).all():
        raise ValueError("Wave48 direct evidence is outside the fixed target")
    if external["accepted_species"].isin(target_species).any():
        raise ValueError("Wave48 external congener entered the fixed target")
    if not direct["evidence_scope"].eq("species_direct").all():
        raise ValueError("Wave48 direct evidence has the wrong scope")
    if not external["evidence_scope"].eq("external_congener_species_direct").all():
        raise ValueError("Wave48 external evidence has the wrong scope")
    if not pd.concat([direct, external])["name_match_method"].eq(
        "strict_wfo_gbif_two_backbone"
    ).all():
        raise ValueError("Wave48 evidence failed the two-backbone identity gate")
    if not pd.concat([direct, external])["axis"].eq("reproductive_assurance").all():
        raise ValueError("Wave48 evidence crossed a strict axis boundary")

    direct_axes = direct[["accepted_species", "axis"]].drop_duplicates()
    direct_status = direct_axes.merge(
        target[["accepted_species", "axis", "quality", "trait_names"]],
        on=["accepted_species", "axis"],
        how="left",
        validate="one_to_one",
    )
    completed = direct_status["quality"].ne("")
    enriched_completed_axes = 0
    if completed.any():
        if not allow_completed_direct_axis_enrichment:
            duplicate = direct_status.loc[
                completed, ["accepted_species", "axis"]
            ].to_dict("records")
            raise ValueError(f"Wave48 tries to reacquire completed cells: {duplicate}")
        packet_traits = (
            direct.groupby(["accepted_species", "axis"])["trait_name"]
            .agg(lambda values: set(map(str, values)))
            .to_dict()
        )
        invalid: list[dict[str, str]] = []
        for row in direct_status.loc[completed].itertuples(index=False):
            # A direct row may replace a Low inference for the same trait.  For
            # an already-direct Medium/High axis, only a genuinely new trait
            # may enrich the species-trait ledger.
            if str(row.quality) == "low":
                continue
            existing = {token for token in str(row.trait_names).split("|") if token}
            incoming = packet_traits[(str(row.accepted_species), str(row.axis))]
            if existing & incoming:
                invalid.append(
                    {
                        "accepted_species": str(row.accepted_species),
                        "axis": str(row.axis),
                        "duplicate_traits": "|".join(sorted(existing & incoming)),
                    }
                )
        if invalid:
            raise ValueError(
                "Wave48 tries to re-ingest completed direct species-traits: "
                f"{invalid}"
            )
        enriched_completed_axes = int(completed.sum())

    evidence = pd.concat([direct, external], ignore_index=True)
    if set(review["record_id"]) != set(evidence["source_record_id"]):
        raise ValueError("Wave48 source review does not cover every evidence row")
    if not review["accepted_correct"].str.casefold().eq("true").all():
        raise ValueError("Wave48 review contains an unaccepted evidence row")
    paired = evidence.merge(
        review[
            [
                "record_id",
                "exact_supporting_quote",
                "source_lineage",
                "content_fingerprint",
            ]
        ],
        left_on="source_record_id",
        right_on="record_id",
        suffixes=("_evidence", "_review"),
        validate="one_to_one",
    )
    if not paired["source_excerpt"].eq(paired["exact_supporting_quote"]).all():
        raise ValueError("Wave48 evidence and review quote differ")
    if not paired["source_lineage_evidence"].eq(paired["source_lineage_review"]).all():
        raise ValueError("Wave48 evidence and review lineage differ")
    if not paired["content_fingerprint_evidence"].eq(
        paired["content_fingerprint_review"]
    ).all():
        raise ValueError("Wave48 evidence and review content fingerprint differ")
    calculated_fingerprints = paired["exact_supporting_quote"].map(_fingerprint)
    if not calculated_fingerprints.eq(paired["content_fingerprint_evidence"]).all():
        raise ValueError("Wave48 content fingerprint does not match its exact quote")

    expected_identity = set(evidence["accepted_species"])
    if set(identity["accepted_species"]) != expected_identity:
        raise ValueError("Wave48 identity audit does not cover every new species")
    fixed_identity = identity["target_universe_status"].eq("fixed_target")
    external_identity = identity["target_universe_status"].eq("external_congener_only")
    identity_ok = (
        identity["name_match_method"].eq("strict_wfo_gbif_two_backbone").all()
        and identity["wfo_match_id"].str.startswith("wfo-").all()
        and identity["wfo_status"].str.casefold().eq("accepted").all()
        and identity["wfo_rank"].str.casefold().eq("species").all()
        and identity["gbif_usage_key"].ne("").all()
        and identity["gbif_status"].eq("ACCEPTED").all()
        and identity["gbif_rank"].eq("SPECIES").all()
        and identity["gbif_match_type"].eq("EXACT").all()
        and identity["gbif_confidence"].astype(int).ge(95).all()
        and identity["family"].eq(identity["gbif_family"]).all()
        and identity.loc[fixed_identity, "accepted_species"].isin(target_species).all()
        and not identity.loc[external_identity, "accepted_species"].isin(target_species).any()
        and int(fixed_identity.sum()) == int(direct["accepted_species"].nunique())
        and int(external_identity.sum()) == int(external["accepted_species"].nunique())
    )
    if not identity_ok:
        raise ValueError("Wave48 strict WFO/GBIF identity gate failed")

    sources = manifest["sources"]
    if len(sources) != int(manifest["selection"]["source_pages_retrieved"]) or sum(
        int(source["reviewed_rows"]) for source in sources
    ) != len(evidence):
        raise ValueError("Wave48 source manifest row accounting changed")
    source_receipts: dict[str, dict[str, Any]] = {}
    manual_receipt_sources: list[str] = []
    for source in sources:
        path = retrieved_source_dir / source["retrieved_filename"]
        if source.get("retrieval_mode") == "manual_browser_pinned":
            manual = source.get("manual_receipt", {})
            source_rows = evidence.loc[
                evidence["source_lineage"].eq(source["source_lineage"])
            ]
            quote_fingerprints = sorted(set(source_rows["content_fingerprint"]))
            declared_fingerprints = sorted(
                str(value) for value in manual.get("reviewed_quote_fingerprints", [])
            )
            digest = str(source.get("retrieved_content_sha256", ""))
            byte_count = int(manual.get("retrieved_bytes", 0))
            page_numbers = manual.get("reviewed_page_numbers", [])
            if (
                not re.fullmatch(r"[0-9a-f]{64}", digest)
                or byte_count <= 0
                or declared_fingerprints != quote_fingerprints
                or not page_numbers
                or not all(isinstance(value, int) and value > 0 for value in page_numbers)
                or not str(manual.get("retrieved_at_utc", "")).endswith("Z")
                or not str(manual.get("official_catalog_url", "")).startswith("https://")
                or manual.get("binary_available_in_ci") is not False
                or path.exists()
            ):
                raise ValueError(
                    f"Wave48 manual browser receipt is incomplete: {source['source_id']}"
                )
            source_receipts[source["source_id"]] = {
                "bytes": byte_count,
                "sha256": digest,
                "retrieval_mode": "manual_browser_pinned",
                "source_binary_present_in_run": False,
                "reviewed_page_numbers": page_numbers,
                "reviewed_quote_fingerprints": declared_fingerprints,
            }
            manual_receipt_sources.append(source["source_id"])
            continue
        if not path.is_file():
            raise ValueError(f"Wave48 retrieved source missing: {path}")
        receipt: dict[str, Any] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
        if expected_sha := source.get("retrieved_content_sha256"):
            if receipt["sha256"] != expected_sha:
                raise ValueError(f"Wave48 source digest changed: {source['source_id']}")
        else:
            source_rows = evidence.loc[
                evidence["source_lineage"].eq(source["source_lineage"])
            ]
            anchors = sorted(set(source_rows["source_excerpt"]))
            if len(anchors) != 1:
                raise ValueError("Wave48 mutable source must have one frozen content anchor")
            if _fingerprint(anchors[0]) != source["content_anchor_sha256"]:
                raise ValueError("Wave48 source anchor fingerprint changed")
            visible = _html_text(path)
            if _normalise_text(anchors[0]) not in visible:
                raise ValueError("Wave48 source page no longer contains the frozen quote")
            receipt["content_anchor_sha256"] = source["content_anchor_sha256"]
        source_receipts[source["source_id"]] = receipt

    policy = manifest["inference_policy"]
    declared = manifest["evidence_counts"]
    if (
        manifest["fixed_target_species"] != expected_species
        or manifest["formal_wave33_baseline"]["run_id"] != FORMAL_WAVE33_RUN_ID
        or manifest["immediate_formal_baseline"]["run_id"] != baseline_formal_run_id
        or declared["direct_rows"] != len(direct)
        or declared["external_rows"] != len(external)
        or declared["review_rows"] != len(review)
        or declared["identity_rows"] != len(identity)
        or declared["rejected_rows"] != len(rejected)
        or policy["join_key"] != "genus x axis x trait_name"
        or policy["minimum_species"] != 3
        or policy["family_inference"]
        or policy["global_fallback"]
        or policy["reproductive_traits_interchangeable"]
        or policy["single_source_lineage_can_unlock_rule"]
    ):
        raise ValueError("Wave48 source or inference contract changed")

    output_dir.mkdir(parents=True, exist_ok=True)
    direct_path = output_dir / f"{packet_label}_reviewed_direct_evidence.csv.gz"
    external_path = output_dir / f"{packet_label}_external_congener_evidence.csv.gz"
    rejected_path = output_dir / f"{packet_label}_rejected_candidates.csv.gz"
    _write_gzip_csv(direct, direct_path)
    _write_gzip_csv(external, external_path)
    _write_gzip_csv(rejected, rejected_path)

    summary: dict[str, Any] = {
        "contract": contract,
        "fixed_target_species": expected_species,
        "formal_wave33_run_id": FORMAL_WAVE33_RUN_ID,
        "baseline_formal_run_id": baseline_formal_run_id,
        "evidence": {
            "new_direct_rows": len(direct),
            "new_direct_species": int(direct["accepted_species"].nunique()),
            "new_direct_species_trait": int(
                direct[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "new_direct_species_axis": int(direct_axes.shape[0]),
            "new_direct_unresolved_species_axis": int((~completed).sum()),
            "completed_direct_axis_enrichments": enriched_completed_axes,
            "new_external_rows": len(external),
            "new_external_species_trait": int(
                external[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
            ),
            "review_rows": len(review),
            "identity_rows": len(identity),
            "rejected_rows": len(rejected),
        },
        "queries": manifest.get(
            "query_accounting",
            {
                "formal_search_api_queries": 0,
                "source_pages_retrieved": len(sources),
                "query_cost_usd": 0,
            },
        ),
        "query_cost_usd": float(
            manifest.get("query_accounting", {}).get("query_cost_usd", 0)
        ),
        "source_receipts": source_receipts,
        "manual_browser_receipt_sources": manual_receipt_sources,
        "source_binary_coverage": {
            "verified_in_run": len(sources) - len(manual_receipt_sources),
            "manual_receipt_only": len(manual_receipt_sources),
        },
        "limitations": {
            "manual_source_binary_not_in_artifact": manual_receipt_sources,
        },
        "checks": {
            "fixed_denominator": True,
            "formal_wave33_baseline_pinned": True,
            baseline_check_label: True,
            "all_new_rows_reviewed": True,
            "retrieved_sources_verified": True,
            "manual_browser_receipts_verified": True,
            "exact_quote_and_provenance_complete": True,
            "content_fingerprints_verified": True,
            "direct_cells_or_novel_axis_traits_only": True,
            "completed_axis_enrichment_trait_specific": bool(
                not enriched_completed_axes
                or allow_completed_direct_axis_enrichment
            ),
            "direct_inside_fixed_target": True,
            "external_outside_fixed_target": True,
            "strict_two_backbone_identity": True,
            "trait_specific_genus_join": True,
            "reproductive_traits_not_interchanged": True,
            "single_lineage_cannot_unlock_rule": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {label: _sha256(path) for label, path in paths.items()},
        "artifact_sha256": {
            direct_path.name: _sha256(direct_path),
            external_path.name: _sha256(external_path),
            rejected_path.name: _sha256(rejected_path),
        },
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet-dir", required=True, type=Path)
    parser.add_argument("--target-coverage-csv", required=True, type=Path)
    parser.add_argument("--retrieved-source-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=EXPECTED_SPECIES)
    args = parser.parse_args()
    summary = validate_packet(
        packet_dir=args.packet_dir,
        target_coverage_csv=args.target_coverage_csv,
        retrieved_source_dir=args.retrieved_source_dir,
        output_dir=args.output_dir,
        output_json=args.output_json,
        expected_species=args.expected_species,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
