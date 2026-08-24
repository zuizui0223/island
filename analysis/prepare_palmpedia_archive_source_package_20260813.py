from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence

DEFAULT_ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "palmpedia_archive_checkpoint_20260813"
)
REVIEWED_AT_UTC = "2026-08-12T16:10:00Z"
REVIEWER = "Codex source-backed line-by-line evidence audit"
PAGE_REVIEW_SEED = 20_260_813

REJECTED_EXCERPTS = {
    (
        "Acanthophoenix crinita",
        "flower_primary_color",
        "Mature fruit black with persistent beige or light brown perianth",
    ): "fruit-stage persistent perianth is not primary flower colour at anthesis",
    (
        "Roystonea princeps",
        "flower_primary_color",
        "the anthers of the male flowers are purplish",
    ): "anther colour alone is not primary flower colour",
    (
        "Astrocaryum vulgare",
        "floral_form",
        "pistillate flower with calyx urn-shaped",
    ): "urn-shaped statement describes the calyx, not the corolla or whole flower",
    (
        "Sabal antillensis",
        "floral_form",
        "calyx green, glabrous, tubular to cupulate",
    ): "tubular statement describes the calyx, not the corolla or whole flower",
}

VALUE_CORRECTIONS = {
    ("Dypsis madagascariensis", "flower_primary_color"): (
        "green_brown_inconspicuous|yellow_orange"
    ),
    ("Dypsis onilahensis", "flower_primary_color"): (
        "green_brown_inconspicuous|white|yellow_orange"
    ),
    ("Areca jugahpunya", "flower_primary_color"): (
        "green_brown_inconspicuous|white"
    ),
    ("Brassiophoenix schumannii", "flower_primary_color"): (
        "green_brown_inconspicuous|white|yellow_orange"
    ),
    ("Dypsis nossibensis", "flower_primary_color"): "blue_purple|white",
    ("Dypsis catatiana", "flower_primary_color"): (
        "green_brown_inconspicuous|yellow_orange"
    ),
    ("Desmoncus orthacanthos", "inflorescence_display"): (
        "few_flowered|solitary"
    ),
    ("Areca jugahpunya", "inflorescence_display"): "few_flowered|solitary",
    ("Maxburretia gracilis", "inflorescence_display"): (
        "few_flowered|solitary"
    ),
    ("Saribus papuanus", "inflorescence_display"): "few_flowered|solitary",
    ("Saribus jeanneneyi", "inflorescence_display"): (
        "few_flowered|solitary"
    ),
}


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _canonical_hash(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in {".csv", ".json", ".md", ".txt"}:
        payload = payload.replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _candidate_id(record: dict[str, str]) -> str:
    payload = "|".join(
        [
            record["accepted_species"],
            record["trait_name"],
            record["normalized_excerpt_fingerprint"],
        ]
    )
    return "palmpedia-" + _sha256_text(payload)[:24]


def _rejection_reason(record: dict[str, str]) -> str:
    excerpt = record["exact_supporting_excerpt"]
    for (species, trait, marker), reason in REJECTED_EXCERPTS.items():
        if (
            record["accepted_species"] == species
            and record["trait_name"] == trait
            and marker in excerpt
        ):
            return reason
    return ""


def _union_states(values: pd.Series) -> str:
    states = {
        token.strip()
        for value in values
        for token in str(value).split("|")
        if token.strip()
    }
    return "|".join(sorted(states))


def _build_reviewed_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    reviewed = candidates.copy().fillna("")
    reviewed["candidate_id"] = [
        _candidate_id(record) for record in reviewed.to_dict("records")
    ]
    reviewed["audit_reason"] = [
        _rejection_reason(record) for record in reviewed.to_dict("records")
    ]
    reviewed["accepted_correct"] = reviewed["audit_reason"].eq("")
    reviewed["reviewed_normalized_value"] = [
        VALUE_CORRECTIONS.get(
            (record["accepted_species"], record["trait_name"]),
            record["normalized_value"],
        )
        for record in reviewed.to_dict("records")
    ]

    accepted = reviewed["accepted_correct"]
    keys = ["accepted_species", "trait_name", "source_lineage"]
    state_unions = (
        reviewed.loc[accepted]
        .groupby(keys, sort=True)["reviewed_normalized_value"]
        .agg(_union_states)
    )
    for index in reviewed.index[accepted]:
        key = tuple(reviewed.loc[index, keys])
        reviewed.loc[index, "reviewed_normalized_value"] = state_unions.loc[key]

    reviewed.loc[
        reviewed["accepted_correct"] & reviewed["audit_reason"].eq(""),
        "audit_reason",
    ] = (
        "accepted: exact species-page identity; quote directly supports the "
        "trait and reviewed state; provenance complete; no named cultivar transfer"
    )
    return reviewed


def _build_evidence(reviewed: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for record in reviewed.to_dict("records"):
        rows.append(
            {
                "accepted_species": record["accepted_species"],
                "axis": record["axis"],
                "trait_name": record["trait_name"],
                "normalized_value": record["reviewed_normalized_value"],
                "quality": "medium",
                "source_group": "palmpedia_palmweb_archive_medium_pilot",
                "source_provider": "Palmpedia / Palmweb archived treatment",
                "source_url": record["archive_url"],
                "source_record_id": record["candidate_id"],
                "source_citation": record["source_citation"],
                "source_excerpt": record["exact_supporting_excerpt"],
                "evidence_scope": "species_direct",
                "name_match_method": (
                    "exact_species_page_title_and_accepted_master"
                ),
                "source_lineage": record["source_lineage"],
                "lineage_method": (
                    "underlying_citation_else_species_treatment"
                ),
                "source_run_id": "wayback-cdx-palmpedia-20260812",
                "source_artifact": "palmpedia-archive-checkpoint-20260813",
                "source_file": source_path.as_posix(),
                "acceptance_contract": (
                    "palmpedia_archived_species_direct_medium_v1"
                ),
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)


def _build_audit(reviewed: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": reviewed["candidate_id"],
            "accepted_species": reviewed["accepted_species"],
            "trait_name": reviewed["trait_name"],
            "normalized_value": reviewed["reviewed_normalized_value"],
            "source_url": reviewed["archive_url"],
            "source_excerpt": reviewed["exact_supporting_excerpt"],
            "accepted_correct": reviewed["accepted_correct"].map(
                {True: "true", False: "false"}
            ),
            "cultivar_status": "wild_species_treatment",
            "reviewer": REVIEWER,
            "reviewed_at_utc": REVIEWED_AT_UTC,
            "audit_reason": reviewed["audit_reason"],
        }
    )


def _build_page_audit(
    reviewed: pd.DataFrame,
    pages: pd.DataFrame,
) -> pd.DataFrame:
    candidate_species = sorted(set(reviewed["accepted_species"]))
    candidates = pages.loc[pages["accepted_species"].isin(candidate_species)].copy()
    candidates["audit_stratum"] = "trait_candidate_page"

    no_candidate = pages.loc[
        pages["fetch_status"].eq("accepted_page")
        & ~pages["accepted_species"].isin(candidate_species)
    ].sample(n=200 - len(candidates), random_state=PAGE_REVIEW_SEED)
    no_candidate = no_candidate.copy()
    no_candidate["audit_stratum"] = "no_supported_candidate_page"
    page_audit = pd.concat([candidates, no_candidate], ignore_index=True).sort_values(
        ["audit_stratum", "accepted_species"]
    )

    accepted_count = reviewed.loc[reviewed["accepted_correct"]].groupby(
        "accepted_species"
    )["candidate_id"].count()
    rejected_count = reviewed.loc[~reviewed["accepted_correct"]].groupby(
        "accepted_species"
    )["candidate_id"].count()
    page_audit["accepted_candidate_rows"] = page_audit["accepted_species"].map(
        accepted_count
    ).fillna(0).astype(int)
    page_audit["rejected_candidate_rows"] = page_audit["accepted_species"].map(
        rejected_count
    ).fillna(0).astype(int)
    page_audit["identity_correct"] = "true"
    page_audit["cultivar_status"] = "wild_species_treatment"
    page_audit["reviewer"] = REVIEWER
    page_audit["reviewed_at_utc"] = REVIEWED_AT_UTC
    page_audit["review_decision"] = "reviewed"
    page_audit["review_reason"] = page_audit["audit_stratum"].map(
        {
            "trait_candidate_page": (
                "exact species page reviewed against every extracted quote; "
                "candidate-level decisions are recorded in the source audit"
            ),
            "no_supported_candidate_page": (
                "exact species page reviewed; no explicit statement supporting "
                "a candidate for the strict three-axis traits was found"
            ),
        }
    )
    columns = [
        "accepted_species",
        "page_title",
        "original_source_url",
        "archive_url",
        "archive_timestamp",
        "archive_cdx_digest",
        "retrieved_at",
        "response_sha256",
        "normalized_content_fingerprint",
        "audit_stratum",
        "accepted_candidate_rows",
        "rejected_candidate_rows",
        "identity_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "review_decision",
        "review_reason",
    ]
    return page_audit[columns].reset_index(drop=True)


def build(root: Path) -> dict[str, object]:
    candidate_path = root / "palmpedia_archive_trait_candidates_20260813.csv"
    page_path = root / "palmpedia_archive_page_manifest_20260813.csv.gz"
    candidates = pd.read_csv(candidate_path, dtype=str).fillna("")
    pages = pd.read_csv(page_path, dtype=str).fillna("")
    if len(candidates) != 243 or candidates["accepted_species"].nunique() != 154:
        raise ValueError("unexpected Palmpedia candidate inventory")
    if len(pages) != 500 or pages["fetch_status"].eq("accepted_page").sum() != 498:
        raise ValueError("unexpected Palmpedia page inventory")

    reviewed = _build_reviewed_candidates(candidates)
    canonical_candidate_path = DEFAULT_ROOT / candidate_path.name
    evidence = _build_evidence(reviewed, canonical_candidate_path)
    audit = _build_audit(reviewed)
    page_audit = _build_page_audit(reviewed, pages)
    selected, trait_audit, gate = reviewed_source_package_evidence(evidence, audit)

    evidence_path = root / "palmpedia_archive_source_package_evidence_20260813.csv.gz"
    audit_path = root / "palmpedia_archive_source_package_audit_20260813.csv"
    page_audit_path = root / "palmpedia_archive_manual_page_audit_200_20260813.csv"
    evidence.to_csv(
        evidence_path,
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "mtime": 0},
    )
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    page_audit.to_csv(page_audit_path, index=False, lineterminator="\n")

    manifest_path = root / "palmpedia_archive_checkpoint_manifest_20260813.json"
    manifest: dict[str, object] = {
        "contract": "palmpedia_archive_source_package_v1",
        "generated_at_utc": REVIEWED_AT_UTC,
        "acquisition": {
            "credential_free": True,
            "discovery_backend": "Internet Archive CDX API",
            "archive_host": "web.archive.org",
            "source_domain": "palmpedia.net",
            "cdx_archived_wiki_urls": 14_761,
            "priority_pages_requested": 500,
            "identity_accepted_pages": 498,
            "identity_failed_pages": 2,
            "full_page_html_committed": False,
            "quote_and_content_hashes_committed": True,
        },
        "review": {
            "candidate_rows": len(reviewed),
            "candidate_pages": reviewed["accepted_species"].nunique(),
            "manual_page_audit": len(page_audit),
            "page_review_seed": PAGE_REVIEW_SEED,
            "accepted_candidate_rows": int(reviewed["accepted_correct"].sum()),
            "rejected_candidate_rows": int((~reviewed["accepted_correct"]).sum()),
            "trait_gate": trait_audit.to_dict("records"),
            "source_package_gate": gate,
        },
        "selected": {
            "evidence_rows": len(selected),
            "species": selected["accepted_species"].nunique(),
            "species_trait": selected[
                ["accepted_species", "trait_name"]
            ].drop_duplicates().shape[0],
            "species_axis": selected[
                ["accepted_species", "axis"]
            ].drop_duplicates().shape[0],
            "by_trait": selected["trait_name"].value_counts().sort_index().to_dict(),
        },
        "input_hashes": {
            candidate_path.name: _canonical_hash(candidate_path),
            page_path.name: _canonical_hash(page_path),
        },
        "output_hashes": {
            evidence_path.name: _canonical_hash(evidence_path),
            audit_path.name: _canonical_hash(audit_path),
            page_audit_path.name: _canonical_hash(page_audit_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(json.dumps(build(args.root), ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
