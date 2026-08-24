"""Acquire visible floral traits from the SANBI e-Flora of South Africa.

The source is one versioned Darwin Core Archive (DwC-A), rather than a set of
per-species searches.  The archive contains species names, exact morphology
descriptions, and the underlying bibliographic citation for every description.
This adapter retains all three and reuses the shared six-trait extractor.

Promotion is deliberately separate.  The command writes an incremental common-
schema evidence package and a deterministic audit queue; a completed review must
pass the repository's source-package gate before the rows enter formal coverage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    COLOUR_PATTERNS,
    _cultivar_or_hybrid_treatment,
    extract_description,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import load_rules

ARCHIVE_VERSION = "1.42"
ARCHIVE_URL = (
    "https://ipt.sanbi.org.za/archive.do?"
    f"r=flora_descriptions&v={ARCHIVE_VERSION}"
)
LANDING_URL = "https://ipt.sanbi.org.za/resource?r=flora_descriptions"
EXPECTED_ARCHIVE_SHA256 = (
    "157fe16ca6a5698df24fca4796aec2d69a8da20fd38a3bd36633eb0d5cdedfe4"
)
SOURCE_RUN_ID = "sanbi-eflora-south-africa-v1.42-20260809"
SOURCE_ARTIFACT = "sanbi-eflora-south-africa-dwca-v1.42"
SOURCE_FILE = "description.txt"
SOURCE_PROVIDER = "sanbi_eflora_south_africa_v1_42"
ACCEPTANCE_CONTRACT = "official_dwca_flora_exact_species_description_quote_v1"
USER_AGENT = "island-trait-research/1.0 (github.com/zuizui0223/island)"
AUDIT_SEED = "sanbi-eflora-south-africa-audit-20260809"
AUDIT_QUOTA = {
    "flower_primary_color": 70,
    "flower_size_class": 45,
    "inflorescence_display": 35,
    "floral_form": 20,
    "tube_depth_class": 20,
    "floral_symmetry": 10,
}

_ANCILLARY_COLOUR_SUBJECT = (
    r"(?:buds?|bracts?|bracteoles?|involucr\w*|sepals?|calyx|stamens?|"
    r"staminodes?|anthers?|pollen|margins?)"
)
_COLOUR_WORD = (
    r"(?:white|whitish|cream|creamy|ivory|yellow|yellowish|orange|orangish|"
    r"gold|golden|ochre|red|reddish|pink|pinkish|rose|roseate|scarlet|"
    r"crimson|magenta|maroon|salmon|coral|blue|bluish|purple|purplish|"
    r"violet|violaceous|lilac|mauve|lavender|green|greenish|brown|brownish|"
    r"grey|greyish|gray|grayish|fuscous|ferruginous|ferrugineous|black|"
    r"blackish)"
)
ANCILLARY_COLOUR_RE = re.compile(
    rf"\b{_ANCILLARY_COLOUR_SUBJECT}\b[^.;:]{{0,45}}\b{_COLOUR_WORD}\b",
    re.IGNORECASE,
)
COLOUR_RISK_RE = re.compile(
    r"\b(?:buds?|bracts?|bracteoles?|involucr\w*|sepals?|calyx|stamens?|"
    r"staminodes?|anthers?|pollen|nectar(?:y|ies)?|spikelets?|glumes?|"
    r"leaves|leaf|fruits?|stems?|hairs?|ovary|peduncles?|margins?)\b|"
    r"\b(?:probably|possibly|perhaps|uncertain)\b|\?",
    re.IGNORECASE,
)

AUDIT_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "provider",
    "source_url",
    "candidate_id",
    "source_citation",
    "exact_supporting_quote",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "source_file",
    "acceptance_contract",
    "accepted_correct",
    "cultivar_status",
    "reviewer",
    "reviewed_at_utc",
    "audit_reason",
]


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_species_name(scientific_name: str) -> str:
    """Return a species binomial, rejecting hybrids and infraspecific rows."""

    value = _text(scientific_name)
    if "×" in value or re.search(
        r"\s(?:subsp\.?|ssp\.?|var\.?|f\.?|forma)\s+[a-z]", value
    ):
        return ""
    match = re.match(r"^([A-Z][A-Za-z.-]+)\s+([a-z][A-Za-z.-]+)(?:\s|$)", value)
    if not match:
        return ""
    return f"{match.group(1)} {match.group(2)}"


def download_archive(cache_path: Path) -> tuple[Path, str]:
    """Download the fixed archive once and fail if its bytes drift."""

    status = "cached"
    if not cache_path.exists():
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = cache_path.with_suffix(cache_path.suffix + ".partial")
        with requests.get(
            ARCHIVE_URL,
            timeout=600,
            stream=True,
            headers={"User-Agent": USER_AGENT},
        ) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for block in response.iter_content(chunk_size=1024 * 1024):
                    if block:
                        handle.write(block)
        temporary.replace(cache_path)
        status = "fetched"
    observed = _file_sha256(cache_path)
    if observed != EXPECTED_ARCHIVE_SHA256:
        raise ValueError(
            "SANBI archive SHA-256 drifted: "
            f"expected {EXPECTED_ARCHIVE_SHA256}, observed {observed}"
        )
    return cache_path, status


def extract_archive(archive_path: Path, output_dir: Path) -> Path:
    """Extract only the four DwC-A tables used by this adapter."""

    required = {"taxon.txt", "description.txt", "reference.txt", "eml.xml"}
    output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as archive:
        missing = required.difference(archive.namelist())
        if missing:
            raise ValueError(f"SANBI DwC-A is missing files: {sorted(missing)}")
        for name in sorted(required):
            destination = output_dir / name
            if not destination.exists():
                with archive.open(name) as source, destination.open("wb") as target:
                    shutil.copyfileobj(source, target)
    return output_dir


def _citation_lineage(citation: str, wfo_id: str, reference_id: str) -> tuple[str, str]:
    normalized = re.sub(
        r"\s+", " ", re.sub(r"\[CC BY\]", "", citation, flags=re.IGNORECASE)
    ).strip().casefold()
    if normalized:
        return f"citation:{_sha(normalized)[:32]}", "normalized_original_citation"
    return (
        f"sanbi-eflora-treatment:{wfo_id}:{reference_id}",
        "wfo_taxon_and_reference_id_fallback",
    )


def _exact_quote_segments(quote: str, description: str) -> str:
    """Keep a contiguous quote or explicitly delimit exact sentence segments."""

    normalized_quote = _text(quote)
    normalized_description = _text(description)
    if normalized_quote in normalized_description:
        return normalized_quote
    segments = [
        segment.strip()
        for segment in re.split(r"(?<=[.!?])\s+(?=[A-Z])", normalized_quote)
        if segment.strip()
    ]
    if len(segments) > 1 and all(segment in normalized_description for segment in segments):
        return " || ".join(segments)
    return ""


def has_ancillary_colour_conflict(quote: str) -> bool:
    """Reject colour rows whose extracted states include a non-target organ."""

    return bool(ANCILLARY_COLOUR_RE.search(_text(quote)))


def is_clean_primary_colour(quote: str, normalized_value: str) -> bool:
    """Require complete target-organ colour mapping with no risky context."""

    text = _text(quote)
    if COLOUR_RISK_RE.search(text) or has_ancillary_colour_conflict(text):
        return False
    detected = {
        state
        for pattern, state in COLOUR_PATTERNS
        if re.search(pattern, text, re.IGNORECASE)
    }
    expected = {value for value in _text(normalized_value).split("|") if value}
    return bool(detected) and detected == expected


def _build_evidence_row(
    *,
    species: str,
    wfo_id: str,
    reference_id: str,
    citation: str,
    trait: str,
    value: str,
    quote: str,
) -> dict[str, str]:
    lineage, lineage_method = _citation_lineage(citation, wfo_id, reference_id)
    record_hash = _sha(f"{wfo_id}|{reference_id}|{trait}|{value}|{quote}")[:24]
    return {
        "accepted_species": species,
        "axis": AXIS[trait],
        "trait_name": trait,
        "normalized_value": value,
        "quality": "high",
        "source_group": "latest_public_web",
        "source_provider": SOURCE_PROVIDER,
        "source_url": LANDING_URL,
        "source_record_id": (
            f"sanbi:{wfo_id}:{reference_id}:{trait}:{record_hash[:12]}"
        ),
        "source_citation": citation or f"e-Flora of South Africa v{ARCHIVE_VERSION}",
        "source_excerpt": quote,
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact_wfo_species_record",
        "source_lineage": lineage,
        "lineage_method": lineage_method,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_file": SOURCE_FILE,
        "acceptance_contract": ACCEPTANCE_CONTRACT,
    }


def read_dwca(dwca_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    taxon = pd.read_csv(dwca_dir / "taxon.txt", sep="\t", dtype=str).fillna("")
    descriptions = pd.read_csv(
        dwca_dir / "description.txt", sep="\t", dtype=str
    ).fillna("")
    references = pd.read_csv(
        dwca_dir / "reference.txt", sep="\t", dtype=str
    ).fillna("")
    return taxon, descriptions, references


def build_evidence(
    *,
    taxon: pd.DataFrame,
    descriptions: pd.DataFrame,
    references: pd.DataFrame,
    accepted_species: set[str],
    current_direct_pairs: set[tuple[str, str]],
    rules: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Extract incremental exact-species evidence with original citations."""

    names = taxon[["id", "scientificName"]].copy()
    names["accepted_species"] = names["scientificName"].map(canonical_species_name)
    names = names.loc[names["accepted_species"].isin(accepted_species)].copy()
    names = names.drop_duplicates("id")
    name_by_id = dict(zip(names["id"], names["accepted_species"], strict=True))

    morphology = descriptions.loc[
        descriptions["id"].isin(name_by_id) & descriptions["type"].eq("Morphology")
    ].copy()
    citation_table = references[["id", "identifier", "bibliographicCitation"]].drop_duplicates(
        ["id", "identifier"]
    )
    citation_by_key = {
        (row.id, row.identifier): _text(row.bibliographicCitation)
        for row in citation_table.itertuples(index=False)
    }

    rows: list[dict[str, str]] = []
    for record in morphology.sort_values(["id", "source"], kind="stable").itertuples(
        index=False
    ):
        wfo_id = str(record.id)
        reference_id = str(record.source)
        species = name_by_id[wfo_id]
        description = _text(record.description)
        if not description or _cultivar_or_hybrid_treatment(species, description):
            continue
        citation = citation_by_key.get((wfo_id, reference_id), "")
        extracted = extract_description(description, rules)
        for trait, item in extracted.items():
            if (species, trait) in current_direct_pairs:
                continue
            quote = _exact_quote_segments(_text(item["quote"]), description)
            if not quote:
                continue
            value = _text(item["value"])
            if trait == "flower_primary_color" and not is_clean_primary_colour(
                quote, value
            ):
                continue
            rows.append(
                _build_evidence_row(
                    species=species,
                    wfo_id=str(wfo_id),
                    reference_id=str(reference_id),
                    citation=citation,
                    trait=trait,
                    value=value,
                    quote=quote,
                )
            )

    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    if not evidence.empty:
        evidence = evidence.drop_duplicates("source_record_id").sort_values(
            ["trait_name", "accepted_species", "source_lineage", "source_record_id"],
            kind="stable",
        )
    inventory = {
        "source_taxon_rows": len(taxon),
        "source_description_rows": len(descriptions),
        "source_morphology_rows": int(descriptions["type"].eq("Morphology").sum()),
        "source_reference_rows": len(references),
        "exact_fixed_species": int(names["accepted_species"].nunique()),
        "exact_fixed_taxon_ids": len(name_by_id),
        "exact_fixed_morphology_taxon_ids": int(morphology["id"].nunique()),
    }
    return evidence.reset_index(drop=True), inventory


def deterministic_audit_sample(evidence: pd.DataFrame) -> pd.DataFrame:
    samples: list[pd.DataFrame] = []
    for trait, quota in AUDIT_QUOTA.items():
        candidates = evidence.loc[evidence["trait_name"].eq(trait)].copy()
        if len(candidates) < quota:
            raise ValueError(f"not enough {trait} evidence for audit: {len(candidates)} < {quota}")
        seed = int(_sha(f"{AUDIT_SEED}|{trait}")[:8], 16)
        samples.append(candidates.sample(n=quota, random_state=seed))
    audit = pd.concat(samples, ignore_index=True).sort_values(
        ["trait_name", "accepted_species", "source_record_id"], kind="stable"
    )
    audit = audit.rename(
        columns={
            "source_provider": "provider",
            "source_record_id": "candidate_id",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = ""
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""
    return audit.reindex(columns=AUDIT_COLUMNS)


def acquire(
    *,
    archive_path: Path,
    dwca_dir: Path,
    master_csv: Path,
    strict_coverage_csv: Path,
    current_direct_ledger_csv: Path,
    measurement_config: Path,
    output_dir: Path,
) -> dict[str, Any]:
    observed_hash = _file_sha256(archive_path)
    if observed_hash != EXPECTED_ARCHIVE_SHA256:
        raise ValueError(
            "SANBI archive SHA-256 drifted: "
            f"expected {EXPECTED_ARCHIVE_SHA256}, observed {observed_hash}"
        )
    extract_archive(archive_path, dwca_dir)
    taxon, descriptions, references = read_dwca(dwca_dir)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    strict = pd.read_csv(
        strict_coverage_csv, dtype=str, usecols=["accepted_species"]
    ).fillna("")
    strict_names = set(strict["accepted_species"])
    if len(strict_names) != 106_295:
        raise ValueError("strict coverage must contain exactly 106,295 accepted species")
    master = master.loc[master["accepted_species"].isin(strict_names)].copy()
    if len(master) != 106_295 or master["accepted_species"].nunique() != 106_295:
        raise ValueError("fixed master does not match the 106,295-species strict universe")
    direct = pd.read_csv(current_direct_ledger_csv, dtype=str).fillna("")
    if "resolution_status" in direct:
        direct = direct.loc[direct["resolution_status"].eq("resolved")]
    current_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=True))
    evidence, inventory = build_evidence(
        taxon=taxon,
        descriptions=descriptions,
        references=references,
        accepted_species=set(master["accepted_species"]),
        current_direct_pairs=current_pairs,
        rules=load_rules(measurement_config),
    )
    audit = deterministic_audit_sample(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "sanbi_eflora_south_africa_evidence.csv.gz"
    audit_path = output_dir / "sanbi_eflora_south_africa_audit_queue_200.csv"
    evidence.to_csv(evidence_path, index=False, compression="gzip")
    audit.to_csv(audit_path, index=False)
    summary: dict[str, Any] = {
        "contract": "sanbi_eflora_south_africa_source_package_v1",
        "archive_version": ARCHIVE_VERSION,
        "archive_url": ARCHIVE_URL,
        "landing_url": LANDING_URL,
        "archive_sha256": observed_hash,
        "archive_bytes": archive_path.stat().st_size,
        "fixed_universe_species": 106_295,
        **inventory,
        "incremental_evidence_rows": len(evidence),
        "incremental_species_trait": int(
            evidence.drop_duplicates(["accepted_species", "trait_name"]).shape[0]
        ),
        "incremental_species_axis": int(
            evidence.drop_duplicates(["accepted_species", "axis"]).shape[0]
        ),
        "incremental_species": int(evidence["accepted_species"].nunique()),
        "trait_rows": {
            str(key): int(value)
            for key, value in evidence["trait_name"].value_counts().sort_index().items()
        },
        "audit_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "audit_quota": AUDIT_QUOTA,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "quality": "high",
        "promotion_performed": False,
        "strict_exclusions": [
            "infraspecific_and_hybrid_names",
            "non_exact_fixed_master_names",
            "non_morphology_description_types",
            "existing_direct_species_trait_pairs",
            "pollen_vector_and_reward_outside_strict_axes",
            "family_and_global_inference",
        ],
        "outputs": {
            evidence_path.name: {
                "sha256": _file_sha256(evidence_path),
                "bytes": evidence_path.stat().st_size,
            },
            audit_path.name: {
                "sha256": _file_sha256(audit_path),
                "bytes": audit_path.stat().st_size,
            },
        },
    }
    summary_path = output_dir / "sanbi_eflora_south_africa_source_package_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-path", type=Path)
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache/sanbi-eflora"))
    parser.add_argument("--dwca-dir", type=Path)
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument("--current-direct-ledger-csv", type=Path, required=True)
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=Path("config/measurement_classification.yml"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.archive_path is None:
        archive_path, _ = download_archive(
            args.cache_dir / f"flora_descriptions-v{ARCHIVE_VERSION}.zip"
        )
    else:
        archive_path = args.archive_path
    dwca_dir = args.dwca_dir or args.cache_dir / f"dwca-v{ARCHIVE_VERSION}"
    summary = acquire(
        archive_path=archive_path,
        dwca_dir=dwca_dir,
        master_csv=args.master_csv,
        strict_coverage_csv=args.strict_coverage_csv,
        current_direct_ledger_csv=args.current_direct_ledger_csv,
        measurement_config=args.measurement_config,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
