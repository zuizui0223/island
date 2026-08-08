"""Prepare a source-lineage-safe BHL floral-trait package from NHM ADEPT.

Only fields that map directly to the strict three-axis ontology are emitted.
Sex system, pollination, and reward remain outside the strict axes.  Promotion
is intentionally separate and requires the generated deterministic audit queue
to be reviewed under the common source-package contract.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

RESOURCE_ID = "cc636680-abfb-4143-9b99-62dcb562f8da"
DATASET_DOI = "https://doi.org/10.5519/p3dm31kc"
SOURCE_RUN_ID = "nhm-adept-bhl-2026-20260808"
SOURCE_ARTIFACT = "nhm-adept-bhl-2026-datastore"
SOURCE_FILE = "bhl_exact_master_records.csv.gz"
ACCEPTANCE_CONTRACT = "nhm_adept_bhl_exact_record_quote_audit_v1"
AUDIT_QUOTA = {
    "flower_primary_color": 80,
    "tube_depth_class": 50,
    "inflorescence_display": 40,
    "floral_form": 20,
    "floral_symmetry": 10,
}

EVIDENCE_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "source_provider",
    "source_url",
    "source_record_id",
    "source_citation",
    "source_excerpt",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "source_file",
    "acceptance_contract",
]
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

COLOUR_PATTERNS = (
    (r"\b(?:white|whitish|cream|creamy|ivory)\b", "white"),
    (
        r"\b(?:yellow|yellowish|orange|orangish|orange-red|red-orange|gold|golden|stramineous|straw|ochre)\b",
        "yellow_orange",
    ),
    (
        r"\b(?:red|reddish|pink|pinkish|rose|roseate|maroon|scarlet|crimson|magenta|fuchsia|peach|cherry|coral)\b",
        "red_pink",
    ),
    (
        r"\b(?:blue|bluish|azure|purple|purplish|violet|violaceous|lilac|mauve|mauvish|lavender|plum)\b",
        "blue_purple",
    ),
    (
        r"\b(?:green|greenish|brown|brownish|grey|gray|greyish|grayish|beige|fuscous|ferruginous|ferrugineous|inconspicuous)\b",
        "green_brown_inconspicuous",
    ),
    (r"\bblack(?:ish)?\b", "other_described"),
)
FORM_PATTERNS = (
    (r"\bcampanulat", "bell_campanulate"),
    (r"\btubular\b|\btubulose\b", "tubular"),
    (r"\bsalver", "salverform"),
    (r"\binfundibuliform\b|\bfunnel\b|\btrumpet\b", "funnel_trumpet"),
    (r"\burn\b|\burceolat", "urn_urceolate"),
    (r"\bpapilion", "papilionaceous"),
    (r"\bbilabiat|\btwo[- ]lipped\b", "bilabiate"),
    (r"\bspur(?:red)?\b", "spurred"),
    (r"\brotate\b|\bcup\b|\bsaucer\b|\bstar\b", "open_radial"),
    (r"\bcalyciform\b|\bligulate\b|\bgibbose\b", "other_described"),
)
INFLORESCENCE_PATTERNS = (
    (r"\bsolitary\b", "solitary"),
    (r"\bfew[- ]flower", "few_flowered"),
    (r"\braceme\b|\bspike\b|\bpanicle\b|\bthyrse\b", "raceme_spike_panicle"),
    (r"\bumbel\b|\bcorymb\b|\bcyme\b|\bcincinnus\b", "umbel_corymb"),
    (r"\bheads?\b|\bcapitul(?:um|a)\b", "composite_display"),
    (r"\bbrush\b|\bpuff\b", "brush_puff_display"),
)


def _states(value: str, patterns: tuple[tuple[str, str], ...]) -> str:
    states = {state for pattern, state in patterns if re.search(pattern, value, re.IGNORECASE)}
    return "|".join(sorted(states))


def _base_record(record: dict[str, str], trait: str, value: str, quote: str) -> dict[str, str]:
    page_id = str(record.get("source_id", "")).strip()
    nhm_id = str(record.get("_id", "")).strip()
    species = str(record.get("taxon", "")).strip()
    axis = "flower_colour" if trait == "flower_primary_color" else "floral_structural_complexity"
    return {
        "accepted_species": species,
        "axis": axis,
        "trait_name": trait,
        "normalized_value": value,
        "quality": "high",
        "source_group": "latest_public_web",
        "source_provider": "nhm_adept_bhl_2026",
        "source_url": f"https://www.biodiversitylibrary.org/page/{page_id}",
        "source_record_id": f"nhm:{nhm_id}:{trait}",
        "source_citation": (
            "Araujo, Hartley & Brummitt (2026), Angiosperm functional traits "
            f"extracted from botanical descriptions, {DATASET_DOI}; BHL page {page_id}"
        ),
        "source_excerpt": quote,
        "evidence_scope": "species_direct",
        "name_match_method": "wfo_accepted_name_exact_to_master",
        "source_lineage": f"bhl:page:{page_id}",
        "lineage_method": "original_bhl_page_id",
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_file": SOURCE_FILE,
        "acceptance_contract": ACCEPTANCE_CONTRACT,
    }


def _extract_record(record: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    colour_fields = ["flower colour", "petal colour", "corolla colour", "labellum colour"]
    colour_parts = [f"{field}: {record[field]}" for field in colour_fields if record.get(field)]
    if colour_parts:
        quote = "; ".join(colour_parts)
        value = _states(quote, COLOUR_PATTERNS)
        if value:
            rows.append(_base_record(record, "flower_primary_color", value, quote))

    symmetry = str(record.get("flower symmetry", ""))
    symmetry_states = []
    if re.search(r"\bactinomorphic\b", symmetry, re.IGNORECASE):
        symmetry_states.append("actinomorphic")
    if re.search(r"\bzygomorphic\b", symmetry, re.IGNORECASE):
        symmetry_states.append("zygomorphic")
    if re.search(r"\basymmetric\b", symmetry, re.IGNORECASE):
        symmetry_states.append("asymmetric")
    if symmetry_states:
        rows.append(
            _base_record(
                record,
                "floral_symmetry",
                "|".join(sorted(set(symmetry_states))),
                f"flower symmetry: {symmetry}",
            )
        )

    shape = str(record.get("flower shape", ""))
    value = _states(shape, FORM_PATTERNS)
    if value:
        rows.append(_base_record(record, "floral_form", value, f"flower shape: {shape}"))

    inflorescence = str(record.get("inflorescence arrangement", ""))
    value = _states(inflorescence, INFLORESCENCE_PATTERNS)
    if value:
        rows.append(
            _base_record(
                record,
                "inflorescence_display",
                value,
                f"inflorescence arrangement: {inflorescence}",
            )
        )

    tube_fields = [
        "corolla tube min_ length [cm]",
        "corolla tube max_ length [cm]",
    ]
    tube_parts = [f"{field}: {record[field]}" for field in tube_fields if record.get(field)]
    measurements: list[float] = []
    for field in tube_fields:
        measurements.extend(
            float(number) * 10
            for number in re.findall(r"\d+(?:\.\d+)?", str(record.get(field, "")))
        )
    if measurements:
        midpoint = sum(measurements) / len(measurements)
        value = "shallow" if midpoint <= 5 else "intermediate" if midpoint <= 20 else "deep"
        rows.append(_base_record(record, "tube_depth_class", value, "; ".join(tube_parts)))
    return rows


def prepare(*, input_csv: Path, output_dir: Path) -> dict[str, object]:
    source = pd.read_csv(input_csv, dtype=str).fillna("")
    rows = [item for record in source.to_dict("records") for item in _extract_record(record)]
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    evidence = evidence.drop_duplicates().sort_values(
        ["accepted_species", "trait_name", "source_record_id"], kind="stable"
    )
    samples = []
    for trait, quota in AUDIT_QUOTA.items():
        candidates = evidence.loc[evidence["trait_name"].eq(trait)].drop_duplicates(
            ["accepted_species", "source_lineage", "trait_name"]
        )
        if len(candidates) < quota:
            raise ValueError(f"not enough {trait} rows for audit: {len(candidates)} < {quota}")
        samples.append(candidates.sample(n=quota, random_state=20260808))
    audit = pd.concat(samples, ignore_index=True).sort_values(
        ["trait_name", "accepted_species", "source_record_id"], kind="stable"
    )
    audit = audit.rename(
        columns={"source_provider": "provider", "source_excerpt": "exact_supporting_quote"}
    )
    audit["candidate_id"] = audit["source_record_id"]
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = ""
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "nhm_adept_bhl_evidence.csv.gz"
    queue_path = output_dir / "nhm_adept_bhl_audit_queue_200.csv"
    evidence.reindex(columns=EVIDENCE_COLUMNS).to_csv(
        evidence_path, index=False, compression="gzip"
    )
    audit.reindex(columns=AUDIT_COLUMNS).to_csv(queue_path, index=False)
    summary = {
        "contract": "nhm_adept_bhl_source_package_v1",
        "resource_id": RESOURCE_ID,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "source_rows": len(source),
        "evidence_rows": len(evidence),
        "unique_species_trait": int(
            evidence.drop_duplicates(["accepted_species", "trait_name"]).shape[0]
        ),
        "unique_species": int(evidence["accepted_species"].nunique()),
        "trait_rows": {
            str(key): int(value)
            for key, value in evidence["trait_name"].value_counts().sort_index().items()
        },
        "audit_queue_rows": len(audit),
        "audit_quota": AUDIT_QUOTA,
        "promotion_performed": False,
        "strict_exclusions": [
            "flower_architecture_not_floral_form",
            "reproduction_system_is_sex_system_not_reproductive_assurance",
            "pollination_kept_outside_strict_axes",
            "reward_kept_outside_strict_axes",
            "bract_calyx_pollen_fruit_colour_excluded",
        ],
    }
    (output_dir / "nhm_adept_bhl_source_package_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(input_csv=args.input_csv, output_dir=args.output_dir),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
