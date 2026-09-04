"""Convert a private TRY v6 request into audited island-v2 direct evidence.

Raw TRY rows are never written to the repository.  The adapter keeps original
references/TRY identifiers, collapses redistributions by original-source
lineage, and keeps corolla fusion separate from ``floral_form``.
"""
from __future__ import annotations

import hashlib
import json
import re
import zipfile
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)
TRY_URL = "https://www.try-db.org/"
SOURCE_NAME = "try_v6_request_51897"
MEMBER = "51897.txt"
TARGET_IDS = {"207", "211", "2935", "2936"}
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:a-z0-9]+", re.I)
SPECIES_RE = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_symmetry": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
}
CANDIDATE_COLUMNS = [
    "evidence_id", "source_name", "source_record_id", "accepted_species", "genus",
    "family", "name_match_method", "trait_name", "raw_value", "standardized_value",
    "candidate_kind", "evidence_scope", "source_type", "source_reliability",
    "source_citation", "source_url", "evidence_excerpt", "supporting_taxa",
    "inference_rule", "confidence", "inference_note", "needs_human_review",
    "review_status", "analysis_role", "source_lineage", "try_dataset_id", "try_dataset",
    "try_observation_ids", "try_obs_data_ids", "try_data_ids", "try_data_names",
    "try_species_names", "try_accepted_species_name",
]
COMMON_COLUMNS = [
    "accepted_species", "axis", "trait_name", "normalized_value", "quality",
    "source_group", "source_provider", "source_url", "source_record_id",
    "source_citation", "source_excerpt", "evidence_scope", "name_match_method",
    "source_lineage", "lineage_method", "source_run_id", "source_artifact",
    "source_file", "acceptance_contract",
]
SIDE_COLUMNS = [
    "accepted_species", "genus", "family", "name_match_method", "trait_name",
    "raw_value", "standardized_value", "source_lineage", "source_citation",
    "source_url", "source_record_id", "evidence_scope", "confidence",
    "try_dataset_id", "try_dataset", "try_observation_ids", "try_obs_data_ids",
    "try_data_ids", "try_data_names", "try_species_names", "try_accepted_species_name",
]
EXCLUSION_COLUMNS = [
    "source_row", "try_trait_id", "try_trait_name", "try_data_id", "try_data_name",
    "try_species_name", "try_accepted_species_name", "raw_value", "accepted_species",
    "reason", "source_lineage", "source_citation",
]


def text(value: Any) -> str:
    return "" if value is None or pd.isna(value) else " ".join(str(value).split())


def stable_id(*parts: object) -> str:
    raw = "|".join(text(part) for part in parts)
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source_lineage(reference: object, dataset_id: object, dataset: object) -> tuple[str, str]:
    citation = text(reference)
    if match := DOI_RE.search(citation):
        return f"doi:{match.group(0).rstrip('.,;)').casefold()}", "doi"
    folded = re.sub(r"[^a-z0-9]+", " ", re.sub(r"https?://\S+", " ", citation.casefold()))
    folded = " ".join(folded.split())
    if folded:
        key = hashlib.sha256(folded.encode()).hexdigest()[:20]
        return f"citation:{key}", "citation_fingerprint"
    key = hashlib.sha256(f"{text(dataset_id)}|{text(dataset)}".casefold().encode()).hexdigest()[:20]
    return f"origin:try-dataset:{key}", "try_dataset_fallback"


def load_synonyms(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    frame = pd.read_csv(path, dtype=str).fillna("")
    for left in ("synonym", "source_name", "matched_source_name", "original_name"):
        if {left, "accepted_species"}.issubset(frame.columns):
            return {
                text(a).casefold(): text(b)
                for a, b in frame[[left, "accepted_species"]].itertuples(index=False, name=None)
                if text(a) and text(b)
            }
    raise ValueError(
        "synonym map needs synonym/source_name/matched_source_name/original_name plus accepted_species"
    )


def master_maps(master: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    if "accepted_species" not in master.columns:
        raise ValueError("master list must contain accepted_species")
    canonical, genus, family = {}, {}, {}
    for row in master.fillna("").to_dict("records"):
        sp = text(row.get("accepted_species"))
        if not sp:
            continue
        canonical[sp.casefold()] = sp
        genus[sp] = text(row.get("genus")) or sp.split()[0]
        family[sp] = text(row.get("family"))
    return canonical, genus, family


def binomial(value: object) -> str:
    tokens = text(value).split()
    if len(tokens) < 2 or any(
        x.casefold().rstrip(".") in {"subsp", "ssp", "var", "forma", "f"}
        for x in tokens[2:]
    ):
        return ""
    candidate = " ".join(tokens[:2])
    return candidate if SPECIES_RE.fullmatch(candidate) else ""


def resolve_name(
    acc: object,
    raw: object,
    canonical: dict[str, str],
    synonyms: dict[str, str],
) -> tuple[str, str, str]:
    for value in (text(acc), text(raw)):
        if value.casefold() in canonical:
            return canonical[value.casefold()], "accepted_name_exact", "species_direct"
    for value in (text(acc), text(raw)):
        accepted = synonyms.get(value.casefold(), "")
        if accepted.casefold() in canonical:
            return canonical[accepted.casefold()], "exact_synonym_to_accepted", "synonym_direct"
    for value in (text(acc), text(raw)):
        candidate = binomial(value)
        if candidate.casefold() in canonical:
            return canonical[candidate.casefold()], "accepted_name_author_stripped", "species_direct"
    return "", "", ""


COLOR_TERMS = {
    "white": (
        r"\bwhite\b", r"\bcream\b", r"\bivory\b", r"weiss", r"weiß", r"blanc[oa]?", r"creme"
    ),
    "yellow_orange": (
        r"\byellow\b", r"\borange\b", r"gelb", r"gold", r"schwefel", r"zitrone", r"amarill[oa]?", r"naranja"
    ),
    "red_pink": (
        r"\bred\b", r"\bpink\b", r"\brosa\b", r"rosad[oa]?", r"\brojo\b", r"\brot", r"röt"
    ),
    "blue_purple": (
        r"\bblue\b", r"\bpurple\b", r"\bviolet\b", r"blau", r"violett", r"purpur", r"\blila\b", r"\bazul\b", r"morado"
    ),
    "green_brown_inconspicuous": (
        r"\bgreen\b", r"\bbrown\b", r"grün", r"grun", r"braun", r"\bverde\b"
    ),
}


def normalize_colour(value: object) -> str:
    raw = text(value).casefold()
    if not raw or raw in {"/", "-", "unknown", "absent", "sin flor"}:
        return ""
    if re.search(r"\b(various|variable|multicolou?r|verschieden)\b", raw):
        return "multicolored_variable"
    states = {
        state
        for state, terms in COLOR_TERMS.items()
        if any(re.search(term, raw, re.I) for term in terms)
    }
    if len(states) == 1:
        return next(iter(states))
    if len(states) > 1:
        return "multicolored_variable"
    return "other_described" if raw in {"black", "schwarz", "negro"} else ""


def normalize_symmetry(record: dict[str, Any]) -> tuple[str, str]:
    raw = text(record.get("OrigValueStr")).casefold()
    name = text(record.get("DataName")).casefold()
    if raw == "radial":
        return "actinomorphic", ""
    if raw == "bilateral":
        return "zygomorphic", ""
    for label, state in (("zygomorphic", "zygomorphic"), ("actinomorphic", "actinomorphic")):
        if label in name:
            if raw == "yes":
                return state, ""
            if raw == "no":
                return "", "negative_boolean_not_opposite_evidence"
            return "", "unsupported_floral_symmetry_state"
    return "", "unsupported_floral_symmetry_state"


def normalize_si(record: dict[str, Any]) -> tuple[str, str]:
    raw = text(record.get("OrigValueStr")).casefold()
    name = text(record.get("DataName")).casefold()
    data_id = text(record.get("DataID"))
    if not raw:
        return "", "missing_raw_value"
    if "genus" in raw or "family" in raw:
        return "", "non_species_si_scope"
    if data_id == "2045" or name == "inbreeding":
        return "", "inbreeding_measure_not_self_incompatibility"
    if "inbreeding depression" in raw or "lethal allele" in raw:
        return "", "related_reproductive_measure_not_si_state"
    if raw == "unknown mechanism":
        return "", "unknown_si_mechanism_without_state"
    if raw.startswith("±") or ("self compatible" in raw and "self incompatible" in raw):
        return "mixed_or_variable", ""
    if "self-compatibel species" in raw or "self-compatible species" in raw or raw == "self compatible":
        return "SC", ""
    if "self-incompatibel species" in raw or "self-incompatible species" in raw or raw == "self incompatible":
        return "SI", ""
    if raw == "present, mechanism unknown":
        return "SI", ""
    mechanisms = (
        "gametophytic", "sporophytic", "heteromorphic self-incompatibility",
        "postcygotic self-incompatibility", "cryptic self-incompatibility",
    )
    if any(x in raw for x in mechanisms) and (
        "self-incompatibility" in name or "incompatibility systems" in name
    ):
        return "SI", ""
    return "", "unsupported_or_non_si_state"


def normalize_corolla(record: dict[str, Any]) -> tuple[str, str]:
    raw = text(record.get("OrigValueStr")).casefold()
    if raw == "free":
        return "free", ""
    if raw == "(partly) united":
        return "partly_or_fully_fused", ""
    return "", "unsupported_corolla_fusion_state"


def reference(record: dict[str, Any]) -> str:
    return text(record.get("Reference")) or (
        f"TRY dataset {text(record.get('DatasetID'))}: {text(record.get('Dataset'))}"
    )


def joined(records: list[dict[str, Any]], key: str) -> str:
    return "|".join(sorted({text(r.get(key)) for r in records if text(r.get(key))}))


def excerpt(records: list[dict[str, Any]]) -> str:
    out = []
    for r in records[:4]:
        out.append(
            "; ".join(
                [
                    f"TRY DatasetID={text(r.get('DatasetID'))}",
                    f"ObservationID={text(r.get('ObservationID'))}",
                    f"DataID={text(r.get('DataID'))}",
                    f"DataName={text(r.get('DataName'))}",
                    f"SpeciesName={text(r.get('SpeciesName'))}",
                    f"AccSpeciesName={text(r.get('AccSpeciesName'))}",
                    f"OrigValueStr={text(r.get('OrigValueStr'))}",
                ]
            )
        )
    return " || ".join(out)


def exclusion(
    record: dict[str, Any], reason: str, accepted: str, lineage: str, citation: str
) -> dict[str, str]:
    return {
        "source_row": text(record.get("source_row")),
        "try_trait_id": text(record.get("TraitID")),
        "try_trait_name": text(record.get("TraitName")),
        "try_data_id": text(record.get("DataID")),
        "try_data_name": text(record.get("DataName")),
        "try_species_name": text(record.get("SpeciesName")),
        "try_accepted_species_name": text(record.get("AccSpeciesName")),
        "raw_value": text(record.get("OrigValueStr")),
        "accepted_species": accepted,
        "reason": reason,
        "source_lineage": lineage,
        "source_citation": citation,
    }


def build(
    source: pd.DataFrame,
    master: pd.DataFrame,
    synonyms: dict[str, str] | None = None,
    source_name: str = SOURCE_NAME,
):
    required = {
        "DatasetID", "Dataset", "SpeciesName", "AccSpeciesName", "ObservationID",
        "ObsDataID", "TraitID", "TraitName", "DataID", "DataName", "OrigValueStr", "Reference",
    }
    if missing := required.difference(source.columns):
        raise ValueError(f"TRY source missing columns: {sorted(missing)}")
    canonical, genus_map, family_map = master_maps(master)
    synonyms = synonyms or {}
    data = source.fillna("").copy()
    data["source_row"] = data.index + 2
    data = data.loc[data["TraitID"].astype(str).isin(TARGET_IDS)].copy()
    prepared: list[dict[str, Any]] = []
    excluded: list[dict[str, str]] = []
    name_rows: list[dict[str, str]] = []

    for acc, raw in data[["AccSpeciesName", "SpeciesName"]].drop_duplicates().itertuples(
        index=False, name=None
    ):
        sp, method, scope = resolve_name(acc, raw, canonical, synonyms)
        name_rows.append(
            {
                "try_accepted_species_name": text(acc),
                "try_species_name": text(raw),
                "accepted_species": sp,
                "name_match_method": method,
                "evidence_scope": scope,
                "matched": str(bool(sp)),
            }
        )

    for record in data.to_dict("records"):
        sp, method, scope = resolve_name(
            record.get("AccSpeciesName"), record.get("SpeciesName"), canonical, synonyms
        )
        citation = reference(record)
        lineage, lineage_method = source_lineage(
            citation, record.get("DatasetID"), record.get("Dataset")
        )
        trait_id = text(record.get("TraitID"))
        raw = text(record.get("OrigValueStr"))
        trait = value = reason = ""
        if not sp:
            reason = "name_not_in_island_master"
        elif trait_id == "207":
            trait, value = "flower_primary_color", normalize_colour(raw)
            reason = "" if value else "unsupported_or_missing_flower_colour"
        elif trait_id == "2935":
            trait = "floral_symmetry"
            value, reason = normalize_symmetry(record)
        elif trait_id == "211":
            trait = "self_incompatibility"
            value, reason = normalize_si(record)
        elif trait_id == "2936":
            trait = "corolla_fusion"
            value, reason = normalize_corolla(record)
        if reason:
            excluded.append(exclusion(record, reason, sp, lineage, citation))
            continue
        prepared.append(
            {
                **record,
                "accepted_species": sp,
                "name_match_method": method,
                "evidence_scope": scope,
                "mapped_trait": trait,
                "mapped_value": value,
                "source_lineage": lineage,
                "lineage_method": lineage_method,
                "source_citation": citation,
            }
        )

    candidates: list[dict[str, str]] = []
    side: list[dict[str, str]] = []
    if prepared:
        prep = pd.DataFrame(prepared)
        for (_, trait, lineage), group in prep.groupby(
            ["accepted_species", "mapped_trait", "source_lineage"], sort=True
        ):
            records = group.to_dict("records")
            sp = text(records[0]["accepted_species"])
            states = sorted(
                {text(r.get("mapped_value")) for r in records if text(r.get("mapped_value"))}
            )
            if trait == "flower_primary_color" and len(states) > 1:
                state = "multicolored_variable"
            elif len(states) == 1:
                state = states[0]
            else:
                excluded.extend(
                    exclusion(
                        r,
                        "within_original_source_lineage_conflict",
                        sp,
                        lineage,
                        text(r.get("source_citation")),
                    )
                    for r in records
                )
                continue
            first = records[0]
            raw_values = " | ".join(
                sorted({text(r.get("OrigValueStr")) for r in records if text(r.get("OrigValueStr"))})
            )
            rid = f"try:{stable_id(source_name, sp, trait, lineage)}"
            provenance = {
                "accepted_species": sp,
                "genus": genus_map.get(sp, sp.split()[0]),
                "family": family_map.get(sp, ""),
                "name_match_method": text(first.get("name_match_method")),
                "trait_name": trait,
                "raw_value": raw_values,
                "standardized_value": state,
                "source_lineage": lineage,
                "source_citation": text(first.get("source_citation")),
                "source_url": TRY_URL,
                "source_record_id": rid,
                "evidence_scope": text(first.get("evidence_scope")),
                "confidence": "medium",
                "try_dataset_id": joined(records, "DatasetID"),
                "try_dataset": joined(records, "Dataset"),
                "try_observation_ids": joined(records, "ObservationID"),
                "try_obs_data_ids": joined(records, "ObsDataID"),
                "try_data_ids": joined(records, "DataID"),
                "try_data_names": joined(records, "DataName"),
                "try_species_names": joined(records, "SpeciesName"),
                "try_accepted_species_name": joined(records, "AccSpeciesName"),
            }
            if trait == "corolla_fusion":
                side.append(provenance)
                continue
            candidates.append(
                {
                    "evidence_id": rid,
                    "source_name": source_name,
                    "source_record_id": rid,
                    **{
                        k: provenance[k]
                        for k in (
                            "accepted_species", "genus", "family", "name_match_method",
                            "trait_name", "raw_value", "standardized_value",
                        )
                    },
                    "candidate_kind": "source_backed",
                    "evidence_scope": provenance["evidence_scope"],
                    "source_type": "institutional_trait_database",
                    "source_reliability": "A_institutional_trait_database",
                    "source_citation": provenance["source_citation"],
                    "source_url": TRY_URL,
                    "evidence_excerpt": excerpt(records),
                    "supporting_taxa": "",
                    "inference_rule": "",
                    "confidence": "medium",
                    "inference_note": (
                        "TRY redistribution mapped conservatively; original-source lineage preserved"
                    ),
                    "needs_human_review": "False",
                    "review_status": "structured_provider_mapping",
                    "analysis_role": "strict_species_direct_candidate",
                    "source_lineage": lineage,
                    **{
                        k: provenance[k]
                        for k in (
                            "try_dataset_id", "try_dataset", "try_observation_ids",
                            "try_obs_data_ids", "try_data_ids", "try_data_names",
                            "try_species_names", "try_accepted_species_name",
                        )
                    },
                }
            )

    c = pd.DataFrame(candidates, columns=CANDIDATE_COLUMNS)
    s = pd.DataFrame(side, columns=SIDE_COLUMNS)
    e = pd.DataFrame(excluded, columns=EXCLUSION_COLUMNS)
    n = pd.DataFrame(name_rows)
    summary = {
        "source_name": source_name,
        "try_rows_total": int(len(source)),
        "target_trait_rows": int(len(data)),
        "master_species": len(canonical),
        "candidate_rows": len(c),
        "candidate_species": int(c.accepted_species.nunique()) if len(c) else 0,
        "corolla_fusion_rows": len(s),
        "corolla_fusion_species": int(s.accepted_species.nunique()) if len(s) else 0,
        "exclusion_rows": len(e),
        "candidates_by_trait": (
            c.groupby("trait_name").size().astype(int).to_dict() if len(c) else {}
        ),
        "candidate_species_by_trait": (
            c.groupby("trait_name").accepted_species.nunique().astype(int).to_dict()
            if len(c)
            else {}
        ),
        "exclusions_by_reason": e.reason.value_counts().astype(int).to_dict() if len(e) else {},
        "matched_try_name_pairs": int((n.matched == "True").sum()) if len(n) else 0,
        "unmatched_try_name_pairs": int((n.matched != "True").sum()) if len(n) else 0,
    }
    return c, s, n, e, summary


def common_evidence(candidates: pd.DataFrame, source_name: str = SOURCE_NAME) -> pd.DataFrame:
    rows = []
    for r in candidates.fillna("").to_dict("records"):
        trait = text(r.get("trait_name"))
        axis = TRAIT_AXIS.get(trait, "")
        if axis:
            rows.append(
                {
                    "accepted_species": text(r.get("accepted_species")),
                    "axis": axis,
                    "trait_name": trait,
                    "normalized_value": text(r.get("standardized_value")),
                    "quality": "medium",
                    "source_group": "try",
                    "source_provider": source_name,
                    "source_url": TRY_URL,
                    "source_record_id": text(r.get("source_record_id")),
                    "source_citation": text(r.get("source_citation")),
                    "source_excerpt": text(r.get("evidence_excerpt")),
                    "evidence_scope": text(r.get("evidence_scope")),
                    "name_match_method": text(r.get("name_match_method")),
                    "source_lineage": text(r.get("source_lineage")),
                    "lineage_method": "try_original_reference_lineage",
                    "source_run_id": source_name,
                    "source_artifact": "private_try_request_local_prepare",
                    "source_file": "trait_candidates.csv.gz",
                    "acceptance_contract": "try_v6_species_direct_source_lineage_v1",
                }
            )
    return pd.DataFrame(rows, columns=COMMON_COLUMNS)


def read_zip(path: Path, member: str) -> pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        if member not in z.namelist():
            txt = [name for name in z.namelist() if name.casefold().endswith(".txt")]
            if len(txt) != 1:
                raise ValueError(f"TRY ZIP text member is ambiguous: {txt}")
            member = txt[0]
        with z.open(member) as f:
            return pd.read_csv(
                f, sep="\t", dtype=str, encoding="cp1252", low_memory=False
            ).fillna("")


def write_gz(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


@app.command("prepare")
def prepare(
    try_zip: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    member: str = MEMBER,
    synonym_map_csv: Annotated[Path | None, typer.Option(exists=True, dir_okay=False)] = None,
    source_name: str = SOURCE_NAME,
) -> None:
    source = read_zip(try_zip, member)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    candidates, side, names, exclusions, summary = build(
        source, master, load_synonyms(synonym_map_csv), source_name
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_gz(candidates, output_dir / "trait_candidates.csv.gz")
    write_gz(
        common_evidence(candidates, source_name),
        output_dir / "try_common_direct_evidence.csv.gz",
    )
    write_gz(side, output_dir / "try_corolla_fusion_sidecar.csv.gz")
    write_gz(names, output_dir / "try_name_audit.csv.gz")
    write_gz(exclusions, output_dir / "try_exclusions.csv.gz")
    summary.update(
        {
            "try_zip_sha256": sha256(try_zip),
            "master_csv_sha256": sha256(master_csv),
            "member": member,
            "raw_try_redistributed": False,
        }
    )
    (output_dir / "try_integration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    app()
