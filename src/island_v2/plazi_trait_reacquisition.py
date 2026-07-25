"""Retrieve one Plazi treatment XML and promote explicit focal floral measurements."""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "plazi_species_direct_treatment_v1"
PLAZI_FETCH_URL = "https://api.plazi.org/v1/Treatments/fetch"
USER_AGENT = "island-floral-traits/0.1 (local Plazi treatment reacquisition)"
MAX_XML_BYTES = 25_000_000
COROLLA_MEASUREMENT = re.compile(
    r"\bcorolla\s+"
    r"(?P<low>\d+(?:\.\d+)?)"
    r"(?:\s*[–—-]\s*(?P<high>\d+(?:\.\d+)?))?"
    r"\s*(?P<unit>mm|cm)\s+long\b",
    flags=re.IGNORECASE,
)
COROLLA_TUBE_MEASUREMENT = re.compile(
    r"\bcorolla\s+tube\s+"
    r"(?P<low>\d+(?:\.\d+)?)"
    r"(?:\s*[–—-]\s*(?P<high>\d+(?:\.\d+)?))?"
    r"\s*(?P<unit>mm|cm)\s+long\b",
    flags=re.IGNORECASE,
)
EVIDENCE_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "evidence_quality",
    "evidence_scope",
    "source_group",
    "source_provider",
    "source_url",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "name_match_method",
    "inference_rule",
    "raw_content_sha256",
    "source_license",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).replace("\xa0", " ").split())


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_plazi_fetch_url(treatment_uuid: str) -> str:
    value = _text(treatment_uuid).upper()
    if not re.fullmatch(r"[0-9A-F]{32}", value):
        raise ValueError("treatment_uuid must be 32 hexadecimal characters")
    return f"{PLAZI_FETCH_URL}?{urllib.parse.urlencode({'UUID': value})}"


def fetch_plazi_treatment_xml(
    treatment_uuid: str,
    *,
    timeout_seconds: float = 45,
    max_attempts: int = 3,
) -> bytes:
    url = build_plazi_fetch_url(treatment_uuid)
    for attempt in range(1, max_attempts + 1):
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/xml,text/xml,*/*",
                "User-Agent": USER_AGENT,
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                xml_bytes = response.read(MAX_XML_BYTES + 1)
                if len(xml_bytes) > MAX_XML_BYTES:
                    raise ValueError("Plazi treatment XML exceeded 25 MB")
                return xml_bytes
        except urllib.error.HTTPError as exc:
            if exc.code not in {408, 425, 429, 500, 502, 503, 504} or attempt == max_attempts:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == max_attempts:
                raise
        time.sleep(min(30.0, 2.0 ** (attempt - 1)))
    raise RuntimeError("unreachable")


def _classify_measurement(
    match: re.Match[str],
    *,
    trait_name: str,
) -> str:
    low = float(match.group("low"))
    high = float(match.group("high") or low)
    midpoint_mm = ((low + high) / 2.0) * (
        10.0 if match.group("unit").casefold() == "cm" else 1.0
    )
    thresholds = (
        (
            (2.0, "very_small"),
            (10.0, "small"),
            (30.0, "medium"),
            (100.0, "large"),
            (float("inf"), "very_large"),
        )
        if trait_name == "flower_size_class"
        else (
            (1.0, "absent_or_open"),
            (5.0, "shallow"),
            (15.0, "intermediate"),
            (float("inf"), "deep"),
        )
    )
    return next(label for maximum, label in thresholds if midpoint_mm <= maximum)


def _excerpt(text: str, start: int, end: int, pad: int = 180) -> str:
    return _text(text[max(0, start - pad) : min(len(text), end + pad)])


def extract_plazi_treatment_evidence(
    xml_bytes: bytes,
    *,
    accepted_species: str,
    treatment_uuid: str,
    source_url: str | None = None,
) -> pd.DataFrame:
    """Extract controlled structural values from the exact focal treatment."""

    if len(xml_bytes) > MAX_XML_BYTES:
        raise ValueError("Plazi treatment XML exceeded 25 MB")
    try:
        root = ElementTree.fromstring(xml_bytes)
    except ElementTree.ParseError as exc:
        raise ValueError(f"invalid Plazi treatment XML: {exc}") from exc
    species = _text(accepted_species)
    treatment = next(
        (
            element
            for element in root.iter()
            if str(element.tag).rsplit("}", 1)[-1] == "treatment"
        ),
        None,
    )
    if treatment is None:
        raise ValueError("Plazi XML contains no treatment")
    treatment_text = _text(" ".join(treatment.itertext()))
    if not re.search(
        rf"(?<![\w-]){re.escape(species)}(?![\w-])",
        treatment_text,
        flags=re.IGNORECASE,
    ):
        raise ValueError("Plazi treatment does not contain the exact accepted species")

    descriptions = [
        element
        for element in treatment.iter()
        if (
            str(element.tag).rsplit("}", 1)[-1] == "subSubSection"
            and _text(element.attrib.get("type")).casefold() == "description"
        )
    ]
    if not descriptions:
        raise ValueError("Plazi treatment contains no description section")
    description = _text(
        " ".join(
            " ".join(element.itertext())
            for element in descriptions
        )
    )
    document_doi = _text(root.attrib.get("ID-DOI")).casefold()
    if document_doi.startswith("https://doi.org/"):
        document_doi = document_doi.removeprefix("https://doi.org/")
    if document_doi.startswith("http://doi.org/"):
        document_doi = document_doi.removeprefix("http://doi.org/")
    if not document_doi:
        document_doi = _text(root.attrib.get("docSource")).casefold()
        document_doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", document_doi)
    lineage = (
        f"doi:{document_doi}"
        if document_doi
        else f"plazi:master-document:{_text(root.attrib.get('masterDocId')).casefold()}"
    )
    xml_sha = _sha256_bytes(xml_bytes)
    citation = " | ".join(
        part
        for part in (
            _text(root.attrib.get("masterDocTitle")),
            f"DOI:{document_doi}" if document_doi else "",
            f"Plazi treatment:{treatment_uuid.upper()}",
            f"raw_xml_sha256={xml_sha}",
        )
        if part
    )
    license_value = _text(root.attrib.get("zenodo-license-treatments"))
    url = source_url or build_plazi_fetch_url(treatment_uuid)
    candidates: list[tuple[str, re.Match[str]]] = []
    candidates.extend(
        ("flower_size_class", match)
        for match in COROLLA_MEASUREMENT.finditer(description)
    )
    candidates.extend(
        ("tube_depth_class", match)
        for match in COROLLA_TUBE_MEASUREMENT.finditer(description)
    )
    rows: list[dict[str, str]] = []
    for trait_name, match in candidates:
        rows.append(
            {
                "accepted_species": species,
                "axis": "floral_structural_complexity",
                "trait_name": trait_name,
                "normalized_value": _classify_measurement(
                    match,
                    trait_name=trait_name,
                ),
                "evidence_quality": "medium",
                "evidence_scope": "species_direct",
                "source_group": "direct_web_treatment",
                "source_provider": "Plazi TreatmentBank",
                "source_url": url,
                "source_citation": citation,
                "source_excerpt": _excerpt(
                    description,
                    match.start(),
                    match.end(),
                ),
                "source_record_id": (
                    f"{treatment_uuid.upper()}:{trait_name}:"
                    f"{hashlib.sha256(match.group(0).encode('utf-8')).hexdigest()}"
                ),
                "source_lineage": lineage,
                "name_match_method": "accepted_exact",
                "inference_rule": "",
                "raw_content_sha256": xml_sha,
                "source_license": license_value,
            }
        )
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).drop_duplicates(
        [
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineage",
        ]
    )


def run_plazi_treatment_pilot(
    *,
    accepted_species: str,
    treatment_uuid: str,
    output_dir: Path,
) -> dict[str, Any]:
    """Fetch, extract, and freeze one exact Plazi treatment."""

    source_url = build_plazi_fetch_url(treatment_uuid)
    xml_bytes = fetch_plazi_treatment_xml(treatment_uuid)
    evidence = extract_plazi_treatment_evidence(
        xml_bytes,
        accepted_species=accepted_species,
        treatment_uuid=treatment_uuid,
        source_url=source_url,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / f"{treatment_uuid.upper()}.xml"
    raw_path.write_bytes(xml_bytes)
    evidence_path = output_dir / "direct_web_plazi_evidence.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression="gzip")
    summary = {
        "contract_version": CONTRACT_VERSION,
        "accepted_species": _text(accepted_species),
        "source_url": source_url,
        "treatment_uuid": treatment_uuid.upper(),
        "n_evidence_rows": len(evidence),
        "n_species_axis_pure_gain_candidate": int(not evidence.empty),
        "family_inference_used": False,
        "global_fallback_used": False,
        "files": {
            raw_path.name: {
                "sha256": _file_sha256(raw_path),
                "bytes": raw_path.stat().st_size,
            },
            evidence_path.name: {
                "sha256": _file_sha256(evidence_path),
                "rows": len(evidence),
            },
        },
    }
    summary_path = output_dir / "direct_web_plazi_summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command("treatment")
def treatment_command(
    accepted_species: str = typer.Option(...),
    treatment_uuid: str = typer.Option(...),
    output_dir: Path = typer.Option(..., file_okay=False),
) -> None:
    """Run one exact Plazi treatment acquisition locally."""

    summary = run_plazi_treatment_pilot(
        accepted_species=accepted_species,
        treatment_uuid=treatment_uuid,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))
