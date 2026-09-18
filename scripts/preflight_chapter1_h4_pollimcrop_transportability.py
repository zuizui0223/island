"""Outcome-blind PolLimCrop transportability preflight for Chapter 1 H4."""

from __future__ import annotations

import hashlib
import io
import json
import re
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h4_prospective_temporal_replication import (
    _h4a_support,
    _h4b_support,
    prepare_frozen_trait_states,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_h4_pollimcrop_transportability_preflight_v1"


def _request_bytes(url: str, *, timeout: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "island-h4-pollimcrop-preflight/1.0"},
    )
    retryable = {429, 500, 502, 503, 504}
    for attempt in range(7):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in retryable or attempt == 6:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = (
                float(retry_after)
                if retry_after and retry_after.isdigit()
                else min(30.0, 2.0**attempt)
            )
            time.sleep(delay)
        except urllib.error.URLError:
            if attempt == 6:
                raise
            time.sleep(min(30.0, 2.0**attempt))
    raise RuntimeError("unreachable PolLimCrop download retry loop")


def _get_json(url: str) -> Any:
    return json.loads(_request_bytes(url, timeout=60).decode("utf-8"))


def _get_bytes(url: str) -> bytes:
    return _request_bytes(url, timeout=120)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected PolLimCrop H4 preflight contract")
    return value


def discover_dataset_file(config: dict[str, Any]) -> dict[str, Any]:
    collection_id = int(config["source"]["figshare_collection_id"])
    articles = _get_json(
        f"https://api.figshare.com/v2/collections/{collection_id}/articles?page_size=1000"
    )
    if not isinstance(articles, list) or not articles:
        raise typer.BadParameter("Figshare collection returned no articles")

    candidates: list[dict[str, Any]] = []
    for item in articles:
        article_id = item.get("id")
        if article_id is None:
            continue
        detail = _get_json(f"https://api.figshare.com/v2/articles/{article_id}")
        for file_info in detail.get("files", []):
            name = str(file_info.get("name", ""))
            download_url = str(file_info.get("download_url", ""))
            if not download_url:
                continue
            lower = name.casefold()
            score = 0
            if "pollimcrop" in lower:
                score += 10
            if "dataset" in lower:
                score += 5
            if lower.endswith(".zip"):
                score += 3
            if lower.endswith(".csv"):
                score += 2
            candidates.append(
                {
                    "article_id": int(article_id),
                    "article_title": str(detail.get("title", "")),
                    "file_id": int(file_info.get("id", 0) or 0),
                    "name": name,
                    "download_url": download_url,
                    "score": score,
                }
            )
    if not candidates:
        raise typer.BadParameter("no downloadable Figshare files found")
    candidates.sort(key=lambda x: (-int(x["score"]), str(x["name"])))
    return candidates[0]


def extract_dataset_csv(payload: bytes, file_name: str) -> tuple[bytes, str]:
    if file_name.casefold().endswith(".csv"):
        return payload, file_name
    if not zipfile.is_zipfile(io.BytesIO(payload)):
        raise typer.BadParameter(f"selected Figshare file is not CSV/ZIP: {file_name}")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        members = [
            name
            for name in archive.namelist()
            if name.casefold().endswith(".csv") and not name.endswith("/")
        ]
        if not members:
            raise typer.BadParameter("Figshare archive contains no CSV")
        members.sort(
            key=lambda name: (
                -int("pollimcrop" in name.casefold()),
                -int("dataset" in name.casefold()),
                name,
            )
        )
        chosen = members[0]
        return archive.read(chosen), chosen


def metadata_only_frame(csv_bytes: bytes, config: dict[str, Any]) -> tuple[pd.DataFrame, list[str]]:
    csv_format = config["source"]["csv_format"]
    delimiter = str(csv_format["delimiter"])
    decimal_mark = str(csv_format["decimal_mark"])
    encoding = str(csv_format["encoding"])
    header = pd.read_csv(
        io.BytesIO(csv_bytes),
        nrows=0,
        sep=delimiter,
        encoding=encoding,
    )
    columns = [str(x) for x in header.columns]
    allowed = {
        str(x)
        for x in config["outcome_blind_preflight"]["allowed_dataset_columns"]
    }
    required = {"species", "article_code"}
    if missing := required - set(columns):
        raise typer.BadParameter(
            f"PolLimCrop schema missing required metadata: {sorted(missing)}"
        )
    selected = [column for column in columns if column in allowed]
    frame = pd.read_csv(
        io.BytesIO(csv_bytes),
        usecols=selected,
        sep=delimiter,
        decimal=decimal_mark,
        encoding=encoding,
        dtype=str,
    ).fillna("")
    return frame, columns


def exact_binomial(value: object) -> str:
    tokens = str(value or "").strip().split()
    if len(tokens) != 2:
        return ""
    genus, epithet = tokens
    if not re.fullmatch(r"[A-Z][A-Za-z-]+", genus):
        return ""
    if not re.fullmatch(r"[a-z][A-Za-z-]+", epithet):
        return ""
    return f"{genus} {epithet}"


def build_trait_overlap(
    metadata: pd.DataFrame,
    traits: pd.DataFrame,
    parent_contract: dict[str, Any],
) -> pd.DataFrame:
    frozen = prepare_frozen_trait_states(traits, parent_contract)
    work = metadata.copy()
    work["accepted_species"] = work["species"].map(exact_binomial)
    work["study_key"] = work["article_code"].astype(str).str.strip()
    work = work.loc[
        work["accepted_species"].ne("") & work["study_key"].ne("")
    ].copy()
    work = work.merge(
        frozen,
        on="accepted_species",
        how="left",
        validate="many_to_one",
    )
    return work


def support_summary(
    matched: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    support_contract = {
        "primary_hypotheses": {
            "H4a_reproductive_assurance": {
                "support_gate": config["support_gates"]["H4a_reproductive_assurance"]
            },
            "H4b_accessibility_generalization": {
                "support_gate": config["support_gates"]["H4b_accessibility_generalization"]
            },
        }
    }
    return {
        "H4a_reproductive_assurance": _h4a_support(matched, support_contract),
        "H4b_accessibility_generalization": _h4b_support(matched, support_contract),
    }


@app.command("run")
def run(
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    parent_contract_path: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    parent = yaml.safe_load(parent_contract_path.read_text(encoding="utf-8"))
    selected = discover_dataset_file(config)
    payload = _get_bytes(selected["download_url"])
    csv_bytes, csv_member = extract_dataset_csv(payload, selected["name"])
    metadata, schema = metadata_only_frame(csv_bytes, config)
    traits = pd.read_csv(trait_states_csv)
    matched = build_trait_overlap(metadata, traits, parent)
    support = support_summary(matched, config)

    species = sorted(set(matched["accepted_species"].astype(str)) - {""})
    result = {
        "contract": CONTRACT,
        "inferential_role": "secondary_external_domain_transportability_preflight",
        "outcomes_read": False,
        "wild_temporal_replication_repaired": False,
        "source_format": {
            "delimiter": str(config["source"]["csv_format"]["delimiter"]),
            "decimal_mark": str(config["source"]["csv_format"]["decimal_mark"]),
            "encoding": str(config["source"]["csv_format"]["encoding"]),
        },
        "figshare": {
            "collection_id": int(config["source"]["figshare_collection_id"]),
            "article_id": selected["article_id"],
            "file_id": selected["file_id"],
            "file_name": selected["name"],
            "file_sha256": _sha256_bytes(payload),
            "csv_member": csv_member,
            "csv_sha256": _sha256_bytes(csv_bytes),
        },
        "schema_columns": schema,
        "metadata_columns_read": list(metadata.columns),
        "n_metadata_rows": int(len(metadata)),
        "n_metadata_species": int(
            metadata["species"].astype(str).replace("", np.nan).nunique(dropna=True)
        ),
        "n_exact_binomial_rows": int(matched["accepted_species"].ne("").sum()),
        "n_exact_binomial_species": int(len(species)),
        "n_studies": int(matched["study_key"].nunique()),
        "support": support,
        "response_mapping_status": "not_frozen_no_outcome_values_read",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output_dir / "pollimcrop_metadata_only.csv", index=False)
    matched[
        [
            "study_key",
            "accepted_species",
            "autonomous_selfing",
            "generalized_form",
            "actinomorphic_symmetry",
            "shallow_open_tube",
            "accessibility_component_count",
            "accessibility_generalization_score",
        ]
    ].drop_duplicates().to_csv(
        output_dir / "pollimcrop_trait_overlap_outcome_blind.csv",
        index=False,
    )
    (output_dir / "PREFLIGHT.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
