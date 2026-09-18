"""Download the frozen PolLimCrop CSV only after a committed support-result lock exists.

This utility never parses outcome values. It resolves the same Figshare source used by
the outcome-blind preflight, verifies the frozen file/member SHA-256 values recorded in
the committed preflight result lock, and writes the exact CSV bytes for the separately
frozen analysis module.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

PREFLIGHT_SCRIPT = Path("scripts/preflight_chapter1_h4_pollimcrop_transportability.py")


def _load_preflight_module():
    spec = importlib.util.spec_from_file_location("pollimcrop_preflight", PREFLIGHT_SCRIPT)
    if spec is None or spec.loader is None:
        raise typer.BadParameter("cannot load PolLimCrop preflight module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


HYPOTHESES = (
    "H4a_reproductive_assurance",
    "H4b_accessibility_generalization",
)


def validate_support_lock(lock: dict) -> list[str]:
    expected_contract = "chapter1_h4_pollimcrop_transportability_preflight_v1"
    expected_lock_contract = (
        "chapter1_h4_pollimcrop_transportability_preflight_result_lock_v1"
    )
    if lock.get("contract") != expected_contract:
        raise typer.BadParameter(
            "unexpected PolLimCrop preflight result-lock contract"
        )
    if lock.get("lock_contract") != expected_lock_contract:
        raise typer.BadParameter(
            "PolLimCrop dataset access requires the committed support-result lock"
        )
    if lock.get("outcomes_read") is not False:
        raise typer.BadParameter(
            "preflight result lock must record outcomes_read=false"
        )

    support = lock.get("support", {})
    admitted = [
        name for name in HYPOTHESES if bool(support.get(name, {}).get("evaluable"))
    ]
    decision = lock.get("decision", {})
    recorded = [str(x) for x in decision.get("admitted_hypotheses", [])]
    if recorded != admitted:
        raise typer.BadParameter(
            "PolLimCrop support-result lock admission list is inconsistent"
        )
    authorized = bool(decision.get("outcome_extraction_authorized"))
    if authorized != bool(admitted):
        raise typer.BadParameter(
            "PolLimCrop support-result lock authorization is inconsistent"
        )
    if not admitted:
        raise typer.BadParameter(
            "no co-primary hypothesis passed frozen PolLimCrop support gate"
        )
    return admitted


@app.command("run")
def run(
    preflight_result_lock: Path = typer.Option(..., exists=True, dir_okay=False),
    preflight_config: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    lock = json.loads(preflight_result_lock.read_text(encoding="utf-8"))
    validate_support_lock(lock)

    expected_contract = "chapter1_h4_pollimcrop_transportability_preflight_v1"
    config = yaml.safe_load(preflight_config.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != expected_contract:
        raise typer.BadParameter("unexpected PolLimCrop preflight config")

    module = _load_preflight_module()
    selected = module.discover_dataset_file(config)
    figshare = lock.get("figshare", {})
    if int(selected.get("article_id", 0)) != int(figshare.get("article_id", -1)):
        raise typer.BadParameter("Figshare article changed since support lock")
    if int(selected.get("file_id", 0)) != int(figshare.get("file_id", -1)):
        raise typer.BadParameter("Figshare file changed since support lock")

    payload = module._get_bytes(selected["download_url"])
    if _sha256_bytes(payload) != str(figshare.get("file_sha256", "")):
        raise typer.BadParameter("PolLimCrop source-file SHA-256 changed since support lock")
    csv_bytes, csv_member = module.extract_dataset_csv(payload, selected["name"])
    if csv_member != str(figshare.get("csv_member", "")):
        raise typer.BadParameter("PolLimCrop CSV member changed since support lock")
    if _sha256_bytes(csv_bytes) != str(figshare.get("csv_sha256", "")):
        raise typer.BadParameter("PolLimCrop CSV SHA-256 changed since support lock")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_csv.write_bytes(csv_bytes)
    typer.echo(
        json.dumps(
            {
                "article_id": selected["article_id"],
                "file_id": selected["file_id"],
                "csv_member": csv_member,
                "csv_sha256": _sha256_bytes(csv_bytes),
                "outcome_values_parsed_by_downloader": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
