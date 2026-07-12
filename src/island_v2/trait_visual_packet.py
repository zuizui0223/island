"""Freeze GBIF occurrence images into auditable visual-trait review packets.

The packet builder downloads a bounded number of already taxon-matched media
rows, hashes the bytes, and creates contact sheets for visual inspection. It
does not infer any trait value. A later human or vision model can score the
frozen images while retaining exact media provenance and image hashes.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Callable
from io import BytesIO
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import typer
from PIL import Image, ImageDraw, ImageFont, ImageOps

app = typer.Typer(add_completion=False, no_args_is_help=True)

IMAGE_COLUMNS = [
    "accepted_species",
    "image_rank",
    "media_identifier",
    "occurrence_key",
    "media_license",
    "local_image",
    "image_sha256",
    "download_status",
    "download_error",
]

Downloader = Callable[[str], bytes]


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "species"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _default_downloader(url: str) -> bytes:
    response = httpx.get(
        url,
        timeout=60.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.1 (visual-trait-packet)"},
    )
    response.raise_for_status()
    return response.content


def _load_image(data: bytes) -> Image.Image:
    image = Image.open(BytesIO(data))
    image.load()
    return image.convert("RGB")


def freeze_images(
    media: pd.DataFrame,
    output_dir: Path,
    *,
    images_per_species: int,
    downloader: Downloader = _default_downloader,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {"accepted_species", "media_identifier"}
    missing = required.difference(media.columns)
    if missing:
        raise typer.BadParameter(f"media manifest missing columns: {sorted(missing)}")

    image_dir = output_dir / "images"
    image_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for species, group in media.fillna("").groupby("accepted_species", sort=True):
        seen: set[str] = set()
        rank = 0
        for record in group.to_dict("records"):
            url = _text(record.get("media_identifier"))
            if not url or url in seen or rank >= images_per_species:
                continue
            seen.add(url)
            rank += 1
            local_name = f"{_slug(str(species))}__{rank:02d}.jpg"
            local_path = image_dir / local_name
            status = "downloaded"
            error = ""
            sha = ""
            try:
                data = downloader(url)
                image = _load_image(data)
                image.save(local_path, format="JPEG", quality=92)
                sha = _sha256(local_path.read_bytes())
            except Exception as exc:  # noqa: BLE001 - preserve per-image failures
                status = "failed"
                error = str(exc)[:1000]
                local_name = ""
            rows.append(
                {
                    "accepted_species": str(species),
                    "image_rank": str(rank),
                    "media_identifier": url,
                    "occurrence_key": _text(record.get("occurrence_key")),
                    "media_license": _text(record.get("media_license")),
                    "local_image": local_name,
                    "image_sha256": sha,
                    "download_status": status,
                    "download_error": error,
                }
            )
    frame = pd.DataFrame(rows, columns=IMAGE_COLUMNS)
    report = {
        "n_species_requested": int(media["accepted_species"].nunique()) if len(media) else 0,
        "n_species_with_downloaded_image": int(
            frame.loc[frame["download_status"].eq("downloaded"), "accepted_species"].nunique()
        )
        if len(frame)
        else 0,
        "n_images_attempted": len(frame),
        "n_images_downloaded": int(frame["download_status"].eq("downloaded").sum())
        if len(frame)
        else 0,
        "n_images_failed": int(frame["download_status"].eq("failed").sum()) if len(frame) else 0,
        "images_per_species": images_per_species,
        "policy": (
            "Frozen images are visual evidence only. Trait values must be extracted in a later "
            "step and remain review candidates until validated."
        ),
    }
    return frame, report


def build_contact_sheets(
    frozen: pd.DataFrame,
    output_dir: Path,
    *,
    species_per_sheet: int = 5,
    thumb_size: tuple[int, int] = (320, 240),
) -> list[Path]:
    successful = frozen.loc[frozen["download_status"].eq("downloaded")].copy()
    if successful.empty:
        return []
    species_order = list(dict.fromkeys(successful["accepted_species"].astype(str)))
    sheets_dir = output_dir / "contact_sheets"
    sheets_dir.mkdir(parents=True, exist_ok=True)
    font = ImageFont.load_default()
    paths: list[Path] = []
    image_dir = output_dir / "images"
    cols = max(1, int(successful.groupby("accepted_species").size().max()))
    label_height = 34
    for page_index, start in enumerate(range(0, len(species_order), species_per_sheet), start=1):
        batch = species_order[start : start + species_per_sheet]
        canvas = Image.new(
            "RGB",
            (cols * thumb_size[0], len(batch) * (thumb_size[1] + label_height)),
            "white",
        )
        draw = ImageDraw.Draw(canvas)
        for row_index, species in enumerate(batch):
            y = row_index * (thumb_size[1] + label_height)
            draw.text((8, y + 8), species, fill="black", font=font)
            subset = successful.loc[successful["accepted_species"].eq(species)]
            for col_index, record in enumerate(subset.to_dict("records")):
                image = Image.open(image_dir / record["local_image"]).convert("RGB")
                thumb = ImageOps.contain(image, thumb_size)
                x = col_index * thumb_size[0] + (thumb_size[0] - thumb.width) // 2
                iy = y + label_height + (thumb_size[1] - thumb.height) // 2
                canvas.paste(thumb, (x, iy))
                draw.text(
                    (col_index * thumb_size[0] + 6, y + label_height + 6),
                    f"#{record['image_rank']}",
                    fill="black",
                    font=font,
                )
        path = sheets_dir / f"visual_packet_{page_index:03d}.jpg"
        canvas.save(path, format="JPEG", quality=90)
        paths.append(path)
    return paths


@app.command()
def build(
    media_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    images_per_species: int = typer.Option(5, min=1, max=10),
    species_per_sheet: int = typer.Option(5, min=1, max=10),
) -> None:
    media = pd.read_csv(media_csv, dtype=str).fillna("")
    output_dir.mkdir(parents=True, exist_ok=True)
    frozen, report = freeze_images(
        media,
        output_dir,
        images_per_species=images_per_species,
    )
    frozen.to_csv(output_dir / "frozen_visual_evidence.csv", index=False)
    sheets = build_contact_sheets(
        frozen,
        output_dir,
        species_per_sheet=species_per_sheet,
    )
    report["n_contact_sheets"] = len(sheets)
    report["contact_sheets"] = [str(path.relative_to(output_dir)) for path in sheets]
    (output_dir / "visual_packet_manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
