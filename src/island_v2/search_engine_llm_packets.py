"""Build reproducible LLM packets from search-engine result context.

The search engine is used only to collect a frozen context window. Trait assignment is
intentionally deferred to an LLM/human extractor using the fixed 9-column schema.
"""

from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx
import pandas as pd
import typer

from island_v2.search_engine_trait_recovery import clean, search_with_fallback

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERY_TEMPLATES = (
    '"{name}" flower color shape pollination',
    '"{name}" pollinator self-compatible self-incompatible selfing outcrossing',
    '"{name}" floral morphology breeding system mating system',
)

PROMPT_VERSION = "v1_style_search_context_v1"

PROMPT = """You are an expert in plant reproductive biology and pollination ecology.

Using only the frozen search-result context in this packet, infer the following traits for the focal species:
- typical flower color(s)
- flower shape / symmetry
- main pollination guild (bees, bumblebees, flies, butterflies, moths, birds, wind, self, mixed, unknown)
- mating system (obligate_outcrossing, mixed_mating, mainly_selfing, obligate_selfing, unknown)
- self-incompatibility (SI, SC, likely_SI, likely_SC, unknown)
- short free-text notes on pollination biology
- evidence type (field_study, review, flora, horticulture, inference)
- confidence (high, medium, low)

Rules:
1. Return exactly one CSV row and no prose.
2. Column order must be:
species,flower_color,flower_shape,pollination_guild,pollination_notes,mating_system,self_incompatibility,evidence_type,confidence
3. Use direct species-level wording when available.
4. If species-level evidence is incomplete, reasonable biological inference from genus/family-level or convergent search-result context is allowed, but mark evidence_type=inference and use likely_SI/likely_SC where appropriate.
5. Do not equate self-compatibility with autonomous selfing.
6. If information is insufficient, use unknown and confidence=low.
7. Fields containing commas must be enclosed in double quotes.
"""


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def collect_species_context(species: str, results_per_query: int) -> dict[str, object]:
    headers = {
        "User-Agent": "island-floral-v2/0.1 search-engine-llm-packets",
        "Accept-Language": "en;q=0.9,*;q=0.5",
    }
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    with httpx.Client(headers=headers, follow_redirects=True) as client:
        for template in QUERY_TEMPLATES:
            query = template.format(name=species)
            for rank, (engine, result) in enumerate(
                search_with_fallback(client, query, results_per_query), start=1
            ):
                url = clean(result.get("url"))
                if not url or url in seen:
                    continue
                seen.add(url)
                rows.append(
                    {
                        "query": query,
                        "engine": engine,
                        "rank": str(rank),
                        "title": clean(result.get("title")),
                        "snippet": clean(result.get("snippet")),
                        "url": url,
                    }
                )
    return {
        "species": species,
        "prompt_version": PROMPT_VERSION,
        "prompt": PROMPT,
        "results": rows,
    }


@app.command()
def build(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=100000),
    results_per_query: int = typer.Option(6, min=1, max=20),
    workers: int = typer.Option(24, min=1, max=64),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species_list = [
        clean(value)
        for value in table["accepted_species"].astype(str)
        if clean(value)
    ][:max_taxa]

    output_dir.mkdir(parents=True, exist_ok=True)
    packet_dir = output_dir / "packets"
    packet_dir.mkdir(parents=True, exist_ok=True)

    packets: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(collect_species_context, species, results_per_query): species
            for species in species_list
        }
        for future in as_completed(futures):
            species = futures[future]
            try:
                packet = future.result()
                packets.append(packet)
            except Exception as exc:
                errors.append({"species": species, "error": str(exc)})
                packets.append(
                    {
                        "species": species,
                        "prompt_version": PROMPT_VERSION,
                        "prompt": PROMPT,
                        "results": [],
                    }
                )

    order = {species: index for index, species in enumerate(species_list)}
    packets.sort(key=lambda row: order[str(row["species"])])

    index_rows: list[dict[str, object]] = []
    jsonl_lines: list[str] = []
    for i, packet in enumerate(packets, start=1):
        species = str(packet["species"])
        text = json.dumps(packet, ensure_ascii=False, sort_keys=True)
        packet_id = f"species_{i:05d}_{sha256_text(species)[:10]}"
        path = packet_dir / f"{packet_id}.json"
        path.write_text(json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        jsonl_lines.append(text)
        index_rows.append(
            {
                "packet_id": packet_id,
                "species": species,
                "n_results": len(packet.get("results", [])),
                "packet_sha256": sha256_text(text),
                "packet_path": str(path),
            }
        )

    (output_dir / "search_context_packets.jsonl").write_text(
        "\n".join(jsonl_lines) + ("\n" if jsonl_lines else ""), encoding="utf-8"
    )
    pd.DataFrame(index_rows).to_csv(output_dir / "packet_index.csv", index=False)
    pd.DataFrame(errors, columns=["species", "error"]).to_csv(
        output_dir / "packet_errors.csv", index=False
    )
    report = {
        "n_species": len(species_list),
        "n_packets": len(packets),
        "n_species_with_any_search_result": sum(
            1 for packet in packets if packet.get("results")
        ),
        "n_total_search_results": sum(
            len(packet.get("results", [])) for packet in packets
        ),
        "n_errors": len(errors),
        "queries_per_species": len(QUERY_TEMPLATES),
        "prompt_version": PROMPT_VERSION,
        "policy": "freeze broad search-result context first; assign the 9-column trait row semantically afterward",
    }
    (output_dir / "search_context_packet_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
