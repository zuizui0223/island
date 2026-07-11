"""Extract reported and likely trait candidates directly from web search results.

Search engines often expose enough species-scoped text in result titles/snippets to
recover traits even when the destination page blocks automated retrieval. This module
therefore treats the frozen search-result record itself as an auditable evidence
carrier. It is intentionally source-agnostic: acceptance depends on focal-species
identity and the wording of the snippet, not on a site whitelist.

Two evidence classes are kept separate:

* ``reported``: the snippet explicitly states the trait for the focal species;
* ``likely``: the snippet supports a biologically reasonable inference but does not
  explicitly state the final ontology value.

No row is curated truth. Every row keeps the query, rank, result URL, snippet and the
inference rule so later review can reproduce or reject the decision.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

OUTPUT_COLUMNS = [
    "accepted_species",
    "trait_name",
    "provisional_candidate_value",
    "candidate_class",
    "inference_status",
    "confidence",
    "evidence_quote",
    "inference_rule",
    "searched_name",
    "search_query",
    "search_engine",
    "search_rank",
    "source_url",
    "source_title",
    "source_type",
    "evidence_scope",
    "review_status",
]

COLOR_TERMS = {
    "white": ("white", "whitish", "cream", "creamy", "ivory"),
    "green_brown_inconspicuous": ("green", "greenish", "brown", "brownish", "ochre"),
    "yellow_orange": ("yellow", "yellowish", "orange", "ochre-yellow", "golden"),
    "red_pink": ("red", "reddish", "pink", "pinkish", "rose", "crimson", "scarlet"),
    "blue_purple": ("blue", "bluish", "purple", "purplish", "violet", "lavender"),
}

FORM_TERMS = {
    "tubular": ("tubular", "corolla tube", "flower tube"),
    "bell_campanulate": ("bell-shaped", "bell shaped", "campanulate"),
    "funnel_trumpet": ("funnel-shaped", "funnel shaped", "trumpet-shaped", "trumpet shaped"),
    "salverform": ("salverform", "salver-shaped", "salver shaped"),
    "papilionaceous": ("papilionaceous", "pea-like flower", "pea-like flowers"),
    "spurred": ("spurred flower", "spurred flowers", "nectar spur"),
    "composite_head": ("capitulum", "capitula", "flower head", "flower heads"),
}

SYMMETRY_TERMS = {
    "actinomorphic": ("actinomorphic", "radially symmetric", "radial symmetry"),
    "zygomorphic": ("zygomorphic", "bilaterally symmetric", "bilateral symmetry", "two-lipped", "two lipped"),
}

DISPLAY_TERMS = {
    "raceme_spike_panicle": ("raceme", "racemes", "spike", "spikes", "panicle", "panicles"),
    "umbel_corymb": ("umbel", "umbels", "corymb", "corymbs"),
    "composite_display": ("capitulum", "capitula", "flower head", "flower heads"),
}

DIRECT_POLLINATION = {
    "abiotic_wind": (
        "wind-pollinated",
        "wind pollinated",
        "pollinated by wind",
        "pollination by wind",
        "pollination is made by wind",
        "anemogamy",
    ),
    "biotic": (
        "insect-pollinated",
        "insect pollinated",
        "pollinated by insects",
        "pollination by insects",
        "pollination is made by insects",
        "entomogamy",
    ),
}

GUILD_TERMS = {
    "bumblebees": ("bumblebee", "bumblebees"),
    "other_bees": ("bee-pollinated", "bee pollinated", "pollinated by bees", "bees"),
    "birds": ("bird-pollinated", "bird pollinated", "pollinated by birds", "hummingbird", "hummingbirds"),
    "butterflies": ("butterfly", "butterflies"),
    "moths": ("moth", "moths"),
    "flies": ("fly-pollinated", "fly pollinated", "pollinated by flies"),
    "bats": ("bat-pollinated", "bat pollinated", "pollinated by bats"),
}

WIND_LIKE_FAMILIES = {"Poaceae", "Cyperaceae", "Juncaceae"}
BIOTIC_LIKELY_FAMILIES = {
    "Orchidaceae",
    "Lamiaceae",
    "Gesneriaceae",
    "Campanulaceae",
    "Apocynaceae",
    "Rubiaceae",
    "Ericaceae",
    "Myrtaceae",
    "Fabaceae",
    "Malvaceae",
}


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def _binomial(name: str) -> str:
    parts = _text(name).casefold().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else " ".join(parts)


def _species_scoped(text: str, accepted: str, searched: str) -> bool:
    low = _text(text).casefold()
    names = {_binomial(accepted), _binomial(searched)}
    return any(name and name in low for name in names)


def _sentences(text: str) -> list[str]:
    return [part.strip() for part in re.split(r"(?<=[.!?;])\s+", _text(text)) if part.strip()]


def _flower_context(sentence: str) -> bool:
    low = sentence.casefold()
    return any(term in low for term in ("flower", "flowers", "floral", "corolla", "bloom", "inflorescence"))


def _find_terms(sentence: str, mapping: dict[str, tuple[str, ...]]) -> list[tuple[str, str]]:
    low = sentence.casefold()
    found: list[tuple[str, str]] = []
    for value, terms in mapping.items():
        for term in terms:
            if re.search(rf"(?<!\w){re.escape(term)}(?!\w)", low):
                found.append((value, term))
                break
    return found


def _row(
    *,
    base: dict[str, str],
    trait_name: str,
    value: str,
    candidate_class: str,
    confidence: str,
    quote: str,
    rule: str,
) -> dict[str, str]:
    return {
        "accepted_species": base["accepted_species"],
        "trait_name": trait_name,
        "provisional_candidate_value": value,
        "candidate_class": candidate_class,
        "inference_status": "reported" if candidate_class == "reported" else "likely",
        "confidence": confidence,
        "evidence_quote": quote,
        "inference_rule": rule,
        "searched_name": base["searched_name"],
        "search_query": base["search_query"],
        "search_engine": base["search_engine"],
        "search_rank": base["search_rank"],
        "source_url": base["source_url"],
        "source_title": base["source_title"],
        "source_type": "search_result_snippet",
        "evidence_scope": "species_direct_search_result",
        "review_status": "unreviewed_search_result_candidate",
    }


def extract_search_result_candidates(
    search_results: pd.DataFrame,
    species_table: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {
        "accepted_species",
        "searched_name",
        "search_query",
        "search_engine",
        "search_rank",
        "result_title",
        "result_url",
        "result_snippet",
    }
    missing = required.difference(search_results.columns)
    if missing:
        raise typer.BadParameter(f"search result table missing columns: {sorted(missing)}")
    if "accepted_species" not in species_table.columns:
        raise typer.BadParameter("species table missing accepted_species")

    family_by_species = {
        _text(row.get("accepted_species")): _text(row.get("family"))
        for row in species_table.fillna("").to_dict("records")
    }

    rows: list[dict[str, str]] = []
    species_with_scoped_result: set[str] = set()
    for record in search_results.fillna("").to_dict("records"):
        accepted = _text(record.get("accepted_species"))
        searched = _text(record.get("searched_name"))
        title = _text(record.get("result_title"))
        snippet = _text(record.get("result_snippet"))
        combined = f"{title}. {snippet}".strip()
        if not accepted or not _species_scoped(combined, accepted, searched):
            continue
        species_with_scoped_result.add(accepted)
        base = {
            "accepted_species": accepted,
            "searched_name": searched,
            "search_query": _text(record.get("search_query")),
            "search_engine": _text(record.get("search_engine")),
            "search_rank": _text(record.get("search_rank")),
            "source_url": _text(record.get("result_url")),
            "source_title": title,
        }

        for sentence in _sentences(combined):
            if not _flower_context(sentence):
                continue
            color_hits = _find_terms(sentence, COLOR_TERMS)
            color_values = {value for value, _ in color_hits}
            if len(color_values) == 1:
                value = next(iter(color_values))
                rows.append(
                    _row(
                        base=base,
                        trait_name="flower_primary_color",
                        value=value,
                        candidate_class="reported",
                        confidence="medium",
                        quote=sentence,
                        rule="explicit_colour_term_in_flower_context",
                    )
                )
            elif len(color_values) > 1:
                rows.append(
                    _row(
                        base=base,
                        trait_name="flower_primary_color",
                        value="multicolored_variable",
                        candidate_class="reported",
                        confidence="medium",
                        quote=sentence,
                        rule="multiple_colour_terms_in_flower_context",
                    )
                )

            for value, _ in _find_terms(sentence, FORM_TERMS):
                rows.append(
                    _row(
                        base=base,
                        trait_name="floral_form",
                        value=value,
                        candidate_class="reported",
                        confidence="medium",
                        quote=sentence,
                        rule="explicit_form_term_in_flower_context",
                    )
                )
            for value, _ in _find_terms(sentence, SYMMETRY_TERMS):
                rows.append(
                    _row(
                        base=base,
                        trait_name="floral_symmetry",
                        value=value,
                        candidate_class="reported",
                        confidence="medium",
                        quote=sentence,
                        rule="explicit_symmetry_term_in_flower_context",
                    )
                )
            for value, _ in _find_terms(sentence, DISPLAY_TERMS):
                rows.append(
                    _row(
                        base=base,
                        trait_name="inflorescence_display",
                        value=value,
                        candidate_class="reported",
                        confidence="medium",
                        quote=sentence,
                        rule="explicit_inflorescence_term_in_flower_context",
                    )
                )

        low = combined.casefold()
        for value, terms in DIRECT_POLLINATION.items():
            if any(term in low for term in terms):
                rows.append(
                    _row(
                        base=base,
                        trait_name="pollen_vector_mode",
                        value=value,
                        candidate_class="reported",
                        confidence="high",
                        quote=combined,
                        rule="explicit_pollination_vector_phrase",
                    )
                )
        for value, terms in GUILD_TERMS.items():
            if any(term in low for term in terms) and "pollinat" in low:
                rows.append(
                    _row(
                        base=base,
                        trait_name="pollination_functional_guild",
                        value=value,
                        candidate_class="reported",
                        confidence="medium",
                        quote=combined,
                        rule="explicit_pollinator_guild_phrase",
                    )
                )
        if "attracts pollinators" in low or "pollinator-friendly" in low or "pollinator friendly" in low:
            rows.append(
                _row(
                    base=base,
                    trait_name="pollen_vector_mode",
                    value="biotic",
                    candidate_class="likely",
                    confidence="low",
                    quote=combined,
                    rule="likely_biotic_from_attracts_pollinators_wording",
                )
            )

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(
            [
                "accepted_species",
                "trait_name",
                "provisional_candidate_value",
                "candidate_class",
                "source_url",
                "evidence_quote",
            ]
        )

    # Family priors are deliberately weak and only fill pollen-vector gaps. They are
    # useful for sensitivity analyses, never for the reported-evidence analysis.
    existing_vector_species = set(
        frame.loc[frame["trait_name"].eq("pollen_vector_mode"), "accepted_species"]
    ) if not frame.empty else set()
    prior_rows: list[dict[str, str]] = []
    for species, family in family_by_species.items():
        if not species or species in existing_vector_species:
            continue
        if family in WIND_LIKE_FAMILIES:
            value = "abiotic_wind"
            rule = "likely_wind_pollination_from_declared_family_prior"
        elif family in BIOTIC_LIKELY_FAMILIES:
            value = "biotic"
            rule = "likely_biotic_pollination_from_declared_family_prior"
        else:
            continue
        prior_rows.append(
            {
                "accepted_species": species,
                "trait_name": "pollen_vector_mode",
                "provisional_candidate_value": value,
                "candidate_class": "proxy",
                "inference_status": "likely",
                "confidence": "low",
                "evidence_quote": f"family={family}",
                "inference_rule": rule,
                "searched_name": species,
                "search_query": "",
                "search_engine": "",
                "search_rank": "",
                "source_url": "",
                "source_title": "",
                "source_type": "declared_family_prior",
                "evidence_scope": "proxy_not_reported",
                "review_status": "unreviewed_proxy",
            }
        )
    if prior_rows:
        frame = pd.concat([frame, pd.DataFrame(prior_rows, columns=OUTPUT_COLUMNS)], ignore_index=True)

    reported = frame.loc[frame["candidate_class"].eq("reported")] if not frame.empty else frame
    likely = frame.loc[frame["inference_status"].eq("likely")] if not frame.empty else frame
    report = {
        "n_search_result_rows": int(len(search_results)),
        "n_species_with_species_scoped_search_result": len(species_with_scoped_result),
        "n_candidate_rows": int(len(frame)),
        "n_species_with_any_candidate": int(frame["accepted_species"].nunique()) if len(frame) else 0,
        "n_reported_rows": int(len(reported)),
        "n_likely_rows": int(len(likely)),
        "n_species_reported": int(reported["accepted_species"].nunique()) if len(reported) else 0,
        "n_species_likely": int(likely["accepted_species"].nunique()) if len(likely) else 0,
        "coverage_by_trait": {
            str(k): int(v)
            for k, v in frame.groupby("trait_name")["accepted_species"].nunique().sort_index().items()
        }
        if len(frame)
        else {},
        "policy": (
            "Search-result snippets are accepted as frozen evidence carriers when the focal binomial "
            "occurs in the title/snippet. Explicit statements are reported candidates; biological "
            "inferences and family priors are marked likely/proxy and kept out of reported-only analyses."
        ),
    }
    return frame, report


@app.command()
def extract(
    search_results_csv: Path = typer.Option(..., exists=True),
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    search_results = pd.read_csv(search_results_csv, dtype=str).fillna("")
    species_table = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = extract_search_result_candidates(search_results, species_table)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "search_result_trait_candidates.csv", index=False)
    (output_dir / "search_result_trait_candidates_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
