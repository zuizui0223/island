"""Prepare a fail-closed, high-yield WFO official-flora source package.

The package combines three already source-pinned acquisition lanes:

* species treatments parsed from the Flora Malesiana and Flora of the Guianas
  FlorML archives;
* new official WFO Darwin Core flora archives; and
* flower-colour rows recovered from the previously acquired WFO global-six
  package after a stricter organ/context gate.

This module does not promote evidence.  It removes already-direct
``accepted_species x trait_name`` cells, rejects known context ambiguities,
and produces a deterministic 200-row trait-stratified audit queue.  Promotion
remains the responsibility of :mod:`island_v2.open_web_finalize`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.integrate_wfo_bulk_round2_20260807 import COLOUR_WORDS, colour_gate
from island_v2.all_evidence_trait_audit import canonical_trait_name
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

SOURCE_RUN_ID = "wfo-high-yield-official-flora-20260810"
SOURCE_ARTIFACT = "wfo-high-yield-official-flora-source-package-20260810"
AUDIT_SIZE = 200
AUDIT_QUOTA = {
    "floral_form": 30,
    "flower_primary_color": 60,
    "flower_size_class": 40,
    "inflorescence_display": 70,
}
AUDIT_REPLACEMENT_SEED = "wfo-high-yield-v3-replacement"
PACKAGE_TRAITS = frozenset(AUDIT_QUOTA)

COLOUR_CONFOUNDER = re.compile(
    r"^\s*Fruit\s*:|"
    r"\b(?:glands?|spots?|trichomes?|hairs?|fruiting|immature|pirple|"
    r"crests?|appendages?)\b|"
    r"\b(?:when dry|drying|dried)\b|"
    r"\bcolour of the\b.*\bcolour of the\b",
    re.IGNORECASE,
)
BROWN_INDUMENTUM = re.compile(
    r"\bbrown\s+(?:pubescent|tomentose|hairy)|\bbrown\s+hairs?\b",
    re.IGNORECASE,
)
DISPLAY_PARTIAL_UNIT = re.compile(
    r"female flowers solitary or subtending|"
    r"partial peduncles?\s+2-flowered|"
    r"flowers solitary or(?: more rarely)? in twos",
    re.IGNORECASE,
)

DISPLAY_STATE_TERMS: dict[str, re.Pattern[str]] = {
    "raceme_spike_panicle": re.compile(
        r"\b(?:rac[eèé]m\w*|panic\w*|spic\w*|spike\w*|catkin\w*|ament\w*)\b",
        re.IGNORECASE,
    ),
    "umbel_corymb": re.compile(
        r"\b(?:cym\w*|umbel\w*|corymb\w*)\b",
        re.IGNORECASE,
    ),
    "solitary": re.compile(
        r"\b(?:flowers?\s+solitary|solitary\s+flowers?|fleurs?\s+solitaires?)\b",
        re.IGNORECASE,
    ),
    "composite_display": re.compile(
        r"\b(?:capitul\w*|flower heads?)\b",
        re.IGNORECASE,
    ),
    "few_flowered": re.compile(
        r"\b(?:few[- ]flowered|pauciflor\w*|[1-5]\s*[-–—]+\s*[1-5]"
        r"(?:\([^)]*\))?\s*[- ]*flowered|1\s+or\s+2(?:\(or\s+3\))?"
        r"\s*[- ]*flowered)\b",
        re.IGNORECASE,
    ),
}


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _canonical_direct_pairs(direct: pd.DataFrame) -> set[tuple[str, str]]:
    frame = direct.copy().fillna("")
    if "resolution_status" in frame:
        frame = frame.loc[frame["resolution_status"].eq("resolved")]
    frame["trait_name"] = frame["trait_name"].map(canonical_trait_name)
    return set(zip(frame["accepted_species"], frame["trait_name"], strict=True))


def _supported_colour_states(quote: str) -> set[str]:
    supported: set[str] = set()
    for state in COLOUR_WORDS:
        probe = pd.Series(
            {
                "exact_supporting_quote": quote,
                "normalized_value": state,
                "raw_value": state,
            }
        )
        if colour_gate(probe)[0]:
            supported.add(state)
    return supported


def safe_colour_row(row: pd.Series) -> bool:
    """Require complete, explicit floral colour support and no confounder."""

    quote = _text(row["source_excerpt"])
    requested = {part for part in _text(row["normalized_value"]).split("|") if part}
    if not requested or not requested.issubset(COLOUR_WORDS):
        return False
    if COLOUR_CONFOUNDER.search(quote):
        return False
    if (
        "green_brown_inconspicuous" in requested
        and BROWN_INDUMENTUM.search(quote)
    ):
        return False
    return _supported_colour_states(quote) == requested


def ambiguous_display_quote(quote: str) -> bool:
    """Reject sex-specific or partial-inflorescence unit counts."""

    if DISPLAY_PARTIAL_UNIT.search(quote):
        return True
    lower = quote.casefold()
    return all(
        marker in lower
        for marker in ("in male", "in female", "many-flowered", "few-flowered")
    )


def safe_display_row(row: pd.Series) -> bool:
    """Require every explicit display state in the quote to be retained.

    The acquisition regexes deliberately search broadly.  This final gate is
    fail-closed: if a quote explicitly names a second inflorescence state but
    the normalized state set omits it, the row stays out of production.  It
    prevents a valid mention such as ``capitula in corymbs`` from being
    flattened to only one of the two display states.
    """

    quote = _text(row["source_excerpt"])
    requested = {part for part in _text(row["normalized_value"]).split("|") if part}
    if not requested or ambiguous_display_quote(quote):
        return False
    supported = {
        state for state, pattern in DISPLAY_STATE_TERMS.items() if pattern.search(quote)
    }
    return bool(supported) and supported.issubset(requested)


def prepare_evidence(
    *,
    florml: pd.DataFrame,
    regional: pd.DataFrame,
    prior_global: pd.DataFrame,
    current_direct: pd.DataFrame,
) -> pd.DataFrame:
    """Combine the three acquisition lanes and apply the final strict gates."""

    direct_pairs = _canonical_direct_pairs(current_direct)
    frames: list[pd.DataFrame] = []
    for source, frame in (
        ("florml", florml),
        ("regional", regional),
        ("prior_global", prior_global),
    ):
        candidate = frame.copy().fillna("")
        missing = set(EVIDENCE_COLUMNS).difference(candidate.columns)
        if missing:
            raise ValueError(f"{source} evidence is missing columns: {sorted(missing)}")
        candidate["trait_name"] = candidate["trait_name"].map(canonical_trait_name)
        if source == "prior_global":
            candidate = candidate.loc[
                candidate["trait_name"].eq("flower_primary_color")
            ]
        frames.append(candidate.reindex(columns=EVIDENCE_COLUMNS))

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.loc[combined["trait_name"].isin(PACKAGE_TRAITS)].copy()
    combined = combined.loc[
        combined.apply(
            lambda row: (row["accepted_species"], row["trait_name"])
            not in direct_pairs,
            axis=1,
        )
    ].copy()

    colour = combined["trait_name"].eq("flower_primary_color")
    combined = combined.loc[
        ~colour | combined.apply(lambda row: safe_colour_row(row), axis=1)
    ].copy()
    display = combined["trait_name"].eq("inflorescence_display")
    combined = combined.loc[
        ~display
        | ~combined["source_excerpt"].map(_text).map(ambiguous_display_quote)
    ].copy()

    combined = combined.drop_duplicates("source_record_id")
    combined = combined.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"],
        keep="first",
    )
    combined = combined.sort_values(
        ["trait_name", "accepted_species", "source_lineage", "source_record_id"],
        kind="stable",
    )
    return combined.reset_index(drop=True)


def deterministic_audit_queue(
    evidence: pd.DataFrame,
    seed_audit: pd.DataFrame,
) -> pd.DataFrame:
    """Retain the frozen random sample and deterministically replace filtered rows."""

    if seed_audit["candidate_id"].duplicated().any():
        raise ValueError("seed audit contains duplicate candidate_id values")
    evidence_ids = set(evidence["source_record_id"])
    seed_ids = set(seed_audit["candidate_id"])
    kept = seed_audit.loc[seed_audit["candidate_id"].isin(evidence_ids)].copy()
    parts = [kept]
    for trait, quota in AUDIT_QUOTA.items():
        have = len(kept.loc[kept["trait_name"].eq(trait)])
        need = quota - have
        if need < 0:
            raise ValueError(f"seed audit exceeds the quota for {trait}")
        if not need:
            continue
        remainder = evidence.loc[
            evidence["trait_name"].eq(trait)
            & ~evidence["source_record_id"].isin(seed_ids)
        ].copy()
        remainder["_audit_order"] = remainder["source_record_id"].map(
            lambda value, trait_name=trait: _sha(
                f"{AUDIT_REPLACEMENT_SEED}|{trait_name}|{value}"
            )
        )
        if len(remainder) < need:
            raise ValueError(f"not enough replacement candidates for {trait}")
        replacement = remainder.sort_values("_audit_order", kind="stable").head(need)
        replacement = replacement.drop(columns="_audit_order").rename(
            columns={
                "source_provider": "provider",
                "source_record_id": "candidate_id",
                "source_excerpt": "exact_supporting_quote",
            }
        )
        for column in (
            "accepted_correct",
            "cultivar_status",
            "reviewer",
            "reviewed_at_utc",
            "audit_reason",
        ):
            replacement[column] = ""
        parts.append(replacement.reindex(columns=seed_audit.columns))

    audit = pd.concat(parts, ignore_index=True).sort_values(
        ["trait_name", "accepted_species", "candidate_id"], kind="stable"
    )
    if len(audit) != AUDIT_SIZE:
        raise ValueError(f"audit must contain {AUDIT_SIZE} rows, got {len(audit)}")
    return audit.reset_index(drop=True)


def run(args: argparse.Namespace) -> dict[str, Any]:
    evidence = prepare_evidence(
        florml=pd.read_csv(args.florml_evidence, dtype=str).fillna(""),
        regional=pd.read_csv(args.regional_evidence, dtype=str).fillna(""),
        prior_global=pd.read_csv(args.prior_global_evidence, dtype=str).fillna(""),
        current_direct=pd.read_csv(args.current_direct_ledger, dtype=str).fillna(""),
    )
    seed = pd.read_csv(args.audit_seed, dtype=str).fillna("")
    audit = deterministic_audit_queue(evidence, seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(
        args.output_dir / "wfo_high_yield_incremental_evidence.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    audit.to_csv(args.output_dir / "wfo_high_yield_audit_queue_200.csv", index=False)
    summary = {
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "candidate_rows": len(evidence),
        "species_traits": len(
            evidence[["accepted_species", "trait_name"]].drop_duplicates()
        ),
        "species_axes": len(
            evidence[["accepted_species", "axis"]].drop_duplicates()
        ),
        "species": evidence["accepted_species"].nunique(),
        "trait_rows": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "audit_rows": len(audit),
        "audit_contract": "accepted_correct / all_reviewed; production by trait",
        "family_inference": False,
        "global_fallback": False,
    }
    (args.output_dir / "wfo_high_yield_source_package_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--florml-evidence", type=Path, required=True)
    parser.add_argument("--regional-evidence", type=Path, required=True)
    parser.add_argument("--prior-global-evidence", type=Path, required=True)
    parser.add_argument("--current-direct-ledger", type=Path, required=True)
    parser.add_argument("--audit-seed", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2))
