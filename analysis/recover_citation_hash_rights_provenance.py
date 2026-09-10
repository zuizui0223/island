"""Recover normalized source families behind Database 1.0 citation lineages.

Rights/provenance audit only. This script never changes scientific trait values and
never grants redistribution rights. It resolves opaque ``citation:*`` provenance
against the immutable reviewed Public-Web evidence packet, then maps each citation
lineage to a stable source family suitable for the existing rights policy.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


def _tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value or "").split("|") if token.strip()]


def _slug(value: object) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", str(value or "").casefold()).strip("_")
    return text or "unknown"


def _family(provider: object, url: object) -> str:
    p = str(provider or "").strip()
    folded = p.casefold()
    u = str(url or "").strip()
    host = urlparse(u).netloc.casefold().removeprefix("www.") if u else ""
    if "dryad" in folded or "dryad" in host or "doi.org/10.5061/dryad" in u.casefold():
        return "dataset:dryad"
    if folded == "flora_of_australia_official":
        return "provider_treatment:flora_of_australia"
    if folded.startswith("sanbi_eflora_south_africa"):
        return "provider_treatment:sanbi_eflora_south_africa"
    if folded == "proteus":
        return "database:proteus"
    if folded == "gift_v3_2_direct":
        return "database:gift_v3_2"
    if host:
        return f"domain:{host}"
    return f"provider:{_slug(p)}"


def _load_reviewed(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"source_lineage", "source_provider", "source_url", "source_citation"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"reviewed Public-Web ledger missing columns: {sorted(missing)}")
    frame = frame.loc[frame["source_lineage"].str.startswith("citation:")].copy()
    frame["source_family"] = [
        _family(provider, url)
        for provider, url in zip(frame["source_provider"], frame["source_url"], strict=True)
    ]
    return frame


def recover(
    species_axis_path: Path,
    public_web_reviewed_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    if "source_lineages" not in coverage.columns:
        raise ValueError("species-axis ledger missing source_lineages")
    target = sorted(
        {
            token
            for value in coverage["source_lineages"]
            for token in _tokens(value)
            if token.startswith("citation:")
        }
    )
    if not target:
        raise ValueError("no citation lineages in Database 1.0")

    reviewed = _load_reviewed(public_web_reviewed_path)
    reviewed = reviewed.loc[reviewed["source_lineage"].isin(target)].copy()

    rows: list[dict[str, object]] = []
    for lineage in target:
        found = reviewed.loc[reviewed["source_lineage"].eq(lineage)].copy()
        families = sorted(set(found["source_family"]))
        providers = sorted(set(found["source_provider"]))
        urls = sorted({x for x in found["source_url"] if x})
        citations = sorted({x for x in found["source_citation"] if x})
        provenance_source = "reviewed_public_web" if len(found) else ""

        if lineage == "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497":
            families = ["publication:chung_etal_2000_heredity"]
            providers = ["reviewed_manual_literature"]
            urls = ["https://bsapubs.onlinelibrary.wiley.com/doi/10.3732/ajb.1000348"]
            citations = ["Chung et al. (2000), Heredity 85:490-497"]
            provenance_source = "committed_reviewed_provenance"
        elif lineage == "citation:Fischer-Rahelivololona-2007-Adansonia-29-269-315":
            families = ["publication:fischer_rahelivololona_2007_adansonia"]
            providers = ["MNHN Science Press"]
            urls = ["https://sciencepress.mnhn.fr/sites/default/files/articles/pdf/a2007n2a8.pdf"]
            citations = ["Fischer & Rahelivololona (2007), Adansonia 29(2):269-315"]
            provenance_source = "committed_reviewed_provenance"

        if len(families) > 1:
            raise ValueError(f"citation lineage maps to conflicting source families: {lineage} -> {families}")
        rows.append(
            {
                "source_lineage": lineage,
                "source_family": families[0] if families else "",
                "source_provenance_recovered": len(families) == 1,
                "source_providers": "|".join(providers),
                "source_urls": "|".join(urls),
                "source_citations": " || ".join(citations),
                "provenance_source": provenance_source,
            }
        )

    result = pd.DataFrame(rows)
    recovered = int(result["source_provenance_recovered"].sum())
    unresolved = result.loc[~result["source_provenance_recovered"]].copy()
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "CITATION_LINEAGE_FAMILY_MAP.csv", index=False)
    unresolved.to_csv(output_dir / "CITATION_UNRESOLVED.csv", index=False)

    provider_summary = (
        result.loc[result["source_provenance_recovered"]]
        .groupby("source_family", dropna=False)
        .agg(
            citation_lineages=("source_lineage", "nunique"),
            example_lineage=("source_lineage", "first"),
            example_provider=("source_providers", "first"),
        )
        .reset_index()
        .sort_values(["citation_lineages", "source_family"], ascending=[False, True])
    )
    provider_summary.to_csv(output_dir / "CITATION_SOURCE_FAMILY_SUMMARY.csv", index=False)

    summary = {
        "contract": "chapter1_database_v1_citation_lineage_provenance_audit_v1",
        "distinct_citation_lineages": int(len(result)),
        "source_provenance_recovered": recovered,
        "source_provenance_recovery_rate": recovered / len(result),
        "unresolved_citation_lineages": int(len(unresolved)),
        "distinct_recovered_source_families": int(provider_summary["source_family"].nunique()),
        "scientific_database_modified": False,
        "rights_granted": False,
    }
    (output_dir / "CITATION_PROVENANCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--public-web-reviewed", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    recover(args.species_axis, args.public_web_reviewed, args.output_dir)


if __name__ == "__main__":
    main()
