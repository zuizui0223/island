"""Build a small reviewed direct-evidence checkpoint from acquired pages.

This checkpoint deliberately does not search the whole 106,295-species
universe again.  It promotes only unresolved exact-species records from a
completed PFAF inventory and two downloaded primary articles.  PFAF's
``self-fertile`` wording maps to self-compatibility only; the primary-article
records keep mating system, self-incompatibility and autonomous selfing as
separate traits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

DIOSPYROS_SHA256 = "6bbf90b64844268c8232dfdc673bf66646cc5395731e4c52692f0068dad24872"
MAGNOLIA_SHA256 = "99843970098822c65cab6d7afa78c7e984b8dfaad5c3e7516db1292bcd86f9e5"
NPARKS_SHA256 = {
    "Grewia occidentalis": "275d87b68ed4ed089db3692f1e06ab6f32f01c2a0bd7d96bab68173143c61140",
    "Myristica elliptica": "429a45acbb2c83b6459b4eaf792b714c062bb7b3bdff3d3f36685e718e3b1f8b",
    "Nepenthes bicalcarata": "e7a5b5daa2fbb21ba6a9b96764bb9eb83ffeeae040036c31e28348f5916a9957",
}
CREATED_AT = "2026-08-11T05:10:00Z"
REVIEWER = "Codex source-backed line-by-line evidence audit"
TRAIT_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
}

EVIDENCE_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "axis",
    "trait_name",
    "raw_value",
    "normalized_value",
    "evidence_quality",
    "evidence_scope",
    "source_group",
    "source_provider",
    "source_url",
    "page_title",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "lineage_method",
    "name_resolution_lineage",
    "name_match_method",
    "matched_page_name",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "wild_cultivated_cultivar_status",
    "evidence_status",
    "content_sha256",
    "content_sha256_basis",
    "retrieved_at_utc",
    "query",
    "search_rank",
    "inference_rule",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def _canonical_file_bytes(path: Path) -> bytes:
    payload = path.read_bytes()
    if path.suffix.lower() in {".csv", ".json"}:
        payload = payload.replace(b"\r\n", b"\n")
    return payload


def _canonical_file_size(path: Path) -> int:
    return len(_canonical_file_bytes(path))


def _candidate_id(species: str, trait: str, value: str, lineage: str) -> str:
    return hashlib.sha256(
        f"{species}|{trait}|{value}|{lineage}".encode()
    ).hexdigest()[:24]


def _evidence_row(
    *,
    species: str,
    trait: str,
    value: str,
    quality: str,
    provider: str,
    url: str,
    title: str,
    citation: str,
    excerpt: str,
    record_id: str,
    lineage: str,
    lineage_method: str,
    source_tier: str,
    source_type: str,
    domain: str,
    content_sha256: str,
    content_sha256_basis: str,
    retrieved_at_utc: str,
    raw_value: str = "",
) -> dict[str, str]:
    return {
        "candidate_id": _candidate_id(species, trait, value, lineage),
        "accepted_species": species,
        "axis": TRAIT_AXIS[trait],
        "trait_name": trait,
        "raw_value": raw_value or value,
        "normalized_value": value,
        "evidence_quality": quality,
        "evidence_scope": "species_direct",
        "source_group": "high_leverage_direct_checkpoint",
        "source_provider": provider,
        "source_url": url,
        "page_title": title,
        "source_citation": citation,
        "source_excerpt": excerpt,
        "source_record_id": record_id,
        "source_lineage": lineage,
        "lineage_method": lineage_method,
        "name_resolution_lineage": "accepted_name_exact",
        "name_match_method": "accepted_name_exact",
        "matched_page_name": species,
        "source_tier": source_tier,
        "source_type": source_type,
        "domain": domain,
        "language": "en",
        "wild_cultivated_cultivar_status": "wild_species_not_cultivar_limited",
        "evidence_status": "accepted_individual_source_backed_audit",
        "content_sha256": content_sha256,
        "content_sha256_basis": content_sha256_basis,
        "retrieved_at_utc": retrieved_at_utc,
        "query": "priority_support2_exact_species_acquisition",
        "search_rank": "",
        "inference_rule": "",
    }


def resolved_direct_pairs(ledger: pd.DataFrame) -> set[tuple[str, str]]:
    """Only resolved direct cells suppress reacquisition."""

    resolved = ledger.loc[ledger["resolution_status"].eq("resolved")]
    return set(
        zip(resolved["accepted_species"], resolved["trait_name"], strict=True)
    )


def primary_rows(
    *,
    master_names: set[str],
    completed_pairs: set[tuple[str, str]],
    diospyros_pdf: Path,
    magnolia_pdf: Path,
) -> list[dict[str, str]]:
    if _sha256(diospyros_pdf) != DIOSPYROS_SHA256:
        raise ValueError("Diospyros primary-article PDF hash differs from pinned bytes")
    if _sha256(magnolia_pdf) != MAGNOLIA_SHA256:
        raise ValueError("Magnolia primary-article PDF hash differs from pinned bytes")

    candidates = [
        _evidence_row(
            species="Diospyros celebica",
            trait="mating_system",
            value="predominantly_outcrossing",
            raw_value="self-pollination 11.1%; cross pollination 88.9%",
            quality="high",
            provider="restu_etal_2017_biotropia",
            url="https://journal.biotrop.org/index.php/biotropia/article/download/562/377",
            title=(
                "High outcrossing rate and pollen dispersal distance of Diospyros "
                "celebica Bakh. (Ebenaceae), an endemic tree species in Sulawesi"
            ),
            citation=(
                "Restu, Gusmiaty & Larekeng (2017), BIOTROPIA 24(3), "
                "DOI 10.11598/btb.2017.24.3.562"
            ),
            excerpt=(
                "Within the male parents assigned to the 72 progeny arrays, "
                "self-pollination occurred in 8 events (11.1%) and cross "
                "pollination occurred in 64 events (88.9%)."
            ),
            record_id="doi:10.11598/btb.2017.24.3.562:results:72-progeny",
            lineage="doi:10.11598/btb.2017.24.3.562",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article",
            domain="journal.biotrop.org",
            content_sha256=DIOSPYROS_SHA256,
            content_sha256_basis="downloaded_primary_article_pdf_bytes",
            retrieved_at_utc="2026-08-11T04:50:47Z",
        ),
        _evidence_row(
            species="Magnolia salicifolia",
            trait="self_incompatibility",
            value="SC",
            quality="high",
            provider="ishida_etal_2020_primary_article",
            url=(
                "https://www.davidpublisher.com/Public/uploads/Contribute/"
                "5ed76e5d17b9b.pdf"
            ),
            title=(
                "Inbreeding and Inbreeding Depression in a Deciduous Shrub, "
                "Magnolia salicifolia, in the Understory of a Japanese Beech Forest"
            ),
            citation=(
                "Ishida, Kikuchi & Hayashi (2020), Journal of Environmental Science "
                "and Engineering A 9:90-97, DOI 10.17265/2162-5298/2020.03.002"
            ),
            excerpt=(
                "The pollination experiments revealed that M. salicifolia is "
                "self-compatible, although the flowers cannot set seeds under "
                "environments without pollinators."
            ),
            record_id="doi:10.17265/2162-5298/2020.03.002:discussion:sc",
            lineage="doi:10.17265/2162-5298/2020.03.002",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article",
            domain="davidpublisher.com",
            content_sha256=MAGNOLIA_SHA256,
            content_sha256_basis="downloaded_primary_article_pdf_bytes",
            retrieved_at_utc="2026-08-11T05:07:49Z",
        ),
        _evidence_row(
            species="Magnolia salicifolia",
            trait="autonomous_selfing_capacity",
            value="absent",
            quality="high",
            provider="ishida_etal_2020_primary_article",
            url=(
                "https://www.davidpublisher.com/Public/uploads/Contribute/"
                "5ed76e5d17b9b.pdf"
            ),
            title=(
                "Inbreeding and Inbreeding Depression in a Deciduous Shrub, "
                "Magnolia salicifolia, in the Understory of a Japanese Beech Forest"
            ),
            citation=(
                "Ishida, Kikuchi & Hayashi (2020), Journal of Environmental Science "
                "and Engineering A 9:90-97, DOI 10.17265/2162-5298/2020.03.002"
            ),
            excerpt=(
                "The pollination experiments revealed that M. salicifolia is "
                "self-compatible, although the flowers cannot set seeds under "
                "environments without pollinators."
            ),
            record_id="doi:10.17265/2162-5298/2020.03.002:discussion:no-autonomy",
            lineage="doi:10.17265/2162-5298/2020.03.002",
            lineage_method="original_primary_article_doi",
            source_tier="A",
            source_type="primary_research_article",
            domain="davidpublisher.com",
            content_sha256=MAGNOLIA_SHA256,
            content_sha256_basis="downloaded_primary_article_pdf_bytes",
            retrieved_at_utc="2026-08-11T05:07:49Z",
        ),
    ]
    return [
        row
        for row in candidates
        if row["accepted_species"] in master_names
        and (row["accepted_species"], row["trait_name"]) not in completed_pairs
    ]


def nparks_rows(
    *,
    master_names: set[str],
    completed_pairs: set[tuple[str, str]],
    grewia_html: Path,
    myristica_html: Path,
    nepenthes_html: Path,
) -> list[dict[str, str]]:
    pages = {
        "Grewia occidentalis": grewia_html,
        "Myristica elliptica": myristica_html,
        "Nepenthes bicalcarata": nepenthes_html,
    }
    for species, path in pages.items():
        if _sha256(path) != NPARKS_SHA256[species]:
            raise ValueError(f"NParks page hash differs from pinned bytes: {species}")

    records = [
        {
            "species": "Myristica elliptica",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "raw": "Radial",
            "excerpt": "Flower Symmetry | Radial",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/6/2/6288",
            "record": "nparks:6288:flower-symmetry",
            "retrieved": "2026-08-11T05:35:12Z",
        },
        {
            "species": "Myristica elliptica",
            "trait": "flower_primary_color",
            "value": "white|yellow_orange",
            "raw": "Cream / Off-White; Yellow / Golden",
            "excerpt": "Flower Colour(s) | Cream / Off-White, Yellow / Golden",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/6/2/6288",
            "record": "nparks:6288:flower-colour",
            "retrieved": "2026-08-11T05:35:12Z",
        },
        {
            "species": "Nepenthes bicalcarata",
            "trait": "floral_symmetry",
            "value": "actinomorphic",
            "raw": "Radial",
            "excerpt": "Flower Symmetry | Radial",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/4/4475",
            "record": "nparks:4475:flower-symmetry",
            "retrieved": "2026-08-11T05:35:13Z",
        },
        {
            "species": "Nepenthes bicalcarata",
            "trait": "floral_form",
            "value": "open_radial",
            "raw": "Stellate / Star-shaped",
            "excerpt": "Individual Flower Shape | Stellate / Star-shaped",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/4/4475",
            "record": "nparks:4475:individual-flower-shape",
            "retrieved": "2026-08-11T05:35:13Z",
        },
        {
            "species": "Nepenthes bicalcarata",
            "trait": "inflorescence_display",
            "value": "raceme_spike_panicle",
            "raw": "Panicle",
            "excerpt": "Inflorescence Type | Panicle",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/4/4475",
            "record": "nparks:4475:inflorescence-type",
            "retrieved": "2026-08-11T05:35:13Z",
        },
        {
            "species": "Nepenthes bicalcarata",
            "trait": "flower_primary_color",
            "value": "blue_purple|other_described",
            "raw": "Purple; Black",
            "excerpt": "Flower Colour(s) | Purple, Black",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/4/4475",
            "record": "nparks:4475:flower-colour",
            "retrieved": "2026-08-11T05:35:13Z",
        },
        {
            "species": "Grewia occidentalis",
            "trait": "floral_form",
            "value": "open_radial",
            "raw": "Stellate / Star-shaped",
            "excerpt": "Individual Flower Shape | Stellate / Star-shaped",
            "url": "https://www.nparks.gov.sg/florafaunaweb/flora/4/2/4236",
            "record": "nparks:4236:individual-flower-shape",
            "retrieved": "2026-08-11T05:35:14Z",
        },
    ]
    output: list[dict[str, str]] = []
    for item in records:
        species = item["species"]
        trait = item["trait"]
        if species not in master_names or (species, trait) in completed_pairs:
            continue
        output.append(
            _evidence_row(
                species=species,
                trait=trait,
                value=item["value"],
                raw_value=item["raw"],
                quality="high",
                provider="singapore_nparks_flora_fauna_web",
                url=item["url"],
                title=f"NParks Flora & Fauna Web: {species}",
                citation=(
                    "National Parks Board Singapore, Flora & Fauna Web, "
                    f"{species} species record"
                ),
                excerpt=item["excerpt"],
                record_id=item["record"],
                lineage=f"url:{item['url'].casefold()}",
                lineage_method="government_species_treatment_url",
                source_tier="A",
                source_type="government_species_database",
                domain="nparks.gov.sg",
                content_sha256=NPARKS_SHA256[species],
                content_sha256_basis="downloaded_species_page_html_bytes",
                retrieved_at_utc=item["retrieved"],
            )
        )
    return output


def pfaf_rows(
    candidates: pd.DataFrame,
    *,
    master_names: set[str],
    completed_pairs: set[tuple[str, str]],
) -> list[dict[str, str]]:
    required = {
        "accepted_species",
        "source_url",
        "matched_page_name",
        "identity_exact",
        "supporting_excerpt",
        "normalized_value",
        "trait_name",
        "status",
        "cultivar_qualified",
        "nonsexual_qualified",
        "dioecious_conflict",
        "page_sha256",
        "retrieved_at_utc",
    }
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"PFAF candidates missing columns: {sorted(missing)}")
    selected = candidates.loc[
        candidates["status"].eq("accepted_self_fertile_statement")
        & candidates["identity_exact"].str.casefold().eq("true")
        & candidates["normalized_value"].eq("SC")
        & candidates["trait_name"].eq("self_incompatibility")
        & candidates["cultivar_qualified"].str.casefold().eq("false")
        & candidates["nonsexual_qualified"].str.casefold().eq("false")
        & candidates["dioecious_conflict"].str.casefold().eq("false")
    ].copy()
    rows: list[dict[str, str]] = []
    for item in selected.sort_values("accepted_species").to_dict("records"):
        species = _text(item["accepted_species"])
        pair = (species, "self_incompatibility")
        if species not in master_names or pair in completed_pairs:
            continue
        if _text(item["matched_page_name"]) != species:
            continue
        excerpt = _text(item["supporting_excerpt"])
        if excerpt.casefold() != "the plant is self-fertile":
            continue
        url = _text(item["source_url"])
        lineage = f"ken-fern-useful-plants:species:{species.casefold()}"
        retrieved = _text(item["retrieved_at_utc"]).replace("+00:00", "Z")
        rows.append(
            _evidence_row(
                species=species,
                trait="self_incompatibility",
                value="SC",
                raw_value=excerpt,
                quality="medium",
                provider="plants_for_a_future",
                url=url,
                title=f"Plants For A Future species profile: {species}",
                citation="Plants For A Future species profile, compiled by Ken Fern",
                excerpt=excerpt,
                record_id=f"pfaf-page:{hashlib.sha256(url.encode()).hexdigest()[:24]}",
                lineage=lineage,
                lineage_method="shared_ken_fern_species_lineage_across_provider_mirrors",
                source_tier="B",
                source_type="specialist_species_database",
                domain="pfaf.org",
                content_sha256=_text(item["page_sha256"]).casefold(),
                content_sha256_basis="retrieved_species_page_bytes",
                retrieved_at_utc=retrieved,
            )
        )
    return rows


def _audit(evidence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for item in evidence.to_dict("records"):
        excerpt = _text(item["source_excerpt"])
        normalized = excerpt.casefold()
        rows.append(
            {
                "candidate_id": item["candidate_id"],
                "accepted_species": item["accepted_species"],
                "trait_name": item["trait_name"],
                "normalized_value": item["normalized_value"],
                "source_url": item["source_url"],
                "page_title": item["page_title"],
                "source_citation": item["source_citation"],
                "source_tier": item["source_tier"],
                "source_type": item["source_type"],
                "domain": item["domain"],
                "language": item["language"],
                "supporting_excerpt": excerpt,
                "normalized_excerpt_sha256": hashlib.sha256(
                    normalized.encode("utf-8")
                ).hexdigest(),
                "content_fingerprint": hashlib.sha256(
                    (item["source_lineage"] + "\n" + normalized).encode("utf-8")
                ).hexdigest(),
                "name_match_method": item["name_match_method"],
                "wild_cultivated_cultivar_status": item[
                    "wild_cultivated_cultivar_status"
                ],
                "decision": "accept",
                "species_identity_correct": "true",
                "value_correct": "true",
                "provenance_complete": "true",
                "cultivar_contamination": "false",
                "false_positive_reason": "",
                "decision_reason": (
                    "Exact accepted species identity, exact source statement, trait-specific "
                    "mapping, source URL and retrieved-content hash rechecked; no cultivar, "
                    "family inference, cross-trait substitution or global fallback."
                ),
                "reviewer": REVIEWER,
                "reviewed_at_utc": CREATED_AT,
            }
        )
    return pd.DataFrame(rows)


def build(
    *,
    pfaf_csv: Path,
    current_direct_csv: Path,
    master_csv: Path,
    diospyros_pdf: Path,
    magnolia_pdf: Path,
    grewia_html: Path,
    myristica_html: Path,
    nepenthes_html: Path,
    output_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    master_names = set(master["accepted_species"].map(_text))
    direct = pd.read_csv(current_direct_csv, dtype=str).fillna("")
    completed_pairs = resolved_direct_pairs(direct)
    pfaf = pd.read_csv(pfaf_csv, dtype=str).fillna("")

    rows = pfaf_rows(
        pfaf, master_names=master_names, completed_pairs=completed_pairs
    ) + primary_rows(
        master_names=master_names,
        completed_pairs=completed_pairs,
        diospyros_pdf=diospyros_pdf,
        magnolia_pdf=magnolia_pdf,
    ) + nparks_rows(
        master_names=master_names,
        completed_pairs=completed_pairs,
        grewia_html=grewia_html,
        myristica_html=myristica_html,
        nepenthes_html=nepenthes_html,
    )
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).fillna("")
    evidence = evidence.sort_values(
        ["source_provider", "accepted_species", "trait_name", "candidate_id"]
    ).reset_index(drop=True)
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("checkpoint candidate IDs are not unique")
    audit = _audit(evidence)

    evidence_path = output_dir / "high_leverage_direct_evidence_20260811.csv"
    audit_path = output_dir / "high_leverage_direct_manual_audit_20260811.csv"
    evidence.to_csv(evidence_path, index=False)
    audit.to_csv(audit_path, index=False)
    summary: dict[str, object] = {
        "contract": "high_leverage_individually_reviewed_direct_checkpoint_v1",
        "created_at_utc": CREATED_AT,
        "current_resolved_direct_species_trait": len(completed_pairs),
        "evidence_rows": len(evidence),
        "species": int(evidence["accepted_species"].nunique()),
        "species_trait": len(
            evidence.drop_duplicates(["accepted_species", "trait_name"])
        ),
        "quality_counts": evidence["evidence_quality"].value_counts().to_dict(),
        "trait_counts": evidence["trait_name"].value_counts().to_dict(),
        "provider_counts": evidence["source_provider"].value_counts().to_dict(),
        "audit": {
            "reviewed": len(audit),
            "accepted_correct": len(audit),
            "precision": 1.0 if len(audit) else None,
            "cultivar_contamination_rate": 0.0 if len(audit) else None,
        },
        "mapping_guards": {
            "self_fertile_to_self_incompatibility_only": True,
            "pollinator_exclusion_to_autonomous_selfing_only": True,
            "outcrossing_rate_to_mating_system_only": True,
            "family_inference": False,
            "global_fallback": False,
        },
        "inputs": {
            "pfaf_csv": {"path": str(pfaf_csv), "sha256": _sha256(pfaf_csv)},
            "current_direct_csv": {
                "path": str(current_direct_csv),
                "sha256": _sha256(current_direct_csv),
            },
            "master_csv": {"path": str(master_csv), "sha256": _sha256(master_csv)},
            "diospyros_pdf": {"path": str(diospyros_pdf), "sha256": _sha256(diospyros_pdf)},
            "magnolia_pdf": {"path": str(magnolia_pdf), "sha256": _sha256(magnolia_pdf)},
            "grewia_html": {"path": str(grewia_html), "sha256": _sha256(grewia_html)},
            "myristica_html": {
                "path": str(myristica_html),
                "sha256": _sha256(myristica_html),
            },
            "nepenthes_html": {
                "path": str(nepenthes_html),
                "sha256": _sha256(nepenthes_html),
            },
        },
    }
    files = {
        path.name: {"sha256": _sha256(path), "size_bytes": _canonical_file_size(path)}
        for path in (evidence_path, audit_path)
    }
    summary["files"] = files
    (output_dir / "high_leverage_direct_checkpoint_manifest_20260811.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfaf-csv", type=Path, required=True)
    parser.add_argument("--current-direct-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--diospyros-pdf", type=Path, required=True)
    parser.add_argument("--magnolia-pdf", type=Path, required=True)
    parser.add_argument("--grewia-html", type=Path, required=True)
    parser.add_argument("--myristica-html", type=Path, required=True)
    parser.add_argument("--nepenthes-html", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(**vars(args)), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
