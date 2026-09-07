"""Wave57: audit prior output, then discover reproductive literature; never promote.

Public Wave52 is a search baseline, not the private-plus-public canonical ledger.
No abstract co-occurrence becomes a trait value, quality grade or genus-rule vote.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

AXIS = "reproductive_assurance"
TRAITS = {"self_incompatibility", "mating_system", "autonomous_selfing_capacity", "cleistogamy"}
BINOMIAL = re.compile(r"[A-Z][a-z]+ [a-z][a-z-]+")
TERMS = re.compile(r"self[ -]?(?:compatib|incompatib|pollinat|fertili)|selfing|autogam|cleistogam|breeding system|mating system|outcross|bagged|bagging", re.I)
BLOCKED = {"Sideroxylon", "Illicium", "Portulaca", "Eugenia", "Liparis", "Melaleuca", "Callicarpa", "Durio"}
# Only documented synonyms; keys on the right remain frozen analysis-master labels.
ALIASES = {"Lindernia micrantha": "Vandellia micrantha", "Torenia micrantha": "Vandellia micrantha", "Lindernia setulosa": "Vandellia setulosa"}
PRIORITY = ["Schoenoplectiella", "Spermacoce", "Dicliptera", "Lindernia", "Cyanotis", "Myriophyllum", "Vandellia", "Torenia", "Bonnaya", "Cyrtandra", "Pilea", "Lithocarpus"]


def read_csv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row)) or ["status"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def is_species(name: str) -> bool:
    return bool(BINOMIAL.fullmatch(name)) and name.split()[1] not in {"sp", "spp", "cf", "aff", "x"}


def audit_prior(rows: list[dict], universe: set[str]) -> list[dict]:
    audited = []
    for index, row in enumerate(rows, 1):
        name = row.get("accepted_species", "")
        if row.get("trait_name") not in TRAITS:
            reason = "off_reproductive_axis"
        elif not is_species(name):
            reason = "invalid_species_name"
        elif name not in universe:
            reason = "outside_public_analysis_universe"
        elif not re.search(r"(?<![A-Za-z-])" + re.escape(name) + r"(?![A-Za-z-])", row.get("excerpt", "")):
            reason = "no_exact_species_in_excerpt"
        else:
            reason = "pending_source_identity_lineage_review"
        # Preserve original machine metadata; explicitly disallow promotion.
        audited.append({"input_row": index, "accepted_species": name,
                        "trait_name": row.get("trait_name", ""),
                        "provider": row.get("provider", ""), "source_url": row.get("source_url", ""),
                        "machine_quality": row.get("evidence_quality", ""),
                        "decision": reason, "promotion_allowed": "false"})
    return audited


def paper_leads(paper: dict, universe: set[str], quality: dict[str, str]) -> list[dict]:
    text = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", paper["text"])))
    if not TERMS.search(text):
        return []
    names = set(BINOMIAL.findall(text))
    leads = []
    for source_name in sorted(names):
        accepted = ALIASES.get(source_name, source_name)
        if not is_species(accepted) or accepted not in universe:
            continue
        leads.append({"accepted_species": accepted, "source_species_name": source_name,
                      "axis": AXIS, "source_lineage": paper["lineage"],
                      "source_url": paper["url"], "provider": paper["provider"],
                      "wave52_quality": quality.get(accepted, ""),
                      "evidence_scope": "unreviewed_document_cooccurrence",
                      "normalized_value": "", "evidence_quality": "unreviewed",
                      "species_claim_verified": "false", "promotion_allowed": "false",
                      "genus_rule_training_allowed": "false",
                      "conflicted_genus_rule_blocked": str(accepted.split()[0] in BLOCKED).lower(),
                      "review_status": "needs_fulltext_species_trait_attribution"})
    return leads


def fetch_json(url: str, cache: Path) -> dict:
    key = cache / (hashlib.sha256(url.encode()).hexdigest() + ".json")
    if key.exists():
        return json.loads(key.read_text())
    for attempt in range(2):
        try:
            request = Request(url, headers={"User-Agent": "island-wave57/1.0 (https://github.com/zuizui0223/island)"})
            with urlopen(request, timeout=20) as response:
                data = json.load(response)
            key.parent.mkdir(parents=True, exist_ok=True)
            key.write_text(json.dumps(data))
            time.sleep(0.5)
            return data
        except (OSError, ValueError):
            if attempt:
                raise
            time.sleep(2)
    raise RuntimeError("unreachable")


def discover(genus: str, provider: str, cache: Path) -> tuple[list[dict], int]:
    if provider == "europe_pmc":
        query = f'TITLE_ABS:"{genus}" AND ("self-compatible" OR "self-incompatible" OR selfing OR "breeding system" OR "mating system" OR cleistogamy OR autogamy OR bagging)'
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": 100})
        data = fetch_json(url, cache)
        items = data.get("resultList", {}).get("result", [])
        papers = []
        for item in items:
            doi = str(item.get("doi", "")).lower()
            pmcid = item.get("pmcid", "")
            identity = f"{item.get('source', '')}:{item.get('id', '')}"
            if identity == ":" and not doi:
                continue
            source_url = "https://doi.org/" + doi if doi else "https://europepmc.org/article/" + ("PMC/" + pmcid if pmcid else identity.replace(":", "/"))
            papers.append({"lineage": "doi:" + doi if doi else "epmc:" + identity,
                           "url": source_url, "provider": provider,
                           "text": item.get("title", "") + " " + item.get("abstractText", "")})
        return papers, int(data.get("hitCount", len(items)))
    url = "https://api.crossref.org/works?" + urlencode({"query.bibliographic": genus + " breeding system self incompatibility mating pollination", "filter": "type:journal-article", "rows": 25})
    data = fetch_json(url, cache).get("message", {})
    papers = []
    for item in data.get("items", []):
        doi = str(item.get("DOI", "")).lower()
        if not doi:
            continue
        papers.append({"lineage": "doi:" + doi, "url": "https://doi.org/" + doi,
                       "provider": provider,
                       "text": " ".join(item.get("title", [])) + " " + item.get("abstract", "")})
    return papers, int(data.get("total-results", len(papers)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-genera", type=int, default=40)
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.max_genera <= 200:
        parser.error("max-genera must be between 1 and 200")
    coverage = [r for r in read_csv(args.coverage) if r["axis"] == AXIS]
    quality = {r["accepted_species"]: r["quality"] for r in coverage}
    if len(quality) != 106295 or len(coverage) != len(quality):
        raise ValueError("expected unique 106295-species public reproductive-axis baseline")
    universe = set(quality)
    old = read_csv(args.previous)
    audited = audit_prior(old, universe)
    write_csv(args.output / "wave56_machine_audit.csv", audited)
    packet_rows = []
    for wave in (56, 57):
        root = Path(f"data/v2/staging/traits/wave{wave}_reproductive_recovery")
        for path in sorted(root.glob("reviewed_direct_evidence*.csv")):
            packet_rows.extend(read_csv(path))
    packet_counts = Counter(r["accepted_species"] for r in packet_rows if r.get("axis") == AXIS)
    public_check = [{"accepted_species": n, "axis": AXIS, "reviewed_trait_rows": count,
                     "public_wave52_quality": quality.get(n, ""),
                     "public_presence": "out_of_scope" if n not in universe else ("already_resolved" if quality[n] else "previously_unresolved"),
                     "canonical_collision_audit_complete": "false", "promotion_allowed": "false"}
                    for n, count in sorted(packet_counts.items())]
    write_csv(args.output / "reviewed_packet_public_presence.csv", public_check)
    tasks = read_csv(args.plan)
    invalid = [r for r in tasks if r.get("axis") != AXIS or not is_species(r["accepted_species"]) or r["accepted_species"] not in universe]
    write_csv(args.output / "wave56_invalid_tasks.csv", invalid)
    # The earlier priority plan selects genera only, never genus-only target cells.
    valid = [r for r in tasks if r not in invalid]
    ordered = list(dict.fromkeys(PRIORITY + [r["genus"] for r in valid]))[:args.max_genera]
    unresolved = Counter(n.split()[0] for n, q in quality.items() if not q and is_species(n))
    write_csv(args.output / "wave57_genus_search_plan.csv", [
        {"genus": g, "public_wave52_unresolved": unresolved[g], "automatic_promotion": "false"}
        for g in ordered])
    summary = {"contract": "wave57_reproductive_literature_discovery_v1",
               "baseline": "public Wave52 only; private TRY plus Wave55 collision audit still required",
               "canonical_net_change_verified": None, "automatic_promotions": 0,
               "reviewed_packet_rows": sum(packet_counts.values()),
               "reviewed_packet_species": len(packet_counts),
               "reviewed_packet_public_presence": dict(Counter(r["public_presence"] for r in public_check)),
               "wave56_raw_rows": len(old), "wave56_audit": dict(Counter(r["decision"] for r in audited)),
               "wave56_invalid_target_rows": len(invalid), "planned_genera": len(ordered),
               "acquisition_executed": args.acquire, "status": "audit_only"}
    logs, leads, papers = [], {}, {}
    failures = Counter()
    stopped = set()
    if args.acquire:
        for genus in ordered:
            for provider in ("europe_pmc", "crossref"):
                if provider in stopped:
                    continue
                try:
                    found, hits = discover(genus, provider, args.output / "private_cache")
                    matched = 0
                    for paper in found:
                        candidates = paper_leads(paper, universe, quality)
                        if not TERMS.search(paper["text"]):
                            continue
                        papers[(provider, paper["lineage"])] = {"provider": provider, "source_lineage": paper["lineage"], "source_url": paper["url"], "review_status": "discovery_only"}
                        for row in candidates:
                            leads[(row["accepted_species"], row["source_lineage"], provider)] = row
                            matched += 1
                    logs.append({"genus": genus, "provider": provider, "status": "success", "hit_count": hits, "fetched": len(found), "lead_rows": matched, "truncated": hits > len(found), "error": ""})
                    failures[provider] = 0
                except (OSError, ValueError, TypeError) as exc:
                    failures[provider] += 1
                    logs.append({"genus": genus, "provider": provider, "status": "error", "error": str(exc), "hit_count": "", "fetched": 0, "lead_rows": 0, "truncated": ""})
                    if failures[provider] >= 3:
                        stopped.add(provider)
                write_csv(args.output / "query_audit.csv", logs)
                write_csv(args.output / "literature_review_queue.csv", list(leads.values()))
                write_csv(args.output / "source_discovery.csv", list(papers.values()))
        summary.update(status="partial_provider_failure" if any(r["status"] == "error" for r in logs) else "completed",
                       successful_queries=sum(r["status"] == "success" for r in logs),
                       failed_queries=sum(r["status"] == "error" for r in logs),
                       providers_stopped_after_three_errors=sorted(stopped),
                       distinct_species_lineage_leads=len({(s, l) for s, l, _ in leads}),
                       distinct_lead_species=len({s for s, _, _ in leads}),
                       unique_document_lineages=len({l for _, l in papers}),
                       acquisition_truncation_policy="first 100 EuropePMC / 25 Crossref hits per genus; not exhaustive")
    summary["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    (args.output / "wave57_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.acquire and stopped:
        raise SystemExit("provider circuit breaker opened; partial artifacts retained")


if __name__ == "__main__":
    main()
