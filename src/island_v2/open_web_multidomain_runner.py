"""Run multi-domain acquisition under the strict manuscript trait contract.

Reward and pollen-vector traits remain separate outcomes.  This entry point
also applies two source-independent refinements discovered during the real
pilot: compound colour words must not be split into contradictory shorter
colours, and institutional flora URLs may supply an exact focal binomial when
the rendered title is vernacular-only.
"""

from __future__ import annotations

import re
from urllib.parse import parse_qs, unquote_plus, urlsplit

from island_v2 import open_web_evidence as evidence
from island_v2 import open_web_multidomain as pilot

_original_term_match = evidence._term_match
_original_subject = pilot._resolved_subject
_GENERIC_JA_COLOURS = {"白色", "赤色", "黄色", "緑色", "青色", "紫色"}
_COMPOUND_JA_COLOUR = re.compile(
    r"(?:緑白色|黄白色|淡黄白色|青紫色|紫赤色|紅紫色|紫紅色|赤橙色|黄緑色)"
)


def _guarded_term_match(sentence: str, term: str) -> bool:
    if not _original_term_match(sentence, term):
        return False
    if term not in _GENERIC_JA_COLOURS:
        return True
    for compound in _COMPOUND_JA_COLOUR.finditer(sentence):
        start = compound.start()
        end = compound.end()
        for match in re.finditer(re.escape(term), sentence):
            if start <= match.start() and match.end() <= end and compound.group(0) != term:
                return False
    return True


def _url_subject(page: evidence.Page, resolver: evidence.StrictNameResolver) -> tuple[str, str]:
    for url in (page.final_url, page.requested_url):
        query = parse_qs(urlsplit(url).query)
        for key in ("name", "taxon", "species"):
            for raw in query.get(key, []):
                candidate = unquote_plus(raw).replace("~", " ").strip()
                resolution = resolver.resolve(candidate)
                if resolution.status == "accepted":
                    return resolution.accepted_species, candidate
    return "", ""


def _enhanced_subject(
    page: evidence.Page,
    resolver: evidence.StrictNameResolver,
) -> tuple[str, str]:
    species, matched = _url_subject(page, resolver)
    if species:
        return species, matched
    return _original_subject(page, resolver)


evidence._term_match = _guarded_term_match
pilot._resolved_subject = _enhanced_subject

if __name__ == "__main__":
    pilot.app()
