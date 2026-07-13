from __future__ import annotations

from dataclasses import dataclass

from island_v2 import bombus_niche_inputs_production as production


@dataclass
class _Response:
    payload: dict

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return self.payload


def test_discovery_uses_camel_case_facet_limit(monkeypatch) -> None:
    seen_params: list[dict] = []
    counts = [{"name": str(index), "count": 100 - index} for index in range(1, 21)]

    def fake_get(url: str, *, params=None, timeout=None, **kwargs):
        if url.endswith("/occurrence/search"):
            seen_params.append(dict(params or {}))
            return _Response({"facets": [{"counts": counts}]})
        species_key = int(url.rsplit("/", 1)[-1])
        return _Response(
            {
                "canonicalName": f"Bombus species{species_key}",
                "rank": "SPECIES",
            }
        )

    monkeypatch.setattr(production.base, "_gbif_genus_key", lambda: 1)
    monkeypatch.setattr(production.base.httpx, "get", fake_get)

    result = production.discover_bombus_species(50, 20)

    assert len(result) == 20
    assert seen_params[0]["facetLimit"] == 500
    assert "facet_limit" not in seen_params[0]
