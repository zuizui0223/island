import httpx

from island_v2.gbif_submit_strict import classify_rejection


def _response(status_code: int, body: str) -> httpx.Response:
    request = httpx.Request("POST", "https://api.gbif.org/v1/occurrence/download/request")
    return httpx.Response(status_code=status_code, text=body, request=request)


def test_classify_rejection_recognizes_capacity_limit():
    body = "Unable to create a new download: too many simultaneous downloads for this user."
    assert classify_rejection(_response(400, body)) == "deferred_capacity"


def test_classify_rejection_recognizes_geometry_error():
    body = "Invalid geometry: the polygon appears to cross the antimeridian."
    assert classify_rejection(_response(400, body)) == "geometry_repair_pending"


def test_classify_rejection_falls_back_to_generic_failure():
    body = "Unauthorized"
    assert classify_rejection(_response(401, body)) == "submit_failed"
