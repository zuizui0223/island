import numpy as np
import pytest
from shapely.geometry import Point, box

from island_v2.spherical_coast_distance import (
    CoastIndex,
    R,
    geometry_segments,
    segment_distance,
    xyz,
)


def test_crossing_endpoint_and_dateline_minima():
    for coords, degrees in [
        ([[0, 0], [10, 0], [5, -5], [5, 5]], 0),
        ([[0, 0], [10, 0], [15, 0], [20, 0]], 5),
        ([[179, 0], [-179, 0], [180, 1], [180, 2]], 1),
    ]:
        a, b, c, d = xyz(coords)
        assert float(segment_distance(a, b, c, d)) == pytest.approx(np.deg2rad(degrees), abs=1e-12)


def test_index_matches_exhaustive_for_long_and_short_arcs():
    rng = np.random.default_rng(91)
    for width in [0.1, 100]:
        ll = rng.uniform([-60, -40], [60, 40], (80, 2))
        a, b = xyz(ll), xyz(ll + rng.uniform(-width, width, (80, 2)) * [1, 0.25])
        index = CoastIndex(a, b)
        for _ in range(12):
            ll = rng.uniform([-60, -40], [60, 40], (2, 2))
            c, d = xyz(ll), xyz(ll + rng.uniform(-0.1, 0.1, (2, 2)))
            brute = min(float(segment_distance(x, y, a, b).min()) for x, y in zip(c, d)) * R
            assert index.distance(c, d) == pytest.approx(brute, abs=1e-7)



def test_exact_great_circle_normal_does_not_collapse_to_zero_distance():
    index = CoastIndex(xyz([[0, 0]]), xyz([[10, 0]]))
    north_pole = np.array([[0.0, 0.0, 1.0]])
    assert index.distance(north_pole, north_pole) == pytest.approx(np.pi * R / 2, abs=1e-8)


def test_nearshore_gap_is_positive_and_mainland_point_zero():
    land = box(0, 0, 1, 1)
    index = CoastIndex(*geometry_segments(land))
    assert index.point_distance(0.5, 0.5, land) == 0
    assert index.point_distance(1.0001, 0.5, land) > 0
    assert index.distance(*geometry_segments(box(1.0001, 0.1, 1.1, 0.2))) > 0


def test_dateline_seam_not_a_physical_coastline():
    a, _b = geometry_segments(box(179, 0, 180, 1))
    assert len(a) == 3


@pytest.mark.parametrize(
    "a,b",
    [
        (np.empty((0, 3)), np.empty((0, 3))),
        (xyz([[0, 0]]), xyz([[180, 0]])),
        (np.array([[float("nan"), 0, 0]]), xyz([[1, 0]])),
    ],
)
def test_invalid_coast_segments_rejected(a, b):
    with pytest.raises(ValueError):
        CoastIndex(a, b)


def test_empty_geometry_rejected():
    with pytest.raises(ValueError):
        geometry_segments(Point())
