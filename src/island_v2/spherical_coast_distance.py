"""Minimum distance between minor great-circle coastline segments on a sphere.

Sphere radius 6371.0088 km; this is not an ellipsoidal WGS84 distance.
Candidate pruning uses spherical triangle inequalities, not a fixed k nearest set.
"""

from itertools import pairwise

import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Point

R = 6371.0088


def xyz(lonlat):
    q = np.deg2rad(np.asarray(lonlat))
    x, y = q[..., 0], q[..., 1]
    return np.stack([np.cos(y) * np.cos(x), np.cos(y) * np.sin(x), np.sin(y)], axis=-1)


def angle(a, b):
    return np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1))


def on_arc(q, a, b, ab):
    return angle(a, q) + angle(q, b) <= ab + 2e-12


def point_arc(p, a, b):
    n = np.cross(a, b)
    norm = np.linalg.norm(n, axis=-1)
    n = n / np.maximum(norm[..., None], 1e-300)
    q = p - n * np.sum(p * n, axis=-1)[..., None]
    q = q / np.maximum(np.linalg.norm(q, axis=-1)[..., None], 1e-300)
    ab = angle(a, b)
    ends = np.minimum(angle(p, a), angle(p, b))
    return np.where((norm > 1e-14) & on_arc(q, a, b, ab), np.minimum(ends, angle(p, q)), ends)


def segment_distance(a, b, c, d):
    """Vectorized angular minimum between segment ab and each segment cd."""
    v = np.minimum.reduce(
        [point_arc(a, c, d), point_arc(b, c, d), point_arc(c, a, b), point_arc(d, a, b)]
    )
    cross = np.cross(np.cross(a, b), np.cross(c, d))
    nn = np.linalg.norm(cross, axis=-1)
    q = cross / np.maximum(nn[..., None], 1e-300)
    for sign in [1, -1]:
        qq = sign * q
        hit = (nn > 1e-14) & on_arc(qq, a, b, angle(a, b)) & on_arc(qq, c, d, angle(c, d))
        v = np.where(hit, 0.0, v)
    return v


def geometry_segments(geom):
    if geom.is_empty or geom.geom_type not in {"Polygon", "MultiPolygon"}:
        raise ValueError("Expected nonempty polygonal coastline")
    polygons = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
    result = []
    for poly in polygons:
        if poly.geom_type != "Polygon":
            raise ValueError(poly.geom_type)
        for ring in [poly.exterior, *poly.interiors]:
            ll = np.asarray(ring.coords)[:, :2]
            # Split-file dateline boundaries are not physical coastlines.
            seam = (np.abs(np.abs(ll[:-1, 0]) - 180) < 1e-9) & (
                np.abs(np.abs(ll[1:, 0]) - 180) < 1e-9
            )
            a, b = xyz(ll[:-1]), xyz(ll[1:])
            length = angle(a, b)
            keep = (~seam) & (length > 1e-14)
            if not (length[keep] < np.pi - 1e-8).all():
                raise ValueError("Antipodal coastline segment")
            result.append((a[keep], b[keep]))
    return np.concatenate([q[0] for q in result]), np.concatenate([q[1] for q in result])


class CoastIndex:
    def __init__(self, a, b):
        a, b = validate_segments(a, b, allow_points=False)
        self.a, self.b = a, b
        self.vertices = cKDTree(np.concatenate([a, b]))
        self.mid = (a + b) / np.linalg.norm(a + b, axis=1)[:, None]
        self.half = angle(a, b) / 2
        self.groups = []
        bounds = np.array([0, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25, 50, 100, 500, 20001]) / R
        for lo, hi in pairwise(bounds):
            ids = np.where((self.half > lo) & (self.half <= hi))[0]
            if len(ids):
                self.groups.append((ids, cKDTree(self.mid[ids]), self.half[ids].max()))

    def distance(self, a, b):
        a, b = validate_segments(a, b, allow_points=True)
        chord, _ = self.vertices.query(np.concatenate([a, b]), workers=1)
        best = float(2 * np.arcsin(np.clip(chord.min() / 2, 0, 1)))
        for aa, bb in zip(a, b):
            mid = (aa + bb) / np.linalg.norm(aa + bb)
            half = float(angle(aa, bb) / 2)
            for ids, tree, mx in self.groups:
                cap = min(np.pi, best + half + mx + 1e-11)
                ix = tree.query_ball_point(mid, 2 * np.sin(cap / 2) + 1e-12)
                if not ix:
                    continue
                candidates = ids[ix]
                lower = angle(mid, self.mid[candidates]) - half - self.half[candidates]
                candidates = candidates[lower <= best + 1e-11]
                if len(candidates):
                    best = min(
                        best,
                        float(
                            segment_distance(aa, bb, self.a[candidates], self.b[candidates]).min()
                        ),
                    )
            if best < 1e-14:
                return 0.0
        return best * R

    def point_distance(self, lon, lat, land):
        """Site exposure: continental interior/boundary is truly zero."""
        if not np.isfinite([lon, lat]).all() or abs(lon) > 180 or abs(lat) > 90:
            raise ValueError("Invalid longitude/latitude")
        if land.covers(Point(lon, lat)):
            return 0.0
        p = xyz([[lon, lat]])
        return self.distance(p, p)


def validate_segments(a, b, allow_points):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.ndim != 2 or a.shape[1:] != (3,) or a.shape != b.shape or len(a) == 0:
        raise ValueError("Expected nonempty paired N x 3 segment arrays")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Nonfinite coastline coordinates")
    if not np.allclose(np.linalg.norm(a, axis=1), 1, atol=1e-12, rtol=0) or not np.allclose(
        np.linalg.norm(b, axis=1), 1, atol=1e-12, rtol=0
    ):
        raise ValueError("Coastline endpoints must be unit vectors")
    lengths = angle(a, b)
    if np.any(lengths >= np.pi - 1e-8) or (not allow_points and np.any(lengths <= 1e-14)):
        raise ValueError("Degenerate or antipodal coastline segment")
    return a, b
