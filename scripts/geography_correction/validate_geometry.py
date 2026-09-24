import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

from island_v2.spherical_coast_distance import CoastIndex, R, angle, segment_distance, xyz


def slerp(a, b, t):
    omega = angle(a, b)
    return (np.sin((1 - t) * omega) * a + np.sin(t * omega) * b) / np.sin(omega)


rng = np.random.default_rng(20260924)
checks = []
for k in range(24):
    center = rng.uniform([-170, -65], [170, 65])
    ll = center + rng.uniform(-3, 3, (4, 2))
    a, b, c, d = xyz(ll)
    exact = float(segment_distance(a, b, c[None, :], d[None, :])[0]) * R
    numerical = differential_evolution(
        lambda t, a=a, b=b, c=c, d=d: float(angle(slerp(a, b, t[0]), slerp(c, d, t[1]))) * R,
        [(0, 1), (0, 1)],
        seed=k,
        tol=1e-10,
        atol=1e-8,
        popsize=12,
        maxiter=300,
        polish=True,
    )
    # The independent numerical solver gives an upper bound; near crossings can
    # stop at a tiny positive value rather than exact zero.
    assert abs(exact - numerical.fun) < 2e-5, (k, exact, numerical.fun)
    mid1 = slerp(a, b, 0.5)
    mid2 = slerp(c, d, 0.5)
    split = CoastIndex(np.array([c, mid2]), np.array([mid2, d])).distance(
        np.array([a, mid1]), np.array([mid1, b])
    )
    assert abs(exact - split) < 1e-7
    checks.append(
        {
            "case": k,
            "analytic_km": exact,
            "independent_numeric_km": float(numerical.fun),
            "split_difference_km": abs(exact - split),
        }
    )
record = {
    "status": "pass",
    "metric": "great-circle boundary distance on mean-radius sphere, radius 6371.0088 km; not WGS84 ellipsoid",
    "analytic_tests": "crossing, endpoint, dateline, spherical-cap pruning vs exhaustive search for long and short arcs",
    "independent_numerical_cases": checks,
    "max_numeric_difference_km": max(
        abs(x["analytic_km"] - x["independent_numeric_km"]) for x in checks
    ),
    "tolerance_km": 2e-5,
    "claim_limit": "Algorithm validation only; real coastline topology and cohort validation are separate gates.",
}
import argparse

p = argparse.ArgumentParser()
p.add_argument("--output", type=Path, required=True)
args = p.parse_args()
if args.output.exists():
    raise ValueError("Refusing to overwrite validation receipt")
args.output.write_text(json.dumps(record, indent=2), encoding="utf8")
print(record["status"], record["max_numeric_difference_km"])
