"""Fixed-island geo5 vs climate10 source-sector sensitivity with detection gates.

Pure statistics on source-rights-protected data frames; returns only regional
aggregates. Does NOT infer biological loss, arrival, effective visitation or
within-species evolution from static GIFT/GBIF species-wide evidence.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd

MORPHOLOGY = (
    "zygomorphic", "tube_or_bell", "open_radial", "deep_tube",
    "bilateral_tube_or_bell_joint",
)
SCOPES = ("all_recorded", "strict_native", "strict_native_nonendemic")
REGIONS = (
    "northern_midlatitude", "tropical", "southern_extratropical",
    "northern_high_latitude", "geography_unresolved",
)
MODES = ("geo_k5", "geo50_climate10")
MIN_DIRECT_ISLAND = 5
MIN_DIRECT_CANDIDATE = 10
MIN_ISLANDS = 30
MIN_BLOCKS = 10

NEEDED = {
    "island_id", "regime", "scope", "sector", "spatial_block",
    "source_signature", "n_island_total", "n_island_direct",
    "island_sector_share", "n_geo5_total", "n_geo5_direct",
    "geo5_sector_share", "n_climate10_total", "n_climate10_direct",
    "climate10_sector_share", "climate10_source_assigned",
}


def _numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def _valid(frame: pd.DataFrame) -> pd.Series:
    return (
        _numeric(frame.n_island_direct).ge(MIN_DIRECT_ISLAND)
        & _numeric(frame.n_geo5_direct).ge(MIN_DIRECT_CANDIDATE)
        & _numeric(frame.n_climate10_direct).ge(MIN_DIRECT_CANDIDATE)
        & frame.climate10_source_assigned.eq(True)
        & _numeric(frame.island_sector_share).between(0, 1)
        & _numeric(frame.geo5_sector_share).between(0, 1)
        & _numeric(frame.climate10_sector_share).between(0, 1)
    )


def validate_paired_frame(frame: pd.DataFrame, *, expected_islands: int | None = None) -> None:
    missing = NEEDED - set(frame)
    if missing:
        raise ValueError(f"missing original source-mode fields {sorted(missing)}")
    if frame.duplicated(["island_id", "scope", "sector"]).any():
        raise ValueError("duplicate same-island/scope/sector")
    if not frame.scope.isin(SCOPES).all() or not frame.sector.isin(MORPHOLOGY).all():
        raise ValueError("original plant scope or genuine floral morphology drifted")
    if not frame.regime.isin(REGIONS).all():
        raise ValueError("unexpected original GSHHG regime")
    if frame[["island_id", "regime", "scope", "sector"]].isna().any().any():
        raise ValueError("missing island, region or plant scope")
    n = frame.island_id.nunique()
    if expected_islands is not None and n != expected_islands:
        raise ValueError("original GIFT island population contracted")
    if len(frame) != n * len(SCOPES) * len(MORPHOLOGY):
        raise ValueError("full original island×scope×morphology roster incomplete")
    a = frame[["island_id", "regime", "spatial_block", "source_signature"]].drop_duplicates()
    if a.island_id.duplicated().any():
        raise ValueError("original island geography or reference donor varies by phenotype")
    for field in ("n_island_total", "n_island_direct",
                  "n_geo5_total", "n_geo5_direct",
                  "n_climate10_total", "n_climate10_direct"):
        values = _numeric(frame[field])
        if values.isna().any() or (values < 0).any() or (
            (values % 1) != 0
        ).any():
            raise ValueError(f"invalid original trait/source opportunity {field}")
    if (frame.n_island_direct > frame.n_island_total).any():
        raise ValueError("more Direct than original island plants")
    for n_direct, n_total in (
        ("n_geo5_direct", "n_geo5_total"),
        ("n_climate10_direct", "n_climate10_total"),
    ):
        if (frame[n_direct] > frame[n_total]).any():
            raise ValueError("more Direct than original source candidate plants")
    for share, carriers in (
        ("island_sector_share", "n_island_direct"),
        ("geo5_sector_share", "n_geo5_direct"),
        ("climate10_sector_share", "n_climate10_direct"),
    ):
        z = _numeric(frame[share])
        if ((~z.between(0, 1)) & z.notna()).any():
            raise ValueError(f"invalid binary sector proportion {share}")
        if (z.isna() != frame[carriers].eq(0)).any():
            raise ValueError(f"original missing trait share mishandled {share}")
    if (frame.loc[~frame.climate10_source_assigned,
                  "n_climate10_total"].ne(0)).any():
        raise ValueError("climate10 source missingness was imputed")


def _block_summary(frame: pd.DataFrame, *, draws: int, suffix: str) -> dict:
    n = len(frame)
    group = frame.dropna(subset=["spatial_block"]).groupby("spatial_block", sort=True)
    G = group.ngroups
    res = {"n_islands": int(n), "n_physical_blocks": int(G)}
    if n < MIN_ISLANDS or G < MIN_BLOCKS:
        res["status"] = "NOT_IDENTIFIABLE_GEOGRAPHIC_SUPPORT"
        return res
    cols = ("delta_geo5", "delta_climate10", "paired_difference_climate_minus_geo5")
    values = group[list(cols)].mean().to_numpy(float)
    rng = np.random.default_rng(int.from_bytes(sha256(suffix.encode()).digest()[:8], "big"))
    ii = rng.integers(G, size=(draws, G))
    boot = values[ii].mean(axis=1)
    mean = values.mean(axis=0)
    for j, key in enumerate(cols):
        res[key + "_block_mean"] = float(mean[j])
        res[key + "_block_bootstrap95"] = [
            float(v) for v in np.quantile(boot[:, j], [.025, .975])
        ]
    res["same_island_mode_sign_concordance"] = float(np.mean(
        np.sign(frame.delta_geo5.to_numpy(float)) ==
        np.sign(frame.delta_climate10.to_numpy(float))
    ))
    res["status"] = "DESCRIPTIVE_SOURCE_DEFINITION_SENSITIVITY"
    return res


def same_island_source_mode_summary(
    frame: pd.DataFrame, *, expected_islands: int | None = None,
    draws: int = 999,
) -> dict:
    validate_paired_frame(frame, expected_islands=expected_islands)
    if draws < 199:
        raise ValueError("insufficient geographic block replicates")
    result = []
    for region in REGIONS:
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                q = frame.loc[
                    frame.regime.eq(region) & frame.scope.eq(scope)
                    & frame.sector.eq(sector)
                ].copy()
                both = q.loc[_valid(q)].copy()
                both["delta_geo5"] = both.island_sector_share - both.geo5_sector_share
                both["delta_climate10"] = both.island_sector_share - both.climate10_sector_share
                both["paired_difference_climate_minus_geo5"] = (
                    both.delta_climate10 - both.delta_geo5
                )
                if len(both) and not np.allclose(
                    both.paired_difference_climate_minus_geo5,
                    both.geo5_sector_share - both.climate10_sector_share,
                    rtol=0, atol=1e-12,
                ):
                    raise AssertionError("same-island source arithmetic mismatch")
                if len(both):
                    cov_island = both.n_island_direct / both.n_island_total
                    cov_geo = both.n_geo5_direct / both.n_geo5_total
                    cov_clim = both.n_climate10_direct / both.n_climate10_total
                    cov_mask = ((cov_island-cov_geo).abs() <= .10) & (
                        (cov_island-cov_clim).abs() <= .10
                    )
                else:
                    cov_mask = pd.Series([], index=both.index, dtype=bool)
                even = both.loc[cov_mask]
                row = {
                    "region": region, "scope": scope, "morphological_sector": sector,
                    "n_original_region_floras": int(len(q)),
                    "n_with_geo5_source_direct10": int(q.n_geo5_direct.ge(10).sum()),
                    "n_with_climate10_source_direct10": int(q.n_climate10_direct.ge(10).sum()),
                    "n_same_islands_eligible_both_source_definitions": int(len(both)),
                    "n_same_islands_Direct_documentation_gap_both_le10pp": int(len(even)),
                    "full_4505_study_population_not_intersection_only": True,
                }
                row["same_island_paired_all"] = _block_summary(
                    both, draws=draws, suffix=f"{region}|{scope}|{sector}|all"
                )
                row["same_island_documentation_gap_both_le10pp"] = _block_summary(
                    even, draws=draws, suffix=f"{region}|{scope}|{sector}|doc10"
                )
                result.append(row)
    return {
        "status": "ORIGINAL_SOURCE_GEO5_CLIMATE10_SAME_ISLAND_MORPHOLOGY",
        "n_original_unique_botanical_islands": int(frame.island_id.nunique()),
        "n_original_island_scope_sector_rows": int(len(frame)),
        "n_scope_sector_regional_panels": len(result),
        "n_original_islands_with_climate10_source_assignment": int(
            frame.loc[frame.scope.eq(SCOPES[0])
                      & frame.sector.eq(MORPHOLOGY[0])
                      & frame.climate10_source_assigned].island_id.nunique()
        ),
        "panels": result,
        "limitations": {
            "paired_source_contrast_is_algebraically_source_geo_minus_climate":
                True,
            "no_independently_repeated_island_evolution_or_pollinator_loss": True,
            "trait_documentation_conditional_not_randomly_missing": True,
            "same_island_original_taxa_for_each_scope_between_modes": True,
            "no_17D_KDE_or_fabricated_joint_trait_completion": True,
        },
    }
