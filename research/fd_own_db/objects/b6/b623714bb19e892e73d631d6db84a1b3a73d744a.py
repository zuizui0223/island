"""Keep true island pollinator evidence separate from GloBI source-genus breadth.

Source data are indexed from pinned island / GloBI checkpoints. Only fully
anonymized region-level evidence and aggregate score summaries may be exported.
Unknown pollinator observations must NEVER become absences; documented
interaction breadth is not realized local pollinator abundance or efficiency.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse

CHANNELS = (
    "bombus", "non_bombus_bees", "lepidoptera",
    "diptera", "flower_visiting_birds",
)
ALLOWED = {"retained", "disrupted", "unresolved"}
COVARIATES = (
    "rao_marginal", "log_completeness", "trait_globi_fraction",
    "mean_log_genus_refs", "log_island_area_km2",
    "log_distance_to_continent_km",
    "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
)
BRIDGE_FIELDS = ("island_id", "analysis_regime", "spatial_block",
                 "color_arch_excess", "mean_source_genus_effective_channels",
                 *COVARIATES)


def channel_support(
    morphological_islands: pd.DataFrame,
    pollinator_evidence: pd.DataFrame,
    *, min_islands_per_region: int = 20,
    min_state_count: int = 10,
    min_state_blocks: int = 10,
) -> tuple[pd.DataFrame, dict]:
    """Count evidence opportunity, not true local guild richness.

    All five pollinator states must be present as rows for EACH original
    morphology-complete island, including 'unresolved' (never fill zero).
    A region-channel contrast is eligible only if both biological/evidence
    classes independently pass the count AND spatial-block gates.
    """
    if min_islands_per_region < 20 or min_state_count < 10 or min_state_blocks < 5:
        raise ValueError("Source-evidence privacy and independence gate too weak")
    missing = {"island_id", "analysis_regime", "spatial_block"} - set(morphological_islands)
    if missing:
        raise ValueError(f"Missing island metadata: {sorted(missing)}")
    req = {"island_id", "channel_id", "channel_state", "source_state",
           "channel_record_count", "background_record_count"}
    if req - set(pollinator_evidence):
        raise ValueError("Missing frozen exact-island channel fields")
    meta = morphological_islands[["island_id", "analysis_regime", "spatial_block"]].copy()
    if meta.isna().any().any() or meta.island_id.duplicated().any():
        raise ValueError("Invalid unique morphology island IDs")
    z = pollinator_evidence[list(req)].copy()
    if z[["island_id", "channel_id", "channel_state", "source_state"]].isna().any().any():
        raise ValueError("Missing channel/source state; cannot infer absence")
    if z.duplicated(["island_id", "channel_id"]).any():
        raise ValueError("Duplicate island x pollinator channel")
    if set(z.island_id) != set(meta.island_id):
        raise ValueError("Island and pollinator cohorts are not identical")
    if not z.channel_id.isin(CHANNELS).all() or not z.channel_state.isin(ALLOWED).all():
        raise ValueError("Unknown pollinator channel or evidence state")
    if not z.source_state.isin(("available", "unresolved", "structurally_absent")).all():
        raise ValueError("Invalid candidate source state")
    for col in ("channel_record_count", "background_record_count"):
        z[col] = pd.to_numeric(z[col], errors="coerce")
        if not np.isfinite(z[col].to_numpy()).all() or (z[col] < 0).any():
            raise ValueError("Missing or negative effort/count data")
    if ((z.channel_state.eq("retained") & z.channel_record_count.le(0)) |
        (z.channel_state.eq("disrupted") & z.channel_record_count.ne(0))).any():
        raise ValueError("Retained/disrupted record evidence inconsistent")
    if ((z.channel_state.ne("unresolved")) & z.source_state.ne("available")).any():
        raise ValueError("Pollinator channel certified without source opportunity")
    five = z.groupby("island_id").channel_id.nunique()
    if not five.eq(len(CHANNELS)).all() or len(z) != len(meta) * len(CHANNELS):
        raise ValueError("Missing complete five-channel observation rows")
    z = z.merge(meta, on="island_id", validate="many_to_one")
    rows = []
    for (region, channel), g in z.groupby(["analysis_regime", "channel_id"], sort=True):
        n = g.island_id.nunique()
        if n < min_islands_per_region:
            continue
        r, dis, un = (g.channel_state.eq(kind) for kind in ("retained", "disrupted", "unresolved"))
        nr, nd, nu = int(r.sum()), int(dis.sum()), int(un.sum())
        br, bd = g.loc[r, "spatial_block"].nunique(), g.loc[dis, "spatial_block"].nunique()
        rows.append({
            "region": region, "channel": channel, "n_islands": n,
            "n_retained": nr, "n_disrupted": nd, "n_unresolved": nu,
            "n_retained_blocks": int(br), "n_disrupted_blocks": int(bd),
            "can_compare_retained_vs_disrupted":
                bool(nr >= min_state_count and nd >= min_state_count and
                     br >= min_state_blocks and bd >= min_state_blocks),
        })
    state = z.pivot(index="island_id", columns="channel_id", values="channel_state")
    known = state.isin(("retained", "disrupted"))
    retained = state.eq("retained")
    disrupted = state.eq("disrupted")
    n_known, n_retained, n_disrupted = known.sum(axis=1), retained.sum(axis=1), disrupted.sum(axis=1)
    receipt = {
        "n_islands": int(len(meta)),
        "n_channel_rows": int(len(z)),
        "n_islands_all_five_evaluable": int(n_known.eq(5).sum()),
        "n_islands_all_five_retained": int((n_retained.eq(5) & n_known.eq(5)).sum()),
        "n_islands_any_confirmed_disrupted": int(n_disrupted.gt(0).sum()),
        "n_islands_with_retained_and_disrupted":
            int((n_disrupted.gt(0) & n_retained.gt(0)).sum()),
        "spearman_recorded_count_vs_evaluable_count":
            float(stats.spearmanr(n_retained, n_known).statistic),
        "all_unresolved_is_not_absent": True,
    }
    return pd.DataFrame(rows), receipt


def adjusted_globi_breadth(
    island_summary: pd.DataFrame, *, min_islands: int = 40,
    min_blocks: int = 10,
) -> dict:
    """Cluster-robust adjusted genus-linked GloBI breadth association.

    OLS cluster sandwich uses one cluster per original spatial block.
    It is post-hoc observational association and not a pollinator
    assembly or phylogenetic mechanism test.
    """
    if min_islands < 40 or min_blocks < 10:
        raise ValueError("Privacy/uncertainty threshold")
    miss = set(BRIDGE_FIELDS)-set(island_summary)
    if miss:
        raise ValueError(f"Missing bridge fields: {sorted(miss)}")
    x = island_summary[list(BRIDGE_FIELDS)].copy()
    if x.island_id.isna().any() or x.island_id.duplicated().any():
        raise ValueError("Missing or repeated island")
    if x[["analysis_regime", "spatial_block"]].isna().any().any():
        raise ValueError("Incomplete island spatial identity")
    nums = ["color_arch_excess", "mean_source_genus_effective_channels", *COVARIATES]
    for col in nums:
        x[col] = pd.to_numeric(x[col], errors="coerce")
    x = x.dropna(subset=nums)
    if not np.isfinite(x[nums].to_numpy()).all():
        raise ValueError("Nonfinite input score")
    n, groups = len(x), x.spatial_block.nunique()
    if n < min_islands or groups < min_blocks:
        raise ValueError("Insufficient complete-case island/block evidence")
    reg = pd.get_dummies(x.analysis_regime, drop_first=True, dtype=float)
    design = np.column_stack((
        np.ones(n), x["mean_source_genus_effective_channels"].to_numpy(float),
        x[list(COVARIATES)].to_numpy(float), reg.to_numpy(float),
    ))
    y = x.color_arch_excess.to_numpy(float)
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Rank-deficient confounder design")
    coefs = np.linalg.lstsq(design, y, rcond=None)[0]
    residuals = y-design @ coefs
    bread = np.linalg.inv(design.T @ design)
    meat = np.zeros((design.shape[1], design.shape[1]))
    for block in x.spatial_block.unique():
        rows = x.spatial_block.eq(block).to_numpy()
        cluster_score = design[rows].T @ residuals[rows]
        meat += np.outer(cluster_score, cluster_score)
    factor = (groups/(groups-1))*((n-1)/(n-design.shape[1]))
    se = float(np.sqrt(max(0, factor * (bread@meat@bread)[1,1])))
    estimate = float(coefs[1])
    return {
        "n_islands": int(n), "n_spatial_blocks": int(groups),
        "partial_guild_breadth_slope": estimate,
        "cluster_se": se,
        "cluster_normal_ci_low": estimate-1.95996398454*se,
        "cluster_normal_ci_high": estimate+1.95996398454*se,
        "nuisance_columns": int(design.shape[1]-2),
        "actual_local_pollinator_abundance_measured": False,
    }
