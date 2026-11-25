"""Private-only heterogeneity audit of marginal-fixed floral source prediction.

Consumes the previously computed internal matched-island score tables
(block, gain_joint_marginal_fixed, etc.) from authorized source checkpoint
analysis. It outputs nonidentifying GROUPED summaries only. Do not upload
input per-island rows or names of geographic blocks to the repository.

The effect is an observational prediction gain, not floral selection,
evolution, independent island replicates or mainland-source ancestry.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd

SOURCE_FIELD = "gain_joint_marginal_fixed"
OUTPUT_FIELDS = (
    "scope", "n_islands", "n_geographic_blocks",
    "mean_joint_marginal_fixed", "median_joint_marginal_fixed",
    "fraction_positive_islands",
    "trimmed10_mean", "trimmed10_block_q025", "trimmed10_block_q975",
    "winsor10_mean", "winsor10_block_q025", "winsor10_block_q975",
    "median_block_q025", "median_block_q975",
    "share_top10pct_of_positive_gains",
    "share_top10pct_of_net_signed_gain", "status",
)


def trim_winsor(values: np.ndarray, trim_fraction: float = .1
                ) -> tuple[float, float, float]:
    """Symmetric 10% trimmed mean, 10% winsor mean, and median."""
    x = np.sort(np.asarray(values, dtype=float))
    if x.ndim != 1 or len(x) < 10 or not np.isfinite(x).all():
        raise ValueError("Need at least ten finite per-island score records")
    if not 0 < trim_fraction < .5:
        raise ValueError("Invalid robust trim fraction")
    n = int(np.floor(len(x)*trim_fraction))
    trimmed = float(x[n:len(x)-n].mean())
    lo, hi = np.quantile(x, [trim_fraction, 1-trim_fraction])
    return trimmed, float(np.clip(x,lo,hi).mean()), float(np.median(x))


def heterogeneity_summary(data: pd.DataFrame, *, scope: str,
                          draws: int = 1999, seed: int = 20261110
                          ) -> dict:
    """Equal-island location with whole-geographic-block bootstrap.

    Bracket intervals are conditional on model fit and observed geography;
    they do not account for dependence among genera/species across blocks.
    """
    required = {"block", SOURCE_FIELD}
    if not required.issubset(data):
        raise ValueError("Missing private group keys or outcome column")
    if draws < 99 or len(data) < 10 or data.block.isna().any():
        raise ValueError("Invalid support for block-level sensitivity")
    x = data[SOURCE_FIELD].to_numpy(dtype=float)
    if not np.isfinite(x).all():
        raise ValueError("Nonfinite per-island source gain")
    grouped = [sub[SOURCE_FIELD].to_numpy(dtype=float)
               for _,sub in data.groupby("block",sort=True)]
    if len(grouped) < 5:
        raise ValueError("Need at least five separate spatial blocks")
    robust = np.array([trim_winsor(x)[0], trim_winsor(x)[1],
                       trim_winsor(x)[2]])
    rng = np.random.default_rng(seed)
    boots = np.empty((draws,3),dtype=float)
    for j in range(draws):
        selected=rng.integers(len(grouped),size=len(grouped))
        samples=np.concatenate([grouped[k] for k in selected])
        boots[j]=trim_winsor(samples)
    bounds=np.quantile(boots,[.025,.975],axis=0)
    top_n=max(1,int(.1*len(x)))
    top_sum=float(np.sort(x)[-top_n:].sum())
    positive_total=float(x[x>0].sum())
    net=float(x.sum())
    return {
        "scope":scope, "n_islands":len(x),
        "n_geographic_blocks":len(grouped),
        "mean_joint_marginal_fixed":float(x.mean()),
        "median_joint_marginal_fixed":float(np.median(x)),
        "fraction_positive_islands":float(np.mean(x>0)),
        "trimmed10_mean":float(robust[0]),
        "trimmed10_block_q025":float(bounds[0,0]),
        "trimmed10_block_q975":float(bounds[1,0]),
        "winsor10_mean":float(robust[1]),
        "winsor10_block_q025":float(bounds[0,1]),
        "winsor10_block_q975":float(bounds[1,1]),
        "median_block_q025":float(bounds[0,2]),
        "median_block_q975":float(bounds[1,2]),
        "share_top10pct_of_positive_gains":top_sum/positive_total if positive_total>0 else float("nan"),
        "share_top10pct_of_net_signed_gain":top_sum/net if abs(net)>1e-10 else float("nan"),
        "status":"posthoc_heterogeneous_matched_source_prediction_not_selection",
    }


def run(all_files: list[Path], direct_files: list[Path],
        output: Path, *, draws: int = 1999, seed: int = 20261110
        ) -> dict:
    if len(all_files) != 5 or len(direct_files) != 5:
        raise ValueError("Exactly five archived per-fold inputs per scope required")
    result=[]
    for scope, files in (("all", all_files),("direct", direct_files)):
        parts=[]
        for file in files:
            # Explicit usecols prevents accidental loading/output of island IDs
            # or other restricted columns in a future internal checkpoint.
            part=pd.read_csv(file,usecols=["block", SOURCE_FIELD])
            parts.append(part)
        result.append(heterogeneity_summary(
            pd.concat(parts,ignore_index=True),scope=scope,
            draws=draws,seed=seed))
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(result,columns=OUTPUT_FIELDS).to_csv(
        output/"coupling_source_heterogeneity.csv",index=False)
    manifest={
        "schema":"fd_private_coupling_heterogeneity_v1",
        "status":"aggregate_only_posthoc_sensitivity",
        "n_folds_per_scope":5, "draws":draws, "seed":seed,
        "privacy":"input block/island rows never persisted, only two group summaries",
        "limitations":[
            "observational source assignment and known test island genera",
            "All and Direct source scopes overlap",
            "spatial block bootstrap does not establish species or phylogenetic independence",
            "top decile share of signed NET gain can exceed one due cancellation by negative source effects",
        ],
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"n_public_rows":len(result),"status":manifest["status"]}))
    return manifest


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all-files",type=Path,nargs=5,required=True)
    p.add_argument("--direct-files",type=Path,nargs=5,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=1999)
    a=p.parse_args()
    run(a.all_files,a.direct_files,a.output,draws=a.draws)


if __name__=="__main__":
    main()
