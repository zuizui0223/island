import json
import os
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_v13_functional_bridge import run_functional_bridge

src = Path(os.environ["ISLAND_CORRECTION_SOURCES"])
out = Path(os.environ["ISLAND_CORRECTION_OUTPUT"])
cfg = yaml.safe_load((src / "chapter1_v13_functional_bridge_v1.yml").read_text(encoding="utf8"))
for mode in ["original", "corrected"]:
    ra = pd.read_csv(
        src / "h4-ra/out/MATCHED_EFFECT_ROWS.csv.gz"
        if mode == "original"
        else out / "h4_ra_corrected_effect_rows.csv.gz"
    )
    arch = pd.read_csv(
        src / "h4-arch/out/MATCHED_EFFECT_ROWS.csv.gz"
        if mode == "original"
        else out / "h4_arch_corrected_effect_rows.csv.gz"
    )
    result, manifest = run_functional_bridge(ra, arch, cfg)
    result.to_csv(out / f"h4_atomic_{mode}.csv", index=False)
    (out / f"h4_atomic_{mode}_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf8"
    )
