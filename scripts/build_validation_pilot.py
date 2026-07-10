"""Build the stratified gold-standard validation pilot from the collected taxa.

Reproducible builder for the Stage-6 trait-evidence validation pilot. It reads the
collected island taxa table, restricts the pool to angiosperms (a floral-syndrome
study excludes ferns, lycophytes, and gymnosperms), seeds the trait-derived
stratification axes as ``unknown`` (fail-closed: those labels are exactly what the
pilot's human review supplies), and runs the stratified selector so the pilot
spreads evenly across families rather than tracking GBIF record count.

Usage:
    python scripts/build_validation_pilot.py \
        --taxa data/v2/staging/gbif/collected/island_taxa.csv \
        --output-dir data/v2/staging/traits/validation_pilot

The output is a candidate species list for human review, not an analysis input.
No trait value, native/establishment status, Bombus applicability, or analysis
inclusion is decided here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from island_v2.trait_pilot_eval import build_stratified_pilot, load_config  # noqa: E402

# Non-angiosperm families excluded from a floral-syndrome pilot (ferns & allies,
# gymnosperms). Coarse, family-level, and review-refinable.
NON_ANGIOSPERM_FAMILIES = {
    # ferns and fern allies
    "Dryopteridaceae", "Thelypteridaceae", "Pteridaceae", "Polypodiaceae", "Aspleniaceae",
    "Athyriaceae", "Blechnaceae", "Cyatheaceae", "Dennstaedtiaceae", "Hymenophyllaceae",
    "Lomariopsidaceae", "Nephrolepidaceae", "Marattiaceae", "Ophioglossaceae", "Osmundaceae",
    "Schizaeaceae", "Selaginellaceae", "Lycopodiaceae", "Equisetaceae", "Gleicheniaceae",
    "Cystopteridaceae", "Woodsiaceae", "Tectariaceae", "Davalliaceae", "Oleandraceae",
    "Marsileaceae", "Salviniaceae", "Azollaceae", "Psilotaceae", "Isoetaceae", "Cibotiaceae",
    "Dicksoniaceae", "Dipteridaceae", "Metaxyaceae", "Saccolomataceae", "Lindsaeaceae",
    "Plagiogyriaceae", "Culcitaceae", "Loxsomataceae", "Matoniaceae", "Anemiaceae",
    "Lygodiaceae", "Onocleaceae", "Diplaziopsidaceae", "Rhachidosoraceae", "Didymochlaenaceae",
    "Hypodematiaceae", "Pteridiaceae", "Grammitidaceae", "Vittariaceae",
    # gymnosperms
    "Pinaceae", "Cupressaceae", "Podocarpaceae", "Araucariaceae", "Cycadaceae", "Zamiaceae",
    "Ginkgoaceae", "Ephedraceae", "Gnetaceae", "Taxaceae", "Cephalotaxaceae", "Sciadopityaceae",
    "Welwitschiaceae", "Phyllocladaceae",
}


def build_candidate_pool(taxa: pd.DataFrame) -> pd.DataFrame:
    """Angiosperm species with family kept and trait axes seeded ``unknown``."""
    df = taxa.copy()
    df["n_islands"] = pd.to_numeric(df.get("n_islands"), errors="coerce").fillna(0).astype(int)
    df["n_records"] = pd.to_numeric(df.get("n_records"), errors="coerce").fillna(0).astype(int)
    ang = df[~df["family"].isin(NON_ANGIOSPERM_FAMILIES)].copy()
    ang["accepted_species"] = ang["accepted_species"].str.strip('"').str.strip()
    ang = ang[~ang["accepted_species"].str.contains("\u00d7", na=False)]
    ang = ang[ang["accepted_species"].str.split().str.len() >= 2]
    return pd.DataFrame(
        {
            "accepted_species": ang["accepted_species"],
            "family": ang["family"],
            # Fail-closed: trait-derived axes are unknown until the pilot's review
            # supplies them; they are never auto-decided here.
            "native_status": "unknown",
            "growth_form": "unknown",
            "flower_conspicuousness": "unknown",
            # Transparency only -- NOT used for selection.
            "n_islands": ang["n_islands"],
            "n_records": ang["n_records"],
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxa", type=Path, default=Path("data/v2/staging/gbif/collected/island_taxa.csv"))
    parser.add_argument("--config", type=Path, default=Path("config/pilot_stratification.yml"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/v2/staging/traits/validation_pilot"))
    args = parser.parse_args()

    taxa = pd.read_csv(args.taxa, dtype=str).fillna("")
    pool = build_candidate_pool(taxa)
    config = load_config(args.config)
    pilot, summary = build_stratified_pilot(pool, config)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pilot.to_csv(args.output_dir / "validation_pilot_species.csv", index=False)
    import json

    (args.output_dir / "validation_pilot_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(
        f"Angiosperm pool: {len(pool)} species / {pool['family'].nunique()} families. "
        f"Selected {summary['n_selected']} pilot species across "
        f"{summary['n_strata_occupied']} strata."
    )


if __name__ == "__main__":
    main()
