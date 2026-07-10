from pathlib import Path
import re

path = Path("tests/test_global_candidate_index.py")
text = path.read_text(encoding="utf-8")
pattern = re.compile(
    r"def _write_ledger\(campaign_dir\) -> None:\n.*?(?=\ndef test_candidate_index_deduplicates_source_rows_and_retains_wave_history)",
    re.DOTALL,
)
replacement = '''def _write_ledger(campaign_dir) -> None:
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "family": "FamA",
                "n_islands": "3",
                "n_records": "12",
                "machine_biotic_candidate": "True",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "processed",
                "floral_access_wikimedia_status": "processed",
                "alternative_guild_wikimedia_status": "processed",
                "alternative_guild_openalex_status": "pending",
            },
            {
                "accepted_species": "Beta two",
                "family": "FamB",
                "n_islands": "1",
                "n_records": "4",
                "machine_biotic_candidate": "True",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "processed",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
            {
                "accepted_species": "Gamma three",
                "family": "FamC",
                "n_islands": "2",
                "n_records": "8",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "pending",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
        ]
    )
    ledger.to_csv(
        campaign_dir / "campaign_ledger.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

'''
text, count = pattern.subn(replacement, text, count=1)
if count != 1:
    raise SystemExit(f"fixture replacement count={count}")
path.write_text(text, encoding="utf-8")
