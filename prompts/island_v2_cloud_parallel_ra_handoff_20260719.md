# Island v2 cloud handoff — parallel reproductive-assurance collection

Copy the prompt below into the cloud Codex task.

---

Continue Island v2 from the uploaded continuation snapshot. Do not redesign the
project, restart the 106,295-species campaign, introduce a new ontology, or use
GitHub Actions. Network acquisition is allowed. Use papers, breeding-system and
self-incompatibility databases, institutional pages, botanical-garden pages,
publisher abstracts, theses, reports, and ordinary web pages, but accept a
record only when the source explicitly describes the exact species' mating or
pollination system. Never infer from genus, family, morphology, geography, or
related taxa.

Repository and code state:

- repository: `zuizui0223/island`
- PR: `#126`
- branch: `codex/direct-evidence-continuation`
- local-collector implementation commit: `2affe39`
- snapshot directory:
  `data/v2/staging/traits/local_ra_continuation_20260719`

Restore the snapshot under `/mnt/data/local_trait_work`:

```bash
git checkout codex/direct-evidence-continuation
S="$PWD/data/v2/staging/traits/local_ra_continuation_20260719"
R=/mnt/data/local_trait_work
mkdir -p "$R"
mkdir -p /tmp/island-ra-handoff
gh release download island-v2-local-ra-continuation-20260719 \
  --repo zuizui0223/island \
  --pattern 'island-v2-local-ra-continuation-20260719-sources.tar.gz' \
  --dir /tmp/island-ra-handoff
gzip -dc "$S/provider_checkpoint.csv.gz" > "$R/provider_checkpoint.csv"
gzip -dc "$S/new_reproductive_assurance_records.csv.gz" > "$R/new_reproductive_assurance_records.csv"
gzip -dc "$S/reproductive_assurance_priority_next.csv.gz" > "$R/reproductive_assurance_priority_next.csv"
gzip -dc "$S/reproductive_assurance_priority_3449.csv.gz" > "$R/reproductive_assurance_priority_3449.csv"
cp "$S/three_axis_checkpoint.csv.gz" "$R/three_axis_checkpoint.csv.gz"
cp "$S/coverage_summary.json" "$R/coverage_summary.json"
tar -xzf /tmp/island-ra-handoff/island-v2-local-ra-continuation-20260719-sources.tar.gz \
  -C /mnt/data
```

Rebase the old absolute source paths in the restored CSV files before running a
materializer:

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

old = "/Users/rachelzhang/Documents/Codex/2026-07-18/local_trait_work"
new = "/mnt/data/local_trait_work"
root = Path(new)
for name in ["new_reproductive_assurance_records.csv"]:
    path = root / name
    frame = pd.read_csv(path, dtype=str).fillna("")
    frame["local_file_path"] = frame["local_file_path"].str.replace(
        old, new, regex=False
    )
    frame.to_csv(path, index=False)
PY
```

Validate the snapshot before collection:

- 106,295 unique species in `three_axis_checkpoint.csv.gz`
- direct coverage: colour 16,774; floral complexity 4,690;
  reproductive assurance 3,475
- all-three direct: 691
- current priority queue: 3,255 species
- cumulative newly completed from the original 3,449 queue: 194 species
- cumulative new direct records: 243
- provider checkpoint rows: 155,329
- running provider rows: 0
- queue SHA256:
  `9fec1ffce08e552957ca32fff19abb0c3f21f683df42dafa1af53e3cdea9eec4`
- axes SHA256:
  `bc76ad8c54b2c4ac8ba63cfc3a9da1412c07b607c451f4d2a9228f6d562780cb`
- uncompressed provider checkpoint SHA256:
  `24fbbcfe0b690367f646509ea67921422dd11f57e61ea3a667449888248929b7`
- uncompressed records SHA256:
  `ee773116b78b6633348bf247fd3da326d60e4309b4a3dc4c8c9ce276b34b18a1`

First finish the already-reviewed Pladias batch. The snapshot includes:

- `pladias_pending_source_evidence.csv`
- `pladias_pending_reviewed_candidates.csv`
- `pladias_pending_review_decisions.json`

It contains 7 species and 9 trait records: Erigeron acris, Eschscholzia
californica, Eurybia macrophylla, Festuca amethystina, Ficaria verna, Ficus
carica, and Frangula alnus. Rebase the candidates' `local_file_path` prefix to
`/mnt/data/local_trait_work`, then run `materialize-reviewed-cache` with the
restored queue/checkpoints/records. Do not requery these completed Pladias
provider keys.

After that, run acquisition in parallel. Parallelism is the priority:

1. Keep one coordinator as the only writer of the authoritative
   `provider_checkpoint.csv`, `three_axis_checkpoint.csv.gz`, records file,
   queue, and summary.
2. Split the current 3,255-species queue into disjoint immutable shards.
3. Run multiple workers simultaneously, with isolated output directories and
   isolated checkpoint copies. Never let workers write the same checkpoint.
4. Assign separate lanes for Pladias, DOAJ, OpenAIRE, Crossref/Semantic
   Scholar, breeding-system databases, self-incompatibility databases, paper
   search, and explicit ordinary-web descriptions.
5. Each worker must return a reviewed candidate package plus provider-key
   status delta. The coordinator merges by `species × trait × provider`, then
   serially materializes accepted packages and regenerates the queue.
6. Rebase/re-filter each pending package against the newest queue before
   materialization; discard species already completed by another worker.

Accepted evidence includes explicit self-compatible, self-incompatible,
autonomous selfing/self-fertilization, delayed selfing, mixed mating,
facultative selfing, cleistogamy, obligate/predominant outcrossing, breeding
system experiments, and pollination biology experiments. Explicit web text is
allowed. Flora or species pages are allowed only when they directly state the
mating/breeding/pollination system. Do not convert vague phrases, floral
structure, bagging-free observations, outcrossing-rate numbers without a
published category, or absence of pollinators into a trait.

Every accepted record must retain exact species, existing trait/value, raw
text, exact excerpt, URL, DOI if available, publication, year, source type,
local cached file path and SHA256, retrieval date, provider, and source run ID.
Cache the fetched source locally before acceptance. No LLM guessing.

Run focused tests before each code change is used:

```bash
uv run pytest -q tests/test_local_*reproductive*.py tests/test_local_targeted_web_reproductive_import.py
uv run ruff check src/island_v2/local_*reproductive*.py src/island_v2/local_targeted_web_reproductive_import.py tests/test_local_*reproductive*.py tests/test_local_targeted_web_reproductive_import.py
```

Continue until the queue is empty, every provider is terminal, or an entire
full provider pass yields no new direct evidence. Do not stop because one shard
or one provider batch yields zero. Keep all acquisition local; do not trigger
workflows or use GitHub Actions/artifacts.

---
