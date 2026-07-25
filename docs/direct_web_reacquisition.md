# Local direct-web trait reacquisition

This lane resolves selected gaps without GitHub Actions. It starts from the fixed
`106,295 species x 3 axes` coverage ledger and produces immutable local artifacts.

## Providers and promotion gates

- WFO 2026-06 resolves only exact accepted names and exact species-rank synonyms.
  Fuzzy matches and infraspecific-to-species promotion are rejected.
- Europe PMC search metadata is followed to official OA `fullTextXML`. A passage
  is species-direct only when the accepted name or WFO-verified synonym and a
  direct reproductive term occur in the same sentence or table row.
- Plazi TreatmentBank XML is accepted only when the focal treatment contains the
  exact accepted species. Floral measurements are converted with the existing
  controlled thresholds.
- DOI is the preferred source lineage. PMCID or the Plazi master document is used
  only when the original DOI is absent.
- No family inference, global fallback, abbreviated-name expansion, or pending
  candidate auto-promotion is allowed.

Primary references:

- [WFO Plant List 2026-06](https://zenodo.org/records/20782718)
- [Europe PMC REST API](https://europepmc.org/RestfulWebService)
- [Plazi TreatmentBank API](https://api.plazi.org/GgSrvApi/v1/gettingStarted)

## Commands

Run WFO synonym expansion and Europe PMC exact-title/abstract plus OA full-text
retrieval over a bounded target list:

```powershell
island-v2-direct-web-reacquisition europe-pmc-fulltext `
  --coverage-csv integrated_species_axis_coverage.csv.gz `
  --target-species-csv target_species.csv `
  --resolve-wfo-synonyms `
  --max-aliases-per-species 10 `
  --min-interval-seconds 1 `
  --output-dir direct-web-europe-pmc
```

Run one exact Plazi treatment:

```powershell
island-v2-plazi-trait-reacquisition treatment `
  --accepted-species "Abarema cochliacarpos" `
  --treatment-uuid 03A5194AFF97FFF5B0E22FDF5E59F94E `
  --output-dir direct-web-plazi-abarema
```

Combine frozen evidence against the fixed ledger:

```powershell
island-v2-direct-web-reacquisition combine `
  --coverage-csv integrated_species_axis_coverage.csv.gz `
  --evidence-csv direct_web_fulltext_evidence.csv.gz `
  --evidence-csv direct_web_plazi_evidence.csv.gz `
  --output-dir direct-web-combined
```

## 2026-07-25 local pilot

The accepted-name-only Europe PMC pass over 100 two-axis-complete species yielded
no new species-axis. The 13-species WFO synonym pass generated 81 exact synonyms
and recovered reproductive-assurance evidence for `Erythranthe guttata` under
`Mimulus guttatus`.

Two independent articles supported mixed mating, self-compatibility, or autogamy:

- DOI `10.3389/fpls.2024.1378568`
- DOI `10.3389/fpls.2024.1411868`

Plazi treatment `03A5194AFF97FFF5B0E22FDF5E59F94E`, from DOI
`10.11646/phytotaxa.601.1.3`, explicitly reported a `4.1-7.8 mm` corolla for
`Abarema cochliacarpos`. This maps to `flower_size_class=small`.

Against integrated run `30107066525`, the combined pure gain was:

- reproductive assurance: `+1`
- floral structural complexity: `+1`
- total: `+2 species-axis`

The combined ledger retained 318,885 rows. Filled cells changed from 64,026 to
64,028 and unresolved cells from 254,859 to 254,857.
