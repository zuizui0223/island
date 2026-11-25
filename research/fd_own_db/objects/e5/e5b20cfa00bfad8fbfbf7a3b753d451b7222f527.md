# 15 November 2026 — Five original NHM botanical description corpora: verified provider-specific trait opportunity

**ACTUAL ORIGINAL DATA RUN SUCCESS**: [island Actions #43683208060](https://github.com/zuizui0223/island/actions/runs/43683208060), anonymous [artifact #14149691188](https://github.com/zuizui0223/island/actions/runs/43683208060/artifacts/14149691188), extracted JSON SHA256 `eed78876e6d3aac2ae5b4d0a8950855474f1ea6b58e25ede6110bdaf2ced1c4a`. Original full **4,505** island botanical floras and 748 GIFT native mainland source regions restored under protected original SHA-256 pins. Five publicly licensed CC BY 4.0 NHM constituent files fetched at immutable independent NHM repository git commit `4843c22bb2fc91018d9828dcddebf067c83b1097`; **each individual original Git blob SHA** checked against the archived repository (BHL `f46e15b9...`, Ecoflora `dbdc6b4...`, Flora of China `ec87a99...`, Flora of North America `3710b08...`, Flora of Pakistan `7f50641...`). Only anonymous source-specific rows, conflict counts and source-only donor availability escape the rights-restricted GIFT processing runner.

The previously source-pinned [NHM combined source study](INDEPENDENT_NHM_PUBLIC_FLORAL_SOURCE_FEASIBILITY.md) found that the original `combined.csv` has **no source/provider column**. This original source-specific study resolves the true *text-corpus* provenance by retrieving source files separately, never treating the merged category as independent multiple witnesses.

## One corpus at a time: source-only native GIFT plant opportunities

Each species must have an exactly matched accepted botanical binomial, a strict unambiguous NHM flower textual category, and appear as a **GIFT native mainland source** species name while being **ABSENT from all 4,505 original island flora species lists**. Counts below are per textual corpus and therefore overlap across sources; they must NOT be summed to estimate unique species or independent evolutionary samples.

| Original NHM published text corpus | Original description rows | Strict single-colour external names overall | GIFT native source-only strict flower colour species | GIFT native source-only colour after separate corolla/petal expansion |
|---|---:|---:|---:|---:|
| Biodiversity Heritage Library (BHL OCR) | **106,433** | 4,944 | **1,377** | **3,359** |
| Flora of China | 23,622 | 1,983 | **1,111** | **4,807** |
| Flora of North America | 9,482 | 857 | 316 | 973 |
| Flora of Pakistan | 3,537 | 609 | 175 | 350 |
| Ecoflora | 1,989 | 251 | 5 | 16 |

The original BHL textual dataset has **200 source species** with more than one interpreted *primary flower colour* category across descriptions of the same accepted binomial; these are omitted from the strict validated-source reference. Other corpora have no within-corpus colour conflict in this processing. This does **not** prove they are error-free or that flower colour is constant within species.

## Cross-corpus consensus and contradictions, not an invented true mainland phenotype

| Trait evidence scope | Unique taxa with at least one uncontradicted interpreted source-corpus category | Source-only GIFT mainland species after excluding ALL 4,505 island species | Of these, supported by TWO or more distinct original source corpora with the same normalized category | Source corpora contradict on this number of taxa |
|---|---:|---:|---:|---:|
| **Primary whole-flower colour** | **8,136** | **2,896** | **43** | **100** |
| Whole-flower / corolla / petal colour sensitivity | **19,823** | **9,027** | **321** | **297** |
| Flower shape | 922 | 417 | **4** | 1 |
| Flower symmetry | 138 | 37 | **1** | 0 |

**Scientific source-provenance finding:** the existence of 2,896 strict externally described one-corpus source-only flower-colour species (or 9,027 in organ-colour sensitivity) demonstrates an expanded *documented source opportunity*. However, only **43 / 2,896** primary GIFT mainland source-only species have two or more distinct NHM text corpora *agreeing*, and only **321 / 9,027** in expanded colour sensitivity. Multiple corpora may still cite the SAME underlying taxonomic work; even two-corpus agreement is NOT two independently sampled mainland plant populations.

This independently verifies that source-only phenotypes are available at more than zero taxa, but **dataset evidence independence is the dominant remaining bottleneck**, especially flower shape with only FOUR source-only multiply documented species.

## Independent published textual label vs original species-Direct category

Among exact original species and external text corpus consensus traits, whole-flower primary colour **2,779 / 2,924 (95.04%)** agree, broader organ colour **6,027 / 6,219 (96.91%)**, flower shape **337 / 362 (93.09%)**, and symmetry **74 / 76 (97.37%)**. These are accepted taxon **category concordances conditional on documentation opportunity**, NOT causal island-mainland within-species stability, source-independent taxonomic evolution, or direct pollen effect.

## Consequences for the global FD study and next tests

The main paper's source universe stays all **4,505 GIFT island floras**, with original individual Direct floral and reproductive trait-specific samples (colour original n≥5 3,132 islands, colour×form original n≥5 2,324). Do not impose complete-case original colour+form+symmetry or shrink to only 43 multiply described external mainland donor species.

The original NHM primary combined source-only **2,751** donor species and reconstructed five-corpus consensus **2,896** differ by source aggregation and text conflict gates; these are related, nonidentical source definitions and **must not be represented as two independent biological datasets or pooled arithmetically**. Their compatible original 4,505 island name-exclusion filters agree in meaning.

**Next discriminating test:** an explicitly observational global **source-only flower-colour composition prediction** trained only on independently published original NHM textual categories for GIFT mainland native species **ABSENT from every one of the original 4,505 island floras**, compared against the external global colour-frequency prior and TWO outcome-blind source-size-matched pseudo-source assignments. This can estimate a geographic recorded-flora association without requiring congeneric coverage or all-three-trait completeness, but cannot prove colonization, selection or pollen fitness. A positive result should be repeated separately using BHL- and eFloras-only donors, because most source traits have only one corpus of documentation.

Source code: [protected original NHM text-corpus replay](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/fd_nhm_original_provider_audit.py), [synthetic source conflict and provenance tests](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/tests/test_fd_nhm_original_provider_audit.py), [original provider-blob SHA-verifying workflow](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/.github/workflows/run-fd-nhm-original-provider-audit.yml). The exact source code branch is `analysis/fd-discrete-functional-volume-20261113` (links above must use that source branch); original island names and plant entities remain unpublished. Keep fd PR #16 and island PR #326 **DRAFT**.
