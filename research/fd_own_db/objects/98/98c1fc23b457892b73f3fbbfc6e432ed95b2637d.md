# Phylogenetic beta vs floral categorical functional beta — frozen original macro design

**Scope and inference boundary.** This is a cross-sectional observational, source-conditional macroecological relationship; the 2,234 source-direct flowering plant species have ONE fixed recorded categorical flower phenotype per accepted name, with no island-population phenotype measurements, visitor decline or effective pollination outcomes. Original 928 GSHHG island flower floras are the eligibility universe, not 928 independent colonizations.

## Comparison and two genuinely distinct beta measures

Each island and geographically potential five-nearest GIFT **native mainland source** is restricted to accepted plant species present as exact, unmodified published dated GBOTB tree tips and possessing original high/medium Direct colour, form, symmetry. Every identical focal island species is removed from the mainland pool before drawing. Each side contributes precisely **15 distinct accepted species** in each of **20 independent taxon rarefaction draws**; the same 15-plus-15 draw supplies all metrics, and repeated draws are **averaged to one row per physical island**.

The primary explanatory measurement is **PhyloSor phylogenetic beta**: `1 - 2*L_shared/(L_island + L_source)`, using *actual original rooted published branch lengths* and each sample's union of ancestral edges. No imputed genus/family species phylogenetic grafts. This differs from prior one-to-one nearest-tree donor distance and does not require source assignment to be an actual ancestor.

The primary response is **functional beta of original categorical floral support**: a Jaccard dissimilarity `1 − |B_island ∩ B_source|/|B_island ∪ B_source|` of the globally fixed reference consisting of **100 actually recorded** flower-colour × form × symmetry triple profiles, covered by a one-categorical-mismatch Hamming ball around each sampled floral profile. No arbitrary continuous 17D KDE volume. The secondary outcomes are exact-profile Jaccard and signed island-minus-source rare-minus-common neighborhood coverage corrected by an original exactly occupied-profile-k reference null.

All three phylogenies are separately source-pinned to published V.PhyloMaker2 GBOTB.extended TPL/LCVP/WP original dated tip data (74,529 / 73,420 / 72,570 published original tips). Their original Direct taxa with strict independent source-family accepted-name agreement are **1,187 / 1,162 / 1,160**, respectively; insertable but absent taxa are **excluded**, never treated as molecular measurements. CoL XR source taxonomy exact accepted-family correspondences remain unreviewed automatic candidates.

## Geography, taxonomy and observational effort controls

Per region and provenance scope, require at least **20 independent physical islands**, **10 physical 10-degree spatial blocks** and **10 distinct first-nearest GIFT geographic source signatures** before reporting a directional estimate. Equal-weight original physical block and donor signature conditional cluster resampling, **999 draws each**. In each original region, freeze one source-conditional island response row before regression.

Main test is phylogenetic beta vs functional beta, expressed as change in functional beta per **one original within-region SD of PhyloSor beta**. Report raw association and a separate partial model controlling the sampled **genus- and candidate-family Jaccard taxonomic turnover**, original source exact-tip plant richness (`log1p`), and original island log area and log mainland isolation. Do not silently delete a source covariate if design becomes rank deficient; mark the coefficient unestimable. Source- and block-equal cluster results must be compared, not counted as independent biological confirmations.

One physical island enters the regression once, regardless of its species or rarefaction opportunities. Original 20 draws are Monte Carlo precision replicates, not island/colonist biological replications. This guards the previous high-dimensional KDE pitfall, as neither PhyloSor shared branches nor functional Jaccard is defined by how many recorded floral profiles a sampling source happened to provide at unstandardized sample sizes.

## Interpretation stop rules

- A positive slope of original PhyloSor beta on functional beta is an **association between phylogenetic and floral composition turnover**, not a historical or causal effect of isolation/flower evolution.
- If a raw positive association loses support after genus/family turnover adjustment, then **coarse taxonomic assembly or source ascertainment** remains an adequate explanation; a phylogeny-independent functional mechanism was NOT identified.
- If the signed rare-minus-common neighborhood gap is related to phylogenetic beta, it remains a documented original category-space pattern; the globally rare 29 recorded profiles are source Direct-frequency-defined, not observed rare pollination functions. The global common profiles are near a geometric saturation ceiling.
- **Mismatched true tree-tip sampling**, native-only sample attrition and five-nearest GIFT candidate source choice limit external generality. At most the source-independent whole-block/source conditional models will support a selected-subcohort macro description.
- The three GBOTB tree versions share much of the original construction and are **not three independent studies**. Exploratory multiple endpoints and regional models are not confirmatory inferential evidence without a fresh independent taxonomic/phylogenetic source.

Code: [`src/fd/phylo_beta_functional_beta_macro.py`](../src/fd/phylo_beta_functional_beta_macro.py), [`tests/test_phylo_beta_functional_beta_macro.py`](../tests/test_phylo_beta_functional_beta_macro.py), protected original [`island/scripts/run_fd_original_phylo_functional_beta.py`](https://github.com/zuizui0223/island/blob/analysis/fd-discrete-functional-volume-20261113/scripts/run_fd_original_phylo_functional_beta.py). Original source-only CI uploads *anonymous regional statistics*; no accepted plant names, GIFT source-island IDs, original family/genera, exact coordinates or original tree RDA file copies.
