# Broad-web Medium acquisition and genus-level inference

Research date: 2026-07-25

## Executive recommendation

The next material reduction in unresolved species-axis cells should come from a
domain-first web acquisition lane, not from another undirected search-engine
wave.

1. Harvest exact-species descriptions exposed through GBIF ChecklistBank and
   source DwC-A archives, then WFO content packages and an authorized POWO data
   route.
2. Add adapters for institutional floras, botanical gardens, herbaria, and
   government species accounts. Use their sitemaps, APIs, CSVs, or dataset
   archives rather than issuing one search-engine query per species.
3. Expand scholarly full text through Europe PMC, CORE, Crossref/OpenAlex plus
   Unpaywall, Zenodo, and Plazi/BLR. Supplements are especially important for
   reproductive assurance.
4. Use general web search only for the residual set and only as discovery.
   Fetch, validate, and cite the publisher page or PDF; never treat a search
   snippet as evidence.
5. Apply botanical-Latin and multilingual extraction to all acquired text.
6. Recompute the repository's masked genus-consensus contract after the new
   direct-evidence wave. Emit only rules newly enabled by the added direct
   records as Validated Low; never reclaim pre-existing rules as a new gain.
   Keep direct High/Medium coverage and direct-plus-Validated-Low coverage
   separately reported. Any more permissive genus model belongs in a separate
   operational sensitivity layer, never in High or Medium.

The highest-confidence immediate opportunity is GBIF's species-description
surface. A live probe found an exact species treatment for *Abarema
cochliacarpos*, including the statement `corolla 4.1–7.8 mm long`, linked to
the original Phytotaxa DOI. The important implementation detail is that the
description was attached to a checklist usage, not to the GBIF Backbone usage.
Acquisition therefore needs to inspect exact accepted-name and synonym usages
across checklist datasets and deduplicate them by original source lineage.

POWO is potentially higher-yield for flower colour and structure, but direct
automation from the current environment encountered a Cloudflare managed
challenge. It should be accessed through an authorized download/API agreement,
or indirectly through licensed WFO/GBIF content, rather than by escalating
scraping.

## 100-species implementation pilot

The implemented pilot selected 100 unique unresolved species across 91 genera,
balanced as 34 structural, 33 colour, and 33 reproductive-assurance targets.
The first direct-URL attempt showed why a source index is necessary: all 100
Pladias requests reached the read timeout, whereas Go Botany returned useful
structured fields.

The accepted method now downloads the
[Go Botany sitemap](https://gobotany.nativeplanttrust.org/sitemap.txt) once,
intersects its exact `/species/{genus}/{epithet}/` paths with accepted
binomials, and fetches only those matches. Search results and snippets are not
used. The corrected pilot completed in 41.6 seconds with:

- 45 exact sitemap matches and 45 fetched pages;
- 75 species-direct Medium evidence rows from 30 species;
- 29 new species-axis fills: 22 structure, 1 colour, and 6 reproductive
  assurance;
- 2 existing Validated Low cells upgraded to Medium;
- zero retrieval errors;
- one ambiguous page value that listed both radial and bilateral symmetry
  rejected rather than forced to either state.

Across the fixed baseline, the sitemap contains 3,520 exact species pages and
2,514 intersect species that still have at least one unresolved axis. This
makes one full Actions run feasible without issuing requests for the other
roughly 104,000 unresolved species.

Re-running the existing masked genus-consensus contract on the pilot's new
direct evidence enabled one additional unanimous rule. It yielded 167
incremental structural Validated Low cells, all reported separately from the
29 direct Medium gains. Baseline rules were evaluated on the same unresolved
set and subtracted, so none were relabelled as a new contribution.

## Scope and evidence contract

The denominator remains 106,295 accepted species × three axes:

- flower colour;
- floral structure;
- reproductive assurance.

The aggregation key is `accepted_species × axis`. Species-direct evidence keeps
the precedence High > Medium > Validated Low. A new web record can enter Medium
only when all of the following are true:

1. **Exact taxon resolution.** The treatment, table row, page, or article
   explicitly concerns the WFO-accepted species or one of its exact synonyms.
   Fuzzy matches, abbreviated binomials without disambiguating context,
   unmatched infraspecific names, hybrids, and cultivar-only pages do not pass.
2. **Explicit trait statement.** The source states the value. Images, common
   names, pollinator syndromes, neighboring species, or generic genus
   descriptions do not substitute for text or structured data.
3. **Correct organ and assertion.** For example, `red fruit`, `purple leaves`,
   and `white bracts` are not flower-colour evidence unless the controlled
   ontology explicitly includes the stated organ. Negation, comparison,
   specimen-host, and "unlike species X" contexts must be rejected or reviewed.
4. **Deterministic ontology mapping.** The extracted phrase maps without
   guessing to the target axis. Ambiguous or incomplete statements remain
   pending.
5. **Species-level editorial provenance.** A peer-reviewed paper, taxonomic
   treatment, flora, herbarium, government account, botanical institution, or
   similarly accountable source may support Medium. Commercial, vendor,
   community, unattributed, and generic blog pages remain pending or excluded.
   Two rehosts of the same source do not improve quality.
6. **Reproducible evidence.** Store the accepted name, matched source name,
   source URL, original citation/DOI/treatment identifier, page or section
   locator, minimal supporting quote, retrieval time, content hash, extraction
   version, language, licence/rights state, and taxon-match reason.
7. **Lineage deduplication.** Deduplicate on the original DOI, ISBN/treatment
   UUID, dataset record, or canonical publication. If no stable identifier is
   available, use a normalized text/page hash plus citation fingerprint.
   Provider names such as GBIF, Plazi, POWO, and WFO are distribution channels,
   not automatically independent evidence lineages.

One authoritative exact-species treatment can support Medium. Multiple weak
pages should not be promoted merely by vote count. Conflicts remain visible
and unresolved unless the quality and assertion rules resolve them
deterministically.

## Acquisition plan

### 1. GBIF ChecklistBank descriptions and source archives

The [GBIF Species API](https://techdocs.gbif.org/en/openapi/v1/species)
documents species search and the `/species/{usageKey}/descriptions` endpoint.
A live exact-name probe on 2026-07-25 showed that checklist usages can expose
full species treatments with source citations even when the Backbone usage has
no description.

Recommended workflow:

1. Resolve the accepted name and all exact synonyms locally from WFO.
2. Search GBIF for each name and retain only exact canonical-name matches with
   rank species and kingdom Plantae.
3. Query descriptions for every qualifying usage, not only the Backbone key.
4. Record dataset key, usage key, description source, original citation/DOI,
   language, and text hash.
5. Once a productive dataset is identified, download its source DwC-A or
   ChecklistBank export once and process it locally. Do not make millions of
   redundant record calls.

The live search also returned unrelated homonyms and non-plant results.
Canonical equality, rank, kingdom, accepted/synonym resolution, and source-name
matching are mandatory.

The exact identifiers from the 2026-07-25 live probe are:

- GBIF taxon usage key:
  [`211579716`](https://api.gbif.org/v1/species/211579716), an accepted
  species-level checklist usage of *Abarema cochliacarpos* with Plazi treatment
  taxon ID `03A5194AFF97FFF5B0E22FDF5E59F94E.taxon` and treatment URL
  `https://treatment.plazi.org/id/03A5194AFF97FFF5B0E22FDF5E59F94E`;
- GBIF dataset key:
  [`b31ca15a-de1e-4e51-bde8-7bc74455535b`](https://api.gbif.org/v1/dataset/b31ca15a-de1e-4e51-bde8-7bc74455535b);
- dataset title: *Circumscription of Abarema (Leguminosae, Caesalpinioideae,
  Mimosoid clade)*;
- dataset DOI: [`10.15468/g97d99`](https://doi.org/10.15468/g97d99), the Plazi
  checklist dataset DOI, which is a provider/dataset identifier rather than the
  original publication lineage;
- description record key:
  [`435641593`](https://api.gbif.org/v1/species/211579716/descriptions), type
  `description`, language `eng`;
- description source/citation as returned by GBIF: Guerra, Ethiéne; Soares,
  Marcos Vinicius Batista; Morim, Marli Pires; Iganci, João Ricardo Vieira
  (2023), “Circumscription of *Abarema* (Leguminosae, Caesalpinioideae,
  Mimosoid clade),” *Phytotaxa* 601(1): 51–60;
- original article DOI and deduplication lineage:
  [`10.11646/phytotaxa.601.1.3`](https://doi.org/10.11646/phytotaxa.601.1.3).

The record text contains the exact assertion `corolla 4.1–7.8 mm long`. The
original article DOI, not the GBIF dataset DOI or Plazi provider name, is the
primary lineage key for deduplicating copies of this treatment.

This lane is the first implementation priority because it is open, structured,
lineage-rich, and already demonstrated to contain exact treatment text useful
for the structure axis. It can also surface Plazi treatments and national flora
datasets without separate site-specific crawlers.

### 2. WFO content and POWO

[World Flora Online contributor guidance](https://about.worldfloraonline.org/images/uploads/documents/WFOGuidelinesforContentDataContributorsV._2.04.pdf)
defines description records with taxon identifiers, language, description type,
text, and source. Use the existing local WFO snapshot for accepted names and
aliases, and obtain content packages or bounded GraphQL/portal records for
descriptions. WFO should retain the original content source as lineage; WFO
name acceptance is not independent trait evidence.

[Plants of the World Online](https://powo.science.kew.org/) reports more than
1.4 million plant names and hundreds of thousands of detailed descriptions.
Its [search help](https://powo.science.kew.org/search-help) documents name,
location, and characteristic searches, and its
[about page](https://powo.science.kew.org/about) explains that descriptive
data are attached to an accepted WCVP name while retaining the original name
and source.

POWO should be pursued through an authorized bulk download, an agreed API
route, or a licensed redistribution already exposed through WFO/GBIF. The
undocumented backend route used by some clients returned HTTP 403 with a
Cloudflare managed challenge in this environment on 2026-07-25. That is an
operational stop signal for blind scraping, not a reason to bypass site
controls. Cite POWO following Kew's
[citation guidance](https://powo.science.kew.org/cite-us), and preserve the
underlying source lineage when present.

### 3. Institutional flora and species-account adapters

After the two broad aggregators, build small domain adapters. Each adapter
should enumerate source records through a sitemap, API, CSV, dataset download,
or stable search endpoint, then run the same exact-name and evidence gates.
Promising targets include:

- [NSW PlantNET](https://www.botanicgardens.org.au/our-science/our-services/nsw-flora-online-plantnet),
  which publishes taxonomic descriptions and images. Its
  [Data NSW record](https://data.nsw.gov.au/data/dataset/plantnet) exposes the
  dataset and licence information.
- [PlantZAfrica](https://pza.sanbi.org/about), with authored institutional
  species articles. Its [copyright terms](https://pza.sanbi.org/copyright)
  restrict copying and server redistribution, so retain short evidence
  excerpts and provenance rather than mirroring full articles.
- [Flora do Brasil 2020](https://www.gbif.org/dataset/aacd816d-662c-49d2-ad1a-97e66e2a2908),
  available through GBIF with a CC BY licence and suitable for source-archive
  processing.
- [NC State Extension Gardener Plant Toolbox](https://plants.ces.ncsu.edu/),
  whose [help page](https://plants.ces.ncsu.edu/help/) documents structured
  flower attributes.
- [Missouri Botanical Garden Plant Finder](https://www.missouribotanicalgarden.org/PlantFinder/PlantFinderSearch.aspx),
  which exposes structured bloom descriptions for thousands of plants.
- [US Forest Service FEIS](https://research.fs.usda.gov/feis) and its
  [species reviews](https://research.fs.usda.gov/feis/species-reviews), which
  provide reviewed ecological syntheses and downloadable metadata.
- USDA PLANTS plant guides, when an exact species PDF contains explicit flower
  or reproduction statements.

Horticultural institutional pages require an extra wild-species check. Exclude
cultivar, hybrid, and garden-selection descriptions unless the asserted value
is explicitly stated for the accepted wild species.

Rather than hand-curating a fixed list forever, use productive general-search
results to discover additional institutional domains. Promote a domain into
the adapter tier only after reviewing its editorial identity, taxon specificity,
terms, robots policy, and record stability.

### 4. Scholarly full text, supplements, and repositories

Reproductive assurance will rely more heavily on papers and supplementary
tables than colour and structure.

- [Europe PMC REST](https://europepmc.org/RestfulWebService) provides article
  search, open-access full-text XML, and supplementary-file services.
- [CORE](https://core.ac.uk/services/api) provides machine access to a large
  corpus of repository metadata and open full text and can recover papers
  outside PubMed Central.
- [Crossref REST](https://www.crossref.org/documentation/retrieve-metadata/rest-api/)
  supplies DOI metadata, abstracts, and publisher links. Its
  [text-and-data-mining documentation](https://www.crossref.org/documentation/retrieve-metadata/rest-api/text-and-data-mining/)
  makes clear that a full-text link does not itself grant access; publisher
  rights and rate limits still apply.
- [OpenAlex](https://developers.openalex.org/api-reference/introduction) adds
  scholarly discovery, concepts, citations, and open-access locations. Current
  [authentication and usage rules](https://developers.openalex.org/api-reference/authentication)
  should be incorporated into the acquisition manifest.
- [Unpaywall](https://unpaywall.org/products/api) resolves DOI-level legal
  open-access locations; it is a resolver, not a full-text search engine.
- [Zenodo REST](https://developers.zenodo.org/) exposes records and files.
  Prefer OAI-PMH or metadata exports for bulk discovery.
- [Plazi TreatmentBank](https://plazi.org/treatmentbank/treatment-data-access/)
  and the [Biodiversity Literature Repository](https://plazi.github.io/blr/index.html)
  expose taxonomic treatments and source DOI lineage.

Query accepted names and exact synonyms. Extract table cells and supplement
rows with the same taxon rigor as prose. For reproductive assurance, distinguish
at least self-incompatibility, autonomous selfing, mating system, apomixis, and
cleistogamy. Do not map one to another simply because they are correlated.

### 5. Residual general-web discovery

The [Brave Web Search API](https://api-dashboard.search.brave.com/app/documentation/web-search/get-started)
supports quoted terms, `site:`, `filetype:pdf`, language/country options, and
additional snippets. It is suitable for the residual set, especially for exact
binomials plus trait terms in several languages. Its
[pricing](https://api-dashboard.search.brave.com/documentation/pricing) makes
domain-first harvesting more economical than multiple queries for every one
of 106,295 species.

Its [terms](https://api-dashboard.search.brave.com/documentation/resources/terms-of-service)
also matter: search results are for transient use and should not become a
stored derivative search database. The pipeline should store only the fetched
publisher URL/content and its own provenance, not cached result pages, ranks,
or snippets. A search snippet is discovery material and can never be the
supporting quote.

Google Custom Search is not a new long-term foundation: the
[official overview](https://developers.google.com/custom-search/v1/overview)
states that it is closed to new customers and gives an end date for existing
JSON API customers.

[Common Crawl indexes](https://index.commoncrawl.org/) are useful for finding
URLs under already known flora, herbarium, and garden domains, or for creating
a permitted offline domain corpus. They are not a practical full-text
species-query engine. Recover the original publisher URL and rights state, and
do not treat a crawl copy as an independent source.

Suggested residual query families include:

```text
"<accepted binomial>" ("flower color" OR corolla OR petals)
"<accepted binomial>" ("flower shape" OR floral symmetry OR corolla tube)
"<accepted binomial>" ("self-compatible" OR "self-incompatible")
"<accepted binomial>" (autogamy OR "autonomous selfing" OR "mating system")
"<exact synonym>" <same trait pack>
site:<institutional-domain> "<binomial>"
filetype:pdf "<binomial>" <trait term>
```

General queries should be run only after domain and literature lanes, and only
for unresolved species-axis pairs.

### 6. Historical literature and public PDFs

The [Biodiversity Heritage Library API](https://www.biodiversitylibrary.org/docs/api3.html)
supports metadata and OCR access. BHL's
[developer and data tools](https://about.biodiversitylibrary.org/tools-and-services/developer-and-data-tools/)
also provide bulk OCR and public-cloud datasets covering tens of millions of
pages. It is particularly valuable for old floras and botanical Latin.

OCR output should produce candidates, not automatic Medium evidence. A
candidate passes only after the species treatment boundary, original page
image, assertion, and correct organ have been verified. Page ID, volume,
bibliographic lineage, OCR text hash, and image locator must be retained.

## Multilingual and botanical-Latin extraction

Run language identification by treatment or paragraph and maintain
language-specific lexicons and negation/organ rules. Initial high-yield
languages are English, Spanish, Portuguese, French, German, Italian, and
botanical Latin.

Examples of Latin morphology terms:

- organs: `flos`, `flores`, `corolla`, `petala`, `tubus corollae`;
- colours: `albus/alba`, `luteus/lutea`, `ruber/rubra`, `roseus/rosea`,
  `caeruleus/caerulea`, `violaceus/violacea`, `viridis`, `brunneus/brunnea`;
- forms: `tubulosus/tubulosa`, `campanulatus/campanulata`,
  `infundibuliformis`, `rotatus/rotata`, `urceolatus/urceolata`,
  `bilabiatus/bilabiata`, `papilionaceus/papilionacea`,
  `zygomorphus/zygomorpha`, `actinomorphus/actinomorpha`.

Adjective agreement and local syntax matter: a colour adjective near the
binomial can still modify fruit, seed, leaf, calyx, or bract. Require
co-occurrence within the floral assertion and an explicit mapped organ.

Reproduction query packs should include language-specific equivalents of
self-compatible/self-incompatible, autogamous/allogamous, autonomous selfing,
apomixis, and cleistogamy. These terms remain separate ontology concepts until
the source explicitly supports the target output value.

## Rights, robots, and respectful acquisition

Each source adapter should have a machine-readable manifest:

```yaml
domain: example.org
enumeration_method: sitemap|api|csv|dwca|oai_pmh|search
robots_checked_at: timestamp
terms_checked_at: timestamp
licence: SPDX-or-source-text
tdm_reservation: allowed|reserved|unknown
rate_limit: requests-per-second
redistribute_source_text: true|false
contact: project-contact
```

Follow the Robots Exclusion Protocol in
[RFC 9309](https://www.rfc-editor.org/info/rfc9309/) and inspect explicit
text-and-data-mining reservations, including the
[W3C TDM Reservation Protocol](https://w3c.github.io/tdm-reservation-protocol/spec/).
Use an identifying user agent with contact information, conditional requests,
bounded concurrency, caching, exponential backoff, and source-specific rate
limits.

Record the applicable [Creative Commons licence](https://creativecommons.org/share-your-work/use-remix/cc-licenses/)
or source terms. If rights are unknown or redistribution is restricted, retain
the URL, retrieval timestamp, source hash, page/section locator, licence state,
and only the minimal evidence quotation needed for audit. Do not publish a
mirrored full-text corpus.

## Genus-level inference: scientific boundary

Genus inference is not species-direct evidence. It must never enter the High,
Medium, or Validated Low counts and must never silently fill the direct
coverage matrix.

Use a separate output tier, for example:

```text
quality = Genus-inferred
evidence_scope = genus_inference
confidence = low
source_backed = false
value = likely_<state>
```

Report two metrics:

- **direct coverage:** High + Medium + Validated Low;
- **operational coverage:** direct coverage plus separately flagged
  genus-inferred cells.

Genus inference can reduce `operational_unresolved`, but not
`direct_unresolved`.

### Traits that should not be genus-inferred

**Mating system and autonomous selfing.** A 105-species synthesis found broad
among-population variation in plant mating systems
([Whitehead et al. 2018](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2018.00038/full)).
Intraspecific selfing variation is itself an important evolutionary phenomenon
([Moeller et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34632564/)).
Neither trait should be inferred to a species from congeners.

**Flower colour.** Colour repeatedly evolves within genera and may be
polymorphic within species. Examples include repeated transitions in
[*Nicotiana*](https://academic.oup.com/aob/article/115/7/1117/173374),
red/yellow transitions within monkeyflowers
([Nature Communications, 2025](https://www.nature.com/articles/s41467-025-57639-3)),
and repeated colour changes in Antirrhineae
([Landis et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4904171/)).
Disable genus-level colour inference by default. Apparent donor unanimity can
easily be a sampling artifact.

### Traits eligible only after diagnostics

**Self-incompatibility/self-compatibility.** This can be phylogenetically
conserved, but repeated loss occurs even within small clades. *Linanthus*
section Leptosiphon contains both systems with several independent transitions
([Goodwillie 1999](https://academic.oup.com/evolut/article/53/5/1387/6758091)).
Repeated loss is also documented in
[*Capsella*](https://pmc.ncbi.nlm.nih.gov/articles/PMC2659713/) and natural
variation occurs in
[*Arabidopsis*](https://pmc.ncbi.nlm.nih.gov/articles/PMC528763/).
Inference is acceptable only as `likely_SC` or `likely_SI`, after unanimous
direct donors and out-of-sample validation. It must not imply mating system or
autonomous selfing.

**Specific floral-structure states.** There is no defensible blanket
genus-level rule for the whole structure axis. One community study found weak
phylogenetic signal for flower shape but stronger signal for tube depth and
flower-cluster traits
([Thompson et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10500814/)).
Symmetry can vary among species and within individuals, as in
[*Saxifraga*](https://link.springer.com/article/10.1007/s00606-023-01842-6/),
whereas other clades such as Malpighiaceae can show strong floral conservatism
([Davis et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2851953/)).
Evaluate the exact subtrait—symmetry, form, tube presence, tube depth, flower
size, or display—independently. Never use evidence for one subtrait to fill
another.

Quantitative traits require a predeclared discretization and stability test.
For example, donor tube-depth values must have a sufficiently tight
distribution that leave-one-out discretized classes remain stable; a genus
mean alone is not evidence for an unmeasured species.

## Required genus-inference diagnostics

The existing conservative rule of at least five unanimous direct congeners is
a useful floor for `likely_SC`/`likely_SI`, but it is not sufficient for new
trait families. Recompute it after the new direct-evidence integration and
apply these diagnostics per genus and exact trait:

1. **Direct donors only.** Donors must be species-direct High, Medium, or
   Validated Low. Never allow inferred donors.
2. **One species, one vote.** Resolve its lineages and conflicts before it can
   vote. Multiple records or rehosts do not add weight.
3. **Minimum support.** Keep `n ≥ 5` as the absolute SI/SC minimum. Prefer
   `n ≥ 8` for any tested structure trait. Flower colour remains disabled.
4. **Genus sampling.** Require a predeclared observed fraction of the accepted
   genus, such as at least 30%, alongside the absolute count. Large genera need
   a higher absolute sample or a shrinkage model; five donors in a
   500-species genus are not representative.
5. **State concentration.** Require 100% donor agreement for SI/SC. A
   structure-trait threshold should be estimated from held-out data; until
   validated, use unanimity rather than an arbitrary majority.
6. **No direct conflict.** Do not infer when the target has any unresolved
   direct evidence or when donors are variable, mixed, or internally
   conflicting.
7. **Leave-one-species-out validation.** Predict every directly observed
   donor from the other donors in its genus. Evaluate by exact trait and genus,
   not only as a pooled headline number.
8. **Uncertainty floor.** Require point accuracy ≥ 0.95 and a predeclared lower
   95% confidence bound, for example Wilson lower bound ≥ 0.80. For sparse
   groups, consider beta-binomial shrinkage rather than treating a small
   unanimous sample as certainty.
9. **Representativeness.** Record donor geography, sections/subgenera when
   available, and flags for hybrids, cultivars, polyploids, infraspecific
   variation, and taxonomic instability.
10. **Target eligibility.** The target must be a stable WFO-accepted species in
    the same accepted genus, not a hybrid, unplaced name, or unresolved complex.
11. **Audit record.** Store donor species IDs, resolved values, evidence
    lineages, donor-set hash, state distribution, accepted genus size,
    observed fraction, entropy, exclusions, leave-one-out prediction and
    observed value, accuracy, and confidence interval.
12. **Direct override.** Any later species-direct evidence replaces the
    genus-inferred operational value. A conflict returns the cell to direct
    review rather than preserving the inference.

The previously reported pooled leave-one-out result for conservative
self-compatibility inference—372/384, or 96.9%—is encouraging but must be
recomputed on the fully integrated corpus and reported by trait and genus. It
does not justify flower-colour, mating-system, or general structure inference.

## Implementation order by expected return

| Priority | Lane | Expected value | Main constraint |
|---|---|---|---|
| 1 | GBIF exact checklist usages, descriptions, then source DwC-A | High; live-validated exact treatment text and strong lineage | Must search beyond Backbone and strictly filter homonyms |
| 2 | WFO content packages and authorized POWO access | Very high for colour and structure | POWO automation currently challenged; permission/download route needed |
| 3 | Institutional flora/garden/herbarium/government adapters | High, geographically complementary, stable editorial provenance | Domain-specific parsers, cultivar filtering, heterogeneous rights |
| 4 | Europe PMC + CORE + Crossref/OpenAlex/Unpaywall + supplements | High for reproduction; moderate for morphology | Full-text rights, table parsing, name disambiguation |
| 5 | Brave residual discovery and productive-domain promotion | Moderate long-tail recovery | Query cost and search-result storage restrictions |
| 6 | BHL OCR and historical/public PDFs | Valuable long tail and Latin descriptions | OCR/treatment-boundary verification burden |
| 7 | Recomputed genus inference for SI/SC and validated exact structure traits | Potentially thousands of operational cells | Must remain separate from direct coverage and pass held-out diagnostics |

## Suggested acceptance tests and outputs

Before scaling, construct stratified gold sets by axis, language, source class,
and taxonomic-name condition. Include positives, negatives, organ confusions,
negation, synonyms, hybrids, cultivar pages, comparison sentences, OCR errors,
tables, and mixed-value species.

For each adapter and extractor, report:

- exact-name precision and rejection reasons;
- assertion precision by axis and language;
- ontology-mapping precision;
- lineage duplicate rate;
- conflict rate;
- pending rate;
- HTTP/status and robots/rights exclusions;
- accepted Medium cells per 1,000 fetched species treatments;
- unique accepted species-axis gain after High > Medium > Validated Low
  precedence;
- direct-Low-to-Medium upgrades;
- source-specific and source-lineage-specific pure gain.

The final artifacts should include:

1. a source manifest with run ID, code SHA, query/adapter version, retrieval
   interval, source rights, and checksums;
2. raw candidate evidence with status and exclusion reason;
3. accepted direct evidence with exact quote and locator;
4. lineage clusters and provider duplicates;
5. the recomputed 106,295 × 3 direct coverage matrix;
6. a separate genus-inference diagnostics table and operational matrix;
7. before/after High, Medium, Validated Low, direct unresolved, operational
   unresolved, axis rates, 0/1/2/3-axis species counts, and per-run pure gain.

## Bottom line

A broad web pass is scientifically defensible when it is broad in **sources**
but narrow in **acceptance**. The scalable strategy is to discover productive
institutional domains and bulk archives, fetch their records once, and apply a
single fail-closed species-direct Medium gate. Search engines are the long-tail
locator, not the evidence store.

Genus inference can be revisited after this direct wave, but it cannot be used
to make the direct coverage number look better. Keep it as a low-confidence,
fully audited operational overlay. Restrict the first production use to
unanimous, held-out-validated `likely_SC`/`likely_SI` and any exact
floral-structure subtrait that independently passes the same diagnostics.
Disable genus inference for flower colour, mating system, and autonomous
selfing.
