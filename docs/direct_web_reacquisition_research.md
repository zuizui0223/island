# Direct local web reacquisition for unresolved three-axis traits

## Recommendation

The first local lane should be:

1. build an exact accepted-name/synonym index from the fixed WFO Plant List
   2026-06 snapshot;
2. query Plazi TreatmentBank with every exact WFO name for species treatments;
3. query Europe PMC for open-access full text and supplementary tables, especially
   for reproductive assurance; and
4. use BHL page OCR only as a rights-checked, image-verified fallback.

This can run from a workstation without GitHub Actions. WFO is a name-resolution
source, not independent trait evidence. No fuzzy match, family inference, global
fallback, or unreviewed auto-promotion is permitted.

## 1. WFO: make synonym recovery local and reproducible

Use the June 2026 WFO Plant List release, DOI
[`10.5281/zenodo.20782718`](https://doi.org/10.5281/zenodo.20782718).
WFO releases a new snapshot every June and December. The release contains a
Catalogue of Life package, JSON, SQL, Darwin Core archives, deprecated-name and
deduplicated-ID lookup tables; the record is CC0. The JSON is explicitly
recommended by WFO for a local read-only copy of the API index
([official release and file inventory](https://zenodo.org/records/20782718),
[license](https://zenodo.org/records/20782718#rights)).

For the alias lane, download once:

```text
https://zenodo.org/api/records/20782718/files/_uber.zip/content
https://zenodo.org/api/records/20782718/files/deduplicated_ids_lookup.csv.gz/content
https://zenodo.org/api/records/20782718/files/deprecated_names_lookup.csv.gz/content
```

`_uber.zip` includes accepted, synonym, and deprecated names. Its release MD5 is
`0d95b8e4ab12e3d7dac8d0aa07747e10`; the acquisition manifest should also record
a locally computed SHA-256. If deprecated names are not wanted, the smaller
`_DwC_backbone_R.zip` contains the non-deprecated backbone.

Join the project `accepted_species` to WFO only by exact canonical name plus rank
and, when present, authorship. Resolve duplicate WFO IDs through the prescribed-ID
lookup. Retain:

```text
accepted_species
search_name
name_role                 # accepted | synonym
observed_wfo_id
prescribed_wfo_id
accepted_wfo_id
versioned_usage_id        # e.g. wfo-0001333299-2026-06
wfo_release
wfo_release_doi
archive_sha256
```

The public API is useful for spot checks, not for 106,295-species traversal:

```text
GET https://list.worldfloraonline.org/matching_rest.php
    ?input_string=Mimulus%20guttatus
    &fuzzy_names=0
    &fuzzy_authors=0
    &check_homonyms=true
    &check_rank=true

POST https://list.worldfloraonline.org/gql.php
Content-Type: application/json

{"query":"query { taxonNameById(nameId: \"wfo-0000450805\") { id stableUri fullNameStringPlain currentPreferredUsage { id hasName { id stableUri fullNameStringPlain } hasSynonym { id stableUri fullNameStringPlain } } } }"}
```

The API requires no key. WFO publishes no numeric request limit, but warns that
throttling or tokens may be introduced if service load requires it and recommends
local data for heavy use. Therefore use the snapshot for production and no more
than one serial request per second for bounded diagnostic checks
([WFO documentation hub](https://plant-list-docs.rbge.info/),
[WFO Plant List service](https://list.worldfloraonline.org/index.php)).

The lineage of a name-resolution assertion is
`wfo:10.5281/zenodo.20782718:<observed_wfo_id>`. It must never be counted as a
trait-evidence lineage.

## 2. Europe PMC: OA XML and supplementary tables

The official REST root is:

```text
https://www.ebi.ac.uk/europepmc/webservices/rest
```

For each accepted species, group the accepted name and exact WFO synonyms into
bounded OR clauses. Run a high-precision title/abstract query first, followed by
an unfielded OA query for full-text-only mentions:

```text
GET /search?query=OPEN_ACCESS:Y AND
  ("Erythranthe guttata" OR "Mimulus guttatus") AND
  (selfing OR outcross* OR autogam* OR "self-incompatib*" OR
   "self-compatib*" OR "breeding system" OR "mating system")
  &resultType=core&pageSize=100&format=json
```

`HAS_SUPPL:Y` can be added for a supplementary-table pass. Use the returned
`nextCursorMark` unchanged when pagination is required. Search metadata should
retain `source`, `id`, `pmid`, `pmcid`, `doi`, `isOpenAccess`, `hasSuppl`,
`license`, `dateOfCreation`, and `dateOfRevision`
([official REST reference](https://dev.europepmc.org/RestfulWebService)).

Fetch only candidates with a PMCID:

```text
GET https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML
GET https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles?includeInlineImage=false
```

`fullTextXML` returns JATS/NLM XML for the open-access subset and
`supplementaryFiles` returns a ZIP
([official endpoint definitions](https://dev.europepmc.org/RestfulWebService)).
Parse paragraphs, captions, tables, and supplementary CSV/TSV/XLSX files. A
candidate is species-direct only when the accepted name or a WFO-verified synonym
occurs in the same treatment, paragraph, caption, or table row as the extracted
trait statement. Merely appearing elsewhere in the article is insufficient.

Europe PMC permits automated retrieval of OA content through REST, SOAP, OAI, and
FTP, but does not permit automated bulk download of other content. OA articles
remain copyrighted and licenses differ by article, so preserve both the search
record's `license` and the JATS `<license>` text
([developer policy](https://europepmc.org/developers),
[OA subset and licensing](https://europepmc.org/downloads/openaccess)).
There is no documented numeric REST request rate. For this pilot use serial
requests at 1 request/second, send a descriptive User-Agent and optional `email`,
honour `Retry-After`, and exponentially back off on 429/5xx.

Evidence lineage is the normalized original DOI, otherwise PMCID:

```text
doi:<normalized-doi>
pmcid:<PMCID>                       # DOI absent
```

Keep XML XPath or table/sheet/row as an evidence locator. A supplement gets a
subrecord identifier such as
`<article-lineage>#supp:<zip-member>:<sheet>:<row>`, but it remains the same
independent source lineage as its parent article.

## 3. Botanical descriptions that support automated retrieval

### Plazi TreatmentBank — first botanical-description source

Plazi provides an official API specifically for automated TreatmentBank access.
Its records include taxonomic treatments, tables, figures, citations, and
bibliographic references
([API overview](https://api.plazi.org/GgSrvApi/v1/aboutThis),
[official tutorial](https://api.plazi.org/GgSrvApi/v1/gettingStarted)).
No API key is required.

Exact species search and XML retrieval are immediately runnable:

```text
GET https://api.plazi.org/v1/Treatments/search
    ?genus=Abarema&species=cochliacarpos&format=JSON&limit=20

GET https://api.plazi.org/v1/Treatments/fetch
    ?UUID=03A5194AFF97FFF5B0E22FDF5E59F94E
```

Search every exact WFO name, but accept a treatment only when its XML taxonomic
name, rank, and accepted-species mapping agree. Extract only from the focal
`<treatment>` and its description/diagnosis sections. Plazi also exposes stable
treatment URIs and XML/RDF/DwC representations
([TreatmentBank access documentation](https://plazi.org/data-apis-tools/treatmentbank-api/)).

Use the original publication DOI from `ID-DOI`/`docSource` as lineage. If absent,
use `plazi:master-document:<masterDocId>`. The treatment UUID is a locator, not a
second independent source. This makes a Plazi copy of a Europe PMC or publisher
article deduplicate against the same DOI.

Plazi describes treatment data as open-access and FAIR, but the reuser remains
responsible for rights in the underlying publication. Preserve the document's
original DOI, source, language, and any treatment/figure license fields. The
official API gives no numeric request limit; use serial 1 request/second access
and cache raw XML.

### Biodiversity Heritage Library — fallback only

BHL API v3 supports canonical-name lookup, page metadata, full-text search, and
page OCR. It requires a free API key
([official API v3 documentation](https://www.biodiversitylibrary.org/docs/api3.html),
[developer tools](https://about.biodiversitylibrary.org/tools-and-services/developer-and-data-tools/)).

```text
GET https://www.biodiversitylibrary.org/api3
    ?op=GetNameMetadata&name={exact_name}&format=json&apikey={key}

GET https://www.biodiversitylibrary.org/api3
    ?op=GetPageMetadata&pageid={page_id}&ocr=t&names=t
    &format=json&apikey={key}
```

Metadata is CC0, but the scanned work and OCR inherit item-level rights. Check and
record `Item/Rights` and `Item/CopyrightStatus`; public-domain works remain public
domain, while some in-copyright BHL-hosted works are CC BY-NC-SA
([BHL copyright and reuse](https://about.biodiversitylibrary.org/help/copyright-and-reuse/)).
OCR text alone cannot be accepted. Promote to Medium only after the page image
confirms the exact species treatment and the excerpt; otherwise retain it as a
pending candidate.

Lineage is the original DOI when supplied, otherwise `bhl:item:<ItemID>`, with
`PageID` as the locator. This also deduplicates BHL pages reused by WFO, Plazi, or
another provider.

## 4. Promotion contract

All downloaded material first enters a candidate ledger. Promotion requires:

- exact `accepted_species` or exact WFO synonym-to-accepted mapping;
- an explicit focal-species mention in the same paragraph/table row/treatment;
- a controlled value already supported by the three-axis ontology;
- source URL, original citation/DOI, source-record ID, excerpt, raw-content hash,
  retrieval time, license, and source lineage;
- no family/genus inference, global fallback, or automatic acceptance of pending
  rows; and
- conflict retention rather than winner-by-provider.

Assign **High** only to explicit species-specific experimental or structured
evidence, such as an outcrossing estimate, self-compatibility experiment, or an
unambiguous numeric/table value, after deterministic checks or review. Assign
**Medium** to an explicit species treatment or prose statement whose controlled
value is unambiguous. BHL OCR cannot be High. WFO supplies no quality level
because it is only the name-resolution bridge.

Provider copies sharing the same original DOI/BHL item/source work count as one
lineage. DOI takes precedence over PMCID, Plazi UUID, BHL page ID, and provider
URL for deduplication.

## 5. Small local pilot

Use `unresolved_after_bulk.csv.gz` from integrated run `30107066525` and select
300 distinct accepted species:

- 150 unresolved for reproductive assurance;
- 100 unresolved for floral structural complexity; and
- 50 unresolved for flower colour.

Prefer species with two or three unresolved axes, but do not select by family or
use family traits. Run:

1. WFO local exact join for all 300; write accepted names and every exact synonym.
2. Plazi exact treatment search for all names; fetch each unique treatment UUID.
3. Europe PMC OA search for the 150 reproductive species and their synonyms;
   fetch each unique PMCID once, including supplements when present.
4. BHL lookup for at most 50 Plazi/Europe-PMC misses after an API key is supplied.
5. Extract candidates, review/promote them under the contract, then compare
   species-axis keys against the fixed integrated ledger.

Freeze:

```text
direct_web_aliases.csv
direct_web_candidates.csv
direct_web_candidate_audit.csv
direct_web_errors.jsonl
direct_web_source_manifest.json
raw/{provider}/{source_record_id}.{xml|json|zip}
```

The manifest records endpoint and query, WFO release DOI and checksums, HTTP
status/attempts, raw SHA-256, license, original DOI/PMCID/treatment UUID/BHL IDs,
and extraction-code version. Report candidates, reviewed High/Medium promotions,
conflicts, lineage duplicates, and pure new species-axis gain separately.

## 6. Live feasibility checks on the current unresolved ledger

- `Erythranthe guttata` is still unresolved for reproductive assurance.
  WFO 2026-06 maps `Mimulus guttatus` (`wfo-0000450805`) as a synonym of
  `Erythranthe guttata` (`wfo-0001333299`). Searching Europe PMC under the
  synonym returned eight OA records, including a progeny-array study explicitly
  estimating outcrossing rates
  ([article](https://doi.org/10.3389/fpls.2024.1411868),
  [OA XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11617154/fullTextXML)).
  This is a concrete High-candidate path that accepted-name-only search misses.
- `Abarema cochliacarpos` is unresolved on all three axes. Exact Plazi search
  returned treatment `03A5194AFF97FFF5B0E22FDF5E59F94E`, sourced from
  [`10.11646/phytotaxa.601.1.3`](https://doi.org/10.11646/phytotaxa.601.1.3).
  Its focal treatment contains species-specific corolla, stamen, and staminal-tube
  measurements, making it a concrete structural-trait candidate after ontology
  mapping and review.

These checks show that the highest-leverage change is not another generic web
wave: it is WFO exact-synonym expansion followed by source-native XML retrieval
from Plazi and Europe PMC.
