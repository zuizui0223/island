# PR138 Saint-Pierre source audit: island occurrence versus archipelago origin

Date checked: 2026-08-26

## Decision

The Saint-Pierre route is **admissible only as a partial, rule-based crosswalk**. It does not license an exhaustive native/introduced partition of the 311 PR138 flora rows.

- An exact occurrence on **Saint-Pierre Island** joined to a taxon-level VASCAN assertion of `introduced` for the whole Saint-Pierre-et-Miquelon region can be admitted as `introduced`, provided the VASCAN record is `present`, the accepted-name reconciliation is unique, and the occurrence is not on another island such as Ile aux Marins, Langlade, or Miquelon.
- The same join with VASCAN `native` **cannot** be admitted as island-native. Regional nativity does not show that the taxon reached Saint-Pierre Island without human agency.
- A Le Hors account can support `native` or `introduced` only where its prose explicitly and unambiguously attaches that origin class to the taxon and independently places the taxon on Saint-Pierre Island.
- The unmarked checklist in Le Hors cannot be treated as native. The manuscript gives archipelago totals but does not mark origin or island locality in the checklist itself.

The best scalable occurrence source found is the live NatureSPM taxon index, which reports island names on taxon pages. It can provide exact Saint-Pierre occurrence, but it does not itself expose a native/introduced field and no reuse licence was found. VASCAN remains the machine-readable origin source, at archipelago rather than island resolution.

## 1. Historical flora: source, access, and semantics

### 1.1 What is actually online

The accessible document is Mathurin Le Hors, *La Flore des iles Saint-Pierre et Miquelon*, an undated manuscript inferred by the transcriber to have been written in 1947-1950. The [GrandColombier source page](https://grandcolombier.com/2008/10/23/1947-1950-la-flore-des-iles-saint-pierre-et-miquelon-par-mathurin-le-hors/) labels it unpublished and links the [direct PDF](https://www.grandcolombier.com/wp-content/uploads/pdf/190109/flore_spm_le_hors_derniere_version.pdf).

This is **not a facsimile of the original manuscript**. Roger Etcheberry's preface says he had a photocopy supplied by the author's son and `"decide de le retaper"` (PDF p. 2). It also states that accents were added, obvious typographical errors were corrected, scientific names were italicized, and the author's manuscript annotations were retained. Therefore it is a close, source-provenanced transcription of primary content, but an editorially mediated copy. Historical spelling and status scope require manual QA.

Access verified on 2026-08-26:

- HTTPS PDF: 53 pages, text-extractable, 212,428 bytes.
- SHA-256: `CE2F660B1C33C577F81E2029F59E4265A6B6987B80BE6847E44081630FD72F78`.
- An identical-byte mirror was found at [guilligomarch.com](https://www.guilligomarch.com/wp-content/uploads/2019/11/mathurinlehorsflore.pdf), but GrandColombier is the linked source page used here.
- No explicit open-data or redistribution licence was found in the PDF or source page. Cite and retrieve from the source URL; do not vendor the PDF into the repository without separate permission review.

### 1.2 Origin convention

The general remarks state an archipelago total of `"594 ... dont 88 sont introduites et 506 indigenes"` (PDF p. 3). That demonstrates that Le Hors distinguished origin classes, but it does **not** define a symbol or residual rule for the later checklist.

The document has two different components:

1. **Narrative accounts, PDF pp. 4-36.** Many taxa receive locality prose, and some receive explicit origin prose. These claims may be extracted taxon by taxon.
2. **Checklist, PDF pp. 37-53.** It is titled as a list of plants of the islands and prints taxon names by family. It contains no visible native/introduced marker and no island-locality column. Its unmarked names cannot reproduce the reported 506/88 split.

The aggregate totals therefore do not create an exhaustive implicit-native convention. Unlike the Block Island appendix, there is no statement that every unmarked checklist entry is indigenous.

### 1.3 What the prose can support

The narrative sometimes supplies both status and target-island occurrence in one taxon account. For example, *Trifolium agrarium* is described as `"introduit en 1928 ... a Saint-Pierre"` (PDF p. 25). On PDF p. 36, three named *Hieracium* immediately preceding `"recemment introduites"` are placed at Cap a l'Aigle or Savoyard, both Saint-Pierre localities. These are direct candidate introduced assertions after modern-name reconciliation.

The prose also contains explicit native groupings. On PDF p. 21, *Ranunculus acris* is separated as naturalized from Europe and the following buttercups are called indigenous; some of those following accounts have Saint-Pierre localities such as Savoyard. Such group-scoped statements can support native candidates only after a deterministic scope parser plus manual review confirms exactly which names the phrase governs.

Important locality constraints:

- Accept explicit `St-Pierre` / `Saint-Pierre` or a named site verified inside the target Saint-Pierre polygon.
- Do not equate the Saint-Pierre commune or archipelago with the target island.
- Exclude exact records on Ile aux Marins, Grand Colombier, Langlade, and Miquelon. Le Hors names these separately.
- `dans les trois iles`, `dans le pays`, and an unqualified checklist presence are archipelago occurrence, not Saint-Pierre occurrence.
- Historical occurrence does not by itself establish present-day persistence. It may qualify the origin join only when the PR138 flora row already has an independently admitted contemporary island occurrence basis.

## 2. VASCAN: scope, fields, endpoints, and licence

### 2.1 Geographic and status semantics

VASCAN's [official About page](https://data.canadensys.net/vascan/about/#distribution) defines:

- `Native`: present through natural processes only, without human agency.
- `Introduced`: established/naturalized outside its original range through deliberate or accidental human activity.
- `Ephemeral`: recurring but not permanently established.
- `Excluded`, `Extirpated`, `Doubtful`, and `Absent`: separate non-presence or uncertainty classes.

The [official IPT metadata](https://data.canadensys.net/ipt/resource.do?r=vascan) says distribution status is reported **per region**, and lists `Saint Pierre and Miquelon` as one covered region. It does not split Saint-Pierre, Miquelon, or Langlade. Its geographic scope is the 242-km2 collectivity, not Saint-Pierre Island.

Consequently:

- `introduced` for the containing region plus exact Saint-Pierre occurrence is a defensible introduced assignment under the PR138 admission rule.
- `native` for the containing region is only archipelago/regional nativity and cannot establish island-native status.
- `ephemeral`, `excluded`, `doubtful`, and `absent` must remain distinct from established binary presence.

### 2.2 Programmatic identifiers: `PM`, not `SM`

The existing inventory's instruction to use `SM` is unsafe for programmatic access.

- The checklist table displays the human-facing column label `SM`.
- The checklist query parameter is `province=PM`, as shown by VASCAN's own [Saint-Pierre-et-Miquelon example](https://data.canadensys.net/vascan/checklist?combination=anyof&habit=tree&province=PM&status=doubtful&status=ephemeral&status=excluded&status=extirpated&status=introduced&status=native&taxon=0).
- The Search API returns `locality: "PM"` and `locationID: "ISO 3166-2:FR-PM"`.
- The Darwin Core Archive uses `locality: "Saint Pierre and Miquelon"`, `countryCode: "FR"`, and `locationID: "ISO3166-2:FR-PM"`.

Any implementation should key on `locationID` after normalizing whitespace in the ISO prefix, not on the UI abbreviation.

### 2.3 Stable machine-readable routes

1. **Search API:** [documentation](https://data.canadensys.net/vascan/api), with versioned endpoint `https://data.canadensys.net/vascan/api/0.1/search.json`. GET accepts one name; POST accepts at most 200 newline-separated names. The response supplies `acceptedNameUsageID`, synonym/accepted status, distribution `locationID`, `locality`, `establishmentMeans`, and `occurrenceStatus`.
2. **Versioned full DwC-A:** [IPT resource](https://data.canadensys.net/ipt/resource.do?r=vascan), verified archive `https://data.canadensys.net/ipt/archive.do?r=vascan&v=37.17`, and versioned metadata `https://data.canadensys.net/ipt/eml.do?r=vascan&v=37.17`. Version 37.17 was published 2026-08-04. Its SHA-256 on 2026-08-26 was `B613DDDA9E790CD7698760546F044420FC98D59BE77BD50F0465242601515334`.
3. **Checklist builder:** supports filtered text and DwC-A creation, but the generated file page states that files remain on the server for only 24 hours. It is reproducible input, not a stable artefact URL.
4. **OpenRefine reconciliation:** `https://data.canadensys.net/vascan/refine/0.1/reconcile` is documented by VASCAN, but exact accepted-name joins should remain the admission path.

In DwC-A v37.17, `meta.xml` maps the distribution extension to `locationID`, `locality`, `countryCode`, `occurrenceStatus`, `establishmentMeans`, `source`, and `occurrenceRemarks`. The Saint-Pierre-et-Miquelon rows inspected in that frozen version were:

| `establishmentMeans` | `occurrenceStatus` in DwC-A | Rows |
|---|---:|---:|
| `native` | `present` | 473 |
| `introduced` | `present` | 191 |
| `introduced` | `irregular` | 12 |
| empty | `excluded` | 69 |
| empty | `doubtful` | 21 |
| empty | `absent` (API: `extirpated`) | 8 |

These total 774 taxon rows (598 species, 114 subspecies, and 62 varieties). For the same non-permanent class, the Search API exposes `occurrenceStatus: "ephemeral"` where the DwC-A uses `irregular`. The eight DwC-A `absent` rows are `extirpated` in the API; they are not ordinary unreported/null absences. Preserve the source interface and version in provenance, and do not silently merge any of these classes with established `present` taxa.

Use `establishmentMeans` for origin. Do not infer origin from the DwC-A `occurrenceStatus`, which is `present`/`irregular` rather than `native`/`introduced`.

The distribution extension's first column is the core identifier: join `distribution.txt:id` to `taxon.txt:id`, not to the separate Darwin Core `taxonID` field. Of the 774 `PM` distribution rows, 691 cite NatureSPM in the row-level `source`; 23 have an empty `source`. The blank-source subset includes 10 binary-native and 4 binary-introduced rows. If the admission contract requires a traceable row-level underlying citation in addition to the VASCAN assertion itself, those 14 binary candidates must remain `source_missing` pending review.

Historical synonyms require a second lookup. For example, querying *Hieracium aurantiacum* returns a synonym assertion with accepted usage *Pilosella aurantiaca* but no distribution on the synonym object; querying the returned `acceptedNameUsageID` supplies the `PM` introduced assertion. Admission requires a unique accepted concept and no genus-only or fuzzy fallback.

### 2.4 Rights

The VASCAN [About page](https://data.canadensys.net/vascan/about/#rights) and IPT metadata release the data under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). The IPT assigns DOI [10.5886/zw3aqw](https://doi.org/10.5886/zw3aqw). The fetched archive should be retained by version and SHA-256 because the live API is continuously updated.

## 3. Exact Saint-Pierre occurrence sources

### 3.1 NatureSPM legacy flora portal

VASCAN v37.17 cites `Nature SPM. 2010+. Portail de la biodiversite des Iles. Plantes vasculaires` as the source for most Saint-Pierre-et-Miquelon distribution assertions. The live [vascular-plant index](https://www.naturespm.com/xxxcl.php?code=flecl&debut=0&mv=2&page=0) is server-rendered HTML with stable-looking taxon links (`xxxobs.php?code=flecl&...&numero=...`). Its last page reports `numtot=795` checklist records.

Taxon pages explicitly separate islands. For example:

- [*Equisetum arvense*](https://www.naturespm.com/xxxobs.php?code=flecl&mv=2&numero=1264622) reports observations at St-Pierre, Miquelon, and Langlade and gives its first archipelago observation at `St-Pierre Savoyard`.
- [*Lycopodium annotinum*](https://www.naturespm.com/xxxobs.php?code=flecl&mv=2&numero=1222422) reports observations on all three islands while giving the first observation specifically at Miquelon. The summary island list, not only the first-observation locality, is therefore the relevant occurrence field.
- Taxon pages can omit Saint-Pierre and list only Miquelon/Langlade, demonstrating that the island summary is discriminatory rather than an archipelago boilerplate.

This route is machine-parseable HTML but is not a documented API or versioned download. No explicit data licence or terms of reuse were found on the index, taxon pages, or home page. Use URLs, retrieval date, raw-page hashes, and parsed field provenance; do not treat it as a redistributable open dataset without permission clarification.

### 3.2 Current Collectivite portal

The Collectivite's current [SPM Boreal portal](https://www.spm-patrimoine-naturel.fr/) exposes first-party JSON endpoints used by its web application, including `https://www.spm-patrimoine-naturel.fr/api/observations/search` and `/api/observations/count`. They returned public JSON on 2026-08-26. This is an observation supplement rather than an exhaustive vascular-flora checklist; it should not replace the NatureSPM taxon index or convert non-detection to absence. Its licence and stable API contract were not verified in this audit.

## 4. Join/admission contract

| Exact Saint-Pierre occurrence | Origin assertion | Binary admission | Reason |
|---|---|---|---|
| Yes | VASCAN `introduced` + `present` at `ISO3166-2:FR-PM` | **Introduced: yes** | The containing region is explicitly outside the taxon's original range; exact target-island occurrence supplies locality. |
| Yes | VASCAN `introduced` + `irregular` / API `ephemeral` | **No, keep qualified** | Not permanently established; preserve as `ephemeral_introduced` unless the analysis contract explicitly admits ephemeral taxa. |
| Yes | VASCAN `native` at `ISO3166-2:FR-PM` | **Native: no** | Native somewhere in the archipelago does not establish natural arrival on Saint-Pierre Island. |
| Yes | Le Hors explicitly says the same taxon is introduced and places it on Saint-Pierre | **Introduced: yes** | Direct historical island-status assertion, subject to name and locality QA. |
| Yes | Le Hors explicitly says the same taxon/group is indigenous and places it on Saint-Pierre | **Native: yes, reviewed subset only** | The source itself links origin and target locality; group scope must be unambiguous and manually audited. |
| Yes | Le Hors checklist name only, or aggregate 506/88 totals | **No** | No taxon-level status mapping or residual convention. |
| No / archipelago only | Any regional status | **No island status** | The target-island occurrence premise is missing. |

Minimum provenance fields for every admitted row:

- original queried name;
- accepted VASCAN taxon ID/name and taxonomic status;
- occurrence source URL, retrieval date, raw hash, and exact island field/locality;
- origin source/version, `locationID`, `establishmentMeans`, `occurrenceStatus`, and source citation;
- `status_basis` such as `pm_introduced_present_plus_exact_saint_pierre_occurrence` or `le_hors_direct_island_native_reviewed`;
- confidence/review flag and exclusion reason for every unresolved candidate.

## 5. Blockers and claim ceiling

1. The accessible Le Hors document is a 2009 transcription, not the original manuscript image.
2. The complete checklist has neither an origin marker nor an island-locality column; only prose subsets are directly usable.
3. VASCAN is archipelago-wide. It can license an introduced subset, not an island-native residual.
4. The existing `SM` acquisition note is programmatically wrong/stale: use `PM` / `ISO3166-2:FR-PM`.
5. Historical synonyms require accepted-concept traversal; synonym hits can lack distribution fields.
6. NatureSPM supplies island occurrence but has no verified versioned export or reuse licence.
7. Ile aux Marins and other named islands must be excluded even when administratively associated with Saint-Pierre.
8. Historical occurrence does not prove contemporary persistence; independent PR138 flora occurrence evidence remains necessary.

**Claim ceiling:** this route can add a defensible exact-island **introduced subset** through NatureSPM/Le Hors occurrence plus VASCAN `introduced`, and a smaller manually reviewed native subset from direct Le Hors prose. It cannot fill all unmatched Saint-Pierre rows as native, cannot infer absence, and cannot by itself produce an exhaustive binary island checklist.
