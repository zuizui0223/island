# Strict open-Web trait acquisition

This lane searches broadly but accepts evidence narrowly. It supplements the
fixed GBIF, Wikimedia, World Flora, OpenAlex, Europe PMC, eFloras, and bulk
providers. It does not replace their immutable artifacts or rescan completed
expanded-acquisition ranges.

## Evidence flow

```text
prioritized unresolved species x trait queue
  -> multilingual exact-name / synonym / vernacular queries
  -> Google Custom Search JSON API or Brave Search API
  -> discovery URL (title/snippet are not evidence)
  -> approved domain registry
  -> robots-aware source-page fetch
  -> exact accepted-name or two-backbone synonym identity gate
  -> rules-only trait extraction with an exact page quote
  -> source-lineage and mirror deduplication
  -> 100-page independent audit
  -> accepted species-direct evidence
  -> strict coverage and all-direct Validated Low rebuild
```

The search-backend seam has three adapters:

- `GoogleCustomSearchBackend`, configured only through
  `GOOGLE_CUSTOM_SEARCH_API_KEY` and `GOOGLE_CUSTOM_SEARCH_CX`;
- `BraveSearchBackend`, configured only through `BRAVE_SEARCH_API_KEY`; and
- `FixtureSearchBackend`, used by tests and reproducible dry runs.

Keys are never written to query rows, logs, manifests, artifacts, or the
repository. The current repository has no Brave/Google search secret, so the
official-search phase records a zero-query blocked status while the approved
domain-registry pilot remains runnable.

## Fixed denominator and evidence hierarchy

Coverage remains fixed at 106,295 accepted angiosperm species and three axes,
or 318,885 species-axis cells. Evidence is selected in this order:

1. species/synonym-direct High;
2. species/synonym-direct Medium;
3. validated genus-level Low.

Family inference and global fallback are excluded. Missing values remain
missing, not biological absence. The 12 traits in
`config/open_web_trait_acquisition.yml` remain separate until the axis ledger
is calculated. In particular:

- `self_fertile` may support `self_incompatibility=SC`, but never
  `autonomous_selfing_capacity`;
- pollen vector, SI/SC, autonomous selfing, mating system, and cleistogamy are
  distinct traits;
- Tier C nursery/horticultural pages can contribute only visible floral
  morphology;
- multicolour statements remain `multicolored_variable` with a state-set field;
- cultivar-only and horticultural-hybrid pages are rejected.

## Domain registry

`config/open_web_domain_registry.yml` stores:

- domain and source name;
- source, region, language, and organization types;
- robots/access notes and allowed page pattern;
- scientific-name and cultivar flags;
- provenance quality and recommended tier;
- allowed traits, pilot success rate, and false-positive audit status.

An unknown domain is discovery-only. It must be independently reviewed and
added to the registry before any page on it can emit a candidate. A domain's
page pattern is also enforced, so an approved organization's search, index, or
blog page cannot silently enter as a species treatment.

The first registry inventory is Flora of Mikawa (`mikawanoyasou.org`). Its
adapter is generic: inventory URLs are followed to a bounded depth and only
URLs matching the registered species-page pattern enter the page pipeline.
There is no source-specific extraction code. A registry-level
`treatment_end_pattern` stops extraction before a page's appended genus
catalogue; other regional floras can declare their own boundary without a
Python adapter.

## Species identity gate

A page passes only when:

- a binomial is present in page text/title;
- it exactly matches the fixed accepted-species denominator, or an exact
  synonym on which at least two named backbone snapshots agree;
- the page subject and resolved rank are species (not `subsp.`, `var.`, or
  `f.`), and the master family is not contradicted;
- the page is not cultivar/hybrid-only;
- the registered stable page pattern matches; and
- the supporting excerpt is present in the fetched source page.

Fuzzy, ambiguous, one-backbone synonym, family-conflict, genus-only, and
cultivar matches are audit rows. WFO/WCVP, POWO/WCVP, and GBIF exports can be
normalized into a synonym CSV with columns
`synonym,accepted_species,family,backbone`; the resolver requires two distinct
`backbone` values to agree.

## Multilingual query and extraction policy

English, botanical Latin, Japanese, Chinese, Spanish, Portuguese, French, and
German query terms are versioned in the acquisition config. Query generation
uses accepted names first, then exact synonyms and vernacular names. A separate
PDF query is available for dissertations, reports, full text, and
supplementary tables.
Synonym inputs use `synonym,accepted_species,backbone`; a query synonym is
emitted only when at least two named backbone snapshots agree. Optional local
names use `accepted_species,vernacular_name,language`. Missing snapshots are
visible in task composition and are never replaced by fuzzy matching.

The rules-only extractor guards against negation and leaf/fruit/seed colour
contamination. Flower and corolla measurements are retained in the exact quote
and deterministically mapped to size/depth classes. LLM extraction is not
needed for the pilot. Any future LLM adapter must receive only fetched page
text, reproduce an exact quote found in that text, use fail-closed ontology
mapping, and freeze prompt/model/hash manifests.

Individual flower form and inflorescence display stay separate. In particular,
`頭花` / `capitulum` maps only to `inflorescence_display`, not
`floral_form`; ray and tubular florets can remain a multistate floral-form
composition. Measurements are accepted only within a bounded window after an
explicit flower, corolla, perianth/petal, or tube term. Dimensions of leaves,
fruits, pedicels, inflorescences, bracts, lemmas, receptacles, and cultivar-only
descriptions are excluded; tube depth is not silently reused as whole-flower
size.

## Pilot

The pilot workflow:

1. downloads the formal integrated-coverage artifact;
2. independently verifies 106,295 x 3 and recalculates High/Medium/Low,
   axes, unresolved cells, providers, domains, languages, source types, and
   exclusion reasons;
3. builds a 300-species information-value queue, prioritizing a
   genus-trait with two agreeing direct species and many unresolved congeners;
4. creates multilingual official-search tasks;
5. uses an official search adapter only when its secret exists;
6. runs a bounded approved-domain inventory pilot;
7. when search credentials are absent, crawls the registered inventories of
   at least five independently approved domains (currently Mikawa, Botanic.jp,
   Tsukuba Botanical Garden, PlantNET NSW, Go Botany, and NZPCN);
8. freezes a reproducibly hashed, domain × trait-stratified,
   200-unique-page audit sheet; and
9. promotes no evidence until reviewed decisions are supplied.

Applying the audit does not contain a second Validated Low implementation.
Reviewed rows are converted to the shared evidence schema, then PR #131's
`all_evidence_trait_audit` lineage deduplication, conflict resolution,
`build_rule_audit`, `apply_genus_rules`, coverage, and acquisition-queue
functions run directly. Rules are applied by `genus × trait_name`; axis is
reporting metadata and is never the inference join. `reward_type` and
`pollen_vector_mode` remain queryable ledgers but never fill the strict flower
structure axis.

Run it with:

```bash
gh workflow run open-web-trait-pilot.yml \
  --repo zuizui0223/island \
  --ref <branch> \
  -f baseline_run_id=30107066525 \
  -f pilot_species=300 \
  -f registry_domain=mikawanoyasou.org \
  -f registry_max_pages=600 \
  -f max_search_queries=1200 \
  -f search_backend=auto
```

Fill `manual_audit_multidomain.csv` with:

- `decision=accept|reject`;
- `species_identity_correct`;
- `value_correct`;
- `provenance_complete`;
- `cultivar_contamination`;
- false-positive reason and reviewer metadata.

Only `accept` rows whose three accuracy fields are true and cultivar field is
false count as `accepted_correct`; precision is exactly
`accepted_correct / all reviewed`. Production is gated on at least 200
reviewed pages from at least five domains, overall precision at least 0.95,
and cultivar contamination no more than 0.02. Within a passing pilot, each
domain × trait scope must itself have at least 10 reviews, precision at least
0.95, and cultivar contamination no more than 0.02. Unreviewed or weak scopes
remain blocked.

## Production and resumption

`open-web-trait-production.yml` refuses to start unless the referenced reviewed
artifact satisfies that gate and the selected official-search credential
exists. Each query has a stable task ID. Before a shard runs, all available
prior checkpoint manifests on the branch are loaded; a different run whose
half-open range overlaps the request is rejected, and every completed task ID
is removed. An exact incomplete range may be resumed only by naming its source
run. A completed range cannot be resumed or searched again.

Without Brave or Google CSE credentials, the six-domain workflow is the honest
credential-free acquisition lane. Its query count and cost remain zero; it
does not label an inventory crawl as an official search.

Production outputs remain review-only until their own audit/promotion step.
They include search tasks and costs, URL discovery, page fetches, species
identity outcomes, evidence candidates, source-lineage duplicates, and
provider/domain/language summaries.

Successive reviewed source-package waves must chain from the immediately
preceding formal Web artifact by setting `prior_public_web_run_id`,
`prior_public_web_artifact_name`, and `prior_public_web_file_name`. This keeps
all previously promoted direct rows in the append-only ledger while the
completed-page manifests prevent reacquisition.

## Analysis contract

- Confirmatory: source-lineage-deduplicated species/synonym-direct evidence,
  primarily Tier A/B, exact quote present, conflicts resolved or represented as
  multistate.
- Secondary robustness: audited Tier C morphology, image-derived evidence, and
  validated/probabilistic genus inference with uncertainty propagation.
- Sensitivity only: family inference, global fallback, weak personal pages,
  and diagnostic two-species genus rules.

No missing cell is converted to zero. A species-level imputation is shared
across islands within an imputation draw.
