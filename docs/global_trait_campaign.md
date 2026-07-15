# Global trait campaign

## Scope

The active v2 acquisition strategy is **global species first**, not a regional
island pilot. The current exact-island GBIF snapshot supplies the unique species
master. A scheduled campaign then advances bounded, family-balanced waves without
using island floral outcomes, Bombus detections, model results, or a preferred
region to choose species.

The Izu records already present in the repository are ordinary members of this
global universe. They are not the campaign's priority stratum and do not define
the method.

## Ordered tasks

1. Retrieve flower colour, form, symmetry, size and access traits from Wikimedia
   species text for every eligible flowering-plant species.
2. Retrieve explicit mating-system, autonomous-selfing and
   self-incompatibility statements from Wikimedia species text.
3. Repeat the reproductive target layer against OpenAlex title/abstract text.
4. Alternative-pollinator functional-guild passes are optional, deferred
   interaction evidence and are not required for primary trait completion.

The flower-trait pass covers the full eligible master before the campaign spends
effort on optional interaction evidence. The active completion fields are flower
colour, flower form/symmetry, mating system and
self-incompatibility. Existing functional-guild evidence is retained, but new
guild retrieval cannot block these fields; flower colour is never used to infer
a bird, bee, moth, or other pollinator.

## Resumption and failure rules

`campaign_ledger.csv.gz` stores one row per global-master species and independent
status/attempt columns for every source task. Each scheduled run selects the next
family-balanced batch, writes a new immutable wave directory, and updates the
ledger. New species entering the master are added without resetting completed
rows.

A successful lookup with no matching statement is a **zero-hit lookup**, not a
biological absence. Transient source failures are retried. After the configured
attempt limit they become `exhausted` audit rows so a few inaccessible pages do
not freeze the entire global campaign.

## Evidence boundary

All output is machine-only candidate evidence with source URL, citation, excerpt,
evidence scope, matched term, and method provenance. It does not:

- accept a trait value;
- decide native, introduced, or cultivated status;
- define an island flora denominator;
- classify Bombus presence or absence;
- enter a confirmatory model without the later evidence-review contract.

## Operation

The production workflow `.github/workflows/run-public-web-trait-shards.yml`
advances 128 resumable all-species shards. The campaign CLI also supports
bounded family-balanced waves. Both routes treat a zero-hit lookup as an audit
result, not as evidence that a trait is absent.


## Per-species streaming scheduler

The scheduler advances bounded family-balanced or hash-stable shard batches.
Flower colour and form/symmetry are attempted for every eligible species;
reproductive enrichment proceeds independently. Optional alternative-guild
work is deferred and never gates primary completion. A forced task selects only
a measurement stage, never a region, island, trait outcome, or Bombus result.


## One-week global primary screen

The primary path is intentionally database-light and resumable. Flower colour,
form and symmetry are screened for the entire 106,295-species angiosperm
master, while OpenAlex remains an optional scholarly-enrichment lane. The
completion target refers to a recorded first-pass attempt for every species,
not to inventing a value when public evidence is absent.

## Free recovery-rate improvement

The primary Wikimedia lane uses multilingual fallback without paid APIs. English
Wikipedia is queried first. If that pass yields no target candidate for a
species, the campaign tries the configured non-English sitelinks from Wikidata
(`es`, `fr`, `de`, `pt`, `it` by default). This keeps the high-throughput
English path cheap while recovering explicit statements that exist only in other
Wikipedia editions. These rows remain machine-only candidates with source URL,
excerpt, language-specific citation, and no curated trait decision.

## Completion hierarchy at 100,000-species scale

The production order is bulk direct datasets first, strict accepted-name and
exact-synonym reconciliation second, then resumable public Flora/plant-site text
retrieval. Source-locked LLM extraction may parse retrieved statements but does
not replace a source. SC/SI gaps may additionally receive a separate genus
inference only when at least five directly supported congeners are unanimous;
the output is low-confidence `likely_SC`/`likely_SI` and never overwrites a
species-level record. Unsupported cells remain `unknown`.

## Landed reproductive evidence (2026-07-15)

Two pinned open datasets now add 3,305 distinct species-trait cells before
deduplication against earlier sources: 2,509 direct SC/SI species from Meyer,
Galloway & Eckert (2026), plus 920 conservative reproductive cells across 585
species from Razanajatovo et al. (2016). After overlap with the previous ledger,
3,218 cells are new direct evidence. A separately labelled unanimous-genus rule
adds 2,353 further new `likely_SC`/`likely_SI` cells; leave-one-out agreement is
372/384 (96.9%). Direct and inferred evidence remain distinguishable in every
export.
