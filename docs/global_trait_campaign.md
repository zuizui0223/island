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

1. Retrieve explicit pollen-vector and reproductive-assurance statements from
   Wikimedia species text.
2. Repeat the same target layer against OpenAlex title/abstract text.
3. Only after both global reproductive tasks terminate, retrieve flower access
   traits for species carrying a direct machine candidate for biotic pollen
   vector.
4. After that layer terminates, make a focused OpenAlex pass for alternative
   pollinator functional guilds.

The order separates the upstream question—whether a species depends on a biotic,
wind, or mixed pollen vector—from self-compatibility and autonomous selfing, and
only then asks how accessible or specialized the flowers are. Functional-guild
claims remain direct text evidence; flower colour is never used to infer a bird,
bee, moth, or other pollinator.

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

The workflow `.github/workflows/drive-global-trait-campaign.yml` runs hourly and
advances one bounded wave. A manual dispatch can override the default batch size
between 1 and 500. The campaign automatically moves between tasks when all
eligible rows in the preceding task are terminal.


## Per-species streaming scheduler

The scheduled workflow advances a bounded bundle every 30 minutes: a broad
family-balanced Wikimedia wave followed by smaller OpenAlex, floral-access,
and alternative-guild waves. Dependencies are evaluated per species, not as
a global phase barrier. A species can therefore reach floral-access and
alternative-pollinator extraction as soon as its own upstream rows are
terminal, while the remaining global master continues through the first
source. The forced task option selects only a measurement stage, never a
region, island, trait outcome, or Bombus result.


## One-week global primary screen

The scheduled primary path is intentionally database-light: 1,000 globally
family-balanced species are screened through Wikimedia every 30 minutes, for a
theoretical capacity of 48,000 species per day. Floral-access and alternative
guild Wikimedia extraction run only for direct machine biotic-vector candidates.
OpenAlex remains an optional scholarly-enrichment lane and is not required for
primary completion. The one-week target refers to complete machine first-pass
coverage, not manual acceptance of every candidate trait.

## Free recovery-rate improvement

The primary Wikimedia lane uses multilingual fallback without paid APIs. English
Wikipedia is queried first. If that pass yields no target candidate for a
species, the campaign tries the configured non-English sitelinks from Wikidata
(`es`, `fr`, `de`, `pt`, `it` by default). This keeps the high-throughput
English path cheap while recovering explicit statements that exist only in other
Wikipedia editions. These rows remain machine-only candidates with source URL,
excerpt, language-specific citation, and no curated trait decision.
