# Minimum endemic-colour manual review checkpoint — 2026-08-23

## Purpose

This checkpoint reviews only the frozen 23-species minimum acquisition queue from run `32614335237`. It does not widen the queue or weaken the evidence contract.

The target remains island model support, not global trait completion:

- northern endemic colour: 31 -> 50 supported islands, maximum queue 19 species;
- tropical endemic colour: 46 -> 50 supported islands, maximum queue 4 species;
- stop a regional campaign immediately once 50 supported islands are reached.

## Automated acquisition checkpoint

`Acquire minimum endemic colour from public Web` completed successfully but produced no reviewable flower-colour candidates from the frozen queue. The zero yield is treated as a provider/retrieval result, not as biological missingness.

Protologue/IPNI+BHL and GBIF specimen-label workflows are wired to the same frozen queue, but their acquisition jobs are manual-dispatch only. PR-triggered runs validate the implementation and intentionally skip the expensive acquisition step.

## Manual and multilingual source review

The complete row-level manual audit is frozen at:

`data/v2/staging/traits/endemic_minimum_acquisition/manual_colour_review_20260823.csv`

A separate PR #132 multilingual checkpoint already contained two independently reviewed Greek species descriptions (`Limonium aphroditae`, `Viola cephalonica`). The current wave24 formal delta unifies those reviewed results with the independently reviewed `Astragalus kawakamii` record so the three additions are rebuilt together from the same prior formal evidence artifact.

### Canonical-safe direct records in the unified delta

1. **Astragalus kawakamii** — Flora of Mikawa gives the species-direct statement `花は淡黄色。` (pale yellow). Canonical value: `yellow_orange`; binary class: non-plain.
2. **Limonium aphroditae** — the Red Data Book of Rare and Threatened Plants of Greece, surfaced as a species description through GBIF, describes the corolla as light violet (`με ανοικτό ιώδες χρώμα`). Canonical value: `blue_purple`; binary class: non-plain.
3. **Viola cephalonica** — independent species-level Greek/English flora descriptions describe the flowers as violet to blue-violet / bluish-violet. Canonical value: `blue_purple`; binary class: non-plain.

All three are individually reviewed species-direct records. The PR #132 wave24 promotion workflow validates candidate identity, ontology value, source tier, provenance, excerpt hash, source-lineage fingerprint and manual-audit agreement before rebuilding the common all-evidence ledger.

The earlier two-record wave24 run (`32637583153`; Astragalus + Viola) completed successfully: both records entered the direct ledger as medium/resolved species-direct evidence and strict flower-colour coverage increased by exactly two cells with no Validated-Low addition/invalidation. This is a workflow validation checkpoint only; the final acquisition result must use the unified three-record run.

### Direct text present but canonical four-colour value unresolved

**Cirsium umezawanum** — Kadota's original species description states that the corolla is pale reddish purple. This is certainly non-plain for the broad contrast, but the current four-colour ontology does not uniquely determine whether the canonical value should be `red_pink` or `blue_purple`. The raw direct statement is retained; no forced canonical category is emitted.

### Direct evidence rejected for the broad binary support target

- **Fritillaria rhodia** — original description: yellowish-green. This crosses `yellow_orange` (non-plain) and `green_brown_inconspicuous` (plain), so it cannot unlock the plain/non-plain outcome without an arbitrary choice.
- **Dracaena fernaldii** — official species descriptions give greenish-yellow / yellowish-green petals or perianth; same cross-class ambiguity.
- **Mammillaria evermanniana** — descriptions combine yellowish-cream/greenish tones with pink-brown markings; this is multicoloured and crosses the binary boundary.
- **Sasa jotanii** — the primary paper describes a reddish-purple *spikelet*, not an explicitly attached flower/perianth colour. Rejected as a subject-attachment mismatch.

### Leads requiring stronger species-direct primary evidence

`Odontarrhena lesbiaca`, `Cirsium tanegashimense`, `Centaurea xylobasis`, `Cirsium ishizuchiense`, and `Erysimum naxense` have secondary or indirect colour indications, but they remain unpromoted until a source satisfies the same direct-evidence contract.

All other queued species remain unresolved rather than inferred from congeners, images alone, or generic floral syndromes.

## Immediate support implication

If the unified three-record wave24 delta completes the formal rebuild unchanged, the frozen queue contract predicts three distinct new northern endemic-colour islands, moving direct support from **31 to 34** and leaving a gap of **16** to the 50-island target.

Tropical support remains **46** because no tropical queue record currently has an ontology-safe plain/non-plain value. This is now clearly a source/evidence-availability bottleneck rather than a reason to widen the queue.

The `Cirsium umezawanum` statement is deliberately not included in the +3 count. A separate binary-only evidence layer could use it in principle, but introducing such a layer solely to increase support would change the evidence contract and is not done here.

## Continue / stop rule

1. Treat the unified three-record formal artifact as the only wave24 result; do not sum independently rebuilt artifacts that share the same prior formal ledger.
2. Rebuild status-stratified direct support after every accepted record.
3. For tropical endemic colour, continue only within the existing four-species queue; do not replace failed targets with globally convenient species unless a new support-gap calculation generates a new frozen queue.
4. For northern endemic colour, prioritize unresolved queue members that extend the direct-supported isolation envelope.
5. Stop at 50 supported islands; do not continue toward global trait completeness.

## Scientific boundary

This acquisition step is a support diagnostic. It cannot rescue the already falsified broad universal syndrome or the unsupported Bombus-hypervolume mechanism. Any new records only improve precision for the endemic-colour stratum under the already frozen status/lineage framework.
