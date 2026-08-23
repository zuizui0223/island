# Minimum endemic-colour manual review checkpoint — 2026-08-23

## Purpose

This checkpoint reviews only the frozen 23-species minimum acquisition queue from run `32614335237`. It does not widen the queue, weaken the evidence contract, or modify the canonical PR #132 direct-evidence ledger.

The target remains island model support, not global trait completion:

- northern endemic colour: 31 -> 50 supported islands, maximum queue 19 species;
- tropical endemic colour: 46 -> 50 supported islands, maximum queue 4 species;
- stop a regional campaign immediately once 50 supported islands are reached.

## Automated acquisition checkpoint

`Acquire minimum endemic colour from public Web` completed successfully but produced no reviewable flower-colour candidates from the frozen queue. The zero yield is treated as a provider/retrieval result, not as biological missingness.

Protologue/IPNI+BHL and GBIF specimen-label workflows are wired to the same frozen queue, but their acquisition jobs are manual-dispatch only. PR-triggered runs validate the implementation and intentionally skip the expensive acquisition step.

## Manual source review

The complete row-level audit is frozen at:

`data/v2/staging/traits/endemic_minimum_acquisition/manual_colour_review_20260823.csv`

Current decisions:

### Canonical-safe direct records

1. **Astragalus kawakamii** — Flora of Mikawa gives the species-direct statement `花は淡黄色。` (pale yellow). Proposed canonical value: `yellow_orange`; binary class: non-plain.
2. **Viola cephalonica** — the official Mt Aenos flora book describes the flowers as bluish-violet. Proposed canonical value: `blue_purple`; binary class: non-plain.

These are review-accepted candidates in this checkpoint only. They are not counted as canonical support until promoted through the PR #132 evidence ledger and rebuilt there.

### Direct text present but canonical four-colour value unresolved

**Cirsium umezawanum** — Kadota's original species description states that the corolla is pale reddish purple. This is certainly non-plain for the broad contrast, but the current four-colour ontology does not uniquely determine whether the canonical value should be `red_pink` or `blue_purple`. The raw direct statement is retained; no forced canonical category is emitted.

### Direct evidence rejected for the broad binary support target

- **Fritillaria rhodia** — original description: yellowish-green. This crosses `yellow_orange` (non-plain) and `green_brown_inconspicuous` (plain), so it cannot unlock the plain/non-plain outcome without an arbitrary choice.
- **Dracaena fernaldii** — perianth described as greenish-yellow / yellowish-green; same cross-class ambiguity.
- **Mammillaria evermanniana** — descriptions combine yellowish-cream/greenish tones with pink-brown markings; this is multicoloured and crosses the binary boundary.
- **Sasa jotanii** — the primary paper describes a reddish-purple *spikelet*, not an explicitly attached flower/perianth colour. Rejected as a subject-attachment mismatch.

### Leads requiring stronger species-direct primary evidence

`Odontarrhena lesbiaca`, `Cirsium tanegashimense`, `Centaurea xylobasis`, `Cirsium ishizuchiense`, and `Erysimum naxense` have secondary or indirect colour indications, but they remain unpromoted until a source satisfies the same direct-evidence contract.

All other queued species remain unresolved rather than inferred from congeners, images alone, or generic floral syndromes.

## Immediate support implication

If the two canonical-safe northern records are independently promoted through PR #132, the frozen queue contract predicts two new northern endemic-colour islands, moving raw support from **31 to 33**, before any other accepted records are considered.

Tropical support remains **46** because no tropical queue record currently has an ontology-safe plain/non-plain value under this review.

The `Cirsium umezawanum` statement is deliberately not included in that +2 count. A separate binary-only evidence layer could in principle use it, but introducing such a layer solely to increase support would change the evidence contract and is not done here.

## Continue / stop rule

1. Promote only independently reviewed species-direct records through the canonical evidence path.
2. Rebuild status-stratified direct support after every accepted record.
3. For tropical endemic colour, continue only within the existing four-species queue; do not replace failed targets with globally convenient species unless a new support-gap calculation generates a new frozen queue.
4. For northern endemic colour, prioritize unresolved queue members that extend the direct-supported isolation envelope.
5. Stop at 50 supported islands; do not continue toward global trait completeness.

## Scientific boundary

This acquisition step is a support diagnostic. It cannot rescue the already falsified broad universal syndrome or the unsupported Bombus-hypervolume mechanism. Any new records only improve precision for the endemic-colour stratum under the already frozen status/lineage framework.
