# Review-only tropical endemic colour packet — 2026-08-22

## Purpose

This packet is **not production evidence**. It contains four independent species-direct colour candidates selected from the corrected endemic-status support audit because tropical corroborated-endemic colour currently has 46 directly supported islands and needs only four additional islands to reach the operational 50-island model-support target.

The four candidates belong to four different currently unsupported tropical islands. If all four pass independent review and the common formal rebuild, the direct endemic-colour island count can move from 46 to 50 without a broad trait campaign.

Production guardrails remain unchanged:

- exact/synonym-direct evidence only;
- source text must explicitly attach colour to the focal species' flowers/corolla/petals;
- taxonomic identity must be safe;
- ambiguous multi-colour states must not be forced into a binary class;
- this packet never directly promotes a row;
- independent manual review is still required.

## Candidate 1 — `Limonium lobinii`

- target island: `gshhg_2.3.7_h_e9684321ab46eed33c40`
- current analysis role: tropical, corroborated endemic, colour-unsupported island
- suggested raw colour: bluish-purple corolla
- suggested repository value: `blue_purple`
- plain-colour contrast side: conspicuous / non-plain
- evidence type: species description / taxonomic paper
- review priority: **A**
- source: Erben & Leyens, *Limonium lobinii ... a new species from the Cape Verde Islands*
- source URL: https://www.researchgate.net/publication/255982604_Limonium_lobinii_Plumbaginaceae_a_new_species_from_the_Cape_Verde_Islands_West_Africa
- short source excerpt: `Corolla bluish-purple, petals 5, free`
- taxonomic identity check: POWO accepts `Limonium lobinii N.Kilian & Leyens` and restricts its native range to Cape Verde (Santiago).
- identity URL: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:989650-1
- caveat: none identified from the source text; reviewer should still verify the excerpt in the original paper/PDF rather than a search snippet.

## Candidate 2 — `Psychotria silhouettae`

- target island: `gshhg_2.3.7_h_4aa93ec0fc19c6708842`
- current analysis role: tropical, corroborated endemic, colour-unsupported island
- suggested raw colour: white
- suggested repository value: `white`
- plain-colour contrast side: plain
- evidence type: conservation action plan quoting the species description
- review priority: **A/B**
- source: Baguette et al., *Species Action Plan – Psychotria silhouettae* (2018)
- source URL: https://www.researchgate.net/publication/327500273_Species_Action_Plan-Psychotria_silhouettae
- short source excerpt: `tiny, white and tubular flowers to 3 mm across`
- independent provenance: BGCI lists the *Psychotria silhouettae* conservation action plan as an associated conservation resource.
- provenance URL: https://www.bgci.org/resources/bgci-tools-and-resources/conservation-action-plans/
- caveat: the action plan explicitly notes that diagnostic characters separating `P. silhouettae` from `P. pervillei` are not fully convincing and calls for further taxonomic work. The colour observation is species-direct, but the taxonomic uncertainty must remain attached to any promoted row.

## Candidate 3 — `Rondeletia anguillensis`

- target island: `gshhg_2.3.7_h_02aab13d659d614914f8`
- current analysis role: tropical, corroborated endemic, colour-unsupported island
- suggested raw colour: pink
- suggested repository value: `red_pink`
- plain-colour contrast side: conspicuous / non-plain
- evidence type: regional botanical account based on the discovery/description history
- review priority: **B**
- source: Don Mitchell, *The Geology and Botany of Anguilla*
- source URL: https://www.aahsanguilla.com/uploads/7/3/7/1/7371196/1._the_geology_and_botany.pdf
- short source excerpt: `its flowers are smaller, and pink rather than white`
- taxonomic identity check: POWO accepts `Rondeletia anguillensis Howard & Kellogg`, first published in *Journal of the Arnold Arboretum* 68:127 (1987).
- identity URL: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:284359-2
- caveat: the colour source is a regional botanical chapter rather than the protologue. Before promotion, search the Howard & Kellogg treatment / herbarium-backed description for a higher-tier species-direct colour statement.

## Candidate 4 — `Scalesia microcephala`

- target island: `gshhg_2.3.7_h_8961f6f3a9eff7528643`
- current analysis role: tropical, corroborated endemic, colour-unsupported island
- suggested raw colour: white
- suggested repository value: `white`
- plain-colour contrast side: plain
- evidence type: species-specific ecological thesis/report citing the Galápagos flora literature
- review priority: **B/C**
- source URL: https://edepot.wur.nl/208411
- short source excerpt: `heads ... contain 15–25 individual white flowers`
- taxonomic identity check: POWO accepts `Scalesia microcephala B.L.Rob.` with native range Galápagos.
- identity URL: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:228722-2
- caveat: source tier is weaker than the first three candidates. Prefer McMullen / Eliasson / another primary Galápagos flora treatment if a species-direct colour statement can be obtained. Do not substitute a genus-level statement that `Scalesia` flowers are white.

## Held / rejected candidates from the same unsupported-island set

### `Dracaena fernaldii`

Species-direct descriptions report greenish-yellow / yellowish-green flowers, but these span repository colour categories that sit on opposite sides of the current plain-versus-conspicuous binary. Keep as raw multistate evidence if independently reviewed; **do not use it to unlock the binary colour island by forcing one state**.

### `Euphorbia origanoides`

Available descriptions include creamy-yellow/white flower structures. The wording is not clean enough for the binary contrast without resolving which floral organ the colour refers to. Hold rather than force a value.

### Generic `Bidens`, `Cyrtandra`, `Mollugo`, or `Scalesia` descriptions

Genus-level flower-colour statements are discovery leads only. They are not species-direct evidence under the PR #132 production contract.

## Expected support effect

The four review candidates above occupy four distinct currently unsupported tropical corrected-endemic islands. The current direct support is 46 islands. Therefore:

```text
accepted candidates 0 -> 46 supported islands
accepted candidates 1 -> at most 47
accepted candidates 2 -> at most 48
accepted candidates 3 -> at most 49
accepted candidates 4 -> at most 50
```

This is an island-support calculation, not a guarantee of formal promotion. Any accepted row still passes the shared evidence rebuild, identity checks, conflict precedence and final support audit.
