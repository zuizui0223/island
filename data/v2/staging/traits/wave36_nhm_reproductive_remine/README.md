# Wave36: NHM descriptions and cached reproductive remine

Wave36 adds reviewed species-direct evidence without refetching completed
public-Web tasks.  It combines two bounded sources:

- Natural History Museum's *Angiosperm functional traits extracted from
  botanical descriptions* (`doi:10.5519/p3dm31kc`), queried across 5,802
  genera that remained unresolved for colour or structure.  The source-specific
  review accepted 28 species x trait cells; four structure candidates remained
  excluded after all-evidence direct conflict/ontology resolution, leaving 24
  direct species x trait and 23 direct species x axis gains.
  Synonym candidates were checked against the pinned WFO Plant List June 2026
  snapshot (`doi:10.5281/zenodo.20782718`, MD5
  `02F989B01B8EB142EC5934BD634B3876`) and GBIF strict species matching.  The
  identity-gate audit records six accepted mappings and nine rejected
  homonym/hybrid mappings.
- A full-text remine of 6,586 already cached Useful Plants/PFAF pages for
  explicit reproductive statements.  Fifteen species x trait candidates were
  manually reviewed, ten were accepted and five apomictic, cultivar-conflicted,
  or otherwise nonsexual claims were excluded.  Provider mirrors and statements
  copied across species from the same cited work share one source lineage.

The all-evidence rebuild keeps `genus x axis x trait_name` throughout and
requires both leave-one-species-out and leave-one-source-lineage-out validation.
Only two new rules become eligible: `Osteomeles x floral_form` and
`Ulex x tube_depth_class`.  They add six Low cells.  No family inference or
global fallback is included, and `reward_type`, `pollen_vector_mode`, and sex
system fields are not collapsed into the strict three axes.

Against formal Wave35 Run `33072865120`, the lossless overlay adds exactly 39
species x axis cells with zero loss: structure +20, colour +9, reproduction
+10.  Of these, 33 are new direct Medium cells and six are trait-specific
Validated Low.  Final coverage is 219,733 / 318,885, leaving 99,152 unresolved.

The formal workflow regenerates both reviewed source packages from the frozen
review tables, checks their hashes, downloads the pinned Wave35 artifact, and
materializes the Wave36 overlay.
