# Harrisia source correction

The 45-cell unpromoted genus frontier included five Harrisia SC candidates.
Source review found that the High statement for H. portoricensis reports
partial self-compatibility, not an unqualified SC state. The original paper
is DOI 10.3732/ajb.0900026, PMID 21622342. Its indexed original abstract was
retrieved successfully through Europe PMC's core metadata API on 2026-09-07.
The publisher body returned 403 and PubMed returned 429; neither is claimed
to have supplied full text. Exact short supporting excerpts, endpoint and
review provenance are recorded in the correction JSON.

## Changes and safeguards

- Replace this species' self_incompatibility SC with mixed_or_variable,
  retaining High quality and a partial-self-compatibility classification.
- Autonomous_selfing_capacity=absent already existed as Medium, through
  Razanajatovo's compilation. Upgrade the existing trait to High using the
  original controlled-pollination report; do not insert another trait row.
- Collapse the PMID URL and the compilation's Rojas-Sandoval / Melendez-Ackerman
  2009 original reference to the same DOI. Do not count them as independent
  experiments. Conflicting lower-tier raw records remain in baseline history.
- Include both now-High traits in this species' axis composition. SC/SI and
  autonomous selfing remain separate individual traits.
- Preserve every other species-axis row exactly, all quality counts, the
  complete denominator and all existing Low cells. No Harrisia reproductive
  Low is present in the baseline. Its colour Low is unrelated and untouched.
- Source records remain disabled for automatic genus retraining. After the
  value correction, Harrisia SC dominance is 6/7, below both 0.95 and 0.90;
  the five unpromoted candidates are held, not counted as Low invalidations.

This changes one existing direct value and upgrades one direct trait. It adds
zero species-trait rows and zero filled axes. Cumulative restart gain remains
+35 relative to Run 34093932899, not +40 or +45. Goal +1,000 remains active.

The optional correction in the restart integration workflow has exact
species, prior-value, prior-lineage, prior-quality and absence-of-duplicate
preconditions. It refuses a second application or changed source receipt.
Generic direct admission remains add-only; this is not a blanket overwrite
permission. Tests cover idempotence rejection, changed lineage, forbidden
genus training and unchanged axis quality. Formal output includes the source
correction receipt and separate value-correction / tier-upgrade counters.
