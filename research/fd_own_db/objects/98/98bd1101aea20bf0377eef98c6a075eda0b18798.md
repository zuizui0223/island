# Shared-cohort pollinator specificity — pre-outcome extension, 2026-11-25

This is a support audit and a proposed next analysis, not a positive ecological result. It does not change the frozen primary bat contrast in PR #68 or the original Bombus result in PR #67.

## Existing support ceiling

PR #64's historical public Direct flower-geometry bridge has 130 tropical bee-linked islands, 32 provisional bird-linked, 58 positive butterfly-record-linked, and only 8 with all three sources. With >=50 floral names the three-source cohort has 4; in the single-state scope it has 8 and 1 respectively. Adding a bat range and diet-completeness condition can only reduce those intersections. An all-four-guild regression cannot be treated as a globally supported ecological test.

## Executable census

`scripts/multiguild_shared_support.py` reads the exact signed public bee/GSHHG bridge and prior bird/butterfly artifacts; it joins newly computed potential bat overlap on the same public geometry IDs. It reports each cohort's count, 10-degree blocks with within-block variation in every relevant predictor, and number of informative islands, for both floral scopes, n>=5/50, bat overlap >=10/50%, area >=20/100 km2.

The pre-outcome gate is >=30 informative blocks and >=200 informative islands. The census does not fit any model even if it passes. It must retain unqualified cohorts. Bird identities remain provisional; butterfly records are positive-only; bat range presence is potential occurrence. No unlisted animal island is coded as a verified zero. The bee geometry match alone is not evidence of a positive native checklist.

## Next estimand: trait specificity rather than total richness

For a sufficiently supported pair of guilds on the SAME cohort, define focal trait endpoints before viewing paired coefficients. Estimate a within-block, group-equal model of each floral category frequency with both guild predictors plus area, isolation, climate and trait documentation. Scale each log1p guild predictor by its common-cohort within-block SD to avoid comparing arbitrary richness units. A meaningful paired specificity contrast is

    S = (b[guild A, trait A] - b[guild A, trait B])
      - (b[guild B, trait A] - b[guild B, trait B]).

For bats vs specialist nectar birds, candidate endpoints are white membership and tubular membership; for Bombus vs other native bees, bilabiate/papilionaceous vs open radial membership. These are coarse recorded state-set proxies, not actual measured pollination syndromes. The same floral carrier set must define both outcomes. Preserve all counterexamples. A paired group bootstrap must resample the same islands/groups across all four coefficients; a plant-genus-by-evidence null is required in protected GIFT before any claim of within-lineage trait matching. Group/evidence limitations may render the contrasts unidentifiable.

Guild complementarity is a separate question: a guild-by-guild interaction or joint category diversity above a fixed-marginal/lineage-preserving null. It does not follow from a positive S, multiple positive single-guild slopes, or total pollinator richness. No such result currently exists.

## Geometry repair

Initial bat run 45537853357 decoded 1,428 map names and diet source, and passed the synthetic statistical test, but stopped on a wide polygon. The old bounding-box rejection was too broad. The correction checks raw ring longitude jumps, unwraps genuine crossing candidates assuming short arcs, aligns holes, clips to longitude bands, and checks degree-area conservation; wide no-jump polygons retain their topology. Synthetic seam/holes, wide and ordinary polygon tests precede data processing. Nonclosing/ambiguous rings fail closed. The repaired simplification cannot recover the original high-resolution geography; the full-resolution product and small-island sensitivity remain future source checks. This correction is outcome-blind and not an ecological result.
