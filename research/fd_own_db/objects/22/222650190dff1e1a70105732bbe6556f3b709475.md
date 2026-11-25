# N1d: fixed 1,000/2,000 km island donor reference for floral-distance SES

**2026-11-09; post-hoc falsification of geographic reach asymmetry.**

The validated N1c nearest-10-islands floral analysis has a major imbalance:
its tenth nearest eligible donor is a median ~387 km away for continental
focal islands but ~1,149 km for oceanic islands (All). A fixed number
of neighbors does not fix the reference's absolute geographic reach.

## Frozen, explicitly exploratory design

The reference is the pinned 496 distinct recorded colour/form state-set
profiles. The distance is ORIGINAL D, equal-weight mean set-Jaccard
colour/form dissimilarity, not a PCoA projection or Gaussian volume.

1. **Absolute donor radii 1,000 km and 2,000 km** are selected after the
   N1c result to bracket the two reach distributions. They are post-hoc
   sensitivities, not preregistered cutpoints. Do not scan additional radii
   to find a preferred geological result.
2. For every geology-typed eligible focal island, select ALL OTHER
   eligible islands in the same original broad analysis regime within
   that radius, by audited centroid great-circle distance. Include
   islands irrespective of their own geological class.
3. Require at least THREE other source islands and enough distinct
   recorded floral profiles among them to draw n distinct profiles.
   Count each profile's occurrence across donor islands, not within-island
   taxon-name carrier abundance. Do not fabricate unrecorded profiles or
   sample with replacement to fill insufficient support.
4. For full support, set n to the focal island's number of recorded
   distinct colour/form profiles. For equal30, uniformly thin the
   observed island to 30 distinct types 99 times and use its average
   observed original-distance score. Compare with 999 fixed-30 reference
   draws. These repeated thinnings are highly correlated.
5. Evaluate two correlated endpoints: mean off-diagonal pairwise
   distance and mean distance to the nearest OTHER recorded profile.
   Use SES=(observed-minus-null-mean)/null sample SD plus the
   unstandardized original-D observed-minus-null distance.
6. Recompute N1b regional island-frequency reference in this run,
   and compare N1b with each absolute-radius N1d ONLY for identical
   focal islands with both references calculable. This eliminates
   apparent contrast changes due purely to different eligible targets.

## Outputs and inference constraints

- fixed_radius_ses_summary.csv: descriptive N1b/N1d group metrics,
  suppress effects for fewer than five islands in a cell
- fixed_radius_ses_support.csv: eligible islands, numbers with >=3 other
  island donors, full and 30-profile support, donor counts and reference
  richness, with small-cell suppression for descriptive quantities
- fixed_radius_matched_reference.csv: paired same-island N1b versus
  N1d SES differences and raw D differences by radius/scope/region/origin
- manifest.json: fixed design and SHA256, no individual island, donor,
  taxon, genus or geography source identifiers

**The fixed radius does not equalize donor density.** A given radius may
include far more recorded neighboring floras for one geological class.
The support output explicitly audits this disparity. Neither N1b nor
N1d conditions on the focal genus composition or mainland native
source opportunity. Nearby donors may share identical taxa with focal
islands; geographic overlap is not independent lineage history.

No new geological-origin hypothesis test, origin p value or multiplicity
claim is made. Two overlapping evidence scopes (All and Direct) and
post-hoc radii cannot establish source-independent evolutionary
adaptation, functional pollination niche shrinkage or actual absences.

## Public-source replay

    python -m fd.fixed_radius_ses \
      --input work/island_public \
      --g0 work/fd_g0 \
      --crosswalk work/fd_crosswalk/island_origin_crosswalk.csv \
      --output results/fixed_radius_ses \
      --draws 999 --thinnings 99

If oceanic proximity approaches zero under one or both identical
geographic radii, the previous negative SES is sensitive to
reference geography. If it remains negative, it is still
observational and needs the private taxon/genus-heldout N2 for
mechanistic source/lineage conditioning.
