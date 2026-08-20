# Targeted support=2 Wave 15 checkpoint

Wave 15 is an acquisition checkpoint, not a separate Validated Low implementation.
It adds exactly 100 individually reviewed `accepted_species × trait_name` rows:

- 83 previously acquired but not formally promoted regional-flora rows from five
  immutable GitHub Actions artifacts;
- 17 newly reviewed primary, flora, university, or professional botanical rows.

The 100 rows target 20 `genus × axis × trait_name` rules that had two agreeing
direct species in Run 31800810018.  No genus value is written by this package.
The shared PR #131 implementation alone recomputes rules and applies dominance,
leave-one-species-out, and leave-one-source-lineage-out validation.

Every row keeps the original URL, exact quote, retrieval date, content hash,
source lineage, name match, and a one-to-one audit decision.  Search snippets,
family inference, global fallback, `n=2` rules, and axis-only joins are excluded.
Pollen-vector and reward traits remain independent and are not mapped to the
floral-structure axis.

Rebuild the package with:

```powershell
python -m island_v2.targeted_support2_wave15_checkpoint
```

The source-artifact snapshot records Run IDs, artifact IDs, artifact names, and
the exact candidate rows used.  Queue potentials are deliberately not reported
as coverage gains; only the formal all-evidence artifact may establish realized
direct and Validated Low increments.
