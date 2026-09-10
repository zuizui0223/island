# Database 1.0 Validated-Low provenance closure

This note records a rights/provenance audit of the immutable Chapter 1 Database 1.0 snapshot. It does **not** change trait values, quality tiers, H1-H5 analyses, or the frozen scientific database.

The final snapshot contains 65,015 resolved Low cells carrying `validated-low:*` derived lineages, representing 3,070 distinct genus × trait × inferred-state rules. The audit first reconstructs current/near-current rules from the hash-verified Wave53 and Run329 artifacts, then completes the append-only historical layer from immutable Wave34, Wave35, Wave39, Wave40, Wave41, Wave42, and Wave43 artifacts.

The release gate requires exact rule-state agreement for semantic recovery and non-empty upstream direct-source provenance. One historical lineage is intentionally retained as a semantic exception: `validated-low:Piptatherum:flower_primary_color:["blue_purple","green_brown_inconspicuous"]`. Its three BaseFlor support lineages are recoverable, but the verified Wave34 frontier and Wave33 secondary material support the historical rule state `[`blue_purple`]`, not the two-state set stored in Database 1.0. The audit therefore preserves source provenance while refusing to relabel the frozen database or claim an exact rule-semantic match.

Expected audit invariants are:

- distinct Validated-Low rule lineages: 3,070;
- exact rule semantics recovered: 3,069 / 3,070;
- upstream source provenance recovered: 3,070 / 3,070;
- semantic mismatches: 1;
- target Low cells with unresolved upstream source provenance: 0;
- target Low cells with zero support lineage: 0;
- scientific database modified: false;
- redistribution rights granted by provenance recovery: false.

Redistribution is evaluated separately. Recovering an upstream source lineage does not itself grant permission to redistribute the derived cell. The Zenodo release gate must evaluate the rights of every upstream source family on which a derived Validated-Low cell depends.
