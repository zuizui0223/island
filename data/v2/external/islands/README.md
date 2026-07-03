# Island polygon source

Place the source vector used to define the island analytical universe here before the first v2 GBIF run.

The initial v2 configuration expects the original `GlbIslands.gdb` source and `BigIslands` layer, because this is the polygon universe used to generate the v1 exploratory sample. v2 does not reuse v1's output: it re-reads the original polygons, repairs and normalizes them, records a source manifest, and queries GBIF using exact geometry.

A replacement source is acceptable only when its release, license, island definition, attributes, and area threshold are explicitly documented in the v2 manifest.
